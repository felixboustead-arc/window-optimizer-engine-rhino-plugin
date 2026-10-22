# -*- coding: utf-8 -*-
"""
glide.py - G-family planners: the window GLIDES along with the thin-beam band instead of hopping between grid poses.

Physics (measured, MODELS.md): at any minute the poses that give a beam of the target width W form a thin, almost
straight band in (ring, aperture) space just above the light edge, r*(a) ~ r0(t) + s(t) a (s ~ 0.06 deg per %).
Inside the band the strip lands in the same place (set by the sun), so the band position along a is free: the planner
uses it to move as little and as smoothly as possible.

G1 band glide, per day (days are independent - the sun sets in between):
  1. band table: every SAMPLE minutes (every minute where the band appears / disappears) and at aperture knots
     A_KNOTS, the ring r*(a) where the main room patch is exactly W wide (exact physics, bisection) + clean flag.
  2. DP over the aperture lattice (1 % steps) per minute: utility = clean beam (1) + QW x quality - target cost;
     ring = r*(a, t) interpolated; transition cost = CV x ((da/100)^2 + (dr/15)^2) (speed^2: a glide is cheap, a hop
     is dear) with |da| <= RATE per minute.
  3. glides (smootherstep) across minutes with no band and across the night.
  4. exact verification of every beam minute; a minute that is not clean is repaired by an exact band search at that
     minute (same aperture).

    python glide.py <model id> <wtag> [json overrides of the module constants]
"""
from __future__ import division

import array
import json
import math
import multiprocessing as mp
import os
import pickle
import sys
import time

HERE = os.path.dirname(os.path.abspath(__file__))
if HERE not in sys.path:
    sys.path.insert(0, HERE)
import labcore as L  # noqa: E402
import physics as PH  # noqa: E402
import predict as PR  # noqa: E402

C = L.C
CFG = {"A_KNOTS": [0, 25, 50, 75, 100], "SAMPLE": 5, "QW": 0.05, "CV": 200.0, "RATE": 3, "NPROC": 6,
       "SMOOTH": 9, "REPAIR_MAX_DR": 0.25, "W_MULT": 1.0, "TARGET": None, "TW": 0.0, "TCAP": 3.0,
       "NIGHT_GLIDE": 60,
       # G1T light target: TARGET_SPEC {depth_m, offset_m, height_m, strength}; WTOL None = width fixed (OFF, default)
       # or [lo, hi] multipliers of the set width (ON, e.g. [0.5, 2.0]); WPEN per octave off the set width
       "WTOL": None, "WPEN": 0.05, "CVM": 2.0, "MRATE": 1}
_PH = None


def _init(ctx_bytes):
    global _PH
    _PH = PH.Physics(pickle.loads(ctx_bytes))


def smootherstep(s):
    s = 0.0 if s < 0 else (1.0 if s > 1 else s)
    return s * s * s * (s * (6 * s - 15) + 10)


def _smooth(xs, win):
    """Centred moving average over win minutes; the window shrinks symmetrically at the ends (end points kept)."""
    if win <= 1 or len(xs) < 3:
        return xs
    h = win // 2
    n = len(xs)
    cs = [0.0]
    for x in xs:
        cs.append(cs[-1] + x)
    out = []
    for k in range(n):
        m = min(h, k, n - 1 - k)
        out.append((cs[k + m + 1] - cs[k - m]) / (2 * m + 1))
    return out


def _light_warp(ph, i, p0, p1, K=16):
    """Cumulative path length of the room-light centroid along the straight pose path p0 -> p1 at minute i
    (normalised 0..1), for retiming a gap glide so the LIGHT moves evenly (fast where the light barely moves)."""
    cs = []
    for k in range(K + 1):
        s = k / float(K)
        e = ph.evaluate(i, (p0[0] + (p1[0] - p0[0]) * s, p0[1] + (p1[1] - p0[1]) * s))
        cs.append(ph.light_centroid(e["faces"])[0])
    cum = [0.0]
    for a, b in zip(cs, cs[1:]):
        d = math.sqrt(sum((u - v) ** 2 for u, v in zip(a, b))) * ph.sc if (a and b) else 0.0
        cum.append(cum[-1] + d + 0.02 / K)         # small floor: the window still moves where the light does not
    return [c / cum[-1] for c in cum]


PARK_CANDS = [(r, a) for r in (9.0, 10.5, 12.0, 13.5, 15.0) for a in (30.0, 40.0, 50.0, 60.0, 70.0, 80.0, 90.0)]


def _park_pose(ph, mins):
    """G2: where to wait while no beam is possible (low sun, day edges): the pose whose sun patch moves slowest
    (m/min, mean over the sample minutes), no false beam, then the smallest lit area."""
    best = None
    for pose in PARK_CANDS:
        sp = area = 0.0
        vis = False
        for m in mins:
            e0, e1 = ph.evaluate(m, pose), ph.evaluate(m + 1, pose)
            c0, a0 = ph.light_centroid(e0["faces"])
            c1, _a1 = ph.light_centroid(e1["faces"])
            if c0 and c1:
                sp += math.sqrt(sum((u - v) ** 2 for u, v in zip(c0, c1))) * ph.sc
            area += a0
            vis = vis or e0["vis"]
        key = (sp / len(mins) + (1.0 if vis else 0.0), area)
        if best is None or key < best[0]:
            best = (key, pose)
    return best[1]


CALM_LAT = [(r, a) for r in (0.0, 3.0, 6.0, 9.0, 12.0, 15.0) for a in (0.0, 20.0, 40.0, 60.0, 80.0, 100.0)]


def _calm(ph, cfg, beam_minutes, out, rows, tpt, stats):
    """G2 CALM PATH: in stretches with no visible beam (low sun at the day's edges, gaps), the sun patch still moves
    - fast when the sun is low.  Re-plan those stretches so the ROOM LIGHT moves slowly: DP over a 6 x 6 pose lattice
    (+ the G1 pose) every CALM_DT minutes; cost = sum of (light-centroid speed)^2 x minutes + a small move cost
    + a penalty for false (visible, unclean) beams; window speed <= CALM_MAXSTEP of a swing per 5 min; fixed ends
    where a beam is visible next to the stretch.  Kept only if the stretch's light really gets calmer."""
    INF = float("inf")
    MINLEN, DT = int(cfg.get("CALM_MIN", 30)), int(cfg.get("CALM_DT", 5))
    VISPEN, CM, MAXSTEP = cfg.get("CALM_VISPEN", 0.02), cfg.get("CALM_CM", 0.2), cfg.get("CALM_MAXSTEP", 0.15)
    # M-012 options: CALM_POW = exponent of the light-speed cost (2 = as M-007; 4 = the fastest minutes count most,
    # like the tester's p99); CALM_ACCEPT "max" = keep a re-plan only if the stretch's FASTEST light minute gets
    # slower, with no extra reversal (tester rule), jerk or hop; CALM_EDGE = only the day's first / last stretch
    POW, ACC, EDGE = cfg.get("CALM_POW", 2), cfg.get("CALM_ACCEPT", "sse"), cfg.get("CALM_EDGE", False)
    bm = sorted(beam_minutes)
    stretches, cur = [], None
    for i in bm:
        if rows[i][3] != 1:
            if cur and i == cur[1] + 1:
                cur[1] = i
            else:
                cur = [i, i]
                stretches.append(cur)
        else:
            cur = None

    def cen(i, p):
        e = ph.evaluate(i, p)
        return ph.light_centroid(e["faces"])[0], (1 if e["vis"] else 0)

    def lcost(c0, c1, dt):
        if c0 is None or c1 is None:
            return 0.0 if c0 is c1 else 0.05
        d = math.sqrt(sum((u - v) ** 2 for u, v in zip(c0, c1))) * ph.sc / dt
        return 0.01 * (d / 0.1) ** POW * dt

    def mcost(p, q, dt):
        dr, da = abs(p[0] - q[0]) / 15.0, abs(p[1] - q[1]) / 100.0
        if max(dr, da) > MAXSTEP * dt / 5.0 + 1e-9:
            return INF
        return CM * (dr * dr + da * da) / dt

    def sse(seq):                                   # sum of squared per-minute light moves of rows seq
        s = 0.0
        for r0, r1 in zip(seq, seq[1:]):
            if r0[14] != "" and r1[14] != "":
                d = math.sqrt(sum((r1[k] - r0[k]) ** 2 for k in (14, 15, 16)))
                s += d * d
        return s

    def vmax(seq):                                  # fastest per-minute light move (m/min) of rows seq
        v = 0.0
        for r0, r1 in zip(seq, seq[1:]):
            if r0[14] != "" and r1[14] != "" and r0[13] > 0.005 and r1[13] > 0.005:
                v = max(v, math.sqrt(sum((r1[k] - r0[k]) ** 2 for k in (14, 15, 16))) * ph.sc)
        return v

    def shape(seq):                                 # (reversals by the tester's rule, max jerk, max step) of rows seq
        R, A = [r[1] for r in seq], [r[2] for r in seq]
        j = s = 0.0
        for k in range(1, len(seq)):
            s = max(s, abs(R[k] - R[k - 1]) / 15.0, abs(A[k] - A[k - 1]) / 100.0)
            if k >= 3:
                j = max(j, abs(R[k] - 3 * R[k - 1] + 3 * R[k - 2] - R[k - 3]) / 15.0,
                        abs(A[k] - 3 * A[k - 1] + 3 * A[k - 2] - A[k - 3]) / 100.0)
        return len(_revs(R, 0.0295)) + len(_revs(A, 0.195)), j, s

    for s0, s1 in stretches:
        if s1 - s0 + 1 < MINLEN:
            continue
        left = out.get(s0 - 1) if (s0 - 1) in rows else None
        right = out.get(s1 + 1) if (s1 + 1) in rows else None
        if EDGE and left is not None and right is not None:
            continue
        ts = list(range(s0, s1 + 1, DT))
        if ts[-1] != s1:
            ts.append(s1)
        S = [[(p,) + cen(t, p) for p in CALM_LAT + [tuple(out[t])]] for t in ts]
        cL = cen(s0 - 1, left)[0] if left else None
        V = []
        for p, c, v in S[0]:
            x = VISPEN * v
            if left:
                x += lcost(cL, c, 1) + mcost(left, p, 1)
            V.append(x)
        bps = []
        for k in range(1, len(ts)):
            dt = ts[k] - ts[k - 1]
            nv, bp = [], []
            for q, cq, vq in S[k]:
                best, bj = INF, -1
                for j, (p, cp, vp) in enumerate(S[k - 1]):
                    if V[j] < INF:
                        x = V[j] + lcost(cp, cq, dt) + mcost(p, q, dt)
                        if x < best:
                            best, bj = x, j
                nv.append(best + VISPEN * vq * dt)
                bp.append(bj)
            V = nv
            bps.append(bp)
        if right:
            cR = cen(s1 + 1, right)[0]
            V = [V[j] + lcost(S[-1][j][1], cR, 1) + mcost(S[-1][j][0], right, 1) if V[j] < INF else INF
                 for j in range(len(V))]
        j = min(range(len(V)), key=lambda x: V[x])
        if V[j] == INF:
            continue
        path = [j]
        for k in range(len(bps) - 1, -1, -1):
            j = bps[k][j]
            path.append(j)
        path.reverse()
        keys = [S[k][path[k]][0] for k in range(len(ts))]
        new = {}
        for k in range(len(ts) - 1):                # linear between samples, then a 9-min average
            a, b = ts[k], ts[k + 1]
            for m in range(a, b):
                f = (m - a) / float(b - a)
                new[m] = tuple(keys[k][z] + (keys[k + 1][z] - keys[k][z]) * f for z in (0, 1))
        new[ts[-1]] = keys[-1]
        mins = list(range(s0, s1 + 1))
        sm = [_smooth([new[m][z] for m in mins], 9) for z in (0, 1)]
        old_rows = [rows[m] for m in mins]
        new_rows, new_out = [], {}
        for n, m in enumerate(mins):
            p = (round(sm[0][n], 2), round(sm[1][n], 2))
            new_out[m] = p
            new_rows.append(ph.pred_row(m, p, tpt))
        ctxL = [rows[s0 - 1]] if (s0 - 1) in rows else []
        ctxR = [rows[s1 + 1]] if (s1 + 1) in rows else []
        if ACC == "max":
            wl = [rows[m] for m in range(s0 - 12, s0) if m in rows]
            wr = [rows[m] for m in range(s1 + 1, s1 + 13) if m in rows]
            b, c = shape(wl + old_rows + wr), shape(wl + new_rows + wr)
            good = vmax(ctxL + new_rows + ctxR) < vmax(ctxL + old_rows + ctxR) - 0.005 and \
                sum(r[3] for r in new_rows) <= sum(r[3] for r in old_rows) and \
                c[0] <= b[0] and c[1] <= max(b[1], 0.01) and c[2] <= max(b[2], 0.04)
        else:
            good = sse(ctxL + new_rows + ctxR) < sse(ctxL + old_rows + ctxR) and \
                sum(r[3] for r in new_rows) <= sum(r[3] for r in old_rows) + 2
        if good:
            for n, m in enumerate(mins):
                out[m] = new_out[m]
                rows[m] = new_rows[n]
            stats["calm_ok"] = stats.get("calm_ok", 0) + 1
        else:
            stats["calm_kept_g1"] = stats.get("calm_kept_g1", 0) + 1


def _warp_inv(cum, u):
    """Pose fraction s at which the light has covered fraction u of its path."""
    K = len(cum) - 1
    for k in range(K):
        if cum[k + 1] >= u:
            f = (u - cum[k]) / max(cum[k + 1] - cum[k], 1e-12)
            return (k + f) / K
    return 1.0


def _dwell(xs, ok=None, eps=0.15, win=12, dwell=13, iters=8):
    """R2 anti-dither: wherever the series turns back within `win` minutes (moves >= eps before and after, opposite
    signs), clip the turn into a flat dwell of >= `dwell` minutes at the turning level (a peak is cut flat from
    below, a valley from above).  ok(k, x) -> False rejects a clipped value (that minute keeps its own value)."""
    xs = list(xs)
    n = len(xs)
    for _it in range(iters):
        changed = False
        last = None
        k = 1
        while k < n:
            d = xs[k] - xs[k - 1]
            if abs(d) < eps:
                k += 1
                continue
            s = 1 if d > 0 else -1
            if last is not None and last[0] != s and k - last[1] <= win:
                sg = last[0]                               # +1: a peak between, -1: a valley
                seg = range(max(0, last[1] - 1), k)
                pk = max(seg, key=lambda j: sg * xs[j])
                a = b = pk
                lev = xs[pk]
                while b - a + 1 < dwell:
                    ca = sg * xs[a - 1] if a > 0 else None
                    cb = sg * xs[b + 1] if b < n - 1 else None
                    if ca is None and cb is None:
                        break
                    if cb is None or (ca is not None and ca >= cb):
                        a -= 1
                        lev = xs[a]
                    else:
                        b += 1
                        lev = xs[b]
                for j in range(a, b + 1):
                    v = min(xs[j], lev) if sg > 0 else max(xs[j], lev)
                    if v != xs[j] and (ok is None or ok(j, v)):
                        xs[j] = v
                        changed = True
                last = None
                k = b + 1
                continue
            last = (s, k)
            k += 1
        if not changed:
            break
    return xs


def _dipfill(xs, H, eps, accept=None):
    """R2 catch-up hold: a move against the local trend that is undone within H minutes (a one-to-few-minute spike or a
    step back that the trend overtakes again) becomes a HOLD at the level before it, so the setting never turns back.
    accept(k, v) -> False keeps minute k as it was.  Returns (new list, changed indices)."""
    xs = list(xs)
    n = len(xs)
    changed = []
    k = 1
    while k < n - 1:
        d = xs[k] - xs[k - 1]
        if abs(d) < eps:
            k += 1
            continue
        s = 1 if d > 0 else -1
        k2, dip = None, False
        for j in range(k + 1, min(n, k + H + 2)):
            if s * (xs[j] - xs[k]) <= -eps:
                dip = True
            elif dip and s * (xs[j] - xs[k]) >= 0:
                k2 = j
                break
        if k2 is None:
            k += 1
            continue
        for m in range(k + 1, k2):
            if s * (xs[m] - xs[k]) < 0 and (accept is None or accept(m, xs[k])):
                xs[m] = xs[k]
                changed.append(m)
        k = k2
    return xs, changed


def _dipfill2(xs, H, eps, accept=None):
    """R3 (M-011) catch-up hold v2: like _dipfill, but the hold is at the RUNNING level (the last level the trend reached
    before it turned back), and a dip counts against that running level.  v1 held at the level where the trend step
    began, which left a smaller turn-back whenever the trend had moved on before the dip (0.33, 0.67, [0.0 -> 0.33],
    1.67): that one-minute aperture dip was most of the E/W jerk and every W reversal.  Returns (new list, changed)."""
    xs = list(xs)
    n = len(xs)
    changed = []
    k = 1
    while k < n - 1:
        d = xs[k] - xs[k - 1]
        if abs(d) < eps:
            k += 1
            continue
        s = 1 if d > 0 else -1
        run, last_raise, dip_at, k2 = xs[k], k, None, None
        j = k + 1
        while j < n:
            v = s * (xs[j] - run)
            if dip_at is None:
                if v <= -eps:
                    dip_at = j
                elif v > 0:
                    run, last_raise = xs[j], j
                elif j - last_raise > H:
                    break
            else:
                if v >= 0:
                    k2 = j
                    break
                if j - dip_at > H:
                    break
            j += 1
        if k2 is None:
            k = (dip_at or k) + 1
            continue
        for m in range(dip_at, k2):
            lvl = xs[m - 1]
            if s * (xs[m] - lvl) < 0 and (accept is None or accept(m, lvl)):
                xs[m] = lvl
                changed.append(m)
        k = k2
    return xs, changed


def _ring_at(rrow, a):
    """Band ring at a fractional aperture a (linear between the 1 % lattice entries)."""
    lo = int(math.floor(a))
    hi = min(100, lo + 1)
    f = a - lo
    if rrow[lo] is not None and (f < 1e-9 or rrow[hi] is None):
        return rrow[lo] if f < 0.5 or rrow[hi] is None else rrow[hi]
    if rrow[lo] is None:
        return rrow[hi] if (rrow[hi] is not None and f >= 0.5) else None
    return rrow[lo] * (1 - f) + rrow[hi] * f

# ------------------------------------------------------------------------------------------------ M-012 ring hold
def _revs(xs, eps, win=10):
    """Reversals by the TESTER's rule (lab\\metrics): a move >= eps opposite to the last move >= eps within win
    minutes.  Returns [(k, j0)]: k = index of the turning move, j0 = index of the last opposite move."""
    out, last = [], None
    for k in range(1, len(xs)):
        d = xs[k] - xs[k - 1]
        if abs(d) < eps:
            continue
        s = 1 if d > 0 else -1
        if last is not None and last[0] != s and k - last[1] <= win:
            out.append((k, last[1]))
        last = (s, k)
    return out


def _hermite(x0, m0, x1, m1, n):
    """The n inner values of a cubic Hermite from x0 (slope m0 per minute) to x1 (slope m1) over n + 1 minutes."""
    T = float(n + 1)
    res = []
    for j in range(1, n + 1):
        t = j / T
        t2, t3 = t * t, t * t * t
        res.append((2 * t3 - 3 * t2 + 1) * x0 + (t3 - 2 * t2 + t) * T * m0 + (3 * t2 - 2 * t3) * x1 + (t3 - t2) * T * m1)
    return res


def _ap_for_ring(rrow, R, a0, maxda):
    """Fractional aperture nearest a0 (within maxda) where this minute's band has ring R (inverse of _ring_at)."""
    best = None
    lo, hi = max(0, int(math.floor(a0 - maxda))), min(99, int(math.ceil(a0 + maxda)))
    for b in range(lo, hi + 1):
        r1, r2 = rrow[b], rrow[b + 1]
        if r1 is None or r2 is None or r1 == r2:
            continue
        if (r1 - R) * (r2 - R) <= 0:
            a = b + (R - r1) / (r2 - r1)
            if best is None or abs(a - a0) < abs(best - a0):
                best = a
    return best


def _td_cands(R, k, j0, D, Gs):
    """Ring-hold candidates for the turn-back (k, j0) of the series R: each holds the ring >= D + 1 minutes between
    the last move one way and the first move back (the tester counts a turn-back only within 10 minutes).
    clip: cut the turn flat at the level where an 11-minute window around it starts and ends (like _dwell);
    delay: hold the level reached at j0, then catch up with the plan (cubic Hermite, G minutes);
    adv: reach the turning level G minutes early (Hermite) and hold it until the turn."""
    n = len(R)
    s_old = 1 if R[j0] - R[j0 - 1] > 0 else -1
    cands = []
    seg = range(max(0, j0 - 1), k)
    pk = max(seg, key=lambda j: s_old * R[j])
    a = b = pk
    lev = R[pk]
    while b - a + 1 < D + 1:
        ca = s_old * R[a - 1] if a > 0 else None
        cb = s_old * R[b + 1] if b < n - 1 else None
        if ca is None and cb is None:
            break
        if cb is None or (ca is not None and ca >= cb):
            a -= 1
            lev = R[a]
        else:
            b += 1
            lev = R[b]
    new = {}
    for j in range(a, b + 1):
        v = min(R[j], lev) if s_old > 0 else max(R[j], lev)
        if v != R[j]:
            new[j] = v
    if new:
        cands.append(("clip", new))
    for G in Gs:
        e = j0 + D + G
        if e + 1 < n:
            new = dict((j, R[j0]) for j in range(j0 + 1, j0 + D + 1))
            for q, v in enumerate(_hermite(R[j0], 0.0, R[e], R[e + 1] - R[e], G - 1)):
                new[j0 + D + 1 + q] = v
            cands.append(("delay%d" % G, new))
    for G in Gs:
        s0 = k - D - 1 - G
        if s0 >= 1:
            L = R[k - 1]
            new = dict((j, L) for j in range(k - D - 1, k))
            for q, v in enumerate(_hermite(R[s0], R[s0] - R[s0 - 1], L, 0.0, G - 1)):
                new[s0 + 1 + q] = v
            cands.append(("adv%d" % G, new))
    return cands


def _turn_dwell(ph, cfg, beam_minutes, out, rows, rrd, tpt, stats, RLO, RHI):
    """M-012 RING HOLD VIA THE APERTURE.  Every ring turn-back the tester would count (both ways within 10 min; on
    tilted and E/W windows at the band's daily turn, at the ring / aperture limits and next to gap glides) is
    replaced by a hold of >= 11 minutes (clip / delay / adv, _td_cands).  Where the minute has a band, the aperture
    moves the pose back onto the band at the held ring, so the beam stays (the band is a line in ring/aperture).
    A change is kept only if, over the day: fewer reversals (ring + aperture), no visible or clean minute lost, and
    in its window no faster light (max m/min, or <= TD_LV), no bigger jerk (or <= TD_J) and no hop."""
    EPS_R, EPS_A = 0.0295, 0.195
    D, Gs = int(cfg.get("TD_D", 10)), tuple(cfg.get("TD_G", (6, 12)))
    LV_OK, J_OK, MAXDA = cfg.get("TD_LV", 0.10), cfg.get("TD_J", 0.01), cfg.get("TD_MAXDA", 25.0)
    runs, cur = [], None
    for i in sorted(beam_minutes):
        if cur and i == cur[-1] + 1:
            cur.append(i)
        else:
            cur = [i]
            runs.append(cur)

    def lit(row):
        return row[13] != "" and row[13] > 0.005 and row[14] != ""

    def lv(r0, r1):
        if not (lit(r0) and lit(r1)):
            return 0.0
        return math.sqrt(sum((r1[z] - r0[z]) ** 2 for z in (14, 15, 16))) * ph.sc

    def win_stats(seq, R, A, rw, lo, hi):
        vmax = jmax = smax = 0.0
        for j in range(max(1, lo), hi + 1):
            vmax = max(vmax, lv(rw(j - 1), rw(j)))
            smax = max(smax, abs(R[j] - R[j - 1]) / 15.0, abs(A[j] - A[j - 1]) / 100.0)
            if j >= 3:
                jmax = max(jmax, abs(R[j] - 3 * R[j - 1] + 3 * R[j - 2] - R[j - 3]) / 15.0,
                           abs(A[j] - 3 * A[j - 1] + 3 * A[j - 2] - A[j - 3]) / 100.0)
        return vmax, jmax, smax

    for seq in runs:
        if len(seq) < 2 * D:
            continue
        n = len(seq)
        R = [out[m][0] for m in seq]
        A = [out[m][1] for m in seq]
        tried = set()
        for _it in range(60):
            evs = [e for e in _revs(R, EPS_R) if seq[e[0]] not in tried]
            if not evs:
                break
            k, j0 = evs[0]
            tried.add(seq[k])
            base_rev = len(_revs(R, EPS_R)) + len(_revs(A, EPS_A))
            done = False
            for name, new in _td_cands(R, k, j0, D, Gs):
                comp_moved = True
                for comp in (True, False):
                    if not comp and not comp_moved:          # same as the compensated try: already refused
                        continue
                    nR, nA, nrow = list(R), list(A), {}
                    for j, v in new.items():
                        v = round(min(RHI, max(RLO, v)), 2)
                        a = A[j]
                        if comp and rrd.get(seq[j]) is not None:
                            aa = _ap_for_ring(rrd[seq[j]], v, a, MAXDA)
                            if aa is not None:
                                a = round(min(100.0, max(0.0, aa)), 2)
                        nR[j], nA[j] = v, a
                    if comp:
                        comp_moved = any(nA[j] != A[j] for j in new)
                    if len(_revs(nR, EPS_R)) + len(_revs(nA, EPS_A)) >= base_rev:
                        continue
                    ok = True
                    vis0 = vis1 = cl0 = cl1 = 0
                    for j in new:
                        if (nR[j], nA[j]) == (R[j], A[j]):
                            continue
                        r = ph.pred_row(seq[j], (nR[j], nA[j]), tpt)
                        nrow[j] = r
                        o = rows[seq[j]]
                        vis0 += o[3]
                        cl0 += o[4]
                        vis1 += r[3]
                        cl1 += r[4]
                    if vis1 < vis0 or cl1 < cl0 or not nrow:
                        continue
                    lo, hi = max(0, min(nrow) - 3), min(n - 1, max(nrow) + 3)
                    b = win_stats(seq, R, A, lambda j: rows[seq[j]], lo, hi)
                    c = win_stats(seq, nR, nA, lambda j: nrow.get(j, rows[seq[j]]), lo, hi)
                    if c[0] > max(b[0], LV_OK) + 1e-9 or c[1] > max(b[1], J_OK) + 1e-9 or \
                            c[2] > max(b[2], 0.04) + 1e-9:
                        ok = False
                    if not ok:
                        continue
                    for j, r in nrow.items():
                        out[seq[j]] = (nR[j], nA[j])
                        rows[seq[j]] = r
                    R, A = nR, nA
                    stats["td_" + name + ("_c" if comp else "")] = stats.get("td_" + name + ("_c" if comp else ""), 0) + 1
                    done = True
                    break
                if done:
                    break
            if not done:
                stats["td_failed"] = stats.get("td_failed", 0) + 1


# ------------------------------------------------------------------------------------------------ width tolerance
def _width_plan(ph, i0, i1, W, tpt, cfg):
    """Width tolerance ON (G1T): per minute width multiplier m in [lo, hi] (cfg WTOL).  Samples the band at NL
    log-spaced widths (one aperture is enough: along the band the landing is fixed), then a smooth DP over a fine
    log-width lattice: value = clean (1) + QW x quality - TW x min(target distance, TCAP)/TCAP - WPEN x |log2 m|,
    move cost CVM x (octaves per minute)^2, at most MRATE lattice steps a minute.  Returns {minute: m}."""
    import bisect
    lo, hi = cfg["WTOL"]
    NL, NF = 5, 17
    l0, l1 = math.log(lo, 2), math.log(hi, 2)
    lv = [l0 + (l1 - l0) * k / (NL - 1) for k in range(NL)]
    lf = [l0 + (l1 - l0) * k / (NF - 1) for k in range(NF)]
    QW, TW, TCAP = cfg["QW"], cfg["TW"], cfg["TCAP"]
    WPEN, CVM, MRATE = cfg.get("WPEN", 0.05), cfg.get("CVM", 2.0), cfg.get("MRATE", 1)
    ts = list(range(i0, i1 + 1, cfg.get("WSAMPLE", 10)))
    if ts[-1] != i1:
        ts.append(i1)
    order = (50, 25, 75, 0, 100)
    warm, samp = {}, {}
    for t in ts:
        row = []
        for k, l in enumerate(lv):
            res = None
            for a in order:
                r, e = ph.band_center(t, float(a), W * 2 ** l, r0=warm.get((k, a)))
                if r is not None:
                    warm[(k, a)] = r
                    ee = ph.evaluate(t, (r, float(a)), polys=True)
                    td = ph.distance_to_lit(tpt, ee["rec"])
                    res = (1.0 if e["clean"] else 0.0, e["q"], TCAP if td is None else min(td, TCAP))
                    break
            row.append(res)
        samp[t] = row
    INF = float("inf")
    mins = list(range(i0, i1 + 1))
    U = []
    for m in mins:
        j = bisect.bisect_left(ts, m)
        s0, s1 = (m, m) if ts[j] == m else (ts[j - 1], ts[j])
        f = 0.0 if s1 == s0 else (m - s0) / float(s1 - s0)
        uk = []
        for k in range(NL):
            a0, a1 = samp[s0][k], samp[s1][k]
            if a0 is None or a1 is None:
                uk.append(None)
                continue
            cl = min(a0[0], a1[0])
            q = a0[1] * (1 - f) + a1[1] * f
            td = a0[2] * (1 - f) + a1[2] * f
            uk.append(cl * (1.0 + QW * q - TW * td / TCAP))
        if all(u is None for u in uk):
            U.append([0.0] * NF)                     # no band at any width: width is free (gap)
            continue
        urow = []
        for l in lf:
            x = (l - l0) / (l1 - l0) * (NL - 1)
            k = min(NL - 2, int(x))
            g = x - k
            ua, ub = uk[k], uk[k + 1]
            if g < 1e-9:
                u = ua
            elif ua is None or ub is None:
                u = None
            else:
                u = ua * (1 - g) + ub * g
            urow.append(-1.0 if u is None else u - WPEN * abs(l))
        U.append(urow)
    step = (l1 - l0) / (NF - 1)
    V = [-u for u in U[0]]
    bps = []
    for n in range(1, len(mins)):
        newV, bp = [INF] * NF, [0] * NF
        for s in range(NF):
            best, bs = INF, s
            for p in range(max(0, s - MRATE), min(NF - 1, s + MRATE) + 1):
                c = V[p] + CVM * ((s - p) * step) ** 2
                if c < best:
                    best, bs = c, p
            newV[s], bp[s] = best - U[n][s], bs
        lo_v = min(newV)
        V = [v - lo_v for v in newV]
        bps.append(bp)
    s = min(range(NF), key=lambda x: V[x])
    path = [s]
    for n in range(len(bps) - 1, -1, -1):
        s = bps[n][s]
        path.append(s)
    path.reverse()
    lpath = _smooth([lf[s] for s in path], int(cfg.get("SMOOTH", 0)))
    return dict((m, 2 ** l) for m, l in zip(mins, lpath))


# ------------------------------------------------------------------------------------------------ worker: one day
def _day(args):
    day, minutes, beam_minutes, cfg = args
    ph = _PH
    W = ph.P["beam_width_m"] * cfg["W_MULT"]
    knots = cfg["A_KNOTS"]
    tpt = tuple(cfg["TARGET"]["point_xyz"]) if cfg.get("TARGET") else None
    stats = {"evals_samples": 0, "repairs": 0, "repair_fail": 0, "sec": 0.0}
    RLO, RHI = float(ph.vals[0][0]), float(ph.vals[0][-1])     # the ring's range, from the GH model's parameters
    t_day = time.time()
    # ---- periods of consecutive candidate minutes
    periods = []
    for i in minutes:
        if periods and i - periods[-1][-1] <= 2:
            periods[-1].append(i)
        else:
            periods.append([i])
    plan = {}                                  # step -> (ring, aperture)
    wmap = {}                                  # step -> width multiplier (width tolerance ON only)
    poss = []
    stats["poss"] = poss
    rrd = {}                                   # minute -> band ring per 1 % aperture (for JOIN_SMOOTH)
    for per in periods:
        i0, i1 = per[0], per[-1]
        # ---- 0. width tolerance ON: choose the beam width per minute first (the landing does not depend on where
        #         along the band the window sits, only on the sun and the width)
        if cfg.get("WTOL") and tpt is not None:
            wmap.update(_width_plan(ph, i0, i1, W, tpt, cfg))
        # ---- 1. band table
        table = {}

        def sample(i, warm):
            row = {}
            for a in knots:
                r, e = ph.band_center(i, float(a), W * wmap.get(i, 1.0), r0=warm.get(a), rmin=RLO, rmax=RHI)
                if r is not None:
                    td = None
                    if tpt is not None and cfg.get("TD_IN_BAND", False):   # landing ~fixed along the band
                        ee = ph.evaluate(i, (r, float(a)), polys=True)
                        td = ph.distance_to_lit(tpt, ee["rec"])
                    row[a] = (r, bool(e["clean"]), e["q"], td)
                    warm[a] = r
            table[i] = row
            return row

        warm = {}
        ts = list(range(i0, i1 + 1, cfg["SAMPLE"]))
        if ts[-1] != i1:
            ts.append(i1)
        prev_t, prev_feas = None, None
        for t in ts:
            row = sample(t, warm)
            if cfg.get("CLEAN_REFINE"):           # G2: any knot gains / loses its band or its clean beam
                feas = tuple((a in row, a in row and row[a][1]) for a in knots)
            else:
                feas = any(v[1] for v in row.values())
            if prev_t is not None and feas != prev_feas and t - prev_t > 1:
                w2 = dict(warm)
                for tt in range(prev_t + 1, t):             # the band appears / disappears: every minute
                    sample(tt, w2)
            prev_t, prev_feas = t, feas
        st = sorted(table)
        # ---- 2. per-minute interpolated band: rr[m][a] (ring) for the a lattice, util[m][a]
        mins = list(range(i0, i1 + 1))
        alat = list(range(0, 101))
        INF = float("inf")
        rr, ut = [], []
        import bisect
        for m in mins:
            k = bisect.bisect_left(st, m)
            if st[k] == m:
                s0 = s1 = m
            else:
                s0, s1 = st[k - 1], st[k]
            f = 0.0 if s1 == s0 else (m - s0) / float(s1 - s0)
            row0, row1 = table[s0], table[s1]
            kv = {}
            for a in knots:
                if a in row0 and a in row1:
                    r = row0[a][0] * (1 - f) + row1[a][0] * f
                    cl = row0[a][1] and row1[a][1]
                    q = row0[a][2] * (1 - f) + row1[a][2] * f
                    td = None
                    if tpt is not None and row0[a][3] is not None and row1[a][3] is not None:
                        td = row0[a][3] * (1 - f) + row1[a][3] * f
                    kv[a] = (r, cl, q, td)
            if kv and any(v[1] for v in kv.values()):
                poss.append(m)                        # a clean band at some knot (planner-side "possible")
            rrow = [None] * 101
            urow = [-INF] * 101
            for a in alat:
                j = min(len(knots) - 2, max(0, bisect.bisect_right(knots, a) - 1))
                ka, kb = knots[j], knots[j + 1]
                if a == ka and ka in kv:
                    va = vb = kv[ka]
                    g = 0.0
                elif a == kb and kb in kv:
                    va = vb = kv[kb]
                    g = 0.0
                elif ka in kv and kb in kv:
                    va, vb = kv[ka], kv[kb]
                    g = (a - ka) / float(kb - ka)
                else:
                    continue
                r = va[0] * (1 - g) + vb[0] * g
                if r < RLO or r > RHI:
                    continue
                rrow[a] = r
                if va[1] and vb[1]:
                    u = 1.0 + cfg["QW"] * (va[2] * (1 - g) + vb[2] * g)
                    if tpt is not None and cfg["TW"] > 0:
                        tdv = None if (va[3] is None or vb[3] is None) else va[3] * (1 - g) + vb[3] * g
                        if tdv is not None:           # (not sampled in the band table by default)
                            u -= cfg["TW"] * min(tdv, cfg["TCAP"]) / cfg["TCAP"]
                    urow[a] = u
                else:
                    urow[a] = 0.0                     # band exists, beam not clean: allowed, worth nothing
            rr.append(rrow)
            ut.append(urow)
            rrd[m] = rrow
        # ---- 3. DP over chains of minutes that have at least one usable state
        CV, RATE = cfg["CV"], cfg["RATE"]
        chains, cur = [], []
        for idx in range(len(mins)):
            if any(x is not None for x in rr[idx]):
                cur.append(idx)
            elif cur:
                chains.append(cur)
                cur = []
        if cur:
            chains.append(cur)
        for ch in chains:
            V = [(-ut[ch[0]][a] if rr[ch[0]][a] is not None else INF) for a in alat]
            bps = []
            for n in range(1, len(ch)):
                pi, ci = ch[n - 1], ch[n]
                if ci != pi + 1:                       # (cannot happen inside a chain)
                    pass
                rp, rc, uc = rr[pi], rr[ci], ut[ci]
                newV = [INF] * 101
                bp = array.array("h", [-1]) * 101
                for a in alat:
                    if rc[a] is None:
                        continue
                    best, ba = INF, -1
                    for b in range(max(0, a - RATE), min(100, a + RATE) + 1):
                        vb = V[b]
                        if vb == INF or rp[b] is None:
                            continue
                        da = (a - b) / 100.0
                        dr = (rc[a] - rp[b]) / 15.0
                        c = vb + CV * (da * da + dr * dr)
                        if c < best:
                            best, ba = c, b
                    if ba >= 0:
                        newV[a] = best - uc[a] if uc[a] > -INF else best
                        bp[a] = ba
                lo = min(newV)
                if lo == INF:                          # dead end (should not happen): restart the chain here
                    newV = [(-uc[a] if rc[a] is not None and uc[a] > -INF else (0.0 if rc[a] is not None else INF))
                            for a in alat]
                    bp = array.array("h", [-2]) * 101
                    lo = min(newV)
                V = [x - lo for x in newV]
                bps.append(bp)
            a = min(alat, key=lambda x: V[x])
            path = [a]
            for n in range(len(bps) - 1, -1, -1):
                b = bps[n][a]
                if b == -2:                            # restart point: keep a
                    b = a
                a = b if b >= 0 else a
                path.append(a)
            path.reverse()
            # integer-% DP path -> smoothed fractional aperture (no 1 % stair steps), ring back on the band
            apath = _smooth([float(x) for x in path], int(cfg.get("SMOOTH", 0)))
            if cfg.get("CORRIDOR"):
                # R4 (M-011): where the band leaves the ring range (E/W windows: the sun crosses the window fast) only
                # a narrowing set of apertures keeps a beam, and the path was forced down / up in one minute (jerk).
                # Smooth INSIDE the feasible corridor [lowest, highest aperture with a band] each minute: alternate
                # clamp + smooth, so the path starts early and meets the limit gently.  No beam is given up.
                lo_b = [next((a for a in alat if rr[idx][a] is not None), None) for idx in ch]
                hi_b = [next((a for a in reversed(alat) if rr[idx][a] is not None), None) for idx in ch]
                if all(v is not None for v in lo_b):
                    win = int(cfg.get("SMOOTH", 9)) or 9
                    for _it in range(int(cfg.get("CORRIDOR_IT", 25))):
                        apath = _smooth([min(max(v, lo_b[n]), hi_b[n]) for n, v in enumerate(apath)], win)
                    apath = [min(max(v, lo_b[n]), hi_b[n]) for n, v in enumerate(apath)]
                    stats["corridor"] = stats.get("corridor", 0) + 1
            if cfg.get("SMOOTH2"):                     # R2: a second pass (triangular kernel): no slope kinks
                apath = _smooth(apath, int(cfg["SMOOTH2"]))
            if cfg.get("DWELL"):                       # R2: no quick back-and-forth of the aperture
                apath = _dwell(apath, ok=lambda j, v: _ring_at(rr[ch[j]], v) is not None,
                               dwell=int(cfg.get("DWELL", 13)))
            r_prev = None
            for n, idx in enumerate(ch):
                aa = apath[n]
                r = _ring_at(rr[idx], aa)
                if cfg.get("EXACT_RING") and r is not None:
                    # R2: the ring solved EXACTLY on the band at this minute and aperture (no 5-min interpolation
                    # kinks); warm start from the interpolated ring; a far-off answer (another branch) is refused
                    re_, ee_ = ph.band_center(mins[idx], aa, W * wmap.get(mins[idx], 1.0),
                                              r0=r_prev if (r_prev is not None and abs(r_prev - r) < 0.3) else r)
                    stats["exact_ring"] = stats.get("exact_ring", 0) + 1
                    if re_ is not None and abs(re_ - r) <= cfg.get("EXACT_MAXDR", 0.5):
                        r = re_
                    else:
                        stats["exact_ring_kept_interp"] = stats.get("exact_ring_kept_interp", 0) + 1
                r_prev = r
                if r is None:
                    # DP restart / dead-end minute: nearest aperture with a band this minute, else leave the
                    # minute unplanned (the gap glide fills it)
                    # R3 FB_SMOOTH (M-011): the band aperture nearest the SMOOTHED value (v1: nearest the raw DP
                    # integer, a one-minute dip of up to ~1 % that DIPFILL then had to repair)
                    c0 = aa if cfg.get("FB_SMOOTH") else path[n]
                    lo_, hi_ = int(math.floor(c0)) - RATE, int(math.ceil(c0)) + RATE
                    near = [b for b in range(max(0, lo_), min(100, hi_) + 1) if rr[idx][b] is not None]
                    if not near:
                        continue
                    b = min(near, key=lambda x: abs(x - c0))
                    aa, r = float(b), rr[idx][b]
                    stats["fallback"] = stats.get("fallback", 0) + 1
                plan[mins[idx]] = (r, aa)
    # ---- 4. glide across minutes of the day without a band (between planned minutes)
    ks = sorted(plan)
    out = dict(plan)
    for x, y in zip(ks, ks[1:]):
        if y - x > 1:
            p0, p1 = plan[x], plan[y]
            warp = _light_warp(ph, (x + y) // 2, p0, p1) if (cfg.get("GAP_RETIME") and y - x >= 4) else None
            for m in range(x + 1, y):
                e = smootherstep((m - x) / float(y - x))
                if warp:                          # G2: move the room light at an even pace, not the window
                    e = _warp_inv(warp, e)
                out[m] = (p0[0] + (p1[0] - p0[0]) * e, p0[1] + (p1[1] - p0[1]) * e)
    # ---- 4z. R4 JOIN_SMOOTH (M-011): smooth the aperture over the WHOLE day, across DP chain joints (E/W windows: a new
    #          chain started a few % away from where the last one ended -> the worst jerk), ring back on that minute's
    #          band; a change that would lose a clean beam is refused
    if cfg.get("JOIN_SMOOTH") and ks:
        seqj = list(range(ks[0], ks[-1] + 1))
        a_old = [out[i][1] for i in seqj]
        a_new = _smooth(a_old, int(cfg["JOIN_SMOOTH"]))
        bsetj = set(beam_minutes)
        for k, i in enumerate(seqj):
            v = a_new[k]
            if abs(v - a_old[k]) < 1e-4 or rrd.get(i) is None:
                continue
            r = _ring_at(rrd[i], v)
            if r is None:
                continue
            if i in bsetj and not ph.evaluate(i, (r, v))["clean"] and ph.evaluate(i, out[i])["clean"]:
                stats["join_refused"] = stats.get("join_refused", 0) + 1
                continue
            out[i] = (r, v)
            stats["join"] = stats.get("join", 0) + 1
    # ---- 4a. R2 DIPFILL: no quick turn-backs (spikes at chain joints / band-branch steps / gap glides): hold instead.
    #          Aperture first (the ring is re-solved on the band at a changed aperture), then the ring; a held pose is
    #          refused where it would lose a clean beam the plan had.
    if cfg.get("DIPFILL") and ks:
        H = int(cfg["DIPFILL"])
        seqm = list(range(ks[0], ks[-1] + 1))
        bset = set(beam_minutes)

        def keeps_beam(i, new):
            if i not in bset:
                return True
            key = "vis" if cfg.get("DIP_KEEP") == "vis" else "clean"     # what a hold may not lose
            e1 = ph.evaluate(i, new)
            if e1[key]:
                return True
            return not ph.evaluate(i, out[i])[key]

        def acc_a(k, v):
            i = seqm[k]
            r0 = out[i][0]
            if cfg.get("EXACT_RING"):
                rb, _e = ph.band_center(i, v, W * wmap.get(i, 1.0), r0=r0)
                if rb is not None and abs(rb - r0) <= cfg.get("EXACT_MAXDR", 0.5):
                    r0 = rb
            if keeps_beam(i, (r0, v)):
                out[i] = (r0, v)
                return True
            return False
        dipf = _dipfill2 if cfg.get("DIPFILL_V2") else _dipfill
        _new, ch_a = dipf([out[i][1] for i in seqm], H, cfg.get("DIP_EPS_A", 0.15), acc_a)

        def acc_r(k, v):
            i = seqm[k]
            if keeps_beam(i, (v, out[i][1])):
                out[i] = (v, out[i][1])
                return True
            return False
        _new, ch_r = dipf([out[i][0] for i in seqm], H, cfg.get("DIP_EPS_R", 0.02), acc_r)
        stats["dip_a"] = stats.get("dip_a", 0) + len(ch_a)
        stats["dip_r"] = stats.get("dip_r", 0) + len(ch_r)
    # ---- 4b. G2 PARK: before the first / after the last beam of the day (low sun, no beam possible) wait at the pose
    #          whose sun patch moves slowest, gliding in / out over PARK_GLIDE minutes (light-even pace)
    if cfg.get("PARK") and ks and beam_minutes:
        G = int(cfg.get("PARK_GLIDE", 20))
        b0, b1 = min(beam_minutes), max(beam_minutes)

        def glide(pa, pb, m0, m1):               # minutes m0..m1 inclusive go pa -> pb (excl. pa, incl. ~pb)
            wp = _light_warp(ph, (m0 + m1) // 2, pa, pb) if cfg.get("GAP_RETIME") else None
            for m in range(m0, m1 + 1):
                e = smootherstep((m - m0 + 1) / float(m1 - m0 + 2))
                e = _warp_inv(wp, e) if wp else e
                out[m] = (pa[0] + (pb[0] - pa[0]) * e, pa[1] + (pb[1] - pa[1]) * e)
        if ks[0] - b0 >= 2 * G:                  # morning
            a, b = b0, ks[0]
            park = _park_pose(ph, [a + (b - a) // 4, (a + b) // 2, b - (b - a) // 4])
            for m in range(a, b - G):
                out[m] = park
            glide(park, plan[b], b - G, b - 1)
            stats["parks"] = stats.get("parks", 0) + 1
        if b1 - ks[-1] >= 2 * G:                 # evening
            a, b = ks[-1], b1
            park = _park_pose(ph, [a + (b - a) // 4, (a + b) // 2, b - (b - a) // 4])
            glide(plan[a], park, a + 1, a + G)
            for m in range(a + G + 1, b + 1):
                out[m] = park
            stats["parks"] = stats.get("parks", 0) + 1
    # ---- 5. exact verification of every beam minute of the day + repair
    rows = {}
    first = ks[0] if ks else None
    for i in beam_minutes:
        pose = out.get(i)
        if pose is None:
            pose = plan[first] if first is not None and i < first else (plan[ks[-1]] if ks else (15.0, 100.0))
            out[i] = pose
        pose = (round(pose[0], 2), round(pose[1], 2))       # the 0.01 slider precision actually played
        out[i] = pose
        row = ph.pred_row(i, pose, tpt)
        if i in plan and not row[4] and cfg.get("REPAIR", True):
            r, e = ph.band_center(i, pose[1], W * wmap.get(i, 1.0), r0=pose[0])
            stats["repairs"] += 1
            # accept only a small ring correction: a big one would be a one-minute spike (a hop and back)
            if r is not None and e["clean"] and abs(r - pose[0]) <= cfg.get("REPAIR_MAX_DR", 0.25):
                out[i] = (round(r, 2), pose[1])
                row = ph.pred_row(i, out[i], tpt)
            else:
                stats["repair_fail"] += 1
                why = "rep_noband" if r is None else ("rep_unclean" if not e["clean"] else "rep_far")
                stats[why] = stats.get(why, 0) + 1
        rows[i] = row
    if cfg.get("CALM") and rows:
        _calm(ph, cfg, beam_minutes, out, rows, tpt, stats)
    if cfg.get("TURN_DWELL") and rows:              # M-012: ring hold via the aperture (no ring turn-backs)
        _turn_dwell(ph, cfg, beam_minutes, out, rows, rrd, tpt, stats, RLO, RHI)
    stats["sec"] = time.time() - t_day
    return day, out, rows, stats


# ------------------------------------------------------------------------------------------------ parent
def candidate_minutes(ev):
    """beam_on minutes whose sun bin lets light onto a room face for at least one grid pose."""
    nk = len(ev.sun.keys)
    lit = bytearray(nk)
    for k in range(nk):
        for p in range(ev.npose):
            r = ev.scan.rec[p][k]
            if r is not None and r[1]:
                lit[k] = 1
                break
    out = []
    for i in range(ev.n):
        if ev.beam_on[i]:
            k = ev.sun.key_of_hour[i]
            if k >= 0 and lit[k]:
                out.append(i)
    return out


def assemble(n, day_poses, night_glide=60):
    """Per-step poses for the whole year from the per-day planned minutes: hold the last pose, and glide to the
    next day's first pose over the night_glide minutes before it (smootherstep)."""
    keys = sorted(day_poses)
    poses = [None] * n
    if not keys:
        return [(15.0, 100.0)] * n
    for i in keys:
        poses[i] = day_poses[i]
    first = keys[0]
    for i in range(first):
        poses[i] = day_poses[first]
    for x, y in zip(keys, keys[1:]):
        if y - x > 1:
            p0, p1 = day_poses[x], day_poses[y]
            g0 = max(x + 1, y - night_glide)
            for m in range(x + 1, y):
                if m < g0:
                    poses[m] = p0
                else:
                    e = smootherstep((m - g0 + 1) / float(y - g0 + 1))
                    poses[m] = tuple(u + (v - u) * e for u, v in zip(p0, p1))
    last = keys[-1]
    for i in range(last + 1, n):
        poses[i] = day_poses[last]
    return poses


def run(model, wtag, cfg_over=None, days=None):
    cfg = dict(CFG)
    if cfg_over:
        cfg.update(cfg_over)
    t0 = time.time()
    scen, blade = cfg.get("SCEN"), cfg.get("BLADE")
    if scen is not None:                        # ANY window (scenario.py): facing, level, scale, width, beam width
        import scenario as SC
        ev, P, jobp = SC.load(scen)
        cand = [i for i in range(ev.n) if ev.beam_on[i]]
    else:
        ev, P, job, _ = L.load(wtag)
        jobp = None
        cand = candidate_minutes(ev)
    ctx = PH.Ctx(ev)
    if blade:
        import blade as B
        ph = B.BladePhysics(ctx, blade)
    else:
        ph = PH.Physics(ctx)
    if cfg.get("TARGET_SPEC"):                  # {"depth_m", "offset_m", "height_m", "strength"}
        ts = cfg["TARGET_SPEC"]
        pt = ph.target_point(ts["depth_m"], ts["offset_m"], ts.get("height_m", 0.0))
        cfg["TARGET"] = dict(ts, point_xyz=list(pt), width_tol=cfg.get("WTOL"))
        cfg["TW"] = float(ts.get("strength", 0.0)) * cfg.get("TW_FULL", 0.9)
    perday = {}
    for i in cand:
        perday.setdefault(i // 1440, []).append(i)
    beamday = {}
    for i in range(ev.n):
        if ev.beam_on[i]:
            beamday.setdefault(i // 1440, []).append(i)
    dlist = sorted(set(perday) | set(beamday))
    if days:
        dlist = [d for d in dlist if d in days]
    jobs = [(d, perday.get(d, []), beamday.get(d, []), cfg) for d in dlist]
    print("%s w%s: %d days, %d candidate minutes, %d beam minutes" % (model, wtag, len(jobs), len(cand),
                                                                        sum(len(v) for v in beamday.values())))
    if scen is None:
        del ev
    blob = pickle.dumps(ctx)
    day_poses, rows, st = {}, [], {"repairs": 0, "repair_fail": 0, "sec": 0.0}
    if blade:
        import blade as B
        pinit, pargs, fn = B._init, (blob, blade), B._bday
    else:
        pinit, pargs, fn = _init, (blob,), _day
    with mp.Pool(cfg["NPROC"], pinit, pargs) as pool:
        for k, (day, out, rws, stats) in enumerate(pool.imap_unordered(fn, jobs, chunksize=1)):
            day_poses.update(out)
            rows.extend(rws.values())
            for kk in stats:
                if kk == "evals_samples":
                    continue
                if isinstance(stats[kk], list):
                    st.setdefault(kk, []).extend(stats[kk])
                else:
                    st[kk] = st.get(kk, 0) + stats[kk]
            if (k + 1) % 30 == 0:
                print("  %d / %d days  %.0f s" % (k + 1, len(jobs), time.time() - t0))
    rows.sort(key=lambda r: r[0])
    poss = sorted(st.pop("poss", []))
    depth_ctrl = st.pop("depth_ctrl", None)
    poses = assemble(ctx.n, day_poses, cfg["NIGHT_GLIDE"])
    runs = L.runs_from_poses(poses, ctx.n, nd=2)
    # predictions at the ROUNDED poses actually played (2 decimals) - recompute rows whose pose was rounded
    pr = L.poses_from_runs(runs, ctx.n)
    fix = [(r[0], pr[r[0]]) for r in rows if (abs(pr[r[0]][0] - r[1]) > 0.0005 or abs(pr[r[0]][1] - r[2]) > 0.0005)]
    if fix and not blade:
        tp = tuple(cfg["TARGET"]["point_xyz"]) if cfg.get("TARGET") else None
        new = dict((r[0], r) for r in PR.predict_poses(ctx, fix, tp, cfg["NPROC"]))
        rows = [new.get(r[0], r) for r in rows]
    st["rounding_fix"] = len(fix)
    tgt = cfg.get("TARGET")
    cfg_rec = dict((k, v) for k, v in cfg.items() if k not in ("TARGET",))
    extra = {"cfg": cfg_rec, "repairs": st}
    if depth_ctrl is not None:
        dc = sorted(depth_ctrl)
        extra["depth_control_m"] = {"n_samples": len(dc), "p50": dc[len(dc) // 2] if dc else None,
                                    "p90": dc[int(0.9 * (len(dc) - 1))] if dc else None,
                                    "share_ge_0.5m": round(sum(1 for x in dc if x >= 0.5) / max(1, len(dc)), 3)}
    if scen is not None:                        # scenario plans: lab\plans\scen\<model>.json, own job + pack
        import scenario as SC
        d = os.path.join(L.LAB, "plans", "scen")
        if not os.path.isdir(d):
            os.makedirs(d)
        packp = os.path.join(SC.OUT, SC.sid(scen), "pack.json")
        SC.make_pack(scen, poss, packp, ev, P)
        path = os.path.join(d, "%s.json" % model)
        data = {"contract": L.CONTRACT, "model": model, "width_m": P["beam_width_m"], "job": jobp, "step_min": 1,
                "n_steps": ctx.n, "params": ["ring", "aperture"], "ranges": [15.0, 100.0], "runs": runs,
                "target": tgt, "pred": None, "scenario": SC.norm(scen), "pack": packp,
                "possible_rule": "planner clean band (lite scenario)"}
        data.update(extra)
        with open(path, "w") as f:
            f.write(json.dumps(data, separators=(",", ":")))
    else:
        path = L.write_plan(model, wtag, runs, target=tgt, extra=extra)
    pred = path[:-5] + "_pred.csv"
    PR.write_pred(pred, rows)
    plan = L.read_plan(path)
    plan["pred"] = os.path.basename(pred)
    if blade:
        import blade as B
        plan["params"] = ["ring", "aperture", "aperture_b"]
        plan["ranges"] = [15.0, 100.0, 100.0]
        plan["blade"] = blade
        plan["blade_spec"] = {"joint": "jb", "new_param": "aperture_b", "gain_factor": B.BLADES[blade],
                              "deg_per_pct": B.JB_DEG_PER_PCT * B.BLADES[blade],
                              "note": "bottom blade on its own drive; today's coupled = aperture / gain_factor"}
    with open(path, "w") as f:
        f.write(json.dumps(plan, separators=(",", ":")))
    print("%s w%s done in %.0f s: %d runs, repairs %s -> %s" % (model, wtag, time.time() - t0, len(runs), st, path))
    return path


if __name__ == "__main__":
    over = json.loads(sys.argv[3]) if len(sys.argv) > 3 else None
    days = None
    if over and "DAYS" in over:
        days = set(over.pop("DAYS"))
    run(sys.argv[1], L.wtag_norm(sys.argv[2]), over, days)
