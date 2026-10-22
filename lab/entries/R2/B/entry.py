# -*- coding: utf-8 -*-
"""
entry.py - MODELER-B round-2 entry (complete: LIGHT + THERMAL + COMBINED, 1-15 windows).  The TESTER calls only this.

    plan_light(scene, weather_epw, params)                    -> {window_id: plan (contract v1, per-minute runs)}
    thermal_model()                                           -> HB-L2 object (make_tparams / assemble / window_gains /
                                                                 run_year), all windows' gains summed in one room
    plan_combined(scene, weather_epw, construction, params, slider) -> {window_id: plan}; slider 0 (comfort) .. 1 (light)
                                                                 or a schedule [[minute, value], ...] (15-min glides)
    param_ranges()                                            -> every parameter: range / choices, default, units

scene: {"windows": [spec, ...]} (or a list of specs, or one spec).  spec = scenario.py window spec {"id", "facing",
"tilt_deg", "level_m", "scale", "width_f", "beam_w"} (+ optional "plane": {"origin": [x,y,z], "normal": [x,y,z]}
used for facing / tilt when "facing" is missing: +Y = north, +Z = up; "area_m2" for the thermal glazing).
LIGHT = the champion G2R3T band glide (declared reuse, frozen copy glide_champion.py, origin A line M-012), planned
ONCE PER WINDOW CLASS (facing, tilt, level, scale, width, beam width; R-003 M1) on the given EPW with the DNI gate
`dni_beam_min` (default 120 W/m2 = the TESTER's gate).  THERMAL = HB-L2 (thermal_hb2.py, B's own).
COMBINED = B's nested-knapsack slider (R-003 S1 / R-002 C1): every (window class, day) item gets a comfort plan
(open pose on heat-need days, shut pose on cool-need days, judged with HB-L2 runs on the summed gains) and a ratio
r = discomfort removed / beam minutes lost; slider s converts the best-r items until (1 - s) of the total benefit is
reached, so the converted sets are NESTED (monotone by construction); converted days glide in / out over 15 min.
"""
from __future__ import division

import array
import copy
import json
import math
import os
import pickle
import sys
import time

HERE = os.path.dirname(os.path.abspath(__file__))
LAB = os.path.dirname(os.path.dirname(os.path.dirname(HERE)))
MODELS = os.path.join(LAB, "models")
for _p in (MODELS, HERE):                         # HERE first: the frozen copies win
    if _p in sys.path:
        sys.path.remove(_p)
    sys.path.insert(0, _p)
import thermal_hb2 as TH  # noqa: E402
import labcore as L  # noqa: E402
import physics as PH  # noqa: E402
import scenario as SC  # noqa: E402
import glide_champion as G  # noqa: E402

C = L.C
VERSION = "B-R2 (MODELER-B, 2026-10-13)"
OUT = os.path.join(LAB, "plans", "B", "R2")
EPW_DEFAULT = os.environ.get("WO_EPW", os.path.join(os.path.dirname(LAB), "data", "weather",
                                                   "AUS_NSW_Newcastle.AP-RAAF.Williamtown.947760_TMYx.2011-2025.epw"))  # not shipped: see README
G2R3T = {"DIPFILL": 6, "DWELL": 13, "DIPFILL_V2": True, "FB_SMOOTH": True, "TURN_DWELL": True}
P_OPEN, P_SHUT = (15.0, 0.0), (0.0, 100.0)

LIGHT_PARAMS = [
    ("beam_w", 0.10, 0.01, 0.60, "m", "year-round beam width (per window spec; this is the default)"),
    ("dni_beam_min", 120.0, 0.0, 400.0, "W/m2", "a beam minute is wanted only at this direct sun or more"),
    ("light_planner", "g2r3t", ["g2r3t"], None, "-", "light planner (champion G2R3T reused)"),
    ("nproc", 1, 1, 12, "processes", "worker processes for the day plans (lab rule: 1 unless the TESTER allows more)"),
    ("night_glide_min", 60, 10, 240, "min", "glide from the evening pose to the next morning's first pose"),
]
COMBINED_PARAMS = [
    ("glide_min", 15, 5, 15, "min", "glide in / out of a comfort day and after a slider change (<= 15)"),
    ("sample_min", 30, 10, 60, "min", "sun / gain sampling step for the comfort plans"),
    ("comfort_category", 2, 1, 3, "-", "EN 16798-1 adaptive category (band +-2 / 3 / 4 K)"),
    ("thermal_case", "A04_felix_open", None, None, "id", "frozen E+ drive (room + weather) the comfort decisions use"),
    ("open_ring", 15.0, 0.0, 15.0, "deg", "comfort pose for heat: ring"),
    ("open_aperture", 0.0, 0.0, 100.0, "%", "comfort pose for heat: aperture close"),
    ("shut_ring", 0.0, 0.0, 15.0, "deg", "comfort pose for cooling: ring"),
    ("shut_aperture", 100.0, 0.0, 100.0, "%", "comfort pose for cooling: aperture close"),
]


def _rows(table, group):
    out = []
    for name, default, lo, hi, units, meaning in table:
        r = {"name": name, "default": default, "units": units, "group": group, "meaning": meaning}
        if isinstance(lo, list):
            r["choices"] = lo
        elif lo is not None:
            r["min"], r["max"] = lo, hi
        out.append(r)
    return out


def param_ranges():
    """Every parameter of the entry: {'light': [...], 'thermal': [...], 'combined': [...]}, each row with name,
    default, min / max or choices, units, meaning (thermal rows also the expected response)."""
    return {"light": _rows(LIGHT_PARAMS, "light"), "thermal": TH.param_table(True),
            "combined": _rows(COMBINED_PARAMS, "combined")}


def thermal_model():
    return TH.ThermalModel()


# ============================================================================================ windows / loading
def windows_of(scene):
    if isinstance(scene, dict) and "windows" in scene:
        wl = scene["windows"]
    elif isinstance(scene, (list, tuple)):
        wl = list(scene)
    else:
        wl = [scene or {}]
    out = []
    for k, w in enumerate(wl):
        w = dict(w or {})
        if w.get("facing") is None and w.get("plane"):
            n = w["plane"].get("normal") or w["plane"].get("z")
            h = math.hypot(n[0], n[1])
            w["facing"] = (math.degrees(math.atan2(n[0], n[1])) % 360.0) if h > 1e-6 else 0.0
            w.setdefault("tilt_deg", max(0.0, min(90.0, math.degrees(math.atan2(n[2], h)))))
        spec = SC.norm(dict((key, w[key]) for key in SC.DEFAULT if key in w))
        if spec["facing"] is None:
            spec["facing"] = 0.0
        out.append((str(w.get("id") or w.get("name") or "w%d" % (k + 1)), spec, w))
    return out


def class_key(spec):
    return tuple(round(float(spec[k]), 4) for k in ("facing", "tilt_deg", "level_m", "scale", "width_f", "beam_w"))


def _epw_tag(epw):
    return os.path.splitext(os.path.basename(epw))[0].split(".")[-2] if epw.count(".") >= 2 else "epw"


_EV = {}


def load_window(spec, epw, dni_min):
    """(ev, P, job path, folder) for one window class on the given weather (scenario.load with the EPW + gate)."""
    key = (class_key(spec), epw, float(dni_min))
    if key in _EV:
        return _EV[key]
    s = SC.norm(spec)
    job = json.load(open(SC.BASE_JOB))
    P0 = job["params"]
    scene = json.load(open(P0["scene_path"]))
    scene = SC.transform_scene(scene, s["level_m"], s["scale"], s["width_f"])
    ov = dict(P0["overrides"])
    ov["beam_width_m"] = float(s["beam_w"])
    ov["dni_beam_min"] = float(dni_min)
    ov["facing_deg"] = float(s["facing"]) % 360.0
    P = C.make_params(ov)
    w = C.prepare_weather(C.read_epw(epw), P["time_step_min"])
    tilt = float(s.get("tilt_deg") or 0.0)
    if tilt:
        C.apply_orientation(scene, float(s["facing"]) % 360.0)
        scene = SC.tilt_scene(scene, tilt)
        sc = C.resolve_scene(scene, dict(P, facing_deg=None))
    else:
        sc = C.resolve_scene(scene, P)
    win = C.split_windows(sc)[0]
    scan = SC.lite_scan(win, w, P)
    ev = C.Evaluator(win, w, scan, P)
    d = os.path.join(OUT, "%s_%s_g%03d" % (SC.sid(s), _epw_tag(epw), int(dni_min)))
    if not os.path.isdir(d):
        os.makedirs(d)
    with open(os.path.join(d, "scene.json"), "w") as f:
        json.dump(scene, f)
    j2 = copy.deepcopy(job)
    j2["params"]["scene_path"] = os.path.join(d, "scene.json")
    j2["params"]["overrides"] = ov
    j2["params"]["epw_path"] = epw
    j2["params"]["cache_dir"] = os.path.join(LAB, "out", "_scan_cache")
    j2["scenario"] = s
    jp = os.path.join(d, "job.json")
    with open(jp, "w") as f:
        json.dump(j2, f, indent=1)
    _EV[key] = (ev, P, jp, d)
    return _EV[key]


# ============================================================================================ LIGHT
def _plan_class(spec, epw, prm):
    ev, P, jp, d = load_window(spec, epw, prm["dni_beam_min"])
    ctx = PH.Ctx(ev)
    cfg = dict(G.CFG)
    cfg.update(G2R3T)
    cfg["NPROC"] = int(prm["nproc"])
    cfg["NIGHT_GLIDE"] = int(prm["night_glide_min"])
    cand = [i for i in range(ev.n) if ev.beam_on[i]]
    perday = {}
    for i in cand:
        perday.setdefault(i // 1440, []).append(i)
    days = sorted(perday)
    if prm.get("days") is not None:
        days = [x for x in days if x in set(prm["days"])]
    jobs = [(x, perday[x], perday[x], cfg) for x in days]
    t0 = time.time()
    day_poses, poss, vis_day = {}, [], {}
    if cfg["NPROC"] > 1 and len(jobs) > 1:
        import multiprocessing as mp
        with mp.Pool(cfg["NPROC"], G._init, (pickle.dumps(ctx),)) as pool:
            res = list(pool.imap_unordered(G._day, jobs, chunksize=1))
    else:
        G._PH = PH.Physics(ctx)
        res = [G._day(j) for j in jobs]
    for day, out, rws, stats in res:
        day_poses.update(out)
        poss.extend(stats.get("poss", []))
        vis_day[day] = sum(1 for r in rws.values() if r[3])
    poses = G.assemble(ctx.n, day_poses, cfg["NIGHT_GLIDE"])
    runs = L.runs_from_poses(poses, ctx.n, nd=2)
    packp = os.path.join(d, "pack.json")
    SC.make_pack(spec, sorted(poss), packp, ev, P)
    plan = {"contract": 1, "model": "B_R2_light", "width_m": P["beam_width_m"], "job": jp, "step_min": 1,
            "n_steps": ctx.n, "params": ["ring", "aperture"], "ranges": [15.0, 100.0], "runs": runs, "target": None,
            "pred": None, "scenario": SC.norm(spec), "pack": packp, "weather": epw, "dni_beam_min": prm["dni_beam_min"],
            "planner": "G2R3T champion (declared reuse, glide_champion.py) per window class", "entry": VERSION,
            "sec": round(time.time() - t0, 1), "visible_planned_min_per_day": vis_day}
    path = os.path.join(d, "plan_light.json")
    with open(path, "w") as f:
        f.write(json.dumps(plan, separators=(",", ":")))
    plan["path"] = path
    return plan, poses, ctx


def _prm(params, table):
    p = dict((r[0], r[1]) for r in table)
    p.update(params or {})
    return p


def plan_light(scene, weather_epw=None, params=None, _keep=None):
    """{window_id: plan}.  Windows of one class (same facing / tilt / level / size / width) share one plan (the band
    does not depend on where the window sits on the wall); each window gets its own copy with its id."""
    prm = _prm(params, LIGHT_PARAMS)
    epw = weather_epw or EPW_DEFAULT
    out, done = {}, {}
    for wid, spec, raw in windows_of(scene):
        if raw.get("beam_w") is None and params and params.get("beam_w"):
            spec["beam_w"] = float(params["beam_w"])
        k = class_key(spec)
        if k not in done:
            done[k] = _plan_class(spec, epw, prm)
        plan = dict(done[k][0])
        plan["window_id"] = wid
        plan["class_key"] = list(k)
        out[wid] = plan
        if _keep is not None:
            _keep[wid] = (k, done[k][1], done[k][2], spec, raw)
    return out


# ============================================================================================ COMBINED
def smootherstep(s):
    s = 0.0 if s < 0 else (1.0 if s > 1 else s)
    return s * s * s * (s * (6 * s - 15) + 10)


def _gain(ph, i, pose):
    """(beam W, diffuse W) before the glass at minute i and pose (exact physics)."""
    cx = ph.cx
    dni, dhi, ghi = cx.w["dni"][i], cx.w["dhi"][i], cx.w["ghi"][i]
    diff = dhi * ph.dsky(pose) + ghi * ph.dgnd(pose)
    if dni <= 0 or cx.alt[i] <= 0:
        return 0.0, diff
    s = ph.sun_vec(i)
    kf = ph.key_frac(s)
    rec = ph.record(pose, s) if kf >= 0.02 else None
    return (dni * rec[0] * ph.a_scale * kf if rec else 0.0), diff


def _series10(ph, poses_fn, n10, sample, days=None):
    """Beam and diffuse W per 10-min step (index i = minutes 10i..10i+9), sampled every `sample` min in daylight."""
    beam = array.array("d", [0.0]) * n10
    dif = array.array("d", [0.0]) * n10
    alt = ph.cx.alt
    step = max(10, int(sample)) // 10
    for i10 in range(0, n10, step):
        if days is not None and (i10 // 144) not in days:
            continue
        m = min(len(alt) - 1, i10 * 10 + 5)
        if alt[m] <= 0:
            continue
        b, d = _gain(ph, m, poses_fn(m))
        for q in range(i10, min(n10, i10 + step)):
            beam[q], dif[q] = b, d
    return beam, dif


def _band(op_ref_te, cat):
    """EN 16798-1 adaptive band per 10-min step from the outdoor running mean."""
    te = op_ref_te
    n = len(te)
    dm = [sum(te[d * 144:(d + 1) * 144]) / 144.0 for d in range(n // 144)]
    trm, r = [], dm[-1]
    for x in dm:
        trm.append(r)
        r = 0.8 * r + 0.2 * x
    w = {1: 2.0, 2: 3.0, 3: 4.0}[int(cat)]
    lo, hi = [], []
    for d in range(n // 144):
        tc = 0.33 * min(30.0, max(10.0, trm[d])) + 18.8
        lo.extend([tc - w - 1.0] * 144)
        hi.extend([tc + w] * 144)
    return lo, hi


def _ddh_days(op, lo, hi):
    nd = len(op) // 144
    cold, warm = [0.0] * nd, [0.0] * nd
    for i in range(nd * 144):
        x = op[i]
        if x < lo[i]:
            cold[i // 144] += (lo[i] - x) / 6.0
        elif x > hi[i]:
            warm[i // 144] += (x - hi[i]) / 6.0
    return cold, warm


def plan_combined(scene, weather_epw=None, construction=None, params=None, slider=1.0):
    """{window_id: plan} at the slider (0 comfort .. 1 light, or [[minute, value], ...])."""
    prm = _prm(params, LIGHT_PARAMS + COMBINED_PARAMS)
    keep = {}
    light = plan_light(scene, weather_epw, params, _keep=keep)
    cons = dict(construction or {})
    import eplus_load as EL
    R = EL.load(cons.get("case") or prm["thermal_case"])
    over = dict(R["over"])
    over.update(dict((k, v) for k, v in cons.items() if k not in ("case", "env")))
    known = set(p[0] for p in TH.PARAMS)
    tp = TH.make_tparams(dict((k, v) for k, v in over.items() if k in known))
    env = copy.deepcopy(cons.get("env") or R["env"])
    wins = list(keep.items())
    env["glazing"] = [{"area_m2": float(raw.get("area_m2") or 10.4 * spec["scale"] ** 2 * spec["width_f"]),
                       "azimuth_deg": float(spec["facing"]), "tilt_deg": 90.0 - float(spec["tilt_deg"]),
                       "name": wid} for wid, (k, poses, ctx, spec, raw) in wins]
    bld = TH.assemble(env, tp)
    n10 = R["n"]
    days = set(prm["days"]) if prm.get("days") is not None else None
    classes = {}
    for wid, (k, poses, ctx, spec, raw) in wins:
        classes.setdefault(k, {"poses": poses, "ctx": ctx, "wids": []})["wids"].append(wid)
    t0 = time.time()
    po = (float(prm["open_ring"]), float(prm["open_aperture"]))
    psh = (float(prm["shut_ring"]), float(prm["shut_aperture"]))
    ser = {}
    for k, c in classes.items():
        ph = PH.Physics(c["ctx"])
        poses = c["poses"]
        ser[k] = {"L": _series10(ph, lambda m, p=poses: p[m], n10, prm["sample_min"], days),
                  "O": _series10(ph, lambda m: po, n10, prm["sample_min"], days),
                  "S": _series10(ph, lambda m: psh, n10, prm["sample_min"], days)}

    def op_for(choice):                              # choice: class key -> "L" / "O" / "S" (all days)
        bl, dl = [], []
        for wid, (k, poses, ctx, spec, raw) in wins:
            b, d = ser[k][choice.get(k, "L")]
            bl.append({"floor": b})
            dl.append(d)
        win = TH.window_gains(bld, bl, dl)
        return TH.run_year(bld, dict(R["wx"]), R["climate"], win=win)["op"]
    lo, hi = _band(R["wx"]["temp"], prm["comfort_category"])
    base = _ddh_days(op_for({}), lo, hi)
    allo = _ddh_days(op_for(dict((k, "O") for k in classes)), lo, hi)
    alls = _ddh_days(op_for(dict((k, "S") for k in classes)), lo, hi)
    nd = len(base[0])
    items = []
    for k, c in classes.items():
        vis = light[c["wids"][0]].get("visible_planned_min_per_day", {})
        bL, _dL = ser[k]["L"]
        bO, _dO = ser[k]["O"]
        bS, _dS = ser[k]["S"]
        for d in range(min(nd, 365)):
            if days is not None and d not in days:
                continue
            gain_o = gain_s = 0.0
            for i in range(d * 144, (d + 1) * 144):
                gain_o += bO[i] - bL[i]
                gain_s += bL[i] - bS[i]
            tot_o = sum(sum(ser[kk]["O"][0][d * 144:(d + 1) * 144]) - sum(ser[kk]["L"][0][d * 144:(d + 1) * 144])
                        for kk in classes) or 1e-9
            tot_s = sum(sum(ser[kk]["L"][0][d * 144:(d + 1) * 144]) - sum(ser[kk]["S"][0][d * 144:(d + 1) * 144])
                        for kk in classes) or 1e-9
            dO = (sum(base[x][d] for x in (0, 1)) - sum(allo[x][d] for x in (0, 1))) * max(0.0, gain_o / tot_o)
            dS = (sum(base[x][d] for x in (0, 1)) - sum(alls[x][d] for x in (0, 1))) * max(0.0, gain_s / tot_s)
            alt, ben = ("O", dO) if dO >= dS else ("S", dS)
            if ben <= 1e-6:
                continue
            lost = float(vis.get(d, vis.get(str(d), 0))) * len(c["wids"])
            items.append((ben / (lost + 1.0), ben, lost, k, d, alt))
    items.sort(key=lambda x: -x[0])
    btot = sum(x[1] for x in items) or 1e-9

    def converted(s):
        s = min(1.0, max(0.0, float(s)))
        need, acc, outc = (1.0 - s) * btot, 0.0, {}
        for it in items:
            if acc >= need - 1e-9:
                break
            outc[(it[3], it[4])] = it[5]
            acc += it[1]
        return outc
    G_MIN = int(min(15, max(1, prm["glide_min"])))

    def poses_at(s, k):
        c = classes[k]
        base_p = c["poses"]
        conv = converted(s)
        alt_ = c["ctx"].alt
        out = list(base_p)
        for (kk, d), a in conv.items():
            if kk != k:
                continue
            tgt = po if a == "O" else psh
            day_m = [m for m in range(d * 1440, (d + 1) * 1440) if alt_[m] > 0]
            if not day_m:
                continue
            m0, m1 = max(0, day_m[0] - 30), min(len(out) - 1, day_m[-1] + 30)
            for m in range(m0, m1 + 1):
                e = min(smootherstep((m - m0) / float(G_MIN)), smootherstep((m1 - m) / float(G_MIN)))
                bp = base_p[m]
                out[m] = (bp[0] + (tgt[0] - bp[0]) * e, bp[1] + (tgt[1] - bp[1]) * e)
        return out
    sched = slider if isinstance(slider, (list, tuple)) else [[0, slider]]
    sched = sorted([int(a), float(b)] for a, b in sched)
    res = {}
    for k, c in classes.items():
        n = len(c["poses"])
        cache = {}
        for _m, v in sched:
            if v not in cache:
                cache[v] = poses_at(v, k)
        out = list(cache[sched[0][1]])
        for (m_prev, v_prev), (m_new, v_new) in zip(sched, sched[1:]):
            pa, pb = cache[v_prev], cache[v_new]
            for m in range(m_new, n):
                e = smootherstep((m - m_new) / float(G_MIN)) if m < m_new + G_MIN else 1.0
                out[m] = (pa[m][0] + (pb[m][0] - pa[m][0]) * e, pa[m][1] + (pb[m][1] - pa[m][1]) * e)
        runs = L.runs_from_poses(out, n, nd=2)
        for wid in c["wids"]:
            p = dict(light[wid])
            p.update({"model": "B_R2_combined", "runs": runs, "slider": slider,
                      "converted_days": sorted([d, a] for (kk, d), a in converted(sched[-1][1]).items() if kk == k),
                      "combined": "nested knapsack over (window class, day); comfort judged by HB-L2 on %s" %
                                  (cons.get("case") or prm["thermal_case"]),
                      "ddh_base_cold_warm": [round(sum(base[0]), 1), round(sum(base[1]), 1)],
                      "ddh_all_open": [round(sum(allo[0]), 1), round(sum(allo[1]), 1)],
                      "ddh_all_shut": [round(sum(alls[0]), 1), round(sum(alls[1]), 1)],
                      "items": len(items), "sec_combined": round(time.time() - t0, 1)})
            p.pop("path", None)
            res[wid] = p
    return res


if __name__ == "__main__":
    import pprint
    pr = param_ranges()
    print("params: light %d thermal %d combined %d" % (len(pr["light"]), len(pr["thermal"]), len(pr["combined"])))
    pprint.pprint(thermal_model().name)
