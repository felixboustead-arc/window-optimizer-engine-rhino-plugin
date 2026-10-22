# -*- coding: utf-8 -*-
"""
physics.py - EXACT beam physics at ANY pose (on or off the 121-pose grid) and at the EXACT sun position of any
minute (no 1-degree sun bins), built only from window_beams_core's own functions (Mechanism.pose_geometry,
clear_pieces, Target.mapping, lit_on_target, beam_metrics) and the same scoring formulas as
Evaluator._beam_eval_full.  Used by the lab planners (off-grid poses) and by predict.py (the per-minute prediction
sidecar of every plan, so every model - the old ones too - is measured the same way).

The only approximation: the diffuse-light factors of a pose (sky / ground view through the window, used for the
room's ambient light in the contrast term) are interpolated bilinearly from the grid poses (they vary slowly).

Light target (window-relative, survives re-orientation):  T = exit_point + depth * r + offset * w + lift * u, with
r = horizontal direction into the room, w = to your RIGHT as you stand in the room looking out of the window,
u = up, lift chosen so T is ``height_m`` above the floor.  target_distance() = distance from T to the nearest
lit point of the beam on the room's surfaces.
"""
from __future__ import division

import math
import os
import sys

CORE_DIR = os.path.abspath(os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", ".."))  # project root
if CORE_DIR not in sys.path:
    sys.path.insert(0, CORE_DIR)
import window_beams_core as C  # noqa: E402


CLEAN_SHARE = 0.7      # a CLEAN beam: visible, and its face carries >= this share of the light landing on the room


class Ctx(object):
    """Everything the exact physics needs, small and picklable (for worker processes): no scan records."""

    def __init__(self, ev):
        import array
        self.P = ev.P
        self.scene = ev.scene
        self.scale = ev.scale
        self.a_scale = ev.a_scale
        self.vals = ev._grid()["vals"]
        self.ranges = ev._ranges()
        self.poses = [tuple(float(v) for v in p) for p in ev.scan.poses]
        self.dsky = list(ev._dsky)
        self.dgnd = list(ev._dgnd)
        self.room_m2 = ev.room_m2
        self.tw = dict(ev._tw)
        self.w = dict((k, array.array("d", ev.w[k])) for k in ("dni", "dhi", "ghi"))
        self.alt = array.array("d", ev.sun.alt)
        self.az = array.array("d", ev.sun.az)
        self.north = ev.sun.north
        self.beam_on = bytearray(1 if b else 0 for b in ev.beam_on)
        self.n = ev.n
        self.step_min = ev.step_min


class Physics(object):
    def __init__(self, ev):
        cx = ev if isinstance(ev, Ctx) else Ctx(ev)
        self.cx = cx
        self.P = cx.P
        self.scene = cx.scene
        self.mech = C.Mechanism(cx.scene)
        self.targets = [C.Target(t) for t in cx.scene["targets"]]
        self.exitp = tuple(cx.scene["exit_point"])
        self.sc = cx.scale
        self.a_scale = cx.a_scale
        self.north = cx.north
        self.vals = cx.vals                        # grid values per parameter
        self.ranges = cx.ranges
        self._index = dict((p, j) for j, p in enumerate(cx.poses))
        self._geo = {}                              # pose -> (mouths, blockers)
        self.ctx = C.build_context(cx.scene, cx.P)
        if self.ctx is not None:
            self._ctx_pts = C.window_front_points(cx.scene)
            self._ctx_skip = cx.P["context_skip_m"] / float(cx.scene.get("units_to_m", 1.0))
        r, w, u = room_axes(cx.scene)
        self.r, self.w, self.u = r, w, u
        self.ext = C.vunit(tuple(cx.scene["exterior_dir"]))
        self.nin = C.vscale(self.ext, -1.0)          # the window's own inward normal (tilts with the window)
        # floor level: the largest up-facing room face
        best = None
        for t in self.targets:
            if C.vdot(t.n, u) > 0.9:
                a = C.poly_area(t.boundary)
                if best is None or a > best[0]:
                    best = (a, C.vdot(t.q0, u))
        self.floor_u = best[1] if best else C.vdot(self.exitp, u) - 1.5 / self.sc
        # room faces (the big ones) vs the window's own reveal / frame faces: a face counts as "room" when some of it
        # is more than 0.6 m in from the window's exit plane, measured along the window's OWN inward normal (= the
        # horizontal depth for a wall window; the drop below the ceiling for a skylight)
        self.room_face = []
        # downstream guard: light can only land on the room side of the exit plane.  Per face, the exit plane as a
        # 2D half-plane (ca*a + cb*b + c0 <= 0 is the room side); None = the whole face is on the room side
        self.half = []
        for t in self.targets:
            pts = [C.vadd(t.q0, C.vadd(C.vscale(t.e1, a), C.vscale(t.e2, b))) for (a, b) in t.boundary]
            depth = max(C.vdot(C.vsub(p, self.exitp), self.nin) for p in pts) * self.sc
            self.room_face.append(depth > 0.6)
            out = max(C.vdot(C.vsub(p, self.exitp), self.ext) for p in pts) * self.sc
            self.half.append(None if out <= 0.002 else
                             (C.vdot(t.e1, self.ext), C.vdot(t.e2, self.ext), C.vdot(C.vsub(t.q0, self.exitp), self.ext)))
        self.n_clipped = 0
        # band search in both directions (strip widening OR narrowing with the ring): on for tilted windows, or forced
        # by P["band_any"]; wall windows keep the measured widening band exactly as before
        self.band_any = bool(cx.scene.get("tilt_deg")) or bool(cx.P.get("band_any"))

    # ---------------------------------------------------------------- geometry
    def geometry(self, pose):
        key = (round(float(pose[0]), 6), round(float(pose[1]), 6))
        g = self._geo.get(key)
        if g is None:
            g = self.mech.pose_geometry(key)
            if len(self._geo) > 4000:
                self._geo.clear()
            self._geo[key] = g
        return g

    def sun_vec(self, i):
        return C.sun_vector(self.cx.alt[i], self.cx.az[i], self.north)

    def key_frac(self, s):
        if self.ctx is None:
            return 1.0
        return self.ctx.fraction(s, self._ctx_pts, self._ctx_skip)

    def record(self, pose, s, polys=False):
        """Like one entry of the scan: (clear area perpendicular to the sun, [(target index, beam_metrics tuple
        [, lit polygons in the face's 2D frame])]) or None when no light passes."""
        mouths, blockers = self.geometry(pose)
        basis = C.sun_basis(s)
        pcs = C.clear_pieces(mouths, blockers, basis)
        if not pcs:
            return None
        a_perp = sum(C.poly_area(pc) for pc in pcs)
        if a_perp < C.EPS_AREA:
            return None
        trecs = []
        for ti, t in enumerate(self.targets):
            mp = t.mapping(s, basis)
            if mp is None:
                continue
            lit = C.lit_on_target(pcs, mp, t)
            if lit and self.half[ti] is not None:
                a0 = sum(C.poly_area(q) for q in lit)
                lit = [q for q in (_clip_half(pc, self.half[ti]) for pc in lit) if q]
                if sum(C.poly_area(q) for q in lit) < a0 - 1e-9:
                    self.n_clipped += 1
            bm = C.beam_metrics(lit, t, mp, self.exitp, s) if lit else None
            if bm is not None:
                trecs.append((ti, bm, lit) if polys else (ti, bm))
        return (a_perp, trecs)

    # ---------------------------------------------------------------- diffuse factors (bilinear on the grid)
    def _bilinear(self, pose, arr):
        (rv, av) = self.vals
        x, y = float(pose[0]), float(pose[1])
        x = min(max(x, rv[0]), rv[-1])
        y = min(max(y, av[0]), av[-1])
        i = max(0, min(len(rv) - 2, int((x - rv[0]) / (rv[1] - rv[0]))))
        j = max(0, min(len(av) - 2, int((y - av[0]) / (av[1] - av[0]))))
        fx = (x - rv[i]) / (rv[i + 1] - rv[i])
        fy = (y - av[j]) / (av[j + 1] - av[j])
        g = lambda a, b: arr[self._index[(rv[a], av[b])]]
        return ((1 - fx) * (1 - fy) * g(i, j) + fx * (1 - fy) * g(i + 1, j) + (1 - fx) * fy * g(i, j + 1)
                + fx * fy * g(i + 1, j + 1))

    def dsky(self, pose):
        return self._bilinear(pose, self.cx.dsky)

    def dgnd(self, pose):
        return self._bilinear(pose, self.cx.dgnd)

    # ---------------------------------------------------------------- scoring (same formulas as _beam_eval_full)
    def score(self, rec, pose, dni, dhi, ghi, kf=1.0):
        """(quality, visible, best detail dict or None, per-face list) for a record from record()."""
        if rec is None:
            return 0.0, False, None, []
        P, sc, cx = self.P, self.sc, self.cx
        v2 = P.get("beam_score_model", "v2") != "v1"
        A = rec[0] * self.a_scale * kf
        beam_w = dni * A
        diff_w = dhi * self.dsky(pose) + ghi * self.dgnd(pose)
        phi_in_lm = diff_w * P["k_diffuse"] + P["rho_target"] * beam_w * P["k_direct"]
        e_amb = phi_in_lm / (cx.room_m2 * (1.0 - P["rho_room"])) + 0.5
        r0, r1 = P["crisp_ratio"]
        tw = cx.tw
        best, bd, present = 0.0, None, False
        faces = []
        vis_room = []
        for item in rec[1]:
            ti, tr = item[0], item[1]
            t = self.targets[ti]
            area, width, length, hero, throw, cosi = tr[:6]
            area_m2 = area * self.a_scale
            faces.append((ti, area_m2, tr))
            if area_m2 < P["beam_min_area_m2"]:
                continue
            W, L, thr = width * sc, length * sc, throw * sc
            if W <= 1e-6:
                continue
            e_beam = dni * P["k_direct"] * cosi
            contrast = e_beam / e_amb
            pen = C.SUN_ANGULAR_DIAMETER * thr
            s_width = C.width_score_v2(W, P) if v2 else C.width_score(W, P)
            s_shape = math.sqrt(C.clip01(L / P["beam_length_min_m"]) *
                                C.clip01((L / W - 1.0) / max(P["beam_aspect_min"] - 1.0, 1e-6)))
            s_contrast = C.clip01(math.log10(max(contrast, 1.0)) / math.log10(max(P["contrast_target"], 1.01)))
            s_crisp = C.crisp_score_v2(pen, W, P) if v2 else C.clip01((r1 - pen / W) / (r1 - r0))
            terms = [(s_width, tw["width"]), (s_shape, tw["shape"]), (s_contrast, tw["contrast"]),
                     (s_crisp, tw["crisp"])]
            if t.heroes:
                terms.append((C.clip01(hero / area), tw["hero"]))
            s = C.wgeo(terms) * t.weight
            vis = ((C.width_score_v2(W, P) if not v2 else s_width) >= P.get("presence_width_min", 0.3) and
                   s_shape >= P.get("presence_shape_min", 0.5) and s_contrast >= P.get("presence_contrast_min", 0.5))
            if vis:
                present = True
                if self.room_face[ti]:
                    vis_room.append(area_m2)
            if s > best:
                best = s
                bd = {"target": t.name, "ti": ti, "width_m": W, "length_m": L, "area_m2": area_m2,
                      "point": tuple(tr[6:9]), "visible": vis, "contrast": contrast}
        room_tot = sum(a for ti, a, tr in faces if self.room_face[ti])
        share = (max(vis_room) / room_tot) if (vis_room and room_tot > 0) else 0.0
        self.last = (present and share >= CLEAN_SHARE, share, room_tot)
        return best, present, bd, faces

    def evaluate(self, i, pose, polys=False):
        """Everything about pose at step i (exact sun of that minute, that minute's weather)."""
        w = self.cx.w
        s = self.sun_vec(i)
        kf = self.key_frac(s)
        rec = self.record(pose, s, polys) if kf >= 0.02 else None
        q, vis, bd, faces = self.score(rec, pose, w["dni"][i], w["dhi"][i], w["ghi"][i], kf)
        clean, share, room_tot = self.last if rec is not None else (False, 0.0, 0.0)
        return {"q": q, "vis": vis, "clean": clean, "share": share, "room_m2": room_tot, "best": bd,
                "faces": faces, "rec": rec}

    # ---------------------------------------------------------------- one prediction row (contract v1 sidecar)
    PRED_COLS = ("step", "ring", "aperture", "visible", "clean", "quality", "width_m", "length_m", "face", "x", "y",
                 "z", "beam_m2", "room_lit_m2", "lx", "ly", "lz", "l_depth", "l_offset", "l_height", "share",
                 "target_d_m", "sun_alt", "sun_az")

    def pred_row(self, i, pose, target_pt=None):
        e = self.evaluate(i, pose, polys=target_pt is not None)
        b = e["best"]
        lc, tot = self.light_centroid(e["faces"])
        wc = self.window_coords(lc) if lc else (None, None, None)
        td = None
        if target_pt is not None and e["rec"] is not None:
            td = self.distance_to_lit(target_pt, e["rec"])
        f = lambda v, n=4: "" if v is None else round(v, n)
        return [i, round(pose[0], 3), round(pose[1], 3), int(e["vis"]), int(e["clean"]), f(e["q"], 3),
                f(b["width_m"] if b else None), f(b["length_m"] if b else None, 3), b["target"] if b else "",
                f(b["point"][0] if b else None, 3), f(b["point"][1] if b else None, 3), f(b["point"][2] if b else None, 3),
                f(b["area_m2"] if b else None), f(tot), f(lc[0] if lc else None, 3), f(lc[1] if lc else None, 3),
                f(lc[2] if lc else None, 3), f(wc[0], 3), f(wc[1], 3), f(wc[2], 3), f(e["share"], 3), f(td, 3),
                round(self.cx.alt[i], 3), round(self.cx.az[i], 3)]

    # ---------------------------------------------------------------- the beam band (see MODELS.md, G1)
    def main_width(self, i, pose):
        """(width m of the largest room patch (0 = no room light), evaluate() dict)."""
        e = self.evaluate(i, pose)
        best = None
        for ti, a, tr in e["faces"]:
            if self.room_face[ti] and (best is None or a > best[0]):
                best = (a, tr[1] * self.sc)
        return (best[1] if best else 0.0), e

    def band_center(self, i, a, W, r0=None, tol=0.005, rmin=0.0, rmax=15.0):
        """Ring r* (at aperture a, step i) where the main room patch is exactly W wide - the centre of the thin-beam
        band, just above the light edge.  Bisection (width grows with ring away from the edge).  Returns
        (r*, evaluate dict) or (None, None).  r0 = warm start (previous minute's answer)."""
        lo = hi = None
        if r0 is not None:
            a_, b_ = max(rmin, r0 - 0.3), min(rmax, r0 + 0.3)
            wa, ea = self.main_width(i, (a_, a))
            if wa < W:
                wb, eb = self.main_width(i, (b_, a))
                if wb >= W:
                    lo, hi, ehi = a_, b_, eb
            elif self.band_any:                            # narrowing branch (tilted windows), M-011
                wb, eb = self.main_width(i, (b_, a))
                if wb < W:
                    return self._bisect_down(i, a, W, a_, b_, ea, tol)
        if lo is None:
            w0, e0 = self.main_width(i, (rmin, a))
            if w0 >= W and not self.band_any:
                return None, None
            w1, e1 = self.main_width(i, (rmax, a))
            if w0 < W <= w1:
                lo, hi, ehi = rmin, rmax, e1
            elif not self.band_any:
                return None, None
            else:
                # M-011 (tilted windows): the strip may NARROW as the ring rises (the rings tilt the tubes away from
                # the sun) or peak inside the range - scan 5 rings, take the first widening bracket, else narrowing
                rs = [rmin + (rmax - rmin) * k / 4.0 for k in range(5)]
                ws = [(w0, e0)] + [self.main_width(i, (r, a)) for r in rs[1:4]] + [(w1, e1)]
                for k in range(4):
                    if ws[k][0] < W <= ws[k + 1][0]:
                        lo, hi, ehi = rs[k], rs[k + 1], ws[k + 1][1]
                        break
                if lo is None:
                    for k in range(4):
                        if ws[k + 1][0] < W <= ws[k][0]:
                            return self._bisect_down(i, a, W, rs[k], rs[k + 1], ws[k][1], tol)
                    return None, None
        while hi - lo > tol:
            m = 0.5 * (lo + hi)
            wm, em = self.main_width(i, (m, a))
            if wm < W:
                lo = m
            else:
                hi, ehi = m, em
        return hi, ehi

    def _bisect_down(self, i, a, W, lo, hi, elo, tol):
        """Narrowing branch: width >= W at lo, < W at hi -> the ring where it is exactly W (the wide side)."""
        while hi - lo > tol:
            m = 0.5 * (lo + hi)
            wm, em = self.main_width(i, (m, a))
            if wm >= W:
                lo, elo = m, em
            else:
                hi = m
        return lo, elo

    # ---------------------------------------------------------------- where the light is
    def light_centroid(self, faces, room_only=True):
        """Area-weighted 3D centroid of the lit patches (room faces only by default) and the lit area (m2)."""
        tot = 0.0
        cx = cy = cz = 0.0
        for ti, a, tr in faces:
            if room_only and not self.room_face[ti]:
                continue
            tot += a
            cx += a * tr[6]
            cy += a * tr[7]
            cz += a * tr[8]
        if tot <= 0:
            return None, 0.0
        return (cx / tot, cy / tot, cz / tot), tot

    def target_point(self, depth_m, offset_m, height_m=0.0):
        k = 1.0 / self.sc
        base = C.vadd(self.exitp, C.vadd(C.vscale(self.r, depth_m * k), C.vscale(self.w, offset_m * k)))
        lift = (self.floor_u + height_m * k) - C.vdot(base, self.u)
        return C.vadd(base, C.vscale(self.u, lift))

    def window_coords(self, p):
        """(depth into the room, offset to the right looking out, height above the floor) in metres of a point."""
        d = C.vsub(p, self.exitp)
        return (C.vdot(d, self.r) * self.sc, C.vdot(d, self.w) * self.sc,
                (C.vdot(p, self.u) - self.floor_u) * self.sc)

    def distance_to_lit(self, T, rec_polys, room_only=True):
        """Distance (m) from point T to the nearest lit point (record() with polys=True), or None if unlit."""
        best = None
        if rec_polys is None:
            return None
        for ti, tr, lit in rec_polys[1]:
            if room_only and not self.room_face[ti]:
                continue
            t = self.targets[ti]
            d = C.vsub(T, t.q0)
            a, b = C.vdot(d, t.e1), C.vdot(d, t.e2)
            hgt = C.vdot(d, t.n)
            for poly in lit:
                dd = _point_poly_dist2d((a, b), poly)
                dist = math.sqrt(dd * dd + hgt * hgt) * self.sc
                if best is None or dist < best:
                    best = dist
        return best


def room_axes(scene):
    """(r, w, u): r = horizontal INTO the room, perpendicular to the tube axis (unchanged by tilt; scene['room_dir']
    for tilted windows, see scenario.tilt_scene), w = across (to the RIGHT looking out), u = world up.  For a wall
    window this is exactly window_beams_core._room_axes."""
    up = (0.0, 0.0, 1.0)
    if scene.get("room_dir"):
        r = C.vunit(tuple(scene["room_dir"]))
    else:
        ext = C.vunit(tuple(scene["exterior_dir"]))
        if math.hypot(ext[0], ext[1]) > 0.2:
            r = C.vunit((-ext[0], -ext[1], 0.0))
        else:                                          # a skylight without room_dir: across the tube axis
            import scenario as _SC
            r = tuple(_SC.room_dir(scene))
    w = C.vunit(C.vcross(up, r))
    return r, w, C.vcross(r, w)


def _clip_half(poly, h):
    """Part of a 2D polygon with ca*a + cb*b + c0 <= 0 (Sutherland-Hodgman, one plane)."""
    ca, cb, c0 = h
    out = []
    n = len(poly)
    for k in range(n):
        p, q = poly[k], poly[(k + 1) % n]
        fp, fq = ca * p[0] + cb * p[1] + c0, ca * q[0] + cb * q[1] + c0
        if fp <= 0:
            out.append(p)
        if (fp < 0 < fq) or (fq < 0 < fp):
            s = fp / (fp - fq)
            out.append((p[0] + (q[0] - p[0]) * s, p[1] + (q[1] - p[1]) * s))
    return out if len(out) >= 3 else None


def _point_poly_dist2d(p, poly):
    """Distance from p to a convex polygon (0 inside; either winding)."""
    pos = neg = False
    n = len(poly)
    best = None
    for k in range(n):
        a, b = poly[k], poly[(k + 1) % n]
        ex, ey = b[0] - a[0], b[1] - a[1]
        cross = ex * (p[1] - a[1]) - ey * (p[0] - a[0])
        if cross > 0:
            pos = True
        elif cross < 0:
            neg = True
        L2 = ex * ex + ey * ey
        t = 0.0 if L2 <= 0 else max(0.0, min(1.0, ((p[0] - a[0]) * ex + (p[1] - a[1]) * ey) / L2))
        qx, qy = a[0] + t * ex, a[1] + t * ey
        d = math.hypot(p[0] - qx, p[1] - qy)
        if best is None or d < best:
            best = d
    if n >= 3 and not (pos and neg):
        return 0.0
    return best if best is not None else 0.0
