# -*- coding: utf-8 -*-
"""
scenario.py - ANY window, not just today's (Felix 2026-10-11: "the model should work for windows at all directions",
"different widths or possibly fully scaled down").  A scenario is

    {"facing": compass bearing the window faces (0 N, 90 E, 180 S, 270 W; None = the scene's own),
     "level_m": window raised (+) / lowered (-) relative to the room floor (raised sill / higher window), m,
     "scale": uniform scale of the WHOLE window assembly (tube, rings, blades, aperture, frame lip) about its wall
              centre (exit point), 1 = today,
     "width_f": extra scale of the window's overall width (along the tube axis), 1 = today,
     "beam_w": the year-round beam width (m),
     "tilt_deg": 0 = today's wall window (facing sideways) ... 90 = skylight in the ceiling (facing straight up)}

TILT (Felix 2026-10-11: "windows on the ceiling or at any angle from facing sideways to facing upwards"): the whole
window assembly turns about its own tube axis so it faces UP by tilt_deg, and sits where its plane touches a quarter
circle of radius R = (window centre below the ceiling) drawn in the wall/ceiling corner.  0 deg = today; 90 deg = a
skylight in the ceiling plane, its centre R in from the wall (level_m moves it: higher window -> skylight nearer the
wall); between = a sloped window across the wall/ceiling corner (a chamfer, roof light).  Tube, rings, blades, lip and frame all
rotate with it (the planner and physics work in the window's own axes); the ROOM stays as it is (floor + walls +
ceiling).  Facing keeps its meaning: the compass direction the window would face if swung back down to vertical - for
a skylight it is just its rotation in the ceiling plane (hinge side towards that direction).


Everything is derived from the GH model's own geometry (scene.json: joints, bodies, exit point, exterior direction,
room faces) - the variant is built IN MEMORY and written only under lab\\out\\scen\\ (Felix's files untouched).
Facing works like Felix's Facing input (the engine's apply_orientation: the sun is rotated, the model stays put).

LITE scan: the sun table, the shading context and the per-pose diffuse factors (sky / ground) - but NOT the
121-pose x sun-bin records (a full scan per window is slow and needs ~1 GB).  The lab physics is exact at any pose and
any minute and does not need them; planner candidates are then all beam-wanted minutes.

    load(spec) -> (ev, P, job_path)        make_pack(spec, possible_minutes) -> tester-harness pack path
"""
from __future__ import division

import copy
import hashlib
import json
import math
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
if HERE not in sys.path:
    sys.path.insert(0, HERE)
import labcore as L  # noqa: E402

C = L.C
BASE_JOB = os.path.join(L.LAB, "jobs", "w005", "job.json")
OUT = os.path.join(L.LAB, "out", "scen")
DEFAULT = {"facing": None, "level_m": 0.0, "scale": 1.0, "width_f": 1.0, "beam_w": 0.05, "tilt_deg": 0.0}


def norm(spec):
    s = dict(DEFAULT)
    s.update(spec or {})
    return s


def sid(spec):
    s = norm(spec)
    f = "own" if s["facing"] is None else "%03d" % int(round(s["facing"]))
    out = "f%s_l%+.2f_s%.2f_wf%.2f_b%03d" % (f, s["level_m"], s["scale"], s["width_f"], int(round(s["beam_w"] * 1000)))
    if s.get("tilt_deg"):                              # ids of vertical windows stay as they were
        out += "_t%02d" % int(round(s["tilt_deg"]))
    return out


def _v(a, b, op):
    return tuple(op(x, y) for x, y in zip(a, b))


def transform_scene(scene, level_m=0.0, scale=1.0, width_f=1.0):
    """Scale the window assembly about its exit point (uniform `scale`, plus `width_f` along the tube axis = the first
    ring joint's axis) and raise it by level_m (model up = +z).  Room faces stay as they are."""
    s = copy.deepcopy(scene)
    if level_m == 0 and scale == 1 and width_f == 1:
        return s
    A = tuple(s["exit_point"])
    j1 = s["joints"][0]["axis"]
    n = math.sqrt(sum(x * x for x in j1))
    ax = tuple(x / n for x in j1)
    up = (0.0, 0.0, 1.0)
    k_al = scale * width_f

    def V(v):
        al = sum(a * b for a, b in zip(v, ax))
        rest = tuple(a - b * al for a, b in zip(v, ax))
        return tuple(b * al * k_al + r * scale for b, r in zip(ax, rest))

    def P(p):
        d = V(_v(p, A, lambda x, y: x - y))
        return tuple(a + b + u * level_m for a, b, u in zip(A, d, up))

    for j in s["joints"]:
        j["origin"] = list(P(j["origin"]))
        j["axis"] = list(V(j["axis"]))
    for b in s["bodies"]:
        for key in ("mouths", "blockers"):
            if b.get(key):
                b[key] = [[list(P(p)) for p in poly] for poly in b[key]]
    s["exit_point"] = list(P(s["exit_point"]))
    return s


UP = (0.0, 0.0, 1.0)


def _unit(v):
    n = math.sqrt(sum(x * x for x in v))
    return tuple(x / n for x in v)


def _dot(a, b):
    return sum(x * y for x, y in zip(a, b))


def _rot_vec(v, k, ang):
    """Rodrigues: v rotated by ang (rad) about the unit axis k."""
    c, s = math.cos(ang), math.sin(ang)
    kx = (k[1] * v[2] - k[2] * v[1], k[2] * v[0] - k[0] * v[2], k[0] * v[1] - k[1] * v[0])
    kd = _dot(k, v) * (1.0 - c)
    return tuple(v[i] * c + kx[i] * s + k[i] * kd for i in range(3))


def room_dir(scene):
    """Horizontal unit vector INTO the room, perpendicular to the tube axis (it does not change when the window
    tilts about that axis).  Taken from the scene's 'room_dir' when present (tilted windows), else from the
    exterior direction (vertical windows), else the horizontal normal of the tube axis pointing at the room."""
    if scene.get("room_dir"):
        return _unit(scene["room_dir"])
    ext = _unit(scene["exterior_dir"])
    h = (-ext[0], -ext[1], 0.0)
    if math.hypot(h[0], h[1]) > 0.2:
        return _unit(h)
    a = _unit(scene["joints"][0]["axis"])
    h = _unit((a[1], -a[0], 0.0))
    pts = [p for t in scene.get("targets", []) for p in [t["origin"]]]
    if pts:
        c = tuple(sum(p[i] for p in pts) / len(pts) for i in range(3))
        if _dot(h, tuple(c[i] - scene["exit_point"][i] for i in range(3))) < 0:
            h = tuple(-x for x in h)
    return h


def ceiling_level(scene):
    """Height (model z) of the ceiling: the largest room face that looks straight down; else the top of the room."""
    best = None
    for t in scene.get("targets", []):
        n = _unit(t["normal"])
        if _dot(n, UP) < -0.9:
            b = t["boundary"]
            a = 0.0
            for k in range(len(b)):
                a += b[k][0] * b[(k + 1) % len(b)][1] - b[(k + 1) % len(b)][0] * b[k][1]
            a = abs(a) / 2.0
            if best is None or a > best[0]:
                best = (a, _dot(t["origin"], UP))
    if best is not None:
        return best[1]
    zs = [_dot(t["origin"], UP) for t in scene.get("targets", [])]
    return max(zs) if zs else _dot(scene["exit_point"], UP) + 1.5 / float(scene.get("units_to_m", 1.0))


def tilt_scene(scene, tilt_deg):
    """Swing the window assembly up by tilt_deg about the hinge line (tube-axis direction, on the window's exit plane
    at ceiling height).  0 = unchanged; 90 = skylight in the ceiling plane.  Room faces stay.  Records
    scene['room_dir'] (horizontal into the room), 'tilt_deg' and 'hinge' for the physics / reports."""
    s = copy.deepcopy(scene)
    rdir = room_dir(scene)
    s["room_dir"] = list(rdir)
    t = float(tilt_deg or 0.0)
    s["tilt_deg"] = t
    if abs(t) < 1e-9:
        return s
    if not (0.0 < t <= 90.0):
        raise ValueError("tilt_deg must be within 0..90 (0 = wall, 90 = ceiling), got %r" % tilt_deg)
    E = tuple(s["exit_point"])
    ext = _unit(s["exterior_dir"])
    if abs(_dot(ext, UP)) > 0.05:
        raise ValueError("tilt needs a vertical (wall) window to start from; this one already looks %.0f deg up"
                         % math.degrees(math.asin(max(-1.0, min(1.0, _dot(ext, UP))))))
    a = _unit(s["joints"][0]["axis"])
    a = _unit((a[0], a[1], 0.0))                          # tilt axis: horizontal, along the tube
    zc = ceiling_level(s)
    R = zc - _dot(E, UP)                                  # exit point below the ceiling = skylight distance from wall
    if R * float(s.get("units_to_m", 1.0)) < 0.05:
        raise ValueError("the window centre is at / above the ceiling (%.2f m) - lower it to tilt it" % R)
    ang = math.radians(t)
    if _dot(_rot_vec(ext, a, ang), UP) < 0:
        ang = -ang                                        # turn so the exterior direction swings UP
    # placement: the window plane is tangent to a quarter circle of radius R about O = E + R*room_dir (R in from the
    # wall, R below the ceiling): t=0 -> E itself (today), t=90 -> in the ceiling plane R in from the wall, between ->
    # a chamfer from the wall to the ceiling with the window at its middle.  (A rigid swing about the wall head cannot
    # do this: swinging the bottom in turns the window to face DOWN, swinging it out puts it outside the room.)
    n_t = _rot_vec(ext, a, ang)
    O = tuple(E[i] + R * rdir[i] for i in range(3))
    Et = tuple(O[i] + R * n_t[i] for i in range(3))
    H = O

    def V(v):
        return list(_rot_vec(tuple(v), a, ang))

    def Pt(p):
        d = _rot_vec(tuple(p[i] - E[i] for i in range(3)), a, ang)
        return [Et[i] + d[i] for i in range(3)]

    for j in s["joints"]:
        j["origin"] = Pt(j["origin"])
        j["axis"] = V(j["axis"])
    for b in s["bodies"]:
        for key in ("mouths", "blockers"):
            if b.get(key):
                b[key] = [[Pt(p) for p in poly] for poly in b[key]]
    s["exit_point"] = Pt(s["exit_point"])
    s["exterior_dir"] = V(s["exterior_dir"])
    # the sloped window cuts the wall/ceiling corner off the room: trim every room face to the room side of the
    # window's plane (drop faces wholly behind it), so EVERY tracer (the tester's too) sees the same room
    trimmed, dropped, kept = 0, 0, []
    for td in s.get("targets", []):
        T = C.Target(td)
        h = (_dot(T.e1, n_t), _dot(T.e2, n_t), _dot(tuple(T.q0[i] - Et[i] for i in range(3)), n_t))
        out = [h[0] * a + h[1] * b + h[2] for (a, b) in T.boundary]
        if max(out) <= 1e-6:
            kept.append(td)
            continue
        if min(out) >= -1e-6:
            dropped += 1
            continue
        import physics as _PH
        poly = _PH._clip_half(list(T.boundary), h)
        if not poly:
            dropped += 1
            continue
        td = dict(td, origin=list(T.q0), e1=list(T.e1), e2=list(T.e2), normal=list(T.n),
                  boundary=[list(p) for p in poly])
        kept.append(td)
        trimmed += 1
    if s.get("targets"):
        s["targets"] = kept
    s["hinge"] = {"faces_trimmed": trimmed, "faces_dropped": dropped,"arc_centre": list(H), "axis": list(a), "radius_m": R * float(s.get("units_to_m", 1.0)),
                  "rule": "window plane tangent to a quarter circle (radius = window centre below the ceiling)"}
    return s


def _geo_sig(scene, P):
    blob = json.dumps({"j": scene["joints"], "b": scene["bodies"], "e": scene["exterior_dir"],
                       "g": P["include_ground"], "params": scene["params"]}, sort_keys=True)
    return hashlib.md5(blob.encode("utf-8")).hexdigest()[:12]


def lite_scan(scene, w, P):
    mech = C.Mechanism(scene)
    ext = C.vunit(tuple(scene["exterior_dir"]))
    north = scene.get("north_deg", 0.0)
    targets = [C.Target(t) for t in scene["targets"]]
    sun = C.SunTable(w, north, ext, P)
    poses = mech.poses(None)
    key_frac = [1.0] * len(sun.keys)
    ctx = C.build_context(scene, P)
    if ctx is not None:
        pts = C.window_front_points(scene)
        skip = P["context_skip_m"] / float(scene.get("units_to_m", 1.0))
        key_frac = [ctx.fraction(s, pts, skip) for s in sun.keys]
    if not os.path.isdir(OUT):
        os.makedirs(OUT)
    dpath = os.path.join(OUT, "diffuse_%s.json" % _geo_sig(scene, P))      # diffuse factors: geometry only
    if os.path.exists(dpath):
        d = json.load(open(dpath))
        dsky, dgnd = d["dsky"], d["dgnd"]
    else:
        patches = C.tregenza_patches()
        sky = [(C.sun_basis(dd), om) for dd, om in patches if C.vdot(dd, ext) > 0.02]
        gnd = [(C.sun_basis((dd[0], dd[1], -dd[2])), om) for dd, om in patches
               if C.vdot((dd[0], dd[1], -dd[2]), ext) > 0.02] if P["include_ground"] else []
        dsky, dgnd = [], []
        for pose in poses:
            mouths, blockers = mech.pose_geometry(pose)
            for lst, out in ((sky, dsky), (gnd, dgnd)):
                tot = 0.0
                for basis, om in lst:
                    pcs = C.clear_pieces(mouths, blockers, basis)
                    if pcs:
                        tot += om * sum(C.poly_area(pc) for pc in pcs)
                out.append(tot)
        with open(dpath, "w") as f:
            json.dump({"dsky": dsky, "dgnd": dgnd}, f)
    res = C.ScanResult()
    res.sun, res.mech, res.targets = sun, mech, targets
    res.poses = [tuple(p) for p in poses]
    res.rec = [[None] * len(sun.keys) for _ in poses]
    res.dsky, res.dgnd, res.key_frac, res.context_info = dsky, dgnd, key_frac, None
    res.signature, res.seconds, res.from_cache = "lite", 0.0, True
    return res


def load(spec, base_job=BASE_JOB, scene_edit=None):
    """-> (ev, P, job_path).  scene_edit(scene) -> scene: an extra in-memory edit (e.g. blade.variant_scene)."""
    s = norm(spec)
    job = json.load(open(base_job))
    for k in ("scene_path", "epw_path", "out_dir", "cache_dir"):   # job paths are relative to the project root
        if job["params"].get(k) and not os.path.isabs(job["params"][k]):
            job["params"][k] = os.path.join(L.CORE_DIR, job["params"][k])
    P0 = job["params"]
    scene = json.load(open(P0["scene_path"]))
    scene = transform_scene(scene, s["level_m"], s["scale"], s["width_f"])
    ov = dict(P0["overrides"])
    ov["beam_width_m"] = float(s["beam_w"])
    if s["facing"] is not None:
        ov["facing_deg"] = float(s["facing"]) % 360.0
    P = C.make_params(ov)
    w = C.prepare_weather(C.read_epw(P0["epw_path"]), P["time_step_min"])
    tilt = float(s.get("tilt_deg") or 0.0)
    if tilt:
        # facing = the compass direction of the window swung back to vertical: set north on the UNtilted window
        # (the engine's apply_orientation ignores a window that looks straight up), then tilt the model
        if s["facing"] is not None:
            C.apply_orientation(scene, float(s["facing"]) % 360.0)
        scene = tilt_scene(scene, tilt)
        sc = C.resolve_scene(scene, dict(P, facing_deg=None))
    else:
        sc = C.resolve_scene(scene, P)
    win = C.split_windows(sc)[0]
    if scene_edit:
        win = scene_edit(win)
    scan = lite_scan(win, w, P)
    ev = C.Evaluator(win, w, scan, P)
    # a job file for this scenario (the tester's tools read width / weather / scene from it)
    d = os.path.join(OUT, sid(s))
    if not os.path.isdir(d):
        os.makedirs(d)
    with open(os.path.join(d, "scene.json"), "w") as f:
        json.dump(scene, f)
    jp = os.path.join(d, "job.json")
    j2 = copy.deepcopy(job)
    j2["params"]["scene_path"] = os.path.join(d, "scene.json")
    j2["params"]["overrides"] = ov
    j2["params"]["cache_dir"] = os.path.join(L.LAB, "out", "_scan_cache")
    j2["scenario"] = s
    with open(jp, "w") as f:
        json.dump(j2, f, indent=1)
    return ev, P, jp


def make_pack(spec, possible_minutes, pack_path, ev=None, P=None, scene_edit=None, extra_param=None):
    """The tester's trace pack (trace.make_pack fields) for a lite scenario.  'possible' = minutes where the
    PLANNER found a clean band (5 aperture knots, exact sun) - a planner-side estimate, NOT the engine's grid search."""
    if ev is None:
        ev, P, _jp = load(spec, scene_edit=scene_edit)
    win = ev.scene
    wanted = [i for i in range(ev.n) if ev.beam_on[i]]
    ps = set(possible_minutes)
    names = list(ev.pose_names)
    ranges = list(ev._ranges())
    pack = {"job": "", "core_dir": L.CORE_DIR, "width_m": P["beam_width_m"], "win": win,
            "P": dict((k, v) for k, v in P.items() if not k.startswith("_") and _jsonable(v)),
            "room_m2": ev.room_m2, "scale": ev.scale, "a_scale": ev.a_scale, "tw": ev._tw,
            "dsky": list(ev._dsky), "dgnd": list(ev._dgnd), "poses": [list(p) for p in ev.scan.poses],
            "names": names, "ranges": ranges, "n_steps": ev.n, "step_min": ev.step_min,
            "wanted": wanted, "possible": [1 if i in ps else 0 for i in wanted],
            "possible_rule": "planner clean band (lite scenario)",
            "dni": [ev.w["dni"][i] for i in wanted], "dhi": [ev.w["dhi"][i] for i in wanted],
            "ghi": [ev.w["ghi"][i] for i in wanted],
            "targets": [t.name for t in ev.scan.targets], "target_normals": [list(t.n) for t in ev.scan.targets]}
    with open(pack_path, "w") as f:
        json.dump(pack, f)
    return pack_path


def _jsonable(v):
    try:
        json.dumps(v)
        return True
    except (TypeError, ValueError):
        return False
