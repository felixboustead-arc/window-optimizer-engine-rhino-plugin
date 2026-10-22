# -*- coding: utf-8 -*-
"""
lab\\engine\\engine.py - THE SHARED R3 ENGINE, one entry point (MODELER-D, D3-001, 2026-10-20).

Assembly only: NO new core algorithm.  Every result comes from the accepted champion parts, called unchanged:
    LIGHT    = G2R3T band glide (A line), as wired in the R2 B entry (real EPW, DNI gate 120, one plan per window
               class)  -> entries\\R2\\B\\entry.py plan_light (frozen folder, imported by path, never edited)
    THERMAL  = HB-L2 (B)                                   -> entries\\R2\\B\\entry.py thermal_model
               (B3-001 is ACCEPTED but lives in the measured-house DRIVER metrics\\tr2_measured.py (Annex 71 inputs),
               not in the HB-L2 core; Felix's room (thermal_case A04_felix_open) is unaffected)
    COMBINED = B's nested knapsack                          -> entries\\R2\\B\\entry.py plan_combined
    (D3-001b, AU-R3-001 L1: the A3-001 light option is REMOVED - J-R3-002 rejected it; only accepted parts remain.)

What D adds (the D-owned area: multi-window, GH integration, per-minute runtime):
    plan(scene, epw, params, slider, mode)  one call for 1-15 windows (light or combined), plans per window class
    Schedule(plans).poses(minute)           per-minute pose lookup for every window (bisect on the plan runs, O(log n))
    minute_of_year(month, day, hour, minute) AEST clock (no DST, non-leap year) -> plan step index
    use_cache_dir(path)                     redirect EVERY disk cache of the engine (cold-cache timing, T-008g lesson)
    code_signature()                        engine VERSION + mtimes of every code file a plan depends on (cache keys)
    dayahead_days(doy)                      [d-1, d, d+1] (mod 365): day d of this plan == day d of the full-year plan
    CLI: python -B engine.py plan <scene.json> <out.json> [--days 171,354] [--slider 0.5] [--nproc 6]
         (the GH component runs this in a background CPython process for full-year plans)

Known (AU-R3-001 L3): the champion registers generic module names (scenario, labcore, physics, glide_champion) in
sys.modules; inside Rhino's shared CPython another plug-in with the same names would collide, and edits to these
files need a Rhino restart (modules are cached).  The GH component shows engine.code_signature() so a stale module
is visible.
Plan contract v1 (unchanged): {"runs": [[step, ring_deg, aperture_pct], ...], "n_steps": 525600, "step_min": 1, ...}.
"""
from __future__ import division

import bisect
import importlib.util
import json
import os
import sys
import time

HERE = os.path.dirname(os.path.abspath(__file__))
LAB = os.path.dirname(HERE)
B_ENTRY = os.path.join(LAB, "entries", "R2", "B", "entry.py")
A3_ENTRY = os.path.join(LAB, "entries", "R3", "A", "entry.py")
PLAN_DIR = os.path.join(LAB, "plans", "engine")          # lab\plans is a junction to E:
VERSION = "ENGINE-R3 D3-001b (LIGHT G2R3T champion, THERMAL HB-L2, COMBINED B nested knapsack)"
CODE_FILES = [os.path.abspath(__file__), B_ENTRY, os.path.join(LAB, "entries", "R2", "B", "glide_champion.py"),
              os.path.join(LAB, "entries", "R2", "B", "thermal_hb2.py"), os.path.join(LAB, "models", "scenario.py"),
              os.path.join(LAB, "models", "physics.py"), os.path.join(LAB, "models", "labcore.py"),
              os.path.join(HERE, "o13.py")]
N_STEPS = 525600
DAYS_BEFORE = [0, 31, 59, 90, 120, 151, 181, 212, 243, 273, 304, 334]

_MOD = {}


def _load(name, path):
    if name not in _MOD:
        spec = importlib.util.spec_from_file_location(name, path)
        m = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(m)
        _MOD[name] = m
    return _MOD[name]


def champion():
    """The R2 B entry module (champion light wrapper + HB-L2 + nested knapsack)."""
    return _load("engine_champion_B_R2", B_ENTRY)


def use_cache_dir(path):
    """Send every disk cache / output of the engine to `path` and drop the in-memory class cache.
    Disk caches found (grep): models\\scenario.py OUT (diffuse factors per geometry) and the B entry OUT (job / plan /
    pack files).  The champion has no plan cache.  Used by the timing protocol (fresh empty folder per repetition).
    AU-R3-001 M9: entry.py:166 writes cache_dir = lab\\out\\_scan_cache into job.json, but the planning path never reads
    it (its only consumer is labcore.load -> C.run_scan, used by glide_champion.run, which the engine does not call;
    the folder does not even exist).  engine_multi.rep records the folder's state before / after every cold rep."""
    B = champion()
    for sub in ("scen", "plans"):
        d = os.path.join(path, sub)
        if not os.path.isdir(d):
            os.makedirs(d)
    sys.modules["scenario"].OUT = os.path.join(path, "scen")
    B.OUT = os.path.join(path, "plans")
    B._EV.clear()
    _MOD["_cache_dir"] = path


# ============================================================================================== planning
def windows(scene):
    """Normalised window list [(id, spec, raw)] exactly as the champion reads it (1-15 windows)."""
    return champion().windows_of(scene)


def classes(scene):
    """{class_key: [window ids]}: windows of one class share one light plan (facing, tilt, level, scale, width, beam)."""
    B = champion()
    out = {}
    for wid, spec, raw in B.windows_of(scene):
        out.setdefault(B.class_key(spec), []).append(wid)
    return out


def plan(scene, weather_epw=None, params=None, slider=1.0, mode=None, construction=None, building=None):
    """{window_id: plan} for 1-15 windows.
    mode "light"    -> LIGHT only (champion G2R3T).
    mode "combined" -> COMBINED nested knapsack at `slider` (0 comfort .. 1 light, or [[minute, value], ...]).
    Default: "light" when slider == 1 (identical plans, no thermal runs), else "combined".
    params: the champion's (param_ranges()) + "days" (day indices, 0 = 1 Jan).
    building: O13 building description (o13.py: zones, doors, vents, MVHR, construction class).  None (default) or
    {"o13": False} = O13 OFF: the R2 champion path, bit-identical to D3-001.  O13 acts on THERMAL / COMBINED only
    (the light planner is unchanged, R3_BRIEF sec. 1); in COMBINED the knapsack's comfort runs use the O13 building."""
    params = dict(params or {})
    light = params.pop("light", "champion")
    if light != "champion":                              # AU-R3-001 L1: only accepted parts (A3-001 rejected)
        raise ValueError("light planner %r is not available: the shared engine runs the accepted champion only" % light)
    p4b = bool(params.pop("light_p4b", False))           # PF-SUBMIT-001c: opt-in LUMEN-P4b span guard, default OFF
    p4b_opt = params.pop("light_p4b_opts", None) or {}   # {"cache_dir": ..., "workers": 1}
    if mode is None:
        mode = "light" if (not isinstance(slider, (list, tuple)) and float(slider) >= 1.0) else "combined"
    if p4b and mode != "light":
        raise ValueError("light_p4b applies to mode 'light' only (LUMEN-P4b is a light-plan layer)")
    t0 = time.time()
    o13_on = bool(building) and building.get("o13", True) and mode == "combined"
    if mode == "light":
        res = champion().plan_light(scene, weather_epw, params)
    elif o13_on:
        O = _load("engine_o13", os.path.join(HERE, "o13.py"))
        res = O.plan_combined(champion(), scene, weather_epw, construction, params, slider, building)
    else:
        res = champion().plan_combined(scene, weather_epw, construction, params, slider)
    if p4b:                                              # lazy import: OFF never touches pathfinder
        P = _load("engine_p4b", os.path.join(LAB, "pathfinder", "light_o1", "p4b_engine.py"))
        p4b_res, by_class = {}, {}
        for wid, p in res.items():                       # one P4b run per window class, own copy per window
            k = tuple(p["class_key"])
            if k not in by_class:
                by_class[k] = P.apply(p, days=params.get("days"), workers=int(p4b_opt.get("workers", 1)),
                                      cache_dir=p4b_opt.get("cache_dir"))
            q = dict(by_class[k])
            q["window_id"] = wid
            p4b_res[wid] = q
        res = p4b_res
    sec = round(time.time() - t0, 1)
    days = params.get("days")
    for p in res.values():
        p["engine"] = VERSION
        p["engine_mode"] = mode
        p["engine_light"] = light
        p["engine_sec"] = sec
        p["engine_days"] = sorted(int(d) for d in days) if days is not None else None   # AU-R3-001 L2
        p["engine_o13"] = bool(o13_on)
        if p4b:
            p["engine_p4b"] = P.LABEL
        if mode == "combined":
            p["engine_comfort_room"] = ("O13 building %s" % (building or {}).get("name", "?")) if o13_on else \
                ((construction or {}).get("case") or "A04_felix_open (default thermal_case)")
    return res


def code_signature():
    """Engine VERSION + mtimes of every code file a plan depends on (AU-R3-001 M6: plan cache keys)."""
    sig = [VERSION]
    for f in CODE_FILES:
        try:
            sig.append("%s:%d" % (os.path.basename(f), int(os.path.getmtime(f))))
        except OSError:
            sig.append("%s:-" % os.path.basename(f))
    return sig


def dayahead_days(doy):
    """Rolling day-ahead plan days for clock day `doy`: [d-1, d, d+1] (mod 365).  G2R3T plans each day on its own,
    and G.assemble only links day d to its planned neighbours (pre-dawn glide from d-1, evening glide to d+1), so every
    minute of day d equals the full-year plan (AU-R3-001 M5; engine_multi.day_ahead_check tests it)."""
    d = int(doy) % 365
    return sorted(set([(d - 1) % 365, d, (d + 1) % 365]))


def thermal_model():
    return champion().thermal_model()


def param_ranges():
    pr = champion().param_ranges()
    pr["engine"] = [
        {"name": "days", "default": None, "units": "day index list",
         "meaning": "plan only these days (0 = 1 Jan); None = full year. Light day plans are independent per day"},
        {"name": "mode", "default": "auto", "choices": ["auto", "light", "combined"], "units": "-",
         "meaning": "auto = light when slider is 1, else combined"},
        {"name": "building", "default": None, "units": "dict",
         "meaning": "O13 building (zones, doors, vents, MVHR, construction class; see o13.py); None = O13 off"}]
    try:
        pr["o13"] = _load("engine_o13", os.path.join(HERE, "o13.py")).param_table()
    except Exception:  # noqa - o13 optional for the R2 path
        pass
    return pr


# ============================================================================================== per-minute runtime
def minute_of_year(month, day, hour=0, minute=0):
    """AEST wall clock (UTC+10, no DST, non-leap year: the EPW / plan timebase, A3-REQ5) -> plan step index.
    Invalid dates / times raise ValueError (AU-R3-001 M8; 29 Feb is not in the non-leap plan year)."""
    month, day, hour, minute = int(month), int(day), int(hour), int(minute)
    mdays = [31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31]
    if not 1 <= month <= 12 or not 1 <= day <= mdays[month - 1] or not 0 <= hour <= 23 or not 0 <= minute <= 59:
        raise ValueError("invalid AEST clock month %d day %d %02d:%02d (non-leap year)" % (month, day, hour, minute))
    d = DAYS_BEFORE[month - 1] + day - 1
    return (d * 1440 + int(hour) * 60 + int(minute)) % N_STEPS


class Schedule(object):
    """Per-minute pose lookup over a set of plans.  pose(window_id, minute) follows labcore.poses_from_runs exactly
    (the run in force at the START of the step; a run starting inside a step counts from that step).
    AU-R3-001 L2: a window without runs is listed in .missing (its pose is the open default, never silently "planned");
    planned(wid, minute) is False outside the plan's engine_days (the pose there is only the held / glided one)."""

    def __init__(self, plans, order=None):
        self.order = list(order or plans.keys())
        self._s, self._v, self._days = {}, {}, {}
        self.missing = []
        for wid in self.order:
            p = plans.get(wid) or {}
            runs = p.get("runs")
            if not runs:
                self.missing.append(wid)
                runs = [[0, 15.0, 100.0]]
            self._s[wid] = [r[0] - 0.5 for r in runs]
            self._v[wid] = [(float(r[1]), float(r[2])) for r in runs]
            dd = p.get("engine_days")
            self._days[wid] = None if dd is None else set(int(x) for x in dd)

    def planned(self, wid, minute):
        if wid in self.missing:
            return False
        dd = self._days.get(wid)
        return dd is None or (int(minute) // 1440) in dd

    def pose(self, wid, minute):
        k = bisect.bisect_right(self._s[wid], minute) - 1
        return self._v[wid][max(0, k)]

    def poses(self, minute):
        """[(ring_deg, aperture_pct)] in window order at this step."""
        return [self.pose(w, minute) for w in self.order]


def save_plans(plans, path):
    d = os.path.dirname(path)
    if d and not os.path.isdir(d):
        os.makedirs(d)
    tmp = path + ".part"
    with open(tmp, "w") as f:
        json.dump(plans, f, separators=(",", ":"), default=str)
    os.replace(tmp, path)
    return path


# ============================================================================================== CLI
def _cli(argv):
    """plan <scene.json> <out.json> [--days a,b] [--slider s] [--nproc n] [--light X] [--epw path] [--cache dir]"""
    if len(argv) < 3 or argv[0] != "plan":
        print(_cli.__doc__)
        return 2
    scene = json.load(open(argv[1]))
    out = argv[2]
    opt = dict(zip(argv[3::2], argv[4::2]))
    params = dict(scene.get("params") or {}) if isinstance(scene, dict) else {}
    if "--days" in opt and opt["--days"] not in ("all", ""):
        params["days"] = [int(x) for x in opt["--days"].split(",")]
    if "--nproc" in opt:
        params["nproc"] = int(opt["--nproc"])
    if "--light" in opt:
        params["light"] = opt["--light"]
    if "--cache" in opt:
        use_cache_dir(opt["--cache"])
    slider = json.loads(opt.get("--slider", "1.0"))
    t0 = time.time()
    status = {"state": "running", "started": time.strftime("%Y-%m-%d %H:%M:%S"), "pid": os.getpid()}
    save_plans(status, out + ".status")
    try:
        bld = scene.get("building") if isinstance(scene, dict) else None
        if "--building" in opt:
            bld = json.load(open(opt["--building"]))
        res = plan(scene, opt.get("--epw"), params, slider, building=bld)
        save_plans(res, out)
        status.update({"state": "done", "sec": round(time.time() - t0, 1), "windows": len(res)})
    except Exception as ex:  # noqa
        import traceback
        status.update({"state": "error", "error": repr(ex)[:400], "tb": traceback.format_exc()[-2000:]})
    save_plans(status, out + ".status")
    print(json.dumps(dict((k, v) for k, v in status.items() if k != "tb")))
    return 0 if status["state"] == "done" else 1


if __name__ == "__main__":
    sys.exit(_cli(sys.argv[1:]))
