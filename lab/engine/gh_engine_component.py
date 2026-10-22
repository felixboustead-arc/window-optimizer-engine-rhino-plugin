#! python 3
# -*- coding: utf-8 -*-
"""
SHARED ENGINE - Grasshopper "Python 3 Script" component (MODELER-D, D3-001, 2026-10-20; D3-001b / D3-002 2026-10-21).  Rhino 8, CPython 3.9.
Paste this file as the component source (or let lab\\test_scenes\\gh_r3_engine\\gh_engine_runner.py add it on a scene
COPY).  It plans 1-15 windows with lab\\engine\\engine.py (champion parts only) and gives the pose of EVERY window
for the clock minute, each solve.  Plans are kept in memory (scriptcontext.sticky) so a clock change costs only the
per-minute lookup (bisect; microseconds), not a re-plan.

INPUTS (all optional; add them on the component with these exact names, type hint in brackets):
  run         [bool]    True = (re)plan when the inputs changed.  False = keep the last plans, only look up poses.
  windows     [str/Plane, list]  one window per item: JSON spec {"id","facing","tilt_deg","level_m","scale",
                                 "width_f","beam_w","area_m2"} or a Plane (normal = outward window normal).
                                 Default: one north window (Felix's).
  epw         [str]     weather file (default Williamtown TMYx = the TESTER's).
  days        [str]     "dayahead" (default: plan days [d-1, d, d+1] around the clock day d; light day plans are
                        independent per day and the night glides only link planned neighbours, so day d equals the
                        full-year plan - checked by engine_multi day_ahead_check, N/E/S, tilts, winter/summer, wrap),
                        "all" (full year: use background=True), or "171,354".
  slider      [float]   1 = light (default), 0 = comfort (B's nested knapsack).  The knapsack judges the whole planned
                        horizon, so use days="all" for COMBINED (day-ahead COMBINED is only an approximation).
  month, day, hour, minute  [int]  AEST clock (no DST).  Default 1 May 10:00 (Felix's time-lapse clock).
  building    [str]     O13 building JSON (or a path to one; see lab\engine\o13.py).  Thermal / COMBINED only.
  background  [bool]    (default True) plan in a separate CPython process (miniconda; full-year plans, nproc workers) and keep
                        solving with the previous plans until the new ones are written; the component re-solves
                        itself every 5 s while waiting; a dead process is relaunched (max 3 tries) and errors are
                        shown in status.  background=False plans on the GH thread (freezes the UI for 0.5-2+ min).
  nproc       [int]     worker processes for background planning (default 4; lab rule: TESTER decides).
OUTPUTS: ring (deg, list), aperture (% close, list), ids (list), status (str), a (JSON summary).
Wire ring / aperture into Felix's chain (Ring Angle -> Rad, Aperture Close -> the two blade expressions) on a COPY.
"""
import hashlib
import json
import os
import subprocess
import sys
import time

ENGINE_DIR = os.environ.get("WO_ENGINE_DIR", r"C:\path\to\window-optimizer-engine\lab\engine")  # EDIT
BG_PY = os.environ.get("WO_PYTHON", "python")  # a CPython 3.9+ for background full-year plans
BG_DIR = os.environ.get("WO_BG_DIR", os.path.join(ENGINE_DIR, "..", "out", "engine_r3", "gh_bg"))
TEST = r"__D3_TEST__"            # the TESTER runner replaces this with a JSON config path (test mode)

try:
    import scriptcontext as _sc
    STICKY = _sc.sticky
except Exception:                # outside Rhino (offline check)
    STICKY = {}


def _g(name, default=None):
    v = globals().get(name)
    if v is None or (isinstance(v, (list, tuple)) and len(v) == 0):
        return default
    return v


def _engine():
    if ENGINE_DIR not in sys.path:
        sys.path.insert(0, ENGINE_DIR)
    import engine as EN
    return EN


def _win_specs(items):
    out = []
    for k, it in enumerate(items or [{"id": "N", "facing": 0.0, "beam_w": 0.10}]):
        if isinstance(it, dict):
            w = dict(it)
        elif isinstance(it, str):
            w = json.loads(it)
        elif hasattr(it, "ZAxis"):                          # Rhino.Geometry.Plane
            z = it.ZAxis
            o = it.Origin
            w = {"plane": {"origin": [o.X, o.Y, o.Z], "normal": [z.X, z.Y, z.Z]}}
        else:
            raise ValueError("window %d: give a JSON spec or a Plane, got %r" % (k, type(it)))
        w.setdefault("id", "w%02d" % (k + 1))
        w.setdefault("beam_w", 0.10)
        out.append(w)
    return out


def _clock(cfg):
    if cfg and os.path.exists(cfg.get("clock_file", "")):
        c = json.load(open(cfg["clock_file"]))
        return int(c["month"]), int(c["day"]), int(c["hour"]), int(c["minute"])
    return (int(_g("month", 5)), int(_g("day", 1)), int(_g("hour", 10)), int(_g("minute", 0)))


def _days(mode, doy, EN=None):
    if mode in (None, "", "dayahead"):                      # [d-1, d, d+1]: day d == the full-year plan (M5)
        return EN.dayahead_days(doy) if EN is not None else sorted(set([(doy - 1) % 365, doy, (doy + 1) % 365]))
    if mode == "all":
        return None
    return [int(x) for x in str(mode).replace(" ", "").split(",") if x != ""]


def _bg_launch(st, tag, scene, params, out, slider, days, epw, building):
    sp = os.path.join(BG_DIR, "scene_%s.json" % tag)
    sc = dict(scene, params={k: v for k, v in params.items() if k != "days"})
    if building:
        sc["building"] = building
    json.dump(sc, open(sp, "w"))
    for f in (out + ".status",):                           # a stale status from an older process is not the truth
        if os.path.exists(f):
            os.remove(f)
    cmd = [BG_PY, "-B", os.path.join(ENGINE_DIR, "engine.py"), "plan", sp, out,
           "--slider", json.dumps(slider), "--nproc", str(int(_g("nproc", 4)))]
    if days is not None:
        cmd += ["--days", ",".join(str(d) for d in days)]
    if epw:
        cmd += ["--epw", epw]
    st["_proc_" + tag] = subprocess.Popen(cmd, creationflags=getattr(subprocess, "CREATE_NO_WINDOW", 0))
    st["_bg_" + tag] = time.time()
    st["_tries_" + tag] = st.get("_tries_" + tag, 0) + 1


def _resolve_later():
    try:
        ghenv.Component.OnPingDocument().ScheduleSolution(  # noqa: F821  (GH global)
            5000, lambda d: ghenv.Component.ExpireSolution(False))  # noqa: F821
    except Exception:
        pass


def main():
    cfg = json.load(open(TEST)) if (TEST != "__D3_" + "TEST__" and os.path.exists(TEST)) else None
    EN = _engine()
    wins = _win_specs(cfg["windows"] if cfg else _g("windows"))
    epw = (cfg or {}).get("epw") or _g("epw")
    building = (cfg or {}).get("building") or _g("building")
    if isinstance(building, str):
        building = json.loads(building) if building.strip().startswith("{") else json.load(open(building))
    mo, da, ho, mi = _clock(cfg)
    minute = EN.minute_of_year(mo, da, ho, mi)              # invalid dates raise a clear ValueError (M8)
    doy = minute // 1440
    dmode = (cfg or {}).get("days") or _g("days", "dayahead")
    days = _days(dmode, doy, EN)
    slider = float((cfg or {}).get("slider", _g("slider", 1.0)))
    run = bool((cfg or {}).get("run", _g("run", True)))
    # M8: planning runs in the background by default (never freezes the GH UI); test mode keeps its own setting
    bg = bool(cfg.get("background", False)) if cfg else bool(_g("background", True))
    params = {"nproc": 1}
    if days is not None:
        params["days"] = days
    scene = {"windows": wins}
    epw_m = int(os.path.getmtime(epw)) if (epw and os.path.exists(epw)) else None
    # M6: the key holds the engine version + code mtimes + EPW mtime, so plans of an older engine are never reused
    key = json.dumps([wins, epw, epw_m, days, slider, building, EN.code_signature()], sort_keys=True, default=str)
    st = STICKY.setdefault("engine_d3", {})
    info = {"minute_of_year": minute, "clock": [mo, da, ho, mi], "days": days if days is None or len(days) < 8 else
            "%d days" % len(days), "slider": slider, "engine": EN.VERSION, "code": EN.code_signature()[1:]}
    warn = []
    if slider < 1.0 and days is not None:
        warn.append("COMBINED on a day subset is an APPROXIMATION (the knapsack judges only the planned days); "
                    "use days=all for the real comfort plan")
    if slider < 1.0:
        info["comfort_room"] = ("O13 building" if building else "A04_felix_open (default thermal_case)")
    if cfg and cfg.get("cold_cache_dir") and not st.get("_cold_done"):
        EN.use_cache_dir(cfg["cold_cache_dir"])             # test mode: every engine disk cache starts empty
        st["_cold_done"] = True
    plans = st.get(key)
    if plans is None and run and not bg:
        t0 = time.time()
        plans = EN.plan(scene, epw, dict(params), slider, building=building)
        info["plan_s"] = round(time.time() - t0, 1)
        st[key] = plans
        st["_last"] = key
    elif plans is None and run and bg:
        if not os.path.isdir(BG_DIR):
            os.makedirs(BG_DIR)
        tag = hashlib.sha1(key.encode()).hexdigest()[:10]
        out = os.path.join(BG_DIR, "plans_%s.json" % tag)
        stf = out + ".status"
        stat = None
        if os.path.exists(stf):
            try:
                stat = json.load(open(stf))
            except Exception:
                stat = None
        if os.path.exists(out) and (stat is None or stat.get("state") == "done"):
            plans = json.load(open(out))
            st[key] = plans
            st["_last"] = key
            info["plan"] = "background plans loaded"
        else:
            proc = st.get("_proc_" + tag)
            dead = proc is not None and proc.poll() is not None
            if stat and stat.get("state") == "error":       # M7: surface the error, no endless 5 s loop
                warn.append("BACKGROUND PLAN FAILED: %s" % stat.get("error", "?"))
                info["plan"] = stat
                if st.get("_tries_" + tag, 0) < 2 and dead:
                    _bg_launch(st, tag, scene, params, out, slider, days, epw, building)   # one automatic retry
                    _resolve_later()
            elif proc is None or (dead and not os.path.exists(out)):
                if proc is not None:                         # M7: process died without a result -> relaunch
                    warn.append("background process exited (code %s) without plans: relaunched" % proc.returncode)
                if st.get("_tries_" + tag, 0) < 3:
                    _bg_launch(st, tag, scene, params, out, slider, days, epw, building)
                    _resolve_later()
                else:
                    warn.append("background planning gave up after 3 tries (see %s)" % stf)
            else:
                info["plan"] = stat or "starting"
                _resolve_later()
            if st.get("_bg_" + tag):
                info["waiting_s"] = round(time.time() - st["_bg_" + tag], 0)
            plans = st.get(st.get("_last")) if st.get("_last") else None
            if plans:
                warn.append("showing the PREVIOUS plans while the new ones are computed")
    elif plans is None:
        plans = st.get(st.get("_last")) if st.get("_last") else None
        info["plan"] = "run is off: showing the last plans" if plans else "run is off and nothing planned yet"
    if warn:
        info["warning"] = warn
    if not plans:
        return [], [], [w["id"] for w in wins], "no plans yet" + ("; " + "; ".join(warn) if warn else ""), \
            json.dumps(info)
    want = [w["id"] for w in wins]
    ids = [i for i in want if i in plans]
    if len(ids) != len(want):                                # L4: an id mismatch is an error, not a silent fallback
        missing = [i for i in want if i not in plans]
        warn.append("windows without a plan: %s (the plans are for %s)" % (missing, sorted(plans)[:15]))
        info["warning"] = warn
        if not ids:
            return [], [], want, "ERROR: " + warn[-1], json.dumps(info)
    t1 = time.perf_counter()
    S = EN.Schedule(plans, ids)
    poses = S.poses(minute)
    info["lookup_ms"] = round(1e3 * (time.perf_counter() - t1), 3)
    ring = [round(p[0], 2) for p in poses]
    ap = [round(p[1], 2) for p in poses]
    unplanned = [i for i in ids if not S.planned(i, minute)]   # L2
    state = "planned" if not unplanned else "NOT PLANNED for this day: %s (held pose)" % unplanned
    status = "%d windows, %s, minute %d: %s" % (len(ids), state, minute,
                                                 "; ".join("%s r%.2f a%.2f" % (i, r, a) for i, r, a in
                                                           zip(ids, ring, ap))[:400])
    if warn:
        status += " | " + " | ".join(warn)
    if cfg:                                                  # test mode: write what the runner / checker compare
        rec = {"info": info, "ids": ids, "ring": ring, "aperture": ap,
               "runs_hash": dict((i, hashlib.sha1(json.dumps(plans[i].get("runs"), separators=(",", ":"),
                                                             default=str).encode()).hexdigest()[:12]) for i in ids)}
        with open(cfg["result_file"], "a") as f:
            f.write(json.dumps(rec) + "\n")
        if cfg.get("plans_file") and not os.path.exists(cfg["plans_file"]):
            json.dump(plans, open(cfg["plans_file"], "w"), separators=(",", ":"))
    return ring, ap, ids, status, json.dumps(info)


try:
    ring, aperture, ids, status, a = main()
except Exception:
    import traceback
    ring, aperture, ids = [], [], []
    status = "ERROR"
    a = traceback.format_exc()[-2000:]
    if TEST != "__D3_" + "TEST__":
        try:
            with open(json.load(open(TEST))["result_file"], "a") as f:
                f.write(json.dumps({"error": a}) + "\n")
        except Exception:
            pass
