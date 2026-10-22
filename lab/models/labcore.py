# -*- coding: utf-8 -*-
"""
labcore.py - shared helpers for the lab planners: load a lab job (scan from the cache -> Evaluator), plan files
(contract in lab\\README.md + board), day / period helpers.
"""
from __future__ import division

import io
import json
import os
import sys
import time

CORE_DIR = os.path.abspath(os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", ".."))  # project root
LAB = os.path.join(CORE_DIR, "lab")
if CORE_DIR not in sys.path:
    sys.path.insert(0, CORE_DIR)
import window_beams_core as C  # noqa: E402

WIDTHS = {"005": 0.05, "010": 0.10, "020": 0.20, "030": 0.30}
CONTRACT = 1


def wtag_norm(w):
    """'5', '05', '005', 'w005', '0.05' -> '005' (beam width in cm, 3 digits - the file / job tag)."""
    s = str(w).strip().lstrip("wW")
    cm = int(round(float(s) * 100)) if "." in s else int(s)
    tag = "%03d" % cm
    if tag not in WIDTHS:
        raise ValueError("unknown beam width tag %r (have %s)" % (w, sorted(WIDTHS)))
    return tag


def job_path(wtag):
    return os.path.join(LAB, "jobs", "w" + wtag_norm(wtag), "job.json")


def load(wtag, select=False, overrides=None, verbose=True):
    """(ev, P, job, sched or None).  select=True also runs the current default planner (C1) for the schedule."""
    t0 = time.time()
    wtag = wtag_norm(wtag)
    job = json.load(open(job_path(wtag)))
    P0 = job["params"]
    scene = json.load(open(P0["scene_path"]))
    wx = C.read_epw(P0["epw_path"])
    ov = dict(P0["overrides"])
    if overrides:
        ov.update(overrides)
    P = C.make_params(ov)
    w = C.prepare_weather(wx, P["time_step_min"])
    sc = C.resolve_scene(scene, P)
    win = C.split_windows(sc)[0]
    scan = C.run_scan(win, w, P, None, None, None, P0["cache_dir"])
    ev = C.Evaluator(win, w, scan, P)
    sched = ev.select() if select else None
    if verbose:
        print("loaded w%s in %.1f s (scan from cache: %s)" % (wtag, time.time() - t0, scan.from_cache))
    return ev, P, job, sched


def write_plan(model, wtag, runs, target=None, pred=None, extra=None, step_min=1, n_steps=525600):
    wtag = wtag_norm(wtag)
    path = os.path.join(LAB, "plans", "%s_w%s.json" % (model, wtag))
    data = {"contract": CONTRACT, "model": model, "width_m": WIDTHS[wtag], "job": job_path(wtag),
            "step_min": step_min, "n_steps": n_steps, "params": ["ring", "aperture"], "ranges": [15.0, 100.0],
            "runs": runs, "target": target, "pred": pred}
    if extra:
        data.update(extra)
    with io.open(path, "w", encoding="utf-8") as f:
        f.write(json.dumps(data, separators=(",", ":")))
    return path


def read_plan(path):
    with io.open(path, "r", encoding="utf-8") as f:
        return json.load(f)


def runs_from_poses(poses, n, step0=0, nd=4):
    """poses = per-step (ring, aperture) list (length n) -> runs [[step, ring, aperture], ...] (one per change)."""
    runs, prev = [], None
    for i in range(n):
        p = tuple(round(float(v), nd) for v in poses[i])
        if p != prev:
            runs.append([i + step0] + list(p))
            prev = p
    return runs


def poses_from_runs(runs, n):
    """runs -> per-step pose list (the run in force at the START of each step; fractional run starts count from the
    step they fall in, i.e. the pose held for most of a minute is the later one when a run starts inside it)."""
    out = [None] * n
    k = 0
    cur = tuple(runs[0][1:])
    for i in range(n):
        while k < len(runs) and runs[k][0] <= i + 0.5:
            cur = tuple(runs[k][1:])
            k += 1
        out[i] = cur
    return out


def sched_poses(ev, sched):
    return [tuple(float(v) for v in ev.scan.poses[p]) for p in sched]
