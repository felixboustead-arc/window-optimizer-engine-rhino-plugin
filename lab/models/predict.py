# -*- coding: utf-8 -*-
"""
predict.py - the per-minute prediction sidecar of ANY plan file (contract v1), with the EXACT physics of
physics.py at the exact sun position of every minute (no 1-degree bins), so every model - the old C1 / v1 too - is
measured the same way.  One row per minute in which a beam is wanted (beam_on), columns = Physics.PRED_COLS.

    python predict.py <plan.json> [<plan.json> ...]        -> writes <plan>_pred.csv next to each plan
    python predict.py --baseline 005                         -> exports C1 (current default) and v1 as plans first

Uses 5 worker processes.
"""
from __future__ import division

import csv
import io
import json
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

C = L.C
NPROC = 5
_PH = None


def _init(ctx_bytes):
    global _PH
    _PH = PH.Physics(pickle.loads(ctx_bytes))


def _rows(args):
    items, target = args
    return [_PH.pred_row(i, pose, target) for i, pose in items]


def predict_poses(ctx, items, target_pt=None, nproc=NPROC):
    """items = [(step, (ring, aperture)), ...] -> rows (same order)."""
    chunks = [items[k:k + 2000] for k in range(0, len(items), 2000)]
    blob = pickle.dumps(ctx)
    with mp.Pool(nproc, _init, (blob,)) as pool:
        parts = pool.map(_rows, [(c, target_pt) for c in chunks])
    return [r for p in parts for r in p]


def write_pred(path, rows):
    with io.open(path, "w", newline="", encoding="utf-8") as f:
        wr = csv.writer(f)
        wr.writerow(PH.Physics.PRED_COLS)
        wr.writerows(rows)
    return path


def predict_plan(plan_path, ctx=None, ev=None):
    plan = L.read_plan(plan_path)
    wtag = os.path.basename(plan_path).rsplit("_w", 1)[1][:3]
    if ctx is None:
        ev, _P, _job, _s = L.load(wtag)
        ctx = PH.Ctx(ev)
    poses = L.poses_from_runs(plan["runs"], ctx.n)
    items = [(i, poses[i]) for i in range(ctx.n) if ctx.beam_on[i]]
    tp = None
    if plan.get("target"):
        tp = tuple(plan["target"]["point_xyz"])
    t0 = time.time()
    rows = predict_poses(ctx, items, tp)
    out = plan_path[:-5] + "_pred.csv"
    write_pred(out, rows)
    plan["pred"] = os.path.basename(out)
    with io.open(plan_path, "w", encoding="utf-8") as f:
        f.write(json.dumps(plan, separators=(",", ":")))
    print("%s: %d rows in %.0f s -> %s" % (os.path.basename(plan_path), len(rows), time.time() - t0, out))
    return out


def export_baseline(wtag):
    """C1 (the current default 'consistency' planner) and v1 as plan files, same contract."""
    ev, P, job, sched = L.load(wtag, select=True)
    paths = []
    runs = C.pose_runs(ev, sched)
    paths.append(L.write_plan("C1", wtag, runs, extra={"note": "current default planner (consistency, C1 tracking), "
                                                                "scheduled runs (before break-up)"}))
    jobp = json.load(open(L.job_path(wtag)))["params"]
    scene = json.load(open(jobp["scene_path"]))
    wx = C.read_epw(jobp["epw_path"])
    ov = dict(jobp["overrides"], optimizer="v1", beam_score_model="v1")
    proj = C.run_project(scene, wx, ov, cache_dir=jobp["cache_dir"])
    paths.append(L.write_plan("v1", wtag, C.pose_runs(proj.ev, proj.scheds[0]),
                              extra={"note": "old optimiser v1 + 25-min median smoothing, scheduled runs"}))
    ctx = PH.Ctx(ev)
    for p in paths:
        predict_plan(p, ctx)
    return paths


if __name__ == "__main__":
    if sys.argv[1] == "--baseline":
        for w in sys.argv[2:]:
            export_baseline(L.wtag_norm(w))
    else:
        for p in sys.argv[1:]:
            predict_plan(p)
