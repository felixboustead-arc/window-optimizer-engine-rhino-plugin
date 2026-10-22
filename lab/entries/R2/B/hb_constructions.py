# -*- coding: utf-8 -*-
"""
hb_constructions.py - MODELER-B.  Layer stacks for the heat-balance model (thermal_hb.py) AND the EnergyPlus reference
runs (eplus_gen.py), built from lab\\research\\constructions.json by ONE rule, so the two always see the same fabric.

Rule (documented for the TESTER):
  * layers are OUTSIDE -> INSIDE as in the library; 'R:x' = a massless resistance (cavity, batt).
  * F02 is listed slab-first in the library although its text and kappa_int (196 = the 85 mm slab) say the EPS is UNDER
    the slab: it is read as [eps, concrete] (outside -> inside).
  * U is matched to the library U (or the u_* override) by changing the LARGEST massless layer's R (framing, ties,
    roof space); if the stack has no massless layer and needs more R, a massless layer is added on the outside.
    Films used for U: Rsi wall 0.13, roof (heat up) 0.10, floor (heat down) 0.17, Rse 0.04 (ground: 0).
  * Ground floors get 0.5 m of soil (lambda 2.0, rho 2000, c 1000) under the slab and a massless deep-ground
    resistance so the U equals the ISO 13370 value; the far side is the earth temperature.
  * A covered (carpeted) floor part gets a massless R 0.3 on top, taken out of the deep-ground R so U stays the same.
Standard library only.
"""
from __future__ import division

import json
import os

LIB_PATH = os.path.abspath(os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "..", "..", "research",
                                        "constructions.json"))
RSI = {"wall": 0.13, "roof": 0.10, "floor": 0.17, "partition": 0.13}
RSE = 0.04
SOIL = {"lambda": 2.0, "rho": 2000.0, "c": 1000.0}
CARPET_R = 0.3
_LIB = None


def load_library(path=LIB_PATH):
    global _LIB
    if _LIB is None or path != LIB_PATH:
        with open(path) as f:
            raw = json.load(f)
        lib = {"materials": raw["materials"], "elements": dict((e["id"], e) for e in raw["elements"]),
               "archetypes": dict((a["id"], a) for a in raw["archetypes"]), "raw": raw}
        if path != LIB_PATH:
            return lib
        _LIB = lib
    return _LIB


def raw_layers(el):
    lay = [list(x) for x in el.get("layers", [])]
    if el["id"] == "F02" and lay and lay[0][0] == "concrete":
        lay = list(reversed(lay))                       # EPS under the slab (see module doc)
    return lay


def _mass_layer(lib, mid, d):
    m = lib["materials"][mid]
    return {"name": mid, "d": float(d), "k": float(m["lambda"]), "rho": float(m["rho"]), "c": float(m["c"]),
            "R": float(d) / float(m["lambda"])}


def _nomass(R, name="R"):
    return {"name": name, "d": 0.0, "k": None, "rho": 0.0, "c": 0.0, "R": max(float(R), 0.001)}


def stack(lib, el, kind, U_target=None, ground=False, covered=False):
    """List of layer dicts (outside -> inside) whose U (with films) equals U_target (library U by default)."""
    U = float(U_target if U_target is not None else el["U"])
    lay = []
    for mid, d in raw_layers(el):
        if mid.startswith("R:"):
            lay.append(_nomass(float(mid[2:]), "cavity_or_batt"))
        else:
            lay.append(_mass_layer(lib, mid, d))
    rsi, rse = RSI[kind], (0.0 if ground else RSE)
    R_target = 1.0 / U
    if ground:
        soil = {"name": "soil", "d": 0.5, "k": SOIL["lambda"], "rho": SOIL["rho"], "c": SOIL["c"],
                "R": 0.5 / SOIL["lambda"]}
        r_now = rsi + sum(l["R"] for l in lay)
        spare = R_target - r_now - (CARPET_R if covered else 0.0)
        if spare < soil["R"] + 0.01:                        # very high U: thinner soil, no deep resistance
            d = max(0.05, (spare - 0.01) * SOIL["lambda"])
            soil["d"], soil["R"] = d, d / SOIL["lambda"]
        lay = [soil] + lay
        if covered:
            lay.append(_nomass(CARPET_R, "carpet"))
        r_now = rsi + sum(l["R"] for l in lay)
        if R_target - r_now > 0.005:
            lay = [_nomass(R_target - r_now, "deep_ground")] + lay
        return lay
    if covered:
        lay.append(_nomass(CARPET_R, "carpet"))
    r_now = rsi + rse + sum(l["R"] for l in lay)
    diff = R_target - r_now
    if abs(diff) > 0.005:
        nm = [l for l in lay if l["k"] is None and l["name"] != "carpet"]
        if nm:
            big = max(nm, key=lambda l: l["R"])
            big["R"] = max(0.01, big["R"] + diff)
        elif diff > 0:
            lay = [_nomass(diff, "adjust")] + lay
    return lay


def u_of(lay, kind, ground=False):
    return 1.0 / (RSI[kind] + (0.0 if ground else RSE) + sum(l["R"] for l in lay))


def kappa_inside(lay, depth=0.1):
    """Room-side areal heat capacity (kJ/m2K) to the first massless layer or `depth` m (for reports only)."""
    k, d = 0.0, 0.0
    for l in reversed(lay):
        if l["k"] is None:
            if l["R"] > 0.1:
                break
            continue
        take = min(l["d"], max(0.0, depth - d))
        k += l["rho"] * l["c"] * take / 1000.0
        d += l["d"]
        if d >= depth:
            break
    return k
