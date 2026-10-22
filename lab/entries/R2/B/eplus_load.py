# -*- coding: utf-8 -*-
"""
eplus_load.py - read one EnergyPlus reference case (E:\\WindowOptimizer_data\\eplus_ref\\<case>) into the inputs BOTH
thermal models take (model A thermal_rc and model B thermal_hb share the interface), plus the E+ answers.

    import eplus_load as EL
    R = EL.load("A04_felix_open")
    tp  = M.make_tparams(R["over"]);  bld = M.assemble(R["env"], tp)
    win = M.window_gains(bld, R["beam_cls"], R["diffuse"])
    out = M.run_year(bld, R["wx"], R["climate"], win=win)          # compare out["op"] with R["ref"]["top"]

Drive = EnergyPlus's OWN boundary data, step by step (10 min, year 2), so the comparison tests the thermal model, not
the solar geometry:  outdoor air, sky temperature (wx['t_sky']), and E+'s incident solar on every outside surface,
pre-filled into the plane-irradiance cache keys both models read: wx[('_plane', az, tilt, 0.2)].
Window: E+ transmitted beam split by where it lands (E+ 'exterior windows incident beam' on the floor / walls /
ceiling), E+ transmitted diffuse + the inward-flowing share of the solar absorbed in the glass (f_in from the film
and glass resistances, 0.30 single .. 0.46 triple; see load()), all divided by g (= SHGC, r_dirt = 1 in 'over')
because the models multiply by g themselves.  R['glass_abs'] = the E+ glass-absorbed solar (W) for models that can
take it on a glass node (B's window_gains(..., glass_abs=...)); NOT used in the common, equal-input comparison.
"""
from __future__ import division

import array
import gzip
import json
import os

ROOT = os.environ.get("WO_EPLUS_REF", os.path.join(os.path.abspath(os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "..", "..", "..")), "data", "eplus_ref"))  # E+ reference cases


def cases(root=ROOT):
    return sorted(d for d in os.listdir(root) if os.path.exists(os.path.join(root, d, "ref.csv.gz")))


def load(case, root=ROOT):
    d = os.path.join(root, case)
    with open(os.path.join(d, "case.json")) as f:
        meta = json.load(f)
    with gzip.open(os.path.join(d, "ref.csv.gz"), "rt") as f:
        head = f.readline().strip().split(",")
        cols = dict((h, array.array("d")) for h in head[1:])
        stamps = []
        order = [cols[h] for h in head[1:]]
        for line in f:
            p = line.rstrip().split(",")
            stamps.append(p[0])
            for a, v in zip(order, p[1:]):
                a.append(float(v))
    n = len(stamps)
    sph = int(meta.get("steps_per_h", 6))
    month = array.array("i", (int(s.strip()[0:2]) for s in stamps))
    import math
    alt = cols["alt"]
    ghi = array.array("d", (cols["dhi"][i] + cols["dni"][i] * max(0.0, math.sin(math.radians(alt[i])))
                            for i in range(n)))
    wx = {"temp": cols["te"], "dni": cols["dni"], "dhi": cols["dhi"], "ghi": ghi, "month": month,
          "t_mid": array.array("d", ((i + 0.5) / sph for i in range(n))), "alt": alt, "az": cols["azi"],
          "step_s": 3600.0 / sph, "n": n, "t_sky": cols["tsky"], "wind": cols["wind"]}
    rho = round(float(meta["over"].get("ground_reflectance", 0.2)), 3)
    for surf, az, tilt in (("WALL_N", 0.0, 90.0), ("WALL_E", 90.0, 90.0), ("WALL_S", 180.0, 90.0),
                           ("WALL_W", 270.0, 90.0), ("ROOF", 0.0, 0.0)):
        wx[("_plane", round(az, 3), round(tilt, 3), rho)] = cols["inc_" + surf]
    g = float(meta["over"]["shgc"]) * float(meta["over"].get("r_dirt", 1.0))
    U = float(meta["over"]["u_glazing"])
    # inward-flowing share of the solar absorbed in the glass (single equivalent E+ layer, absorbed mid-glass):
    # outside path = film (h ~ 15 W/m2K) + half the glass resistance; inside path = film (TARP ~1.5 + radiation 5)
    r_g = max(0.0, 1.0 / U - 0.17)
    f_in = (0.067 + 0.5 * r_g) / (0.067 + r_g + 0.154)
    trb, trd, absg = cols["trb_WIN"], cols["trd_WIN"], cols["absg_WIN"]
    zero = array.array("d", [0.0]) * n
    grp = {"floor": [c for c in cols if c.startswith("bmi_FLOOR")],
           "wall": [c for c in cols if c.startswith("bmi_WALL")],
           "ceiling": [c for c in cols if c == "bmi_ROOF"]}
    beam_cls = dict((k, array.array("d", zero)) for k in grp)
    for i in range(n):
        if trb[i] <= 0:
            continue
        parts = dict((k, sum(cols[c][i] for c in v)) for k, v in grp.items())
        tot = sum(parts.values())
        if tot <= 0:
            beam_cls["floor"][i] = trb[i] / g
            continue
        for k in grp:
            beam_cls[k][i] = trb[i] * parts[k] / tot / g
    diffuse = array.array("d", ((trd[i] + f_in * absg[i]) / g for i in range(n)))
    ref = {"tair": cols["tair"], "top": cols["top"], "tmrt": cols["tmrt"], "stamps": stamps}
    for c in cols:
        if c.startswith("ti_"):
            ref[c] = cols[c]
    return {"case": case, "meta": meta, "env": meta["env"], "over": dict(meta["over"]), "climate": meta["climate"],
            "wx": wx, "beam_cls": beam_cls, "diffuse": diffuse, "ref": ref, "n": n, "cols": cols, "f_in": f_in,
            "glass_abs": absg}
