# -*- coding: utf-8 -*-
"""
lab\\engine\\o13.py - O13 "building designs, doors and ventilation" layer of the shared R3 engine (MODELER-D, D3-002,
2026-10-21).  R3_BRIEF sec. 1; OBJECTIVES O13; ORCH-040.

SWITCH: engine.plan(..., building=None) or building={"o13": False} -> O13 OFF: the R2 champion path, untouched (the
frozen HB-L2 thermal_hb2.run_year, bit-identical).  O13 acts in mode COMBINED only (engine.plan switch; the
THERMAL-only path is not hooked); the light planner is unchanged.

VALIDATION STATUS (D3-002c, AU-R3-002 R4 labels):
  doors        E+ ZoneRefrigerationDoorMixing (B3): mechanism agrees, but that reference is Gosney-Olama (~Cd 0.73
               equivalent) vs Brown-Solvason Cd 0.62 here: an ~18 % flow spread that the 0.06-0.16 K residual cannot
               discriminate.  Door Cd spread 0.62-0.73 OPEN.
  vents        B1_vent = "vs the ASHRAE form in E+" (E+ ZoneVentilation:WindandStackOpenArea uses the SAME correlation):
               checks implementation + control timing, not physics.  The modulating vent (room held at
               max(tin_min, Te + dt_min)) is a modelling assumption (time average of E+'s on/off cycling).
  cross-vent   "IMPLEMENTATION-CHECKED vs E+ AFN, VALIDATED WITH LIMITS" (J-R3-007); accuracy vs a real house UNTESTED.
               Limits: same Cp in E+ and here; E+ opens at Tz > Tout with no 1 K margin; cracks constant 0.6 ach (not
               wind-driven); cross daily error 0.05-0.10 K (gap may be larger at 0.74 ach vent-hour infiltration, run
               not done, AU-R3-003 note 4); single-sided (B4_single) over-cools vs AFN by ~0.1 K (recorded, not tuned);
               door Cd 0.62 vs ~0.73 (see doors); vent modulation is a modelling assumption (see vents).
  fan          UNVALIDATED (no reference; comfort-band change only, applied where the upper limit is >= 25 C).
               (AU-R3-003 R3) The credit (ASHRAE 55 elevated-air-speed band offset) assumes occupants have control of
               the fan / speed; no occupant-control cap is modelled.  Hard step at 25 C: an upper limit of 24.9 C gets
               no credit, 25.0 C gets the full offset.
  shared input (AU-R3-003 note 1)  o13_score.py zone_wins (~:50-58) feeds E+'s own trd / absg columns to base, O13 and
               B4 alike: symmetric, so base-vs-O13 and mechanism comparisons are fair, but the absolute RMSE is
               same-drive (not independent) and glass absorptance is untested.
  drive choice (AU-R3-003 note 2)  not blind: the ORCH-029 drive was chosen after seeing the Williamtown O13 rows vs
               E+, and Nobbys had been scored under both drives in T-R3-009, so the held-out Nobbys 7/7 (T-R3-013) is a
               reproduction on a once-seen case, not a fresh blind test.

WHAT O13 ADDS (one node per zone; zones linked only by partition conduction and door / leak airflow, never merged):
  zones        each zone = the frozen HB-L2 core built by thermal_hb2.assemble(env, tp) (its own layers, solar,
               internal gains, base infiltration n_rest, MVHR n_system / heat_recovery, night purge n_purge ...)
  partitions   interzone walls: a layered element in each zone whose far side sees the OTHER zone's air (lagged one
               10-min step; conduction through a wall is slow)
  doors        bidirectional buoyancy flow (Brown and Solvason 1962; the form of E+ AFN DetailedOpening and Modelica
               Buildings.Airflow.Multizone.DoorDiscretizedOpen):  Q = Cd/3 * A_eff * sqrt(g H dT / T_mean),
               A_eff = W H f_open(t) + leak_m2 (closed-door leakage), Cd 0.62 (brief 0.60-0.65)
  vents        operable windows / vents, single-sided wind + stack exactly in the E+ ZoneVentilation:WindandStackOpenArea
               form (ASHRAE HoF 2009 ch.16 eq. 37-39): Qw = Cw A f V_loc, Qs = Cd A f sqrt(2 g dH_npl |Tz-To| / Tz),
               Q = sqrt(Qw^2 + Qs^2); Cw / Cd 'auto' = E+ autocalculate; V_loc = E+ wind profile at the zone centroid;
               schedule (hourly 0..1, months) + E+ temperature controls (min / max indoor, delta T, min / max outdoor,
               max wind)
  cross-vent   two open vents on facades >= 90 deg apart, in one zone or in two zones joined by an open door:
               Q = Cd A_eff V_loc sqrt(|Cp_up - Cp_down|), A_eff = (sum 1/A_i^2)^-1/2 over the path, Cp = Swami and
               Chandra (1988) low-rise (as E+ AFN SurfaceAverageCalculation); outdoor air enters the upwind zone and
               passes to the downwind zone.  Implementation-checked vs E+ AFN (B4); real-house accuracy untested.
  MVHR         balanced supply ach with sensible effectiveness hr_eff, optional summer bypass (Tout < Tin, Tin > bypass_tin)
  fan          ceiling-fan air speed: comfort only (adaptive band upper limit + ASHRAE 55 elevated air speed offset,
               only where that limit is >= 25 C operative), no heat-balance change.  UNVALIDATED.  Assumes occupant control of
               the fan / speed (no cap modelled); hard step at 25 C (24.9 C: no credit, 25.0 C: full offset).
  class        "light" / "medium" / "heavy" -> archetype A01 / A04 / A15 (the E+ reference stacks) unless a zone sets it

All airflow terms enter each zone's air node implicitly (the zone temperatures are solved together each step, so a
large door flow cannot oscillate); flow COEFFICIENTS use the previous step's temperatures (as E+ ZoneVentilation does).

BUILDING JSON (all fields optional except zones):
  {"name": "...", "o13": true, "class": "medium", "terrain": "suburbs",
   "zones": [{"id": "cell", "env": <thermal env, eplus_ref case.json format>, "over": {<thermal_hb2 params>},
              "windows": true, "z_centroid_m": 1.35}, ...],
   "partitions": [{"a": "cell", "b": "back", "area_m2": 21.6, "type": "stud" | "masonry", "azimuth_a": 180}],
   "doors": [{"a": "cell", "b": "back", "w": 0.9, "h": 2.1, "sched": "always" | "never" | [24 x 0..1] |
              {"open": [7, 22]}, "leak_m2": 0.015, "cd": 0.62}],
   "vents": [{"zone": "cell", "area_m2": 1.0, "height_m": 1.0, "azimuth_deg": 0, "sched": ..., "months": [..],
              "tin_min": 23, "tin_max": 100, "dt_min": 1, "tout_min": -100, "tout_max": 100, "wind_max": 40,
              "cw": "auto", "cd": "auto", "npl_m": null, "cross": true}],
   "mv": [{"zone": "cell", "ach": 0.5, "hr_eff": 0.75, "bypass": false, "bypass_tin": 23}],
   "fan": {"air_speed": 0.6}}
"""
from __future__ import division

import array
import copy
import math
import os
import sys
import warnings

HERE = os.path.dirname(os.path.abspath(__file__))
LAB = os.path.dirname(HERE)
B_DIR = os.path.join(LAB, "entries", "R2", "B")
VERSION = "O13 D3-002c (multi-zone HB-L2: partitions, doors, vents, cross-vent, MVHR, fan>=25C unvalidated)"
G = 9.81
CLASS_ARCH = {"light": "A01", "medium": "A04", "heavy": "A15"}
TERRAIN = {"country": (0.14, 270.0), "suburbs": (0.22, 370.0), "city": (0.33, 460.0), "ocean": (0.10, 210.0)}
PART_STACKS = {
    "stud": [{"name": "plasterboard", "d": 0.010, "k": 0.25, "rho": 900.0, "c": 1000.0, "R": 0.04},
             {"name": "stud_cavity", "d": 0.0, "k": None, "rho": 0.0, "c": 0.0, "R": 0.17},
             {"name": "plasterboard", "d": 0.010, "k": 0.25, "rho": 900.0, "c": 1000.0, "R": 0.04}],
    "masonry": [{"name": "brick_clay", "d": 0.110, "k": 0.70, "rho": 1850.0, "c": 840.0, "R": 0.110 / 0.70}]}
R_FILM_FAR = 0.12                 # far-side surface film of a partition (m2K/W: convection + radiation ~ 8.3 W/m2K)

PARAMS = [
    ("door_cd", 0.62, 0.60, 0.65, "-", "door discharge coefficient (Brown-Solvason bidirectional flow)"),
    ("door_leak_m2", 0.015, 0.0, 0.05, "m2", "closed-door leakage area"),
    ("vent_cw", "auto", 0.30, 0.55, "-", "vent opening effectiveness (auto = E+ angle rule)"),
    ("vent_cd", "auto", 0.40, 0.70, "-", "vent discharge coefficient (auto = 0.40 + 0.0045 dT)"),
    ("vent_area_m2", 0.0, 0.0, 4.0, "m2", "operable vent / window opening area"),
    ("vent_height_m", 1.0, 0.3, 2.0, "m", "opening height (neutral plane = height / 2)"),
    ("cross_cd", 0.60, 0.50, 0.70, "-", "cross-ventilation path discharge coefficient"),
    ("mv_ach", 0.0, 0.0, 1.0, "1/h", "MVHR balanced air changes"),
    ("hr_eff", 0.0, 0.0, 0.85, "-", "MVHR sensible effectiveness"),
    ("air_speed", 0.0, 0.0, 1.2, "m/s", "ceiling fan air speed (comfort band only)"),
    ("class", "medium", ["light", "medium", "heavy"], None, "-", "construction class -> archetype A01 / A04 / A15"),
    ("terrain", "suburbs", ["country", "suburbs", "city", "ocean"], None, "-", "wind profile at the zones (E+ terrain)"),
]


def param_table():
    out = []
    for name, default, lo, hi, units, meaning in PARAMS:
        r = {"name": name, "default": default, "units": units, "group": "o13", "meaning": meaning}
        if isinstance(lo, list):
            r["choices"] = lo
        else:
            r["min"], r["max"] = lo, hi
        out.append(r)
    return out


def _th():
    """The FROZEN HB-L2 module the champion uses (entries\\R2\\B\\thermal_hb2.py; same sys.modules entry)."""
    if "thermal_hb2" in sys.modules:
        return sys.modules["thermal_hb2"]
    for p in (os.path.join(LAB, "models"), B_DIR):
        if p in sys.path:
            sys.path.remove(p)
        sys.path.insert(0, p)
    import thermal_hb2
    return thermal_hb2


def is_on(building):
    return bool(building) and bool(building.get("o13", True))


# ============================================================================================== schedules
def _hourly(s, default=1.0):
    """24 hourly values 0..1 from 'always' / 'never' / [24 values] / {"open": [h0, h1]} (h0 > h1 wraps midnight)."""
    if s is None or s == "always":
        return [float(default)] * 24
    if s == "never":
        return [0.0] * 24
    if isinstance(s, dict):
        h0, h1 = s["open"]
        return [1.0 if ((h0 <= h < h1) if h0 <= h1 else (h >= h0 or h < h1)) else 0.0 for h in range(24)]
    v = [float(x) for x in s]
    if len(v) != 24:
        raise ValueError("an hourly schedule needs 24 values")
    return v


def _blockvals(wx, k, nb, hourly, months=None):
    """Per-block schedule value at the block's first step (same sampling as thermal_hb2's purge mask)."""
    n = len(wx["temp"])
    tm, mo = wx["t_mid"], wx["month"]
    ms = set(int(m) for m in months) if months else None
    out = [0.0] * nb
    for b in range(nb):
        i = min(n - 1, b * k)
        if ms is not None and mo[i] not in ms:
            continue
        out[b] = hourly[int(tm[i] % 24.0)]
    return out


# ============================================================================================== airflow physics
def wind_local(v_met, z, terrain="suburbs"):
    """E+ local wind speed at height z (m) from the met-station wind (10 m, open country)."""
    a, d = TERRAIN.get(terrain, TERRAIN["suburbs"])
    return v_met * (270.0 / 10.0) ** 0.14 * (max(z, 0.1) / d) ** a


def cw_auto(eff_angle, wind_dir):
    dd = abs(float(eff_angle) - float(wind_dir)) % 360.0
    if dd > 180.0:
        dd = 360.0 - dd
    return 0.55 - dd / 180.0 * 0.25


def vent_flow(area, f, v_loc, ta, te, wind_dir, eff_angle, npl, cw="auto", cd="auto"):
    """(Q m3/s, Qw, Qs) - E+ ZoneVentilation:WindandStackOpenArea."""
    if f <= 0.0 or area <= 0.0:
        return 0.0, 0.0, 0.0
    c_w = cw_auto(eff_angle, wind_dir) if cw == "auto" else float(cw)
    dT = abs(ta - te)
    c_d = (0.40 + 0.0045 * dT) if cd == "auto" else float(cd)
    qw = c_w * area * f * v_loc
    qs = c_d * area * f * math.sqrt(2.0 * G * npl * dT / (ta + 273.15))
    return math.sqrt(qw * qw + qs * qs), qw, qs


def door_flow(area_eff, h, ta, tb, cd=0.62):
    """Brown-Solvason bidirectional buoyancy flow through a large vertical opening (m3/s each way)."""
    if area_eff <= 0.0:
        return 0.0
    dT = abs(ta - tb)
    tm = 0.5 * (ta + tb) + 273.15
    return cd / 3.0 * area_eff * math.sqrt(G * h * dT / tm)


def cp_lowrise(theta_deg, side_ratio=1.0):
    """Swami and Chandra (1988) surface-average Cp of a low-rise wall; theta = wind incidence angle to the normal."""
    th = math.radians(abs(((theta_deg + 180.0) % 360.0) - 180.0))
    gl = math.log(max(1e-6, side_ratio))
    x = (1.248 - 0.703 * math.sin(th / 2) - 1.175 * math.sin(th) ** 2 + 0.131 * math.sin(2 * th * gl) ** 3 +
         0.769 * math.cos(th / 2) + 0.07 * gl * gl * math.sin(th / 2) ** 2 + 0.717 * math.cos(th / 2) ** 2)
    return 0.6 * math.log(max(1e-9, x))


# ============================================================================================== zones
def zone_tparams(TH, zone, building):
    over = dict(zone.get("over") or {})
    cls = zone.get("class") or building.get("class")
    if cls and "archetype" not in over:
        over["archetype"] = CLASS_ARCH[cls]
    known = set(p[0] for p in TH.PARAMS)
    return TH.make_tparams(dict((k, v) for k, v in over.items() if k in known))


def _finish(TH, bld, tp):
    """Recompute the radiant star, diffuse / radiant split and beam landing after extra elements were appended
    (same rules as thermal_hb2.assemble, lines 574-614)."""
    els = bld["els"]
    F = [1.0] * len(els)
    for _ in range(30):
        S = sum(e.A * f for e, f in zip(els, F))
        F = [1.0 / max(1e-6, 1.0 - e.A * f / S) for e, f in zip(els, F)]
    for e, f in zip(els, F):
        e.Gr = tp["h_r"] * e.A * f
        e.T = None
    shgc = bld["shgc"]
    a_in = bld["a_in"]
    w_dif = []
    for e in els:
        if e.cls == "glass":
            w_dif.append(e.A * ((e.shgc or shgc) if tp["diffuse_loss"] else 0.0) if e.glass_k is not None else 0.0)
        else:
            w_dif.append(e.A * a_in[e.cls])
    sd = sum(w_dif)
    opq = [e for e in els if e.cls != "glass"]
    A_opq = sum(e.A for e in opq)
    bld["w_dif"] = [w / sd for w in w_dif]
    bld["w_rad"] = [(e.A / A_opq if e.cls != "glass" else 0.0) for e in els]
    win_az = bld["win_az"]
    land = {}
    for cls in ("floor", "wall", "ceiling"):
        cand = [(i, e) for i, e in enumerate(els) if e.cls == cls and
                not (cls == "wall" and e.name.startswith("wall_") and
                     abs(((e.az - win_az + 180) % 360) - 180) < 45)]
        if not cand:
            cand = [(i, e) for i, e in enumerate(els) if e.cls == cls] or [(i, e) for i, e in enumerate(els)
                                                                             if e.cls != "glass"]
        tot = sum(e.A for _i, e in cand)
        land[cls] = [(i, e.A / tot) for i, e in cand]
    bld["land"] = land
    bld["C_tot_kj_m2"] = sum(sum(e.ch["C"]) * e.A for e in els) / 1000.0 / bld["A_f"]


def add_partition(TH, bld, tp, name, area, kind="stud", azimuth=0.0):
    """Append an interzone wall to a zone: far side = film R_FILM_FAR to the other zone's air (drive set per step)."""
    lay = [{"name": "far_film", "d": 0.0, "k": None, "rho": 0.0, "c": 0.0, "R": R_FILM_FAR}] + \
        copy.deepcopy(PART_STACKS[kind])
    dt = float(tp["dt_internal_min"]) * 60.0
    el = TH._make_el(name, float(area), "earth", "wall", lay, tp, dt, az=float(azimuth), tilt=90.0)
    bld["els"].append(el)
    bld.setdefault("iz", []).append(len(bld["els"]) - 1)
    _finish(TH, bld, tp)
    return el


def build(building, TH=None, bld_given=None):
    """[(zone dict, tp, bld)] for every zone, partitions added.  bld_given = {zone_id: (bld, tp)} (the COMBINED window
    zone is assembled by the champion with the planned windows; it is used as is, plus its partitions)."""
    TH = TH or _th()
    out = []
    ids = [z["id"] for z in building["zones"]]
    for z in building["zones"]:
        if bld_given and z["id"] in bld_given:
            bld, tp = bld_given[z["id"]]
        else:
            tp = zone_tparams(TH, z, building)
            bld = TH.assemble(z["env"], tp)
        out.append([z, tp, bld])
    by = dict((z["id"], i) for i, (z, tp, bld) in enumerate(out))
    for p in building.get("partitions") or []:
        ia, ib = by[p["a"]], by[p["b"]]
        az = float(p.get("azimuth_a", 0.0))
        for i, other, a in ((ia, p["b"], az), (ib, p["a"], (az + 180.0) % 360.0)):
            z, tp, bld = out[i]
            el = add_partition(TH, bld, tp, "iz_" + other, p["area_m2"], p.get("type", "stud"), a)
            bld.setdefault("iz_peer", {})[len(bld["els"]) - 1] = by[other]
    if len(set(ids)) != len(ids):
        raise ValueError("zone ids must be unique")
    return out


# ============================================================================================== zone pre-processing
def _prep(TH, bld, wx, climate, win, phi_int, T0):
    """Everything thermal_hb2.simulate computes before its time loop, for one zone (verbatim arithmetic)."""
    tp = bld["tp"]
    els = bld["els"]
    n = len(wx["temp"])
    step_s = float(wx["step_s"])
    k = max(1, int(round(bld["dt"] / step_s)))
    dt = k * step_s
    if abs(dt - bld["dt"]) > 1e-6:
        raise ValueError("O13: the drive step must divide the internal step (got %.0f s vs %.0f s)" % (step_s, bld["dt"]))
    nb = (n + k - 1) // k
    blk = lambda a: TH._block(a, k, nb, n)
    te = blk(wx["temp"])
    tsky_src = wx.get("t_sky")
    tsky = blk(tsky_src) if tsky_src is not None else [t - tp["sky_dtheta_k"] for t in te]
    tg = blk(TH.earth_series(wx, climate["t_earth_c"]))
    qint = blk(phi_int if phi_int is not None else TH.internal_series(wx, tp, bld["A_f"]))
    eta = tp["eta"]
    frad = tp["int_radiant_fraction"]
    h_co = tp["h_out_conv"]
    h_ro = 4.0 * TH.SIGMA * tp["eps_out"] * (283.15 ** 3)
    ne = len(els)
    iz = set(bld.get("iz", []))
    drive = []
    for j, e in enumerate(els):
        if j in iz:
            drive.append([float(T0)] * nb)                    # filled per step with the other zone's air
        elif e.kind == "sun":
            if e.subs:
                At = sum(s[3] for s in e.subs)
                acc = [0.0] * nb
                for (az, tilt, al, Ai) in e.subs:
                    I = blk(TH.plane_series(wx, az, tilt, tp["ground_reflectance"]))
                    sa = TH._solair(e, I, te, tsky, h_co, h_ro)
                    w = Ai / At
                    for b in range(nb):
                        acc[b] += w * sa[b]
                drive.append(acc)
            else:
                I = blk(TH.plane_series(wx, e.az, e.tilt, tp["ground_reflectance"]))
                drive.append(TH._solair(e, I, te, tsky, h_co, h_ro))
        elif e.kind == "glass":
            hs = h_co + h_ro
            fs = e.fsky
            drive.append([(h_co * te[b] + h_ro * (fs * tsky[b] + (1 - fs) * te[b])) / hs for b in range(nb)])
        elif e.kind == "earth":
            drive.append(list(tg))
        elif e.kind == "air":
            drive.append(list(te))
        else:
            drive.append(None)
    src = [[0.0] * nb for _ in range(ne)]
    w_rad = bld["w_rad"]
    for j in range(ne):
        if w_rad[j] > 0:
            f = eta * frad * w_rad[j]
            sj = src[j]
            for b in range(nb):
                sj[b] += f * qint[b]
    conv_int = [eta * (1.0 - frad) * q for q in qint]
    if win is not None:
        w_dif = bld["w_dif"]
        pool = blk(win["diffuse"])
        a_land = bld["a_in"]
        if win.get("alpha_landing") is not None:
            a_land = dict(a_land, floor=win["alpha_landing"])
        for cls, acls in (("floor", "floor"), ("wall", "wall"), ("ceiling", "ceiling"), ("reveal", "wall")):
            if cls not in win:
                continue
            bm = blk(win[cls])
            al = a_land[acls]
            for j, fr in bld["land"][acls]:
                f = eta * al * fr
                sj = src[j]
                for b in range(nb):
                    sj[b] += f * bm[b]
            for b in range(nb):
                pool[b] += (1.0 - al) * bm[b]
        for j in range(ne):
            if w_dif[j] > 0 and els[j].cls != "glass":
                f = eta * w_dif[j]
                sj = src[j]
                for b in range(nb):
                    sj[b] += f * pool[b]
        ga_all = win.get("glass_abs")
        if isinstance(ga_all, dict):
            items = list(ga_all.items())
        elif ga_all is not None:
            items = [(0, ga_all)]
        else:
            items = []
        for kk, ga in items:
            if not bld["glass_idx"]:
                break
            j = bld["glass_idx"][min(int(kk), len(bld["glass_idx"]) - 1)]
            e = els[j]
            gab = blk(ga)
            sj = src[j]
            dv = drive[j]
            ro_per = e.r_o / e.A
            for b in range(nb):
                half = 0.5 * eta * gab[b]
                sj[b] += half
                dv[b] += half * ro_per
    pmask = None
    dH = 0.0
    if tp["n_purge"] > 0:
        months = set(int(m) for m in tp["purge_months"])
        s0, e0 = tp["purge_start_h"], tp["purge_end_h"]
        inwin = (lambda h: s0 <= h or h < e0) if s0 > e0 else (lambda h: s0 <= h < e0)
        tm, mo = wx["t_mid"], wx["month"]
        pmask = [1 if (mo[min(n - 1, b * k)] in months and inwin(tm[min(n - 1, b * k)] % 24.0)) else 0
                 for b in range(nb)]
        dH = tp["rho_cp_air_wh"] * bld["V"] * tp["n_purge"]
    tarp = tp["convection"] == "tarp"
    hmult = tp["h_c_mult"]
    hfix = {"wall": tp["h_c_wall"], "floor": tp["h_c_floor"], "ceiling": tp["h_c_ceiling"], "glass": tp["h_c_wall"]}
    A_all = sum(e.A for e in els)
    wA = [e.A / A_all for e in els]
    mk = [e.A * e.kap for e in els]
    smk = sum(mk) or 1.0
    wM = [m / smk for m in mk]
    D = []
    for j, e in enumerate(els):
        if e.cls == "floor":
            hp, hn = 1.520, 0.760
        elif e.cls == "ceiling":
            hp, hn = 0.760, 1.520
        else:
            hp, hn = 1.31, 1.31
        D.append([e.N, e.cdt, e.a, e.inv, e.e, e.Bconst, e.G_out, None, [0.0] * e.N, drive[j], src[j], e.Gr,
                  e.A * hp * hmult, e.A * hn * hmult, e.A * 0.1 * hmult, hmult * hfix[e.cls] * e.A, float(T0),
                  0.0, 1.0, 0.0, wA[j], wM[j]])
    return {"k": k, "nb": nb, "n": n, "dt": dt, "te": te, "conv_int": conv_int, "D": D, "els": els, "tp": tp,
            "pmask": pmask, "dH": dH, "p_above": tp["purge_above_c"], "p_delta": tp["purge_delta_c"], "tarp": tarp,
            "Ca": bld["C_air"] / dt, "H_ve0": bld["H_ve"], "drive": drive, "V": bld["V"],
            "rcp": tp["rho_cp_air_wh"], "iz_peer": dict(bld.get("iz_peer", {})), "ta": float(T0), "ts": float(T0)}


def _steady_init(Z):
    """thermal_hb2.simulate.steady_init for one zone."""
    nb, te, D, els, dt = Z["nb"], Z["te"], Z["D"], Z["els"], Z["dt"]
    nb0 = min(nb, int(round(86400.0 * 7 / dt)))
    tmean = sum(te[nb - nb0:]) / nb0 if nb0 else 20.0
    troom = tmean + 3.0
    for j, e in enumerate(els):
        dv = Z["drive"][j]
        dmean = (sum(dv[nb - nb0:]) / nb0) if (dv is not None and nb0) else troom
        if j in Z["iz_peer"]:
            dmean = troom
        N = e.N
        if N == 1 or dv is None:
            D[j][7] = [troom] * N
            D[j][16] = troom
            continue
        rs = [1.0 / e.G_out if e.G_out > 0 else 0.0] + [1.0 / g for g in e.a[1:]]
        rin = 1.0 / (2.0 * e.A)
        tot = sum(rs) + rin
        T = []
        acc = 0.0
        for i in range(N):
            acc += rs[i]
            T.append(dmean + (troom - dmean) * acc / tot)
        D[j][7] = T
        D[j][16] = T[-1]
    Z["ta"] = troom
    Z["ts"] = troom


# ============================================================================================== the coupled step
def _reduce(Z, b, Hve, src_out):
    """Eliminate every surface of zone Z at block b (thermal_hb2 run_blocks arithmetic).  Returns a11, a12, a22, r1, r2
    with the outdoor exchange Hve (W/K) at te[b] + extra heat src_out (W) on the air node."""
    ta = Z["ta"]
    tstar = Z["ts"]
    Ca = Z["Ca"]
    teb = Z["te"][b]
    a11 = Ca + Hve
    a12 = 0.0
    a22 = 0.0
    r1 = Ca * ta + Hve * teb + Z["conv_int"][b] + src_out
    r2 = 0.0
    cube = 1.0 / 3.0
    tarp = Z["tarp"]
    for d in Z["D"]:
        if tarp:
            dd = d[16] - ta
            if dd > 0:
                Gc = d[12] * dd ** cube
            else:
                Gc = d[13] * (-dd) ** cube
            if Gc < d[14]:
                Gc = d[14]
        else:
            Gc = d[15]
        N = d[0]
        T = d[7]
        cdt = d[1]
        dv = d[9]
        if N == 1:
            R = cdt[0] * T[0] + d[10][b] + (d[6] * dv[b] if dv is not None else 0.0)
        else:
            inv, a, dps = d[3], d[2], d[8]
            dp = (cdt[0] * T[0] + (d[6] * dv[b] if dv is not None else 0.0)) * inv[0]
            dps[0] = dp
            for i in range(1, N - 1):
                dp = (cdt[i] * T[i] + a[i] * dp) * inv[i]
                dps[i] = dp
            R = cdt[N - 1] * T[N - 1] + d[10][b] + a[N - 1] * dp
        Gr = d[11]
        B = d[5] + Gc + Gr
        d[17] = R
        d[18] = B
        d[19] = Gc
        iB = 1.0 / B
        gB = Gc * iB
        rB = Gr * iB
        a11 += Gc - Gc * gB
        a12 -= Gr * gB
        a22 += Gr - Gr * rB
        r1 += R * gB
        r2 += R * rB
    return a11, a12, a22, r1, r2


def _surfaces(Z, ta, tstar):
    sm_s = 0.0
    sm_m = 0.0
    for d in Z["D"]:
        tsurf = (d[17] + d[19] * ta + d[11] * tstar) / d[18]
        d[16] = tsurf
        T = d[7]
        N = d[0]
        T[N - 1] = tsurf
        if N > 1:
            dps, ee = d[8], d[4]
            nxt = tsurf
            for i in range(N - 2, -1, -1):
                nxt = dps[i] + ee[i] * nxt
                T[i] = nxt
        sm_s += d[20] * tsurf
        sm_m += d[21] * tsurf
    Z["ta"] = ta
    Z["ts"] = tstar
    return sm_s, sm_m


def _solve(A, r):
    """Gaussian elimination with partial pivoting (Z <= 15)."""
    n = len(r)
    M = [list(A[i]) + [r[i]] for i in range(n)]
    for c in range(n):
        p = max(range(c, n), key=lambda i: abs(M[i][c]))
        M[c], M[p] = M[p], M[c]
        piv = M[c][c]
        for i in range(c + 1, n):
            f = M[i][c] / piv
            if f != 0.0:
                Mi, Mc = M[i], M[c]
                for j in range(c, n + 1):
                    Mi[j] -= f * Mc[j]
    x = [0.0] * n
    for i in range(n - 1, -1, -1):
        s = M[i][n] - sum(M[i][j] * x[j] for j in range(i + 1, n))
        x[i] = s / M[i][i]
    return x


_WDIR_WARNED = [False]


class _Net(object):
    """Doors, vents, MVHR and cross-ventilation of a building, evaluated per block with lagged temperatures."""

    def __init__(self, building, zones, wx, k, nb):
        self.zid = [z["id"] for z, tp, bld in zones]
        by = dict((zid, i) for i, zid in enumerate(self.zid))
        n = len(wx["temp"])
        idx = [min(n - 1, b * k) for b in range(nb)]
        self.wind = [float(wx["wind"][i]) for i in idx] if wx.get("wind") is not None else [0.0] * nb
        wd = wx.get("wind_dir")
        self.wdir = [float(wd[i]) for i in idx] if wd is not None else [0.0] * nb
        self.has_wdir = wd is not None
        if wd is None and building.get("vents") and not _WDIR_WARNED[0]:   # D3-002e: was a silent wind_dir = 0 fallback
            _WDIR_WARNED[0] = True
            warnings.warn("o13: weather has no wind_dir; vents use wind_dir = 0 (north) and cross-ventilation is off "
                          "(numerics unchanged)", RuntimeWarning, stacklevel=2)
        terrain = building.get("terrain", "suburbs")
        self.doors = []
        for d in building.get("doors") or []:
            sch = _blockvals(wx, k, nb, _hourly(d.get("sched", "always")), d.get("months"))
            self.doors.append({"a": by[d["a"]], "b": by[d["b"]], "w": float(d.get("w", 0.9)),
                               "h": float(d.get("h", 2.1)), "cd": float(d.get("cd", 0.62)),
                               "leak": float(d.get("leak_m2", 0.015)), "f": sch})
        self.vents = []
        for v in building.get("vents") or []:
            z, tp, bld = zones[by[v["zone"]]]
            zc = float(z.get("z_centroid_m") or 0.5 * float(z["env"].get("height_m", 2.7)))
            self.vents.append({"z": by[v["zone"]], "A": float(v.get("area_m2", 0.0)),
                               "npl": float(v["npl_m"]) if v.get("npl_m") is not None else
                               0.5 * float(v.get("height_m", 1.0)),
                               "az": float(v.get("azimuth_deg", 0.0)), "cw": v.get("cw", "auto"),
                               "cd": v.get("cd", "auto"), "cross": bool(v.get("cross", True)),
                               "f": _blockvals(wx, k, nb, _hourly(v.get("sched", "always")), v.get("months")),
                               "tin_min": float(v.get("tin_min", -100.0)), "tin_max": float(v.get("tin_max", 100.0)),
                               "dt_min": float(v.get("dt_min", -100.0)), "tout_min": float(v.get("tout_min", -100.0)),
                               "tout_max": float(v.get("tout_max", 100.0)), "wind_max": float(v.get("wind_max", 40.0)),
                               "vloc": [wind_local(w, zc, terrain) for w in self.wind]})
        self.mv = []
        for m in building.get("mv") or []:
            z, tp, bld = zones[by[m["zone"]]]
            self.mv.append({"z": by[m["zone"]], "H": tp["rho_cp_air_wh"] * bld["V"] * float(m.get("ach", 0.0)),
                            "eff": float(m.get("hr_eff", 0.0)), "bypass": bool(m.get("bypass", False)),
                            "btin": float(m.get("bypass_tin", 23.0)),
                            "f": _blockvals(wx, k, nb, _hourly(m.get("sched", "always")), m.get("months"))})
        self.vhx = []                       # D3-002e opt-in sensitivity (default off): extra infiltration in vent hours
        vx = building.get("vh_extra")
        if vx:
            fx = _blockvals(wx, k, nb, _hourly(vx.get("sched", "always")), vx.get("months"))
            for zid_, ach_ in (vx.get("zones") or {}).items():
                z, tp, bld = zones[by[zid_]]
                self.vhx.append((by[zid_], tp["rho_cp_air_wh"] * bld["V"] * float(ach_), fx))
        self.cross_cd = float(building.get("cross_cd", 0.60))
        self.side_ratio = float(building.get("side_ratio", 1.0))
        self.stats = {"door_m3h": [0.0] * len(self.doors), "vent_m3h": [0.0] * len(self.zid),
                      "cross_m3h": 0.0, "n": 0, "vent_open_n": [0] * len(self.vents)}

    @staticmethod
    def _route(a, b, door_open, max_doors=3):
        """Shortest zone chain a .. b through doors open this step (None if not connected within max_doors)."""
        if a == b:
            return [a]
        prev = {a: None}
        frontier = [a]
        for _ in range(max_doors):
            nxt = []
            for z in frontier:
                for (p, q), area in door_open.items():
                    if p == z and area > 0.0 and q not in prev:
                        prev[q] = z
                        nxt.append(q)
            if b in prev:
                out = [b]
                while prev[out[-1]] is not None:
                    out.append(prev[out[-1]])
                return out[::-1]
            frontier = nxt
        return None

    def step(self, b, ta, te, rcp, record):
        """(Hout[z] W/K at te, Hin[z][w] W/K from zone w into z) for block b, lagged zone temperatures ta."""
        Z = len(ta)
        Hout = [0.0] * Z
        Hin = [[0.0] * Z for _ in range(Z)]
        door_open = {}
        for k, d in enumerate(self.doors):
            f = d["f"][b]
            q = door_flow(d["w"] * d["h"] * f + d["leak"], d["h"], ta[d["a"]], ta[d["b"]], d["cd"])
            H = rcp * 3600.0 * q
            Hin[d["a"]][d["b"]] += H
            Hin[d["b"]][d["a"]] += H
            door_open[(d["a"], d["b"])] = door_open[(d["b"], d["a"])] = d["w"] * d["h"] * f
            if record:
                self.stats["door_m3h"][k] += q * 3600.0
        open_v = []
        for k, v in enumerate(self.vents):
            f = v["f"][b]
            t_in = ta[v["z"]]
            if f <= 0 or not (v["tin_min"] <= t_in <= v["tin_max"]) or (t_in - te) < v["dt_min"] or \
                    not (v["tout_min"] <= te <= v["tout_max"]) or self.wind[b] > v["wind_max"]:
                continue
            open_v.append((k, f))
            if record:
                self.stats["vent_open_n"][k] += 1             # diagnostics only (D3-002c R1 timing count)
        cross_used = set()
        if self.has_wdir and len(open_v) > 1:
            wdir = self.wdir[b]
            for x in range(len(open_v)):
                for y in range(x + 1, len(open_v)):
                    i, fi = open_v[x]
                    j, fj = open_v[y]
                    vi, vj = self.vents[i], self.vents[j]
                    if not (vi["cross"] and vj["cross"]) or i in cross_used or j in cross_used:
                        continue
                    if abs(((vi["az"] - vj["az"] + 180.0) % 360.0) - 180.0) < 90.0:
                        continue
                    path = [vi["A"] * fi, vj["A"] * fj]
                    route = self._route(vi["z"], vj["z"], door_open)   # zones vi -> vj through open doors
                    if route is None:
                        continue
                    path += [door_open[(p, q_)] for p, q_ in zip(route, route[1:])]
                    aeff = 1.0 / math.sqrt(sum(1.0 / (a * a) for a in path if a > 0))
                    cpi = cp_lowrise(wdir - vi["az"], self.side_ratio)
                    cpj = cp_lowrise(wdir - vj["az"], self.side_ratio)
                    v_loc = 0.5 * (vi["vloc"][b] + vj["vloc"][b])
                    q = self.cross_cd * aeff * v_loc * math.sqrt(abs(cpi - cpj))
                    if cpi < cpj:
                        route = route[::-1]                   # outdoor air enters at the upwind vent
                    H = rcp * 3600.0 * q
                    Hout[route[0]] += H
                    for p, q_ in zip(route, route[1:]):
                        Hin[q_][p] += H                       # directed: each zone's air into the next one
                    cross_used.add(i)
                    cross_used.add(j)
                    if record:
                        self.stats["cross_m3h"] += q * 3600.0
        vctl = {}                                             # zone -> [H of its single-sided vents, T_close]
        for i, f in open_v:
            v = self.vents[i]
            t_in = ta[v["z"]]
            if i in cross_used:                               # its wind part is in the cross flow: stack only
                q, _qw, _qs = vent_flow(v["A"], f, 0.0, t_in, te, self.wdir[b], v["az"], v["npl"], 0.0, v["cd"])
            else:
                q, _qw, _qs = vent_flow(v["A"], f, v["vloc"][b], t_in, te, self.wdir[b], v["az"], v["npl"],
                                        v["cw"], v["cd"])
            H = rcp * 3600.0 * q
            Hout[v["z"]] += H
            c = vctl.setdefault(v["z"], [0.0, -1e9])
            c[0] += H
            c[1] = max(c[1], v["tin_min"], te + v["dt_min"])  # below this the controls close the vent
        for m in self.mv:
            f = m["f"][b]
            if f <= 0:
                continue
            t_in = ta[m["z"]]
            byp = m["bypass"] and t_in > m["btin"] and te < t_in
            Hout[m["z"]] += m["H"] * f * (1.0 if byp else (1.0 - m["eff"]))
        for z_, H_, fx_ in self.vhx:
            Hout[z_] += H_ * fx_[b]
        if record:
            self.stats["n"] += 1
        return Hout, Hin, vctl


# ============================================================================================== simulation
def simulate(building, wx, climate, wins=None, phis=None, keep=("air", "s", "m", "op"), T0=20.0, TH=None,
             zones=None, stats=None):
    """Free-running year of an O13 building.  wins / phis = {zone_id: window_gains dict / internal gain series}.
    Returns {zone_id: {"air", "s", "m", "op", "T_m_end"}} (+ "_flows" diagnostics)."""
    TH = TH or _th()
    zones = zones or build(building, TH)
    wins = wins or {}
    phis = phis or {}
    ZZ = [_prep(TH, bld, wx, climate, wins.get(z["id"]), phis.get(z["id"]), T0) for z, tp, bld in zones]
    nz = len(ZZ)
    nb, k = ZZ[0]["nb"], ZZ[0]["k"]
    if any(Z["nb"] != nb for Z in ZZ):
        raise ValueError("O13: all zones need the same internal step")
    net = _Net(building, zones, wx, k, nb)
    te = ZZ[0]["te"]
    rcp = ZZ[0]["rcp"]
    simple = nz == 1 and not net.doors and not net.vents and not net.mv
    out = [{"air": [0.0] * nb, "s": [0.0] * nb, "m": [0.0] * nb} for _ in range(nz)]

    def run_blocks(order, record):
        for b in order:
            teb = te[b]
            ta_prev = [Z["ta"] for Z in ZZ]
            for Z in ZZ:                                   # partitions: far side = the other zone's lagged air
                for j, w in Z["iz_peer"].items():
                    Z["drive"][j][b] = ta_prev[w]
            Hout, Hin, vctl = (None, None, {}) if simple else net.step(b, ta_prev, teb, rcp, record)
            red = []
            for z, Z in enumerate(ZZ):
                Hve = Z["H_ve0"]
                if Z["pmask"] is not None and Z["pmask"][b] and Z["ta"] > Z["p_above"] and \
                        Z["ta"] > teb + Z["p_delta"]:
                    Hve = Z["H_ve0"] + Z["dH"]
                if Hout is not None:
                    Hve += Hout[z]
                red.append(list(_reduce(Z, b, Hve, 0.0)))

            def solve(red):
                if nz == 1 and (Hin is None or Hin[0][0] == 0.0):
                    a11, a12, a22, r1, r2 = red[0]
                    det = a11 * a22 - a12 * a12
                    return [(r1 * a22 - a12 * r2) / det], [(a11 * r2 - a12 * r1) / det]
                A = [[0.0] * nz for _ in range(nz)]
                rhs = [0.0] * nz
                for z in range(nz):
                    a11, a12, a22, r1, r2 = red[z]
                    A[z][z] = a11 - a12 * a12 / a22 + sum(Hin[z])
                    for w in range(nz):
                        if Hin[z][w]:
                            A[z][w] -= Hin[z][w]
                    rhs[z] = r1 - a12 * r2 / a22
                ta = _solve(A, rhs)
                return ta, [(red[z][4] - red[z][1] * ta[z]) / red[z][2] for z in range(nz)]
            ta_n, ts_n = solve(red)
            # Vent controls act within the step (E+ checks them every system time step, so an open vent that would
            # take the room below its closing temperature cycles): the time-average is a vent that MODULATES to
            # hold the closing temperature.  Solve the vent conductance that gives ta = T_close (others as solved).
            changed = False
            for z, (Hv, Tc) in vctl.items():
                if Hv <= 0.0 or ta_n[z] >= Tc or Tc - teb <= 1e-6:
                    continue
                a11, a12, a22, r1, r2 = red[z]
                inflow = sum(Hin[z][w] * ta_n[w] for w in range(nz)) if Hin is not None else 0.0
                A0 = a11 - Hv - a12 * a12 / a22 + (sum(Hin[z]) if Hin is not None else 0.0)
                B0 = r1 - Hv * teb - a12 * r2 / a22 + inflow
                Hn = min(Hv, max(0.0, (B0 - A0 * Tc) / (Tc - teb)))
                red[z][0] += Hn - Hv
                red[z][3] += (Hn - Hv) * teb
                vctl[z][0] = Hn
                changed = True
            if changed:
                ta_n, ts_n = solve(red)
            if record and vctl:
                for z, (Hv, Tc) in vctl.items():
                    net.stats["vent_m3h"][z] += Hv / (rcp * 3600.0) * 3600.0
            for z, Z in enumerate(ZZ):
                sm_s, sm_m = _surfaces(Z, ta_n[z], ts_n[z])
                if record:
                    o = out[z]
                    o["air"][b] = ta_n[z]
                    o["s"][b] = sm_s
                    o["m"][b] = sm_m

    tp0 = ZZ[0]["tp"]
    if tp0["warmup"] == "year":
        for Z in ZZ:
            for j, e in enumerate(Z["els"]):
                Z["D"][j][7] = [float(T0)] * e.N
        for _ in range(int(tp0["spinup_years"])):
            run_blocks(range(nb), False)
    else:
        for Z in ZZ:
            _steady_init(Z)
        bpd = int(round(86400.0 / ZZ[0]["dt"]))
        wb = min(nb, max(1, int(tp0["warmup_days"]) * bpd))
        order = list(range(nb - wb, nb))
        prev = None
        for p in range(int(tp0["warmup_max"])):
            run_blocks(order, False)
            cur = [x for Z in ZZ for d in Z["D"] for x in d[7]] + [Z["ta"] for Z in ZZ]
            if prev is not None:
                dmax = max(abs(x - y) for x, y in zip(cur, prev))
                if dmax < tp0["warmup_tol_k"] and p >= 1:
                    break
            prev = cur
    run_blocks(range(nb), True)

    def expand(arr):
        if k == 1:
            return array.array("d", arr)
        o = array.array("d")
        prevv = arr[-1]
        for b in range(nb):
            cur = arr[b]
            for q in range(k):
                if b * k + q >= ZZ[0]["n"]:
                    break
                o.append(prevv + (cur - prevv) * (q + 1.0) / k)
            prevv = cur
        return o
    res = {}
    for z, (zz, tp, bld) in enumerate(zones):
        o = out[z]
        r = {}
        if "air" in keep:
            r["air"] = expand(o["air"])
        if "s" in keep:
            r["s"] = expand(o["s"])
        if "m" in keep:
            r["m"] = expand(o["m"])
        if "op" in keep:
            r["op"] = expand([0.5 * (x + y) for x, y in zip(o["air"], o["s"])])
        r["T_m_end"] = o["m"][-1]
        res[zz["id"]] = r
    nrec = max(1, net.stats["n"])
    res["_flows"] = {"door_mean_m3h": [round(x / nrec, 1) for x in net.stats["door_m3h"]], "zones": net.zid,
                     "vent_mean_m3h": [round(x / nrec, 1) for x in net.stats["vent_m3h"]],
                     "cross_mean_m3h": round(net.stats["cross_m3h"] / nrec, 1),
                     "vent_open_blocks": list(net.stats["vent_open_n"]), "blocks": nrec}
    if stats is not None:
        stats.update(res["_flows"])
    return res


# ============================================================================================== COMBINED hook
class _THProxy(object):
    """Stands in for thermal_hb2 inside the champion's plan_combined: run_year runs the O13 building, with the
    champion's own window-zone bld (planned windows) as the window zone, and returns that zone's result."""

    def __init__(self, TH, building, wz):
        self._TH, self._b, self._wz = TH, building, wz
        self.calls = 0

    def __getattr__(self, name):
        return getattr(self._TH, name)

    def run_year(self, bld, wx, climate, win=None, band=None, keep=("air", "s", "m", "op"), T0=20.0):
        if band is not None:
            raise NotImplementedError("O13: ideal-load bands are not supported in the multi-zone path")
        self.calls += 1
        if getattr(self, "_zbld", None) is not bld:          # partitions are appended once per champion bld
            self._zones = build(self._b, self._TH, bld_given={self._wz: (bld, bld["tp"])})
            self._zbld = bld
        zones = self._zones
        if "wind_dir" not in wx and self._b.get("_wind_dir") is not None:
            wx = dict(wx, wind_dir=self._b["_wind_dir"])
        res = simulate(self._b, wx, climate, wins={self._wz: win}, keep=keep, T0=T0, TH=self._TH, zones=zones)
        return res[self._wz]


FAN_MIN_OP_C = 25.0
FAN_LABEL = ("UNVALIDATED: ceiling fan = ASHRAE 55 elevated-air-speed comfort-band offset only, applied where the "
             "upper limit is >= 25 C operative; no reference, no heat-balance effect; assumes occupants control the "
             "fan / speed (no occupant-control cap modelled); hard step at 25 C: an upper limit of 24.9 C gets no "
             "credit, 25.0 C gets the full offset")


def fan_offset(air_speed):
    """ASHRAE 55-2020 adaptive model, elevated air speed: +1.2 K at 0.6 m/s, +1.8 K at 0.9, +2.2 K at >= 1.2."""
    v = float(air_speed or 0.0)
    pts = [(0.3, 0.0), (0.6, 1.2), (0.9, 1.8), (1.2, 2.2)]
    if v <= pts[0][0]:
        return 0.0
    for (v0, d0), (v1, d1) in zip(pts, pts[1:]):
        if v <= v1:
            return d0 + (d1 - d0) * (v - v0) / (v1 - v0)
    return pts[-1][1]


def window_zone(building):
    wz = [z["id"] for z in building["zones"] if z.get("windows")]
    return wz[0] if wz else building["zones"][0]["id"]


def plan_combined(B, scene, weather_epw, construction, params, slider, building):
    """B's nested knapsack with the O13 building as the comfort judge.  The planned windows sit in the window zone
    (building zone with "windows": true); the champion assembles that zone from its env / over, the proxy adds the
    other zones, doors, vents, MVHR and runs them together."""
    TH = B.TH
    wz = window_zone(building)
    zw = [z for z in building["zones"] if z["id"] == wz][0]
    cons = dict(construction or {})
    cons.setdefault("case", building.get("drive_case", "A04_felix_open"))
    cons["env"] = copy.deepcopy(zw["env"])
    over = dict(zw.get("over") or {})
    cls = zw.get("class") or building.get("class")
    if cls and "archetype" not in over:
        over["archetype"] = CLASS_ARCH[cls]
    for k_, v_ in over.items():
        cons.setdefault(k_, v_)
    b2 = dict(building)
    if b2.get("_wind_dir") is None and building.get("epw_wind_dir", True):
        b2["_wind_dir"] = epw_wind_dir(weather_epw or B.EPW_DEFAULT)
    proxy = _THProxy(TH, b2, wz)
    band0 = B._band
    off = fan_offset((building.get("fan") or {}).get("air_speed", 0.0))

    def band_fan(op_ref_te, cat):
        # D3-002c (AU-R3-002 R3): ASHRAE 55 elevated air speed credit only where the operative T is >= 25 C, i.e. only
        # an upper limit already at >= FAN_MIN_OP_C is raised (conservative: no fan credit below 25 C).  UNVALIDATED.
        lo, hi = band0(op_ref_te, cat)
        return lo, [h + off if h >= FAN_MIN_OP_C else h for h in hi] if off else hi
    B.TH = proxy
    B._band = band_fan
    try:
        res = B.plan_combined(scene, weather_epw, cons, params, slider)
    finally:
        B.TH = TH
        B._band = band0
    for p in res.values():
        p["o13"] = {"version": VERSION, "building": building.get("name"), "window_zone": wz,
                    "zones": [z["id"] for z in building["zones"]], "thermal_runs": proxy.calls,
                    "fan_band_offset_k": off, "fan_min_op_c": FAN_MIN_OP_C, "fan_label": FAN_LABEL}
    return res


def epw_wind_dir(epw, sph=6):
    """Wind direction per 10-min step (deg from N) from the EPW (hour-ending values; linear by vector between
    hours, as E+ interpolates weather)."""
    d = []
    with open(epw) as f:
        for i, ln in enumerate(f):
            if i < 8:
                continue
            p = ln.split(",")
            if len(p) > 21:
                d.append(float(p[20]))
    out = array.array("d")
    n = len(d)
    for h in range(n):
        a0, a1 = math.radians(d[h - 1] if h > 0 else d[-1]), math.radians(d[h])
        for s in range(sph):
            w = (s + 1.0) / sph
            x = (1 - w) * math.sin(a0) + w * math.sin(a1)
            y = (1 - w) * math.cos(a0) + w * math.cos(a1)
            out.append(math.degrees(math.atan2(x, y)) % 360.0)
    return out
