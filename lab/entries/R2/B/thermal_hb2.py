# -*- coding: utf-8 -*-
"""
thermal_hb2.py - MODELER-B "HB-L2": the round-2 version of the THERMAL champion HB-L (thermal_hb.py, R1, kept frozen).
Same physics family (layered heat balance: finite-volume nodes through every element's real layers, air node + TARP
convection + Walton radiant star, sun absorbed where the beam lands, implicit internal step), standard library only.

What is new in R2 (SHORTFALLS_R1 thermal items, in order):
  1. SPEED (target <= 0.6 s per case-year, R1 2.6 s):
     * converge-and-stop warm-up (EnergyPlus practice): the LAST `warmup_days` of the year (the days right before
       1 Jan in a cyclic year) are repeated from a steady-state start until no node moves more than `warmup_tol_k`
       between passes (min 2, max `warmup_max` passes) - instead of simulating a whole extra year;
     * reduced wall interiors: walls of the same construction and outside boundary are ONE chain driven by the
       area-weighted sol-air temperature of their facings (linear superposition: exact for the room except TARP's
       cube-root and Walton's view factors, both negligible), soil under ground floors on a stretched grid (fine at
       the slab, coarse at depth), stretched grid inside thick massive layers;
     * a flat, attribute-free inner loop.
  2. GLASS NODE: absorbed solar on the glass (W, `glass_abs`, e.g. EnergyPlus's own) is put half on each face of the
     glass layer (the EnergyPlus window heat balance for one equivalent layer), the outer half through a Thevenin
     sol-air term, so the inward-flowing share comes out of the glass heat balance instead of a fixed f_in.
     glass_model "node" (application mode): incident sun BEFORE the glass is split by the glass's own tau / alpha
     from (U, SHGC) with EnergyPlus's SimpleGlazingSystem correlations (tau lands where the beam lands, alpha heats the
     glass node).  Default "shgc" = R1 semantics (g x incident enters the room), so the R1 drive is unchanged.
  3. U OVERRIDE ON UNINSULATED WALLS: when a stack has no (or not enough) massless insulation to take a higher U, the
     conductivity of its massive layers is scaled (rho*c kept), so U is matched exactly; capped at the films' limit
     (warning in tp['_warnings']).
  4. MULTIPLE WINDOWS: env['glazing'] may list 1-15 windows (own azimuth, tilt, area, optional own u / shgc); each is
     its own glass element (own sky view and sol-air), and window_gains() takes one landing record per window (lists)
     and sums them into the room.

INTERFACE = the R1 one (THERMAL_PARAMS: make_tparams / assemble / window_gains / run_year), plus param_table().
    tp  = make_tparams(over, building)
    bld = assemble(env, tp)
    win = window_gains(bld, beam_cls, diffuse, glass_abs=None)     # dicts / arrays, or LISTS of them (one per window)
    out = run_year(bld, wx, climate, win=win, band=None)           # air, s (MRT), m (mass), op [, hc]
"""
from __future__ import division

import array
import math
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
if HERE not in sys.path:
    sys.path.insert(0, HERE)
CORE_DIR = os.path.abspath(os.path.join(HERE, "..", "..", "..", ".."))  # project root
import hb_constructions as HC  # noqa: E402

VERSION = "HB-L2 (R2, MODELER-B)"
SIGMA = 5.670e-8
ARCH_IDS = ["A%02d" % i for i in range(1, 21)]
GLAZING_IDS = ["G%02d" % i for i in range(0, 10)]
WALL_IDS = ["W%02d" % i for i in range(1, 19)]
ROOF_IDS = ["R%02d" % i for i in range(1, 7)]
FLOOR_IDS = ["F%02d" % i for i in range(1, 7)]

# name, default, min, max (or choices), units, group, meaning, expected response of the mean operative temperature
# when the parameter RISES ('+' warmer, '-' cooler, '0' none / numerics, '~' depends on climate / season, 'c' choice).
# 'A-compat' rows are accepted so one parameter set drives A and B, but do nothing in B; they are NOT advertised as
# adjustable (param_table(advertised_only=True) leaves them out).
PARAMS = [
    ("archetype", "A19", ARCH_IDS, None, "id", "construction", "building type from constructions.json (A01-A20)", "c"),
    ("wall", None, WALL_IDS, None, "id", "construction", "wall construction instead of the archetype's", "c"),
    ("roof", None, ROOF_IDS, None, "id", "construction", "roof / ceiling construction instead of the archetype's", "c"),
    ("floor", None, FLOOR_IDS, None, "id", "construction", "floor construction instead of the archetype's", "c"),
    ("glazing", None, GLAZING_IDS, None, "id", "construction", "glazing product (U + SHGC) of the windows", "c"),
    ("u_wall", None, 0.10, 6.0, "W/m2K", "construction", "wall U (insulation resized; uninsulated: layer k scaled)", "~"),
    ("u_roof", None, 0.08, 4.0, "W/m2K", "construction", "roof / ceiling U", "~"),
    ("u_floor", None, 0.10, 2.5, "W/m2K", "construction", "floor U (slab on ground: ISO 13370 value)", "~"),
    ("u_glazing", None, 0.5, 7.9, "W/m2K", "construction", "window U (whole window)", "~"),
    ("shgc", None, 0.10, 0.85, "-", "construction", "window solar heat gain coefficient (g)", "+"),
    ("kappa_wall", None, 0.0, 400.0, "kJ/m2K", "mass", "room-side heat capacity of the walls (layers inside the "
                                                         "insulation replaced by that much masonry)", "~"),
    ("kappa_roof", None, 0.0, 400.0, "kJ/m2K", "mass", "room-side heat capacity of the ceiling (as above)", "~"),
    ("kappa_floor", None, 0.0, 400.0, "kJ/m2K", "mass", "room-side heat capacity of the exposed floor (as above)", "~"),
    ("kappa_floor_covered", 40.0, 0.0, 400.0, "kJ/m2K", "mass", "A-compat (B models the cover as R 0.3)", "A"),
    ("floor_exposed_fraction", None, 0.0, 1.0, "-", "mass", "share of the floor left exposed (1 = polished slab)", "~"),
    ("floor_ground", None, [True, False], None, "bool", "construction", "floor on the ground or suspended", "c"),
    ("mass_class_cm", None, 50.0, 450.0, "kJ/m2K floor", "mass", "ISO class heat capacity per m2 floor (reached by "
                                                                    "the exposed slab share + masonry partitions)", "~"),
    ("kappa_mode", "rule100", ["rule100", "iso13786", "library"], None, "-", "mass", "A-compat", "A"),
    ("partition_m2_per_floor", 0.0, 0.0, 2.5, "m2/m2", "mass", "internal partition wall area per m2 floor", "~"),
    ("kappa_partition", 9.0, 0.0, 300.0, "kJ/m2K", "mass", "partition heat capacity per side (stud 9, brick 155)", "~"),
    ("a_t_mode", "geometry", ["geometry", "iso"], None, "-", "advanced", "A-compat", "A"),
    ("h_ms", 9.1, 6.0, 12.0, "W/m2K", "advanced", "A-compat (5R1C only)", "A"),
    ("h_is", 3.45, 2.5, 5.0, "W/m2K", "advanced", "A-compat (5R1C only)", "A"),
    ("n_rest", 5.0, 0.0, 10.0, "1/h", "ventilation", "air changes: infiltration + opened windows", "~"),
    ("n_system", 0.0, 0.0, 10.0, "1/h", "ventilation", "mechanical ventilation air changes", "~"),
    ("heat_recovery", 0.0, 0.0, 0.95, "-", "ventilation", "heat recovery efficiency", "~"),
    ("n_purge", 0.0, 0.0, 20.0, "1/h", "ventilation", "extra NIGHT air changes in summer when warmer than outside", "-"),
    ("purge_above_c", 23.0, 18.0, 30.0, "degC", "ventilation", "purge only when the room air is above this", "+"),
    ("purge_delta_c", 2.0, 0.5, 6.0, "K", "ventilation", "purge only when the room is this much warmer than outside", "+"),
    ("purge_start_h", 20.0, 0.0, 24.0, "h", "ventilation", "purge window start", "~"),
    ("purge_end_h", 7.0, 0.0, 24.0, "h", "ventilation", "purge window end", "~"),
    ("purge_months", [11, 12, 1, 2, 3], None, None, "months", "ventilation", "months in which the purge may run", "c"),
    ("alpha_wall", 0.5, 0.1, 0.95, "-", "solar", "outside solar absorptance of the walls", "+"),
    ("alpha_roof", 0.4, 0.1, 0.95, "-", "solar", "outside solar absorptance of the roof", "+"),
    ("sky_dtheta_k", 11.0, 0.0, 20.0, "K", "solar", "sky below air temperature (when the weather has no sky temp)", "-"),
    ("ground_reflectance", 0.2, 0.0, 0.8, "-", "solar", "ground albedo for sun on walls", "+"),
    ("r_dirt", 0.95, 0.7, 1.0, "-", "solar", "dirt / frame factor on the window's solar gain", "+"),
    ("alpha_landing", 0.6, 0.2, 0.95, "-", "solar", "solar absorptance of the floor where the beam lands", "+"),
    ("beam_to_mass", True, [True, False], None, "bool", "solar", "True: beam absorbed where it lands; False: spread "
                                                                 "like diffuse light", "c"),
    ("kappa_mass_full", 100.0, 30.0, 300.0, "kJ/m2K", "solar", "A-compat", "A"),
    ("internal_profile", "dwelling", ["dwelling", "flat", "workbook"], None, "-", "gains", "daily internal-gain shape",
     "c"),
    ("occupants", 4.0, 0.0, 20.0, "people", "gains", "people living in the space", "+"),
    ("w_per_person", 80.0, 50.0, 150.0, "W", "gains", "sensible heat per person", "+"),
    ("home_day_fraction", 0.25, 0.0, 1.0, "-", "gains", "share of the people at home 9-17 h", "+"),
    ("appliance_w_m2", 1.0, 0.0, 10.0, "W/m2", "gains", "appliances + lighting, daily mean per m2 floor", "+"),
    ("heating_days", 9.0, 0.0, 31.0, "days", "gains", "workbook profile only", "+"),
    ("loads", "free", ["free", "ideal"], None, "-", "loads", "free-running or ideal heating/cooling", "c"),
    ("heat_setpoint_c", 20.0, 12.0, 24.0, "degC", "loads", "ideal heating set-point (air)", "+"),
    ("cool_setpoint_c", 26.0, 22.0, 32.0, "degC", "loads", "ideal cooling set-point (air)", "+"),
    ("eta", 1.0, 0.5, 1.0, "-", "numerics", "gain utilisation multiplier on all gains (1 = physical)", "+"),
    ("scheme", "cn", ["cn", "euler"], None, "-", "numerics", "A-compat (B is implicit)", "A"),
    ("spinup_years", 1, 0, 3, "years", "numerics", "used only with warmup 'year' (R1 behaviour)", "0"),
    # ---- B's own
    ("warmup", "converge", ["converge", "year"], None, "-", "numerics", "B: converge-and-stop warm-up on the last days "
                                                                       "of the year, or R1's spin-up years", "0"),
    ("warmup_days", 14, 7, 60, "days", "numerics", "B: length of the warm-up period (days before 1 Jan)", "0"),
    ("warmup_tol_k", 0.05, 0.005, 0.5, "K", "numerics", "B: warm-up converged when no node moves more than this", "0"),
    ("warmup_max", 8, 2, 25, "passes", "numerics", "B: most warm-up passes", "0"),
    ("merge_walls", True, [True, False], None, "bool", "numerics", "B: one chain for walls of the same construction "
                                                                   "(area-weighted sol-air of their facings)", "0"),
    ("grid_stretch", 1.5, 1.0, 2.5, "-", "numerics", "B: node growth factor inside thick layers / soil (1 = even)", "0"),
    ("codegen", True, [True, False], None, "bool", "numerics", "B: run the step as generated straight-line code "
                                                               "(same arithmetic, ~3x faster)", "0"),
    ("glass_model", "shgc", ["shgc", "node"], None, "-", "solar", "B: 'shgc' = g x incident enters (R1); 'node' = tau "
                                                                  "lands, alpha heats the glass node (E+ simple glazing)",
     "c"),
    ("glass_films", "model", ["model", "rating"], None, "-", "solar", "B: glass = its layer resistance (E+ simple "
                                                                       "glazing) + this model's outside film ('model', "
                                                                       "as E+), or the rated U less 0.13 (R1)", "c"),
    ("dt_internal_min", 10.0, 1.0, 60.0, "min", "numerics", "B: internal time step", "0"),
    ("dx_max_m", 0.04, 0.01, 0.10, "m", "numerics", "B: largest node spacing near the surfaces of massive layers", "0"),
    ("convection", "tarp", ["tarp", "fixed"], None, "-", "room", "B: inside convection: TARP (E+ default) or fixed", "c"),
    ("h_c_wall", 2.0, 0.5, 8.0, "W/m2K", "room", "B: fixed-mode convection, walls", "~"),
    ("h_c_floor", 2.0, 0.5, 8.0, "W/m2K", "room", "B: fixed-mode convection, floor", "~"),
    ("h_c_ceiling", 1.0, 0.3, 8.0, "W/m2K", "room", "B: fixed-mode convection, ceiling", "~"),
    ("h_c_mult", 1.0, 0.3, 3.0, "-", "room", "B: multiplier on inside convection (air movement, fans)", "~"),
    ("h_r", 5.0, 3.0, 6.5, "W/m2K", "room", "B: inside long-wave radiation coefficient", "~"),
    ("h_out_conv", 10.0, 3.0, 40.0, "W/m2K", "room", "B: outside convection coefficient (wind-averaged)", "~"),
    ("eps_out", 0.9, 0.1, 0.95, "-", "room", "B: outside long-wave emissivity (sky exchange)", "-"),
    ("alpha_in_wall", 0.5, 0.1, 0.95, "-", "room", "B: inside solar absorptance of walls", "+"),
    ("alpha_in_ceiling", 0.3, 0.1, 0.95, "-", "room", "B: inside solar absorptance of the ceiling", "+"),
    ("diffuse_loss", True, [True, False], None, "bool", "room", "B: part of the reflected / diffuse sun leaves "
                                                                 "through the glass again", "c"),
    ("int_radiant_fraction", 0.5, 0.0, 1.0, "-", "room", "B: radiant share of the internal gains", "~"),
    ("air_cap_mult", 1.0, 1.0, 20.0, "-", "room", "B: air heat capacity multiplier (furniture)", "~"),
    ("rho_cp_air_wh", 0.333, 0.30, 0.36, "Wh/m3K", "room", "B: volumetric heat capacity of air", "~"),
    ("soil_depth_m", 0.5, 0.1, 3.0, "m", "room", "B: soil modelled under a ground floor", "~"),
]
_P = dict((p[0], p) for p in PARAMS)


def make_tparams(over=None, building=None):
    """Defaults <- building inputs (A's names) <- ``over``; numbers clamped (listed in tp['_warnings']).  Unknown names
    raise KeyError, except the comfort / slider names a combined caller may pass along (ignored here)."""
    tp = dict((p[0], p[1] if not isinstance(p[1], list) else list(p[1])) for p in PARAMS)
    warn = []
    b = building or {}
    for src, dst in (("n_rest", "n_rest"), ("n_system", "n_system"), ("heat_recovery", "heat_recovery"),
                     ("occupants", "occupants"), ("w_per_person", "w_per_person"),
                     ("roof_absorptivity", "alpha_roof"), ("heating_days", "heating_days")):
        if b.get(src) is not None:
            tp[dst] = b[src]
    for k, v in (over or {}).items():
        if k not in _P:
            if k in _IGNORED:
                continue
            raise KeyError("unknown thermal parameter %r" % k)
        tp[k] = v
    for p in PARAMS:
        name, default, lo, hi, units = p[0], p[1], p[2], p[3], p[4]
        v = tp[name]
        if v is None:
            continue
        if isinstance(lo, list):
            if v not in lo:
                raise ValueError("%s = %r: choose one of %s" % (name, v, lo))
        elif lo is not None and isinstance(v, (int, float)) and not isinstance(v, bool):
            c = min(max(float(v), lo), hi)
            if c != v:
                warn.append("%s %r clamped to %g (range %g-%g %s)" % (name, v, c, lo, hi, units))
            tp[name] = c if not isinstance(default, int) or isinstance(default, bool) else int(round(c))
    tp["_warnings"] = warn
    return tp


# comfort / slider names that belong to the combined planner (accepted and ignored by the thermal model)
_IGNORED = set(["beam_priority", "comfort_model", "comfort_category", "comfort_lo_c", "comfort_hi_c",
                "occupied_start_h", "occupied_end_h"])


def param_table(advertised_only=True):
    """[{name, default, min, max | choices, units, group, meaning, response}] - every adjustable parameter."""
    out = []
    for name, default, lo, hi, units, group, meaning, resp in PARAMS:
        if advertised_only and resp == "A":
            continue
        row = {"name": name, "default": default, "units": units, "group": group, "meaning": meaning,
               "response": {"+": "warmer when raised", "-": "cooler when raised", "0": "none (numerics)",
                            "~": "depends on season / climate", "c": "choice"}.get(resp, resp)}
        if isinstance(lo, list):
            row["choices"] = lo
        elif lo is None:
            row["choices"] = None
        else:
            row["min"], row["max"] = lo, hi
        out.append(row)
    return out


def write_params_md(path=os.path.join(HERE, "THERMAL_PARAMS_B2.md")):
    lines = ["# Thermal model B (HB-L2, round 2) parameters - generated from thermal_hb2.PARAMS", "",
             "Same names as model A's where they overlap.  'A-compat' rows are accepted but do nothing in B and are not",
             "advertised as adjustable.  response = what the mean operative temperature does when the value rises.", "",
             "| group | name | default | range / choices | units | response | meaning |", "|---|---|---|---|---|---|---|"]
    for name, default, lo, hi, units, group, meaning, resp in PARAMS:
        rng = (", ".join(str(x) for x in lo) if isinstance(lo, list) else ("-" if lo is None else "%g .. %g" % (lo, hi)))
        lines.append("| %s | `%s` | %s | %s | %s | %s | %s |" % (group, name, default, rng, units, resp, meaning))
    with open(path, "w") as f:
        f.write("\n".join(lines) + "\n")
    return path


# ============================================================================================ glazing optics
def simple_glazing(U, shgc):
    """EnergyPlus WindowMaterial:SimpleGlazingSystem correlations (Engineering Reference, 'Simple Window Model'):
    -> dict(tau=solar transmittance, alpha=solar absorptance, f_in=inward-flowing fraction of the absorbed solar at
    the rating conditions, r_layer=glass layer resistance m2K/W).  tau + f_in * alpha = SHGC."""
    U = max(0.3, min(7.9, float(U)))
    g = max(0.01, min(0.95, float(shgc)))
    ri_w = 1.0 / (0.359073 * math.log(U) + 6.949915) if U < 5.85 else 1.0 / (1.788041 * U - 2.886625)
    ro_w = 1.0 / (0.025342 * U + 29.163853)
    r_l = max(1e-4, 1.0 / U - ri_w - ro_w)

    def t_hi(s):                                   # U > 4.5
        return 0.939998 * s * s + 0.20332 * s if s < 0.7206 else 1.30415 * s - 0.30515

    def t_lo(s):                                   # U < 3.4
        return 0.4104 * s if s <= 0.15 else 0.085775 * s * s + 0.963954 * s - 0.084958

    def films(s, d):
        if s == "hi":
            return (1.0 / (29.436546 * d ** 3 - 21.943415 * d ** 2 + 9.945872 * d + 7.426151),
                    1.0 / (2.225824 * d + 20.577080))
        return (1.0 / (199.8208128 * d ** 3 - 90.639733 * d ** 2 + 19.737055 * d + 6.766575),
                1.0 / (5.763355 * d + 20.541528))
    if U >= 4.5:
        f = 1.0
    elif U <= 3.4:
        f = 0.0
    else:
        f = (U - 3.4) / 1.1
    tau = f * t_hi(g) + (1 - f) * t_lo(g)
    tau = max(0.0, min(g, tau))
    d = g - tau
    rih, roh = films("hi", d)
    ril, rol = films("lo", d)
    ri_s, ro_s = f * rih + (1 - f) * ril, f * roh + (1 - f) * rol
    f_in = (ro_s + 0.5 * r_l) / (ro_s + r_l + ri_s)
    alpha = min(1.0 - tau, d / f_in) if f_in > 0 else 0.0
    return {"tau": tau, "alpha": alpha, "f_in": f_in, "r_layer": r_l, "ri_w": ri_w, "ro_w": ro_w}


# ============================================================================================ element chains
_GENERIC = {"masonry": {"k": 1.0, "rho": 2000.0, "c": 1000.0}, "board": {"k": 0.25, "rho": 900.0, "c": 1000.0}}


def _set_kappa(lay, kappa):
    """Replace the layers inside the main insulation by `kappa` kJ/m2K of generic masonry (board if light)."""
    nm = [i for i, l in enumerate(lay) if l["k"] is None and l["name"] not in ("carpet", "deep_ground")]
    cut = (max(nm, key=lambda i: lay[i]["R"]) + 1) if nm else 0
    keep = lay[:cut]
    carpet = [l for l in lay[cut:] if l["name"] == "carpet"]
    g = _GENERIC["masonry" if kappa > 30 else "board"]
    d = max(0.002, kappa * 1000.0 / (g["rho"] * g["c"]))
    return keep + [{"name": "generic_mass", "d": d, "k": g["k"], "rho": g["rho"], "c": g["c"], "R": d / g["k"]}] + carpet


def _match_u(lay, kind, U, ground, warn=None, tag=""):
    """Make the stack's U equal U: resize the largest massless (insulation / cavity) layer; when that cannot go low
    enough (no insulation, uninsulated wall with a higher U asked for) scale the conductivity of the massive layers
    (rho*c kept, so the heat capacity is unchanged).  Capped at the films' limit (warning)."""
    lay = [dict(l) for l in lay]
    films = 1.0 / HC.u_of([], kind, ground)                     # Rsi (+ Rse)
    r_target = 1.0 / U
    diff = r_target - (1.0 / HC.u_of(lay, kind, ground))
    if abs(diff) <= 0.005:
        return lay
    nm = [l for l in lay if l["k"] is None and l["name"] not in ("carpet", "deep_ground", "soil")]
    if diff > 0:
        if nm:
            big = max(nm, key=lambda l: l["R"])
            big["R"] += diff
        else:
            lay.insert(0, {"name": "adjust", "d": 0.0, "k": None, "rho": 0.0, "c": 0.0, "R": diff})
        return lay
    need = -diff                                                 # resistance to remove
    for l in sorted(nm, key=lambda l: -l["R"]):
        take = min(need, max(0.0, l["R"] - 0.01))
        l["R"] -= take
        need -= take
        if need <= 1e-6:
            return lay
    ms = [l for l in lay if l["k"] is not None and l["name"] not in ("soil",)]
    r_m = sum(l["R"] for l in ms)
    if ms and r_m > 1e-6:
        keep_min = 0.02 * r_m                                    # never below 2 % of the layers' own resistance
        new_rm = max(keep_min, r_m - need)
        if r_m - need < keep_min and warn is not None:
            warn.append("%s U %.2f above what the layers allow; capped at U %.2f" % (
                tag or kind, U, 1.0 / (films + sum(l["R"] for l in lay if l["k"] is None) + keep_min)))
        f = r_m / new_rm                                         # conductivity multiplier
        for l in ms:
            l["k"] *= f
            l["R"] = l["d"] / l["k"]
    elif warn is not None:
        warn.append("%s U %.2f cannot be reached (no layers to change)" % (tag or kind, U))
    return lay


def _chain(lay, dx_max, stretch=1.5):
    """Finite-volume nodes for a layer stack (outside -> inside).  Leading massless layers are merged into the outside
    boundary resistance (R_lead), the rest become nodes.  Massive layers: node size <= dx_max near both faces, growing
    by `stretch` towards the middle; soil: fine at the slab, growing with depth.  -> dict(C=[J/m2K], G=[W/m2K], R_lead)."""
    i0 = 0
    R_lead = 0.0
    while i0 < len(lay) and lay[i0]["k"] is None:
        R_lead += lay[i0]["R"]
        i0 += 1
    C = [0.0]
    G = []
    pend_R = 0.0
    for l in lay[i0:]:
        if l["k"] is None:
            pend_R += l["R"]
            continue
        if pend_R > 0:
            G.append(1.0 / pend_R)
            C.append(0.0)
            pend_R = 0.0
        alpha = l["k"] / (l["rho"] * l["c"])
        dx = min(dx_max, max(0.006, 0.8 * math.sqrt(alpha * 3600.0)))
        hs = _cells(l["d"], dx, stretch, soil=(l["name"] == "soil"))
        for h in hs:
            cap = l["rho"] * l["c"] * h
            C[-1] += 0.5 * cap
            G.append(l["k"] / h)
            C.append(0.5 * cap)
    if pend_R > 0:
        G.append(1.0 / pend_R)
        C.append(0.0)
    return {"C": C, "G": G, "R_lead": R_lead}


def _cells(d, dx, q, soil=False):
    """Cell sizes across a layer of thickness d: even (q = 1, at most 12 cells as in R1) or stretched (first cell dx
    at the inside face for soil, at both faces otherwise, growing by q, capped at 4 dx)."""
    if q <= 1.0001 or d <= 3.0 * dx:
        n = max(1, min(12, int(math.ceil(d / dx - 1e-9))))
        return [d / n] * n
    if soil:                                         # outside -> inside: coarse (deep) first, fine at the slab
        hs = []
        h, s = dx, 0.0
        while s + h < d - 1e-9:
            hs.append(h)
            s += h
            h = min(4.0 * dx, h * q)
        rem = d - s
        if hs and rem < 0.5 * hs[-1]:
            hs[-1] += rem
        else:
            hs.append(rem)
        return list(reversed(hs))
    half = []
    h, s = dx, 0.0
    while s + h < d / 2.0 - 1e-9:
        half.append(h)
        s += h
        h = min(4.0 * dx, h * q)
    rem = d / 2.0 - s
    if half and rem < 0.5 * half[-1]:
        half[-1] += rem
        return half + list(reversed(half))
    return half + [2.0 * rem] + list(reversed(half))


class _El(object):
    """One heat-flow path: area A, node chain, outside boundary kind ('sun' / 'earth' / 'air' / 'adiabatic' /
    'glass'), inside class ('wall' / 'floor' / 'ceiling' / 'glass').  subs = [(azimuth, tilt, alpha, area)] of the
    facings merged into a 'sun' element (their area-weighted sol-air drives it)."""
    __slots__ = ("name", "A", "kind", "cls", "az", "tilt", "alpha", "ch", "N", "cdt", "a", "c", "inv", "e",
                 "Bconst", "G_out", "T", "Gr", "Gc", "kap", "fsky", "src", "subs", "r_o", "glass_k", "u_glass",
                 "shgc")


def _make_el(name, A, kind, cls, lay, tp, dt, az=0.0, tilt=90.0, alpha=0.0, u_glass=None, subs=None):
    el = _El()
    el.name, el.A, el.kind, el.cls, el.az, el.tilt, el.alpha = name, A, kind, cls, az, tilt, alpha
    h_ro = 4.0 * SIGMA * tp["eps_out"] * (283.15 ** 3)
    el.r_o = 0.0
    el.u_glass = u_glass
    if kind == "glass":
        ch = {"C": [0.0], "G": [], "R_lead": 0.0}
        r_film = 1.0 / (tp["h_out_conv"] + h_ro)
        if tp.get("glass_films", "model") == "model":
            # R2: the glass LAYER resistance from the rating (EnergyPlus simple glazing: 1/U - rating films) + this
            # model's own outside film (as EnergyPlus does), not the rated U with a 0.13 inside film taken off
            r_out = simple_glazing(u_glass, 0.5)["r_layer"] + r_film
            el.r_o = r_film
        else:                                                   # R1 rule
            r_out = max(1.0 / u_glass - HC.RSI["wall"], 0.005)
            el.r_o = min(r_out, r_film)
        G_out = 1.0 / r_out
    else:
        ch = _chain(lay, tp["dx_max_m"], tp["grid_stretch"])
        if kind == "sun" or kind == "air":
            G_out = 1.0 / (1.0 / (tp["h_out_conv"] + h_ro) + ch["R_lead"])
        elif kind == "earth":
            G_out = 1.0 / max(ch["R_lead"], 1e-4)
        else:
            G_out = 0.0
    el.ch = ch
    N = len(ch["C"])
    el.N = N
    G_out *= A
    el.G_out = G_out
    cdt = [A * c / dt for c in ch["C"]]
    a = [0.0] + [A * g for g in ch["G"]]
    c = [A * g for g in ch["G"]] + [0.0]
    inv, e = [0.0] * N, [0.0] * N
    for i in range(N - 1):
        b = cdt[i] + a[i] + c[i] + (G_out if i == 0 else 0.0)
        den = b - (a[i] * e[i - 1] if i > 0 else 0.0)
        inv[i] = 1.0 / den
        e[i] = c[i] * inv[i]
    el.Bconst = (cdt[N - 1] + a[N - 1] * (1.0 - e[N - 2])) if N > 1 else (cdt[0] + G_out)
    el.cdt, el.a, el.c, el.inv, el.e = cdt, a, c, inv, e
    el.kap = sum(ch["C"][-max(1, min(N, 4)):]) / 1000.0
    fs = (1.0 + math.cos(math.radians(tilt))) / 2.0
    el.fsky = fs * math.sqrt(fs)
    el.subs = subs
    el.src, el.T, el.Gc, el.Gr = None, None, 0.0, 0.0
    el.glass_k = None
    el.shgc = None
    return el


def _archetype(tp, lib):
    a = lib["archetypes"][tp["archetype"]]
    E = lib["elements"]
    uo = a.get("U_overrides") or {}
    wall = tp.get("wall") or a.get("wall") or a.get("wall_mass_like") or "W02"
    roof = tp.get("roof") or a.get("roof") or "R02"
    floor = tp.get("floor") or a.get("floor") or "F01"
    glz = tp.get("glazing") or a.get("glazing") or "G00"
    s = {"wall": E[wall], "roof": E[roof], "floor": E[floor], "glazing": E[glz],
         "u_wall": uo.get("u_wall"), "u_roof": uo.get("u_roof"), "u_floor": uo.get("u_floor"),
         "u_glazing": uo.get("u_glazing"), "shgc": uo.get("g_value"),
         "fexp": a.get("floor_exposed_fraction", 1.0), "cm": None}
    if "wall" not in a and "wall_mass_like" not in a and a.get("cm_per_floor_est"):
        s["cm"] = float(a["cm_per_floor_est"])
    return s


def _stack(lib, el, kind, U, ground=False, covered=False, warn=None):
    """hb_constructions' stack (the rule shared with the E+ runs), then U matched exactly (R2: also when the stack has
    no insulation to resize)."""
    st = HC.stack(lib, el, kind, U, ground=ground, covered=covered)
    if abs(HC.u_of(st, kind, ground) - U) > 0.01 * U and not ground:
        st = _match_u(st, kind, U, ground, warn, el.get("id", kind))
    return st


def assemble(env, tp, lib=None):
    """The heat-balance building for the room ``env`` (window_beams_core.building_envelope format; env['glazing'] may
    list several windows) and parameters ``tp``."""
    lib = lib or HC.load_library()
    s = _archetype(tp, lib)
    warn = tp.setdefault("_warnings", [])
    dt = float(tp["dt_internal_min"]) * 60.0
    A_f, V = float(env["floor_m2"]), float(env["volume_m3"])
    ground = bool(s["floor"].get("ground_coupled")) if tp.get("floor_ground") is None else bool(tp["floor_ground"])
    u_wall = tp.get("u_wall") or s["u_wall"] or s["wall"]["U"]
    u_roof = tp.get("u_roof") or s["u_roof"] or s["roof"]["U"]
    u_floor = tp.get("u_floor") or s["u_floor"] or s["floor"]["U"]
    u_glz = tp.get("u_glazing") or s["u_glazing"] or s["glazing"]["U"]
    shgc = tp.get("shgc") or s["shgc"] or s["glazing"]["SHGC"]
    fexp = tp["floor_exposed_fraction"] if tp.get("floor_exposed_fraction") is not None else float(s["fexp"])
    st_w = _stack(lib, s["wall"], "wall", u_wall, warn=warn)
    st_r = _stack(lib, s["roof"], "roof", u_roof, warn=warn)
    sx = HC.stack(lib, s["floor"], "floor", u_floor, ground=ground)
    sc = HC.stack(lib, s["floor"], "floor", u_floor, ground=ground, covered=True)
    if not ground:
        sx = _match_u(sx, "floor", u_floor, False, warn, "floor")
        sc = _match_u(sc, "floor", u_floor, False, warn, "floor")
    if tp.get("kappa_wall") is not None:
        st_w = _match_u(_set_kappa(st_w, tp["kappa_wall"]), "wall", u_wall, False, warn, "wall")
    if tp.get("kappa_roof") is not None:
        st_r = _match_u(_set_kappa(st_r, tp["kappa_roof"]), "roof", u_roof, False, warn, "roof")
    if tp.get("kappa_floor") is not None:
        sx = _match_u(_set_kappa(sx, tp["kappa_floor"]), "floor", u_floor, ground, warn, "floor")
        sc = _match_u(_set_kappa(sc, tp["kappa_floor"]), "floor", u_floor, ground, warn, "floor")
    if ground and tp["soil_depth_m"] != 0.5:
        for stck in (sx, sc):
            for l in stck:
                if l["name"] == "soil":
                    l["d"] = tp["soil_depth_m"]
                    l["R"] = l["d"] / l["k"]
            dg = [l for l in stck if l["name"] == "deep_ground"]
            r_now = 1.0 / HC.u_of(stck, "floor", True)
            if dg:
                dg[0]["R"] = max(0.001, dg[0]["R"] + (1.0 / u_floor - r_now))
    cm_target = tp.get("mass_class_cm") or s["cm"]
    A_part = tp["partition_m2_per_floor"] * A_f
    k_part = tp["kappa_partition"]
    els = []
    sun_walls = [w for w in env["walls"] if w.get("air_m2", 0) > 1e-6]
    if tp["merge_walls"] and len(sun_walls) > 1:
        At = sum(float(w["air_m2"]) for w in sun_walls)
        subs = [(float(w["azimuth_deg"]), 90.0, tp["alpha_wall"], float(w["air_m2"])) for w in sun_walls]
        els.append(_make_el("walls", At, "sun", "wall", st_w, tp, dt, az=0.0, tilt=90.0, alpha=tp["alpha_wall"],
                            subs=subs))
    else:
        for w in sun_walls:
            els.append(_make_el("wall_" + w["facing"], w["air_m2"], "sun", "wall", st_w, tp, dt,
                                az=float(w["azimuth_deg"]), tilt=90.0, alpha=tp["alpha_wall"]))
    earth = [w for w in env["walls"] if w.get("earth_m2", 0) > 1e-6]
    if earth:
        Ae = sum(float(w["earth_m2"]) for w in earth)
        els.append(_make_el("walls_earth", Ae, "earth", "wall",
                            [{"name": "deep_ground", "d": 0.0, "k": None, "rho": 0, "c": 0, "R": 1.0}] + st_w, tp, dt))
    A_roof = float(env.get("roof_m2") or 0.0)
    if A_roof > 0:
        els.append(_make_el("roof", A_roof, "sun", "ceiling", st_r, tp, dt, az=0.0, tilt=0.0, alpha=tp["alpha_roof"]))
    fk = "earth" if ground else "air"
    if cm_target:
        base = sum(e.A * e.kap for e in els)
        kx = HC.kappa_inside(sx) or 1.0
        kc = HC.kappa_inside(sc)
        fexp = min(1.0, max(0.0, (cm_target * A_f - base - kc * A_f) / max((kx - kc) * A_f, 1e-6)))
        extra = cm_target * A_f - base - (fexp * kx + (1 - fexp) * kc) * A_f
        if extra > 1.0:
            A_part += extra / (2.0 * 155.0)
            k_part = 155.0
    if A_f * fexp > 1e-6:
        els.append(_make_el("floor_x", A_f * fexp, fk, "floor", sx, tp, dt, tilt=180.0))
    if A_f * (1.0 - fexp) > 1e-6:
        els.append(_make_el("floor_c", A_f * (1.0 - fexp), fk, "floor", sc, tp, dt, tilt=180.0))
    if A_part > 1e-6:
        g = _GENERIC["masonry" if k_part > 30 else "board"]
        d = max(0.002, k_part * 1000.0 / (g["rho"] * g["c"]))
        els.append(_make_el("partitions", 2.0 * A_part, "adiabatic", "wall",
                            [{"name": "part", "d": d, "k": g["k"], "rho": g["rho"], "c": g["c"], "R": d / g["k"]}],
                            tp, dt))
    glz = [gg for gg in (env.get("glazing") or []) if float(gg.get("area_m2", 0)) > 1e-6]
    glass_idx = []
    for k, gg in enumerate(glz):
        ug = float(gg.get("u") or gg.get("u_glazing") or u_glz)
        sg = float(gg.get("shgc") or shgc)
        el = _make_el("glass" if len(glz) == 1 else "glass_%d" % (k + 1), float(gg["area_m2"]), "glass", "glass",
                      None, tp, dt, az=float(gg.get("azimuth_deg", 0.0)), tilt=float(gg.get("tilt_deg", 90.0)),
                      u_glass=ug)
        el.glass_k = k
        el.shgc = sg
        glass_idx.append(len(els))
        els.append(el)
    A_door = float(env.get("doors_m2") or 0.0)
    if A_door > 0:
        els.append(_make_el("doors", A_door, "glass", "glass", None, tp, dt, u_glass=2.0))
    F = [1.0] * len(els)
    for _ in range(30):
        S = sum(e.A * f for e, f in zip(els, F))
        F = [1.0 / max(1e-6, 1.0 - e.A * f / S) for e, f in zip(els, F)]
    for e, f in zip(els, F):
        e.Gr = tp["h_r"] * e.A * f
        e.T = None
    n_eff = tp["n_rest"] + tp["n_system"] * (1.0 - tp["heat_recovery"])
    H_ve = tp["rho_cp_air_wh"] * V * n_eff
    a_in = {"wall": tp["alpha_in_wall"], "ceiling": tp["alpha_in_ceiling"], "floor": tp["alpha_landing"]}
    w_dif = []
    for e in els:
        if e.cls == "glass":
            w_dif.append(e.A * ((e.shgc or shgc) if tp["diffuse_loss"] else 0.0) if e.glass_k is not None else 0.0)
        else:
            w_dif.append(e.A * a_in[e.cls])
    sd = sum(w_dif)
    opq = [e for e in els if e.cls != "glass"]
    A_opq = sum(e.A for e in opq)
    A_glz = sum(els[j].A for j in glass_idx)
    win_az = float(glz[0].get("azimuth_deg", 0.0)) if glz else 0.0
    optics = [simple_glazing(els[j].u_glass, els[j].shgc) for j in glass_idx]
    bld = {"els": els, "A_f": A_f, "V": V, "H_ve": H_ve, "n_eff": n_eff, "dt": dt, "tp": tp, "shgc": shgc,
           "u_glazing": u_glz, "g_eff": shgc * tp["r_dirt"], "A_glz": A_glz, "ground": ground, "fexp": fexp,
           "C_air": tp["rho_cp_air_wh"] * 3600.0 * V * tp["air_cap_mult"], "a_in": a_in,
           "w_dif": [w / sd for w in w_dif], "w_rad": [(e.A / A_opq if e.cls != "glass" else 0.0) for e in els],
           "win_az": win_az, "glass_idx": glass_idx, "optics": optics, "n_windows": len(glass_idx),
           "stacks": {"wall": st_w, "roof": st_r, "floor_x": sx, "floor_c": sc},
           "ids": {"wall": s["wall"]["id"], "roof": s["roof"]["id"], "floor": s["floor"]["id"],
                   "glazing": s["glazing"]["id"]}}
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
    bld["_env"] = env
    bld["C_tot_kj_m2"] = sum(sum(e.ch["C"]) * e.A for e in els) / 1000.0 / A_f
    return bld


# ============================================================================================ drives
def plane_series(wx, azimuth_deg, tilt_deg, rho):
    """W/m2 on a plane per step; uses a pre-filled wx[('_plane', az, tilt, rho)] (e.g. EnergyPlus's own incident
    solar) when present - the same cache key as model A and R1."""
    key = ("_plane", round(azimuth_deg, 3), round(tilt_deg, 3), round(rho, 3))
    if key not in wx:
        if CORE_DIR not in sys.path:
            sys.path.insert(0, CORE_DIR)
        import building_energy as BE
        pi = BE.plane_irradiance
        alt, az, dni, dhi, ghi = wx["alt"], wx["az"], wx["dni"], wx["dhi"], wx["ghi"]
        wx[key] = array.array("d", (pi(alt[i], az[i], dni[i], dhi[i], ghi[i], azimuth_deg, tilt_deg, rho)
                                    for i in range(len(dni))))
    return wx[key]


_MID = [15.5, 45.0, 74.5, 105.0, 135.5, 166.0, 196.5, 227.5, 258.0, 288.5, 319.0, 349.5]


def earth_series(wx, t_earth12):
    """Earth temperature per step (monthly values interpolated between month middles; same values as R1, computed
    once per day)."""
    n = len(wx["temp"])
    spd = 86400.0 / wx["step_s"]
    per_day = []
    for dd in range(365):
        d = dd + 0.5
        j = 0
        while j < 12 and _MID[j] <= d:
            j += 1
        m0, m1 = (j - 1) % 12, j % 12
        d0 = _MID[m0] if j > 0 else _MID[11] - 365.0
        d1 = _MID[m1] if j < 12 else _MID[0] + 365.0
        f = (d - d0) / (d1 - d0)
        per_day.append(t_earth12[m0] * (1 - f) + t_earth12[m1] * f)
    return array.array("d", (per_day[int((i // spd) % 365)] for i in range(n)))


def internal_series(wx, tp, A_f):
    n = len(wx["temp"])
    prof = tp["internal_profile"]
    P, W = tp["occupants"], tp["w_per_person"]
    if prof == "workbook":
        if CORE_DIR not in sys.path:
            sys.path.insert(0, CORE_DIR)
        import building_energy as BE
        days = BE.seasonal_values(tp["heating_days"]) or [0.0] * 12
        q = [P * W * min(1.0, d / BE.DAYS[m]) for m, d in enumerate(days)]
        return array.array("d", (q[wx["month"][i] - 1] for i in range(n)))
    app = tp["appliance_w_m2"] * A_f
    if prof == "flat":
        return array.array("d", [P * W + app]) * n
    hd = tp["home_day_fraction"]
    occ = [0.8] * 7 + [1.0, 0.5] + [hd] * 8 + [1.0] * 6 + [0.8]
    shape = [0.6] * 24
    for h in (6, 7):
        shape[h] = 1.0
    for h in range(17, 22):
        shape[h] = 2.2
    sm = sum(shape) / 24.0
    hourly = [P * W * occ[h] + app * shape[h] / sm for h in range(24)]
    tm = wx["t_mid"]
    return array.array("d", (hourly[int(tm[i]) % 24] for i in range(n)))


def _as_list(x):
    return list(x) if isinstance(x, (list, tuple)) and (not x or not isinstance(x[0], (int, float))) else [x]


def window_gains(bld, beam_cls, diffuse, alpha_landing=None, glass_abs=None):
    """The windows' solar gains from the light engine's landing data (same call as R1 / model A):
    ``beam_cls`` = {class: array of beam W landing on that class ('floor', 'wall', 'ceiling', 'reveal'), BEFORE the
    glass}, ``diffuse`` = array of sky + ground W before the glass.  For several windows pass LISTS (one dict / array
    per window, in env['glazing'] order).  ``glass_abs`` (W absorbed in the glass, per window, e.g. EnergyPlus's) is
    put on that window's glass node; then ``diffuse`` must hold only the TRANSMITTED diffuse (divided by g, as R1).
    glass_model 'node': the before-glass sun is split by the glass's own tau / alpha (E+ simple glazing)."""
    tp = bld["tp"]
    bl = _as_list(beam_cls) if isinstance(beam_cls, (list, tuple)) else [beam_cls]
    dl = _as_list(diffuse) if (isinstance(diffuse, (list, tuple)) and isinstance(beam_cls, (list, tuple))) \
        else [diffuse]
    if glass_abs is None:
        gl = [None] * len(bl)
    elif isinstance(glass_abs, (list, tuple)) and isinstance(beam_cls, (list, tuple)):
        gl = list(glass_abs)
    else:
        gl = [glass_abs]
    nwin = max(1, bld["n_windows"])
    n = len(dl[0])
    rec = {"n": n, "diffuse": array.array("d", [0.0]) * n, "glass_abs": {}}
    node = tp["glass_model"] == "node"
    for k, (bc, dif) in enumerate(zip(bl, dl)):
        kk = min(k, nwin - 1)
        if node and bld["optics"]:
            op = bld["optics"][kk]
            gt = op["tau"] * tp["r_dirt"]
            ga = op["alpha"] * tp["r_dirt"]
        else:
            gt = bld["g_eff"] if bld["glass_idx"] == [] else \
                ((bld["els"][bld["glass_idx"][kk]].shgc or bld["shgc"]) * tp["r_dirt"])
            ga = 0.0
        d = rec["diffuse"]
        for i in range(n):
            d[i] += gt * dif[i]
        absb = array.array("d", [0.0]) * n if ga > 0 else None
        if absb is not None:
            for i in range(n):
                absb[i] = ga * dif[i]
        for c in ("floor", "wall", "ceiling", "reveal"):
            if c in bc:
                arr = bc[c]
                tgt = rec.get(c)
                if tgt is None:
                    tgt = rec[c] = array.array("d", [0.0]) * n
                for i in range(n):
                    tgt[i] += gt * arr[i]
                if absb is not None:
                    for i in range(n):
                        absb[i] += ga * arr[i]
        if gl[k] is not None:
            ab = gl[k]
            if absb is None:
                absb = array.array("d", ab)
            else:
                for i in range(n):
                    absb[i] += ab[i]
        if absb is not None:
            prev = rec["glass_abs"].get(kk)
            if prev is None:
                rec["glass_abs"][kk] = absb
            else:
                for i in range(n):
                    prev[i] += absb[i]
    if alpha_landing is not None:
        rec["alpha_landing"] = alpha_landing
    if not tp["beam_to_mass"]:
        for c in ("floor", "wall", "ceiling", "reveal"):
            if c in rec:
                arr = rec.pop(c)
                d = rec["diffuse"]
                for i in range(n):
                    d[i] += arr[i]
    return rec


# ============================================================================================ the integrator
def _block(arr, k, nb, n):
    if k == 1:
        return list(arr)
    out = []
    for b in range(nb):
        s = b * k
        e = min(n, s + k)
        out.append(sum(arr[s:e]) / float(e - s))
    return out


def _solair(e, I, te, tsky, h_co, h_ro):
    hs = h_co + h_ro
    fs = e.fsky
    return [(e.alpha * I[b] + h_co * te[b] + h_ro * (fs * tsky[b] + (1 - fs) * te[b])) / hs for b in range(len(te))]


def simulate(bld, wx, climate, win=None, band=None, keep=("air", "s", "m", "op"), T0=20.0, spinup=None,
             phi_int=None, stats=None):
    tp = bld["tp"]
    els = bld["els"]
    n = len(wx["temp"])
    step_s = float(wx["step_s"])
    k = max(1, int(round(bld["dt"] / step_s)))
    dt = k * step_s
    if abs(dt - bld["dt"]) > 1e-6:
        tp2 = dict(tp)
        tp2["dt_internal_min"] = dt / 60.0
        bld = assemble(bld["_env"], tp2)
        els = bld["els"]
    nb = (n + k - 1) // k
    blk = lambda a: _block(a, k, nb, n)
    te = blk(wx["temp"])
    tsky_src = wx.get("t_sky")
    tsky = blk(tsky_src) if tsky_src is not None else [t - tp["sky_dtheta_k"] for t in te]
    tg = blk(earth_series(wx, climate["t_earth_c"]))
    qint = blk(phi_int if phi_int is not None else internal_series(wx, tp, bld["A_f"]))
    eta = tp["eta"]
    frad = tp["int_radiant_fraction"]
    h_co = tp["h_out_conv"]
    h_ro = 4.0 * SIGMA * tp["eps_out"] * (283.15 ** 3)
    ne = len(els)
    # ---- outside driving temperature per element (W/K-weighted: G_out * drive enters node 0)
    drive = []
    for e in els:
        if e.kind == "sun":
            if e.subs:
                At = sum(s[3] for s in e.subs)
                acc = [0.0] * nb
                for (az, tilt, al, Ai) in e.subs:
                    I = blk(plane_series(wx, az, tilt, tp["ground_reflectance"]))
                    sa = _solair(e, I, te, tsky, h_co, h_ro)
                    w = Ai / At
                    for b in range(nb):
                        acc[b] += w * sa[b]
                drive.append(acc)
            else:
                I = blk(plane_series(wx, e.az, e.tilt, tp["ground_reflectance"]))
                drive.append(_solair(e, I, te, tsky, h_co, h_ro))
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
    # ---- sources on the inside surfaces (W) per block
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
            ro_per = e.r_o / e.A                                # K per W (film of the whole pane)
            for b in range(nb):
                half = 0.5 * eta * gab[b]
                sj[b] += half                                   # inside face
                # outer face: source S_o behind the outside film r_o (Thevenin) -> drive += S_o * r_o / A
                dv[b] += half * ro_per
    # ---- purge mask per block
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
    if band is None and tp["loads"] == "ideal":
        band = (tp["heat_setpoint_c"], tp["cool_setpoint_c"])
    lo_b = hi_b = None
    if band is not None:
        lo_b = blk(band[0]) if hasattr(band[0], "__len__") else [band[0]] * nb
        hi_b = blk(band[1]) if hasattr(band[1], "__len__") else [band[1]] * nb
    tarp = tp["convection"] == "tarp"
    hmult = tp["h_c_mult"]
    hfix = {"wall": tp["h_c_wall"], "floor": tp["h_c_floor"], "ceiling": tp["h_c_ceiling"], "glass": tp["h_c_wall"]}
    Ca = bld["C_air"] / dt
    H_ve0 = bld["H_ve"]
    A_all = sum(e.A for e in els)
    wA = [e.A / A_all for e in els]
    mk = [e.A * e.kap for e in els]
    smk = sum(mk) or 1.0
    wM = [m / smk for m in mk]
    p_above, p_delta = tp["purge_above_c"], tp["purge_delta_c"]
    # ---- flat per-element records (one list each, slots below) - the inner loop touches nothing else
    # 0 N, 1 cdt, 2 a, 3 inv, 4 ee, 5 Bconst, 6 G_out, 7 T, 8 dps, 9 drive, 10 src, 11 Gr, 12 A*h_pos, 13 A*h_neg,
    # 14 A*h_min, 15 fixed Gc, 16 tsi, 17 R, 18 B, 19 Gc, 20 wA, 21 wM
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
    out_air = [0.0] * nb
    out_s = [0.0] * nb
    out_m = [0.0] * nb
    out_hc = [0.0] * nb
    surf_out = [[0.0] * nb for _ in range(ne)] if "surf" in keep else None
    st = {"ta": float(T0), "tstar": float(T0)}
    cube = 1.0 / 3.0

    def steady_init():
        """Linear steady-state start: every chain at the conductance-weighted mix of its outside drive mean and a
        room guess (outdoor mean + 3 K)."""
        nb0 = min(nb, int(round(86400.0 * 7 / dt)))
        tmean = sum(te[nb - nb0:]) / nb0 if nb0 else T0
        troom = tmean + 3.0
        for j, e in enumerate(els):
            dv = drive[j]
            dmean = (sum(dv[nb - nb0:]) / nb0) if (dv is not None and nb0) else troom
            N = e.N
            if N == 1 or dv is None:
                D[j][7] = [troom] * N
                D[j][16] = troom
                continue
            # resistances along the chain (per element, K/W): outside link then node links
            rs = [1.0 / e.G_out if e.G_out > 0 else 0.0] + [1.0 / g for g in e.a[1:]]
            rin = 1.0 / (2.0 * e.A)                          # ~ inside film
            tot = sum(rs) + rin
            T = []
            acc = 0.0
            for i in range(N):
                acc += rs[i]
                T.append(dmean + (troom - dmean) * acc / tot)
            D[j][7] = T
            D[j][16] = T[-1]
        st["ta"] = troom
        st["tstar"] = troom

    def run_blocks(order, record):
        ta = st["ta"]
        tstar = st["tstar"]
        for b in order:
            teb = te[b]
            Hve = H_ve0
            if pmask is not None and pmask[b] and ta > p_above and ta > teb + p_delta:
                Hve = H_ve0 + dH
            a11 = Ca + Hve
            a12 = 0.0
            a22 = 0.0
            r1 = Ca * ta + Hve * teb + conv_int[b]
            r2 = 0.0
            for d in D:
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
            det = a11 * a22 - a12 * a12
            ta_n = (r1 * a22 - a12 * r2) / det
            ts_n = (a11 * r2 - a12 * r1) / det
            phc = 0.0
            if lo_b is not None:
                lo, hi = lo_b[b], hi_b[b]
                if ta_n < lo or ta_n > hi:
                    tset = lo if ta_n < lo else hi
                    s_air = a22 / det
                    phc = (tset - ta_n) / s_air
                    ta_n = tset
                    ts_n = ts_n + (-a12 / det) * phc
            ta, tstar = ta_n, ts_n
            sm_s = 0.0
            sm_m = 0.0
            for d in D:
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
            if record:
                out_air[b] = ta
                out_s[b] = sm_s
                out_m[b] = sm_m
                out_hc[b] = phc
                if surf_out is not None:
                    for j in range(ne):
                        surf_out[j][b] = D[j][16]
        st["ta"] = ta
        st["tstar"] = tstar

    if tp.get("codegen", True) and surf_out is None:
        gen = _gen_runner(D, tarp, Ca, H_ve0, pmask is not None, lo_b is not None)
        loop_blocks = run_blocks

        def run_blocks(order, record):                  # noqa: F811 - the generated straight-line version
            S = [st["ta"], st["tstar"]] + [x for d in D for x in d[7]]
            gen(order, record, S, te, conv_int, out_air, out_s, out_m, out_hc, pmask, dH, p_above, p_delta,
                lo_b, hi_b, *([d[9] for d in D] + [d[10] for d in D]))
            st["ta"], st["tstar"] = S[0], S[1]
            p = 2
            for d in D:
                N = d[0]
                d[7] = S[p:p + N]
                d[16] = d[7][-1]
                p += N
        run_blocks.loop = loop_blocks

    passes_done = 0
    mode = tp["warmup"]
    if mode == "year":
        for j, e in enumerate(els):
            D[j][7] = [float(T0)] * e.N
        sp = int(tp["spinup_years"] if spinup is None else spinup)
        for _ in range(sp):
            run_blocks(range(nb), False)
            passes_done += 1
    else:
        steady_init()
        bpd = int(round(86400.0 / dt))
        wb = min(nb, max(1, int(tp["warmup_days"]) * bpd))
        order = list(range(nb - wb, nb))
        prev = None
        for p in range(int(tp["warmup_max"])):
            run_blocks(order, False)
            passes_done += 1
            cur = [x for d in D for x in d[7]] + [st["ta"]]
            if prev is not None:
                dmax = max(abs(x - y) for x, y in zip(cur, prev))
                if dmax < tp["warmup_tol_k"] and p >= 1:
                    break
            prev = cur
    run_blocks(range(nb), True)
    if stats is not None:
        stats["warmup_passes"] = passes_done
        stats["nodes"] = sum(e.N for e in els)

    res = {}

    def expand(arr):
        if k == 1:
            return array.array("d", arr)
        o = array.array("d")
        prev = arr[-1]
        for b in range(nb):
            cur = arr[b]
            for q in range(k):
                if b * k + q >= n:
                    break
                f = (q + 1.0) / k
                o.append(prev + (cur - prev) * f)
            prev = cur
        return o
    if "air" in keep:
        res["air"] = expand(out_air)
    if "s" in keep:
        res["s"] = expand(out_s)
    if "m" in keep:
        res["m"] = expand(out_m)
    if "op" in keep:
        res["op"] = expand([0.5 * (x + y) for x, y in zip(out_air, out_s)])
    if band is not None:
        res["hc"] = expand(out_hc)
    res["T_m_end"] = out_m[-1]
    if surf_out is not None:
        res["surf"] = dict((els[j].name, expand(surf_out[j])) for j in range(ne))
    return res


_GEN_CACHE = {}


def _gen_runner(D, tarp, Ca, H_ve0, has_purge, has_band):
    """Straight-line Python for one building's step (every node a local variable, every coefficient a literal):
    the same arithmetic as simulate.run_blocks, ~3x faster.  Cached by the generated source."""
    f = repr
    L = []
    ne = len(D)
    names = [["t%d_%d" % (j, i) for i in range(D[j][0])] for j in range(ne)]
    args = ", ".join(["DV%d" % j for j in range(ne)] + ["SR%d" % j for j in range(ne)])
    L.append("def _run(order, record, S, te, ci, oa, osx, om, oh, pm, dH, p_above, p_delta, lo_b, hi_b, %s):" % args)
    L.append("    ta = S[0]; ts = S[1]")
    p = 2
    for j in range(ne):
        for nm in names[j]:
            L.append("    %s = S[%d]" % (nm, p))
            p += 1
    L.append("    C3 = 1.0 / 3.0")
    L.append("    for b in order:")
    I = "        "
    L.append(I + "teb = te[b]")
    L.append(I + "Hve = %s" % f(float(H_ve0)))
    if has_purge:
        L.append(I + "if pm[b] and ta > p_above and ta > teb + p_delta: Hve = %s + dH" % f(float(H_ve0)))
    L.append(I + "a11 = %s + Hve; a12 = 0.0; a22 = 0.0; r1 = %s * ta + Hve * teb + ci[b]; r2 = 0.0" % (
        f(float(Ca)), f(float(Ca))))
    for j, d in enumerate(D):
        N, cdt, a, inv, ee, Bc, Gout, dv = d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[9]
        Gr = d[11]
        tn = names[j]
        if tarp:
            L.append(I + "dd = %s - ta" % tn[-1])
            L.append(I + "if dd > 0.0: g%d = %s * dd ** C3" % (j, f(d[12])))
            L.append(I + "else: g%d = %s * (-dd) ** C3" % (j, f(d[13])))
            L.append(I + "if g%d < %s: g%d = %s" % (j, f(d[14]), j, f(d[14])))
        else:
            L.append(I + "g%d = %s" % (j, f(d[15])))
        drv = (" + %s * DV%d[b]" % (f(Gout), j)) if (dv is not None and Gout != 0.0) else ""
        if N == 1:
            c0 = ("%s * %s + " % (f(cdt[0]), tn[0])) if cdt[0] != 0.0 else ""
            L.append(I + "R%d = %sSR%d[b]%s" % (j, c0, j, drv))
        else:
            L.append(I + "p%d_0 = (%s * %s%s) * %s" % (j, f(cdt[0]), tn[0], drv, f(inv[0])))
            for i in range(1, N - 1):
                L.append(I + "p%d_%d = (%s * %s + %s * p%d_%d) * %s" % (j, i, f(cdt[i]), tn[i], f(a[i]), j, i - 1,
                                                                    f(inv[i])))
            L.append(I + "R%d = %s * %s + SR%d[b] + %s * p%d_%d" % (j, f(cdt[N - 1]), tn[N - 1], j, f(a[N - 1]), j,
                                                                 N - 2))
        L.append(I + "ib%d = 1.0 / (%s + g%d)" % (j, f(Bc + Gr), j))
        L.append(I + "q = g%d * ib%d; r = %s * ib%d" % (j, j, f(Gr), j))
        L.append(I + "a11 += g%d - g%d * q; a12 -= %s * q; a22 += %s - %s * r; r1 += R%d * q; r2 += R%d * r" % (
            j, j, f(Gr), f(Gr), f(Gr), j, j))
    L.append(I + "det = a11 * a22 - a12 * a12")
    L.append(I + "ta_n = (r1 * a22 - a12 * r2) / det")
    L.append(I + "ts_n = (a11 * r2 - a12 * r1) / det")
    L.append(I + "phc = 0.0")
    if has_band:
        L.append(I + "lo = lo_b[b]; hi = hi_b[b]")
        L.append(I + "if ta_n < lo or ta_n > hi:")
        L.append(I + "    tset = lo if ta_n < lo else hi")
        L.append(I + "    phc = (tset - ta_n) / (a22 / det)")
        L.append(I + "    ta_n = tset")
        L.append(I + "    ts_n = ts_n + (-a12 / det) * phc")
    L.append(I + "ta = ta_n; ts = ts_n")
    for j, d in enumerate(D):
        N, ee, Gr = d[0], d[4], d[11]
        tn = names[j]
        L.append(I + "%s = (R%d + g%d * ta + %s * ts) * ib%d" % (tn[-1], j, j, f(Gr), j))
        for i in range(N - 2, -1, -1):
            L.append(I + "%s = p%d_%d + %s * %s" % (tn[i], j, i, f(ee[i]), tn[i + 1]))
    L.append(I + "if record:")
    L.append(I + "    oa[b] = ta; oh[b] = phc")
    L.append(I + "    osx[b] = " + " + ".join("%s * %s" % (f(D[j][20]), names[j][-1]) for j in range(ne)))
    L.append(I + "    om[b] = " + " + ".join("%s * %s" % (f(D[j][21]), names[j][-1]) for j in range(ne)))
    L.append("    S[0] = ta; S[1] = ts")
    p = 2
    for j in range(ne):
        for nm in names[j]:
            L.append("    S[%d] = %s" % (p, nm))
            p += 1
    src = "\n".join(L) + "\n"
    fn = _GEN_CACHE.get(src)
    if fn is None:
        ns = {}
        exec(compile(src, "<hb2_step>", "exec"), ns)
        fn = ns["_run"]
        if len(_GEN_CACHE) > 64:
            _GEN_CACHE.clear()
        _GEN_CACHE[src] = fn
    return fn


def run_year(bld, wx, climate, win=None, band=None, keep=("air", "s", "m", "op"), T0=20.0):
    """Free-running (or ideal-load) year; same call as R1 / model A.  ``win`` = window_gains(...)."""
    return simulate(bld, wx, climate, win=win, band=band, keep=keep, T0=T0)


def summary(bld):
    els = bld["els"]
    return ("B2 %s/%s/%s/%s  C/A_f %.0f kJ/m2K  nodes %d  H_ve %.0f W/K  fexp %.2f  windows %d  %s" % (
        bld["ids"]["wall"], bld["ids"]["roof"], bld["ids"]["floor"], bld["ids"]["glazing"], bld["C_tot_kj_m2"],
        sum(e.N for e in els), bld["H_ve"], bld["fexp"], bld["n_windows"],
        " ".join("%s:%d" % (e.name, e.N) for e in els)))


class ThermalModel(object):
    """The THERMAL_PARAMS interface as an object (entry.thermal_model() returns one)."""
    name = VERSION
    PARAMS = PARAMS

    def make_tparams(self, over=None, building=None):
        return make_tparams(over, building)

    def assemble(self, env, tp, lib=None):
        return assemble(env, tp, lib)

    def window_gains(self, bld, beam_cls, diffuse, alpha_landing=None, glass_abs=None):
        return window_gains(bld, beam_cls, diffuse, alpha_landing, glass_abs)

    def run_year(self, bld, wx, climate, win=None, band=None, keep=("air", "s", "m", "op"), T0=20.0):
        return run_year(bld, wx, climate, win, band, keep, T0)

    def internal_series(self, wx, tp, A_f):
        return internal_series(wx, tp, A_f)

    def param_table(self, advertised_only=True):
        return param_table(advertised_only)

    simple_glazing = staticmethod(simple_glazing)


if __name__ == "__main__":
    print(write_params_md())
