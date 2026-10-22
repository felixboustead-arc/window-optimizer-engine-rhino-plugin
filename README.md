# Window Optimizer engine (R3)

Planning engine for an articulating window: a ring that rotates and an aperture that closes. For every minute of the
year it chooses a window pose (ring angle and aperture closing) that keeps a crisp, steady **beam of sunlight** in the
room. It can also trade that beam against **thermal comfort**. It runs from the command line or live inside
**Rhino 8 / Grasshopper**.

This is the engine as accepted on 21 Oct 2026. It contains only the parts that passed the project's independent
test and judge gate. Nothing experimental is switched on by default.

## What is inside

| Part | What it does | Status |
|---|---|---|
| **LIGHT**: G2R3T band glide (`lab/entries/R2/B/glide_champion.py`) | Beam-first pose plan per window class, real weather, DNI gate 120 W/m² | Accepted (R2 champion) |
| **THERMAL**: HB-L2 (`lab/entries/R2/B/thermal_hb2.py`) | Room heat balance with layered constructions (`lab/research/constructions.json`) | Accepted |
| **COMBINED**: nested knapsack (`lab/entries/R2/B/entry.py`) | Slider 0 (comfort) to 1 (light). Converts beam days into comfort where it pays most. The sets stay nested. | Accepted |
| **Engine** (`lab/engine/engine.py`) | One call for 1-15 windows, per-minute pose lookup, day-ahead planning, cache control, CLI | Accepted final |
| **O13** buildings (`lab/engine/o13.py`) | Optional multi-zone comfort room: partitions, doors, vents, cross-ventilation, MVHR, ceiling-fan comfort band | Accepted final, **opt-in** (off by default; off is bit-identical to the plain engine) |
| **Grasshopper component** (`lab/engine/gh_engine_component.py`) | Plans in a background CPython process and outputs ring / aperture for the clock minute on every solve | Accepted |

**Not included:** the LUMEN-P4b smoothing option. Its `light_p4b` switch exists in `engine.py` but its module is not
shipped, because its multi-window and full-year checks were still running at release. Leave the switch off. Also not
included: the `blade` variant (an optional code path, unused by the engine), the validation tooling (E+ reference
builders, scorers) and the lab record.

## Validation summary

- **Engine identity:**
  - Plans equal the frozen champion's, bit for bit.
  - Day-ahead planning (days d-1, d, d+1) gives the same day d as a full-year plan in 20 of 20 test cases (N/E/S, tilts, winter/summer, year wrap).
- **O13 against EnergyPlus**, using daily-mean operative temperature with the X7 test (RMSE ≤ 1.0 K, |bias| ≤ 0.5 K):
  - Public cases: 51 of 51 pass, against 40 of 51 without O13.
  - Held-out Nobbys weather: 7 of 7 pass, against 3 of 7. This is a reproduction on a case seen once, not a blind test.
  - Light and heavy classes for vents and MVHR: all pass.
- **O13 limits:** they are stated in the docstring of `o13.py`, which you should read before using O13 for design decisions:
  - Cross-ventilation is implementation-checked against the E+ AirflowNetwork model, not against a real house. Its error is 0.05-0.10 K.
  - Door discharge coefficient is 0.62, against about 0.73 in the reference.
  - Light-class vent cases over-cool the hottest-day peak by 0.8-1.8 K.
  - The ceiling fan is **unvalidated**: a comfort-band offset only, which assumes occupants control the fan.
- **Plan time** (CPython, 2 days, idle PC, cold cache):

  | Layout | Time |
  |---|---|
  | 1 window | 67 s |
  | 3 windows | 159 s |
  | 12 windows | 394 s |

  The per-minute pose lookup takes microseconds.

## Requirements

- **Python 3.9+**, standard library only. Inside Rhino 8, the GH component uses Rhino's CPython for the lookup and an
  external CPython for background planning.
- **Weather (not shipped):** the Williamtown (Newcastle NSW) TMYx EPW from climate.onebuilding.org,
  `AUS_NSW_Newcastle.AP-RAAF.Williamtown.947760_TMYx.2011-2025.epw`. Put it in `data/weather/` or point `WO_EPW` at it.
  `SYNTHETIC_TEST_ONLY.epw` is only a placeholder that the scan job template loads. It is not real weather.
- **COMBINED / THERMAL only:** the EnergyPlus reference case for the default comfort room (`A04_felix_open`:
  `case.json` + `ref.csv.gz`). It is not shipped. Put it in `data/eplus_ref/` or point `WO_EPLUS_REF` at the folder.
  The LIGHT mode does not need it.

## Use

Command line, light plan for the windows in a scene file:

```
python -B lab/engine/engine.py plan scene.json plans.json --days 171,354
python -B lab/engine/engine.py plan scene.json plans.json --slider 0.5      # COMBINED (needs the E+ reference)
```

A scene is `{"windows": [{"id": "w00", "facing": 0.0, "tilt_deg": 0.0, "level_m": 0.0, "scale": 1.0,
"width_f": 1.0, "beam_w": 0.10}]}`. `facing` is in degrees from north, `beam_w` in metres. Plans follow contract v1:
`{"runs": [[step, ring_deg, aperture_pct], ...], "n_steps": 525600, "step_min": 1}`.

From Python:

```python
import sys; sys.path.insert(0, "lab/engine")
import engine as EN
plans = EN.plan(scene, None, {"days": [171]}, 1.0)            # light
pose = EN.Schedule(plans).poses(EN.minute_of_year(5, 1, 10, 0))  # 1 May 10:00 AEST
```

In Grasshopper (Rhino 8):
1. Paste `lab/engine/gh_engine_component.py` into a *Python 3 Script* component.
2. Set `WO_ENGINE_DIR` (or edit `ENGINE_DIR` at the top) to this repo's `lab/engine`.
3. Set `WO_PYTHON` to a CPython 3.9+ for background full-year plans.

The inputs and outputs are listed in the component's docstring. `rhino/Windownew.3dm` and `rhino/windowparam.gh` are
the window model and its parametric definition.

O13 (opt-in): pass `building=` to `engine.plan` in COMBINED mode. The format is in the `o13.py` docstring.

## Layout

The folder names (`lab/entries/R2/B`, `lab/models`) are kept from the research lab so that the frozen champion files
run unchanged. All paths are relative to the repository. Runtime caches go to `lab/plans/`, `lab/out/` and `out/_cache/`
(ignored by git).
