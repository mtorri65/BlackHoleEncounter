# Orbital Climate Model — User Manual

A 1-D (latitude-resolved) seasonal **energy-balance model (EBM)** of Earth's
climate under an arbitrary orbit, in the North / Budyko / Sellers lineage. It
computes daily-mean insolation from orbital elements, evolves zonal-mean surface
temperature with meridional heat transport and an ice-albedo feedback, and
supports "sudden orbital change" transient experiments and parallel parameter
sweeps.

The package pairs naturally with the REBOUND black-hole-flyby simulation in the
parent repository: REBOUND tells you the post-flyby orbital elements
`(a, e, λ_p, …)`; this model takes those as an instantaneous step and returns
the climate transient.

---

## 1. Installation

The model is pure Python. From the repository root install the dependencies:

```bash
python -m pip install numpy scipy pandas matplotlib pyyaml
python -m pip install pyarrow      # optional: enables parquet sweep output (CSV fallback otherwise)
python -m pip install pytest       # optional: to run the test suite
```

No build step is required. Run everything from the **repository root** (the
directory that contains the `orbital_climate/` folder) so the package is
importable and the CLI module path resolves.

---

## 2. Package layout

| Module | Responsibility |
|---|---|
| `config.py` | `Config` dataclass + `load_config()` (YAML + overrides). All tunable parameters live here. |
| `kepler.py` | Solves Kepler's equation; true anomaly and heliocentric radius. |
| `insolation.py` | Daily-mean top-of-atmosphere insolation vs. latitude and season; analytic annual mean. |
| `ebm.py` | The `EBM` class: grid, diffusion operator, coalbedo, semi-implicit stepper, spin-up. |
| `experiment.py` | `run_equilibrium`, `run_perturbation`, and output writers (CSV/PNG/summary). |
| `sweep.py` | `run_sweep`: parallel Cartesian-product parameter sweeps → parquet/CSV table. |
| `cli.py` | Command-line entry point with `insolation`, `ebm`, `sweep` subcommands. |
| `tests/` | Physics-validation and plumbing tests (`pytest`). |

---

## 3. Configuration

All parameters are fields of the `Config` dataclass. Values may be set three
ways, in increasing precedence:

1. **Defaults** baked into `Config` (tuned to present-day Earth, ~288 K).
2. **YAML file** — see `input_climate.yaml`. Unknown keys are ignored, so one
   file can carry parameters for every stage plus a `sweep:` block.
3. **Explicit overrides** — `Config(ecc=0.117)` in Python, or `--ecc 0.117` on
   the CLI. These win over the YAML file.

### Parameter reference

**Orbital / stellar forcing**

| Field | Default | Meaning |
|---|---|---|
| `S0` | 1361.0 | Solar constant at 1 AU [W/m²] |
| `a_au` | 1.0 | Semi-major axis [AU] |
| `ecc` | 0.0167 | Orbital eccentricity |
| `obliquity_deg` | 23.44 | Axial tilt ε [deg] |
| `lon_perihelion_deg` | 283.0 | Longitude of perihelion λ_p [deg] (~283° ≈ perihelion on Jan 3) |
| `days_per_year` | 365.256363 | Sidereal year [days] |

**Energy-balance physics**

| Field | Default | Meaning |
|---|---|---|
| `olr_A` | 203.3 | Outgoing-longwave intercept A [W/m²], `OLR = A + B·T` (T in °C) |
| `olr_B` | 2.09 | Outgoing-longwave slope B [W/m²/°C] |
| `diffusion_D` | 0.58 | Meridional heat-transport coefficient [W/m²/°C] |
| `coalbedo_a0` | 0.676 | Mean ice-free coalbedo (tuned to 288 K) |
| `coalbedo_a2` | -0.200 | Latitude structure of coalbedo (× P₂(x)) |
| `coalbedo_ice` | 0.38 | Coalbedo over ice (albedo 0.62) |
| `T_ice_degC` | -10.0 | Ice-formation temperature threshold [°C] |
| `heat_capacity` | 1.05e8 | Mixed-layer heat capacity C [J/m²/°C]; τ = C/B ≈ 1.6 yr (single-surface mode only) |
| `n_lat` | 180 | Number of latitude cells (equal-area in x = sin φ) |
| `diag_lat_deg` | 65.0 | Latitude tracked for seasonal peak/trough + peak insolation (may be negative) |

**Two-surface land/ocean mode** (North & Coakley 1979) — see §9

| Field | Default | Meaning |
|---|---|---|
| `two_surface` | `false` | Enable separate land and ocean temperatures per latitude |
| `heat_capacity_land` | 1.2e6 | C_land [J/m²/°C] (~1 m soil, τ ≈ 7 days) |
| `heat_capacity_ocean` | 2.1e8 | C_ocean [J/m²/°C] (50 m mixed layer, τ ≈ 3.2 yr) |
| `land_ocean_coupling` | 3.5 | ν [W/m²/°C] zonal land↔ocean exchange |
| `land_fraction_override` | `null` | Force a uniform land fraction instead of Earth's zonal profile |

**Time stepping**

| Field | Default | Meaning |
|---|---|---|
| `dt_days` | 2.0 | Timestep [days]; semi-implicit scheme is stable at large dt |
| `spinup_max_years` | 200 | Cap on spin-up length |
| `spinup_tol_degC` | 1e-4 | Year-over-year global-mean change that defines equilibrium |

---

## 4. Command-line interface

General form (run from the repo root):

```bash
python -m orbital_climate.cli <subcommand> [options]
```

All three subcommands accept the shared configuration options:
`--config <path.yaml>`, `--S0`, `--a-au`, `--ecc`, `--obliquity`,
`--lon-perihelion`, `--diag-lat`.

### 4.1 `insolation` — orbital forcing diagnostics

Reports the global-and-annual mean insolation (vs. its analytic value) and the
peak daily-mean insolation at the diagnostic latitude. Optionally writes a
latitude × season insolation map.

| Option | Default | Meaning |
|---|---|---|
| `--peak-lat` | `diag_lat_deg` | Latitude for the peak diagnostic [deg] |
| `--n-time` | 2000 | Orbital samples for time averages / peak search |
| `--plot <path>` | — | Write a latitude × season PNG |

### 4.2 `ebm` — equilibrium and sudden-perturbation transient

Spins the configuration up to equilibrium and prints diagnostics. If
`--perturb-ecc` is given, it then switches the orbit at t=0 and integrates the
transient, writing outputs to a timestamped folder.

| Option | Default | Meaning |
|---|---|---|
| `--perturb-ecc <e>` | — | Run a transient with this perturbed eccentricity |
| `--perturb-lon-perihelion <deg>` | — | Optional perturbed λ_p |
| `--years <n>` | 60 | Transient length [years] |
| `--out-dir <path>` | `climate_runs` | Parent directory for outputs |

**Transient outputs** (under `<out-dir>/<timestamp>/`):
- `transient.csv` — per-year diagnostics (global/NH/SH mean T, ice-edge latitude, `Tdiag_summer_max`, `Tdiag_winter_min`)
- `transient.png` — 4-panel figure (mean-T transient, diagnostic-latitude summer peak, ice-edge migration, temperature profiles)
- `summary.txt` — before/after summary

### 4.3 `sweep` — parallel parameter sweep

Runs the Cartesian product of a `sweep:` block in the YAML config, spinning each
combination up to equilibrium and collecting scalar diagnostics into one table.

| Option | Default | Meaning |
|---|---|---|
| `--workers <n>` | 1 | Worker processes (`>1` runs combos in parallel) |
| `--out-dir <path>` | `climate_sweeps` | Parent directory for the results table |

The `sweep:` block maps `Config` field names to lists of values, e.g.:

```yaml
sweep:
  ecc:         [0.0167, 0.06, 0.117]
  diffusion_D: [0.4, 0.58, 0.8]
  lon_perihelion_deg: [103.0, 283.0]
```

Output: `<out-dir>/<timestamp>/sweep_results.parquet` (or `.csv` without pyarrow),
one row per combination with the swept parameters plus `spinup_years`,
`global_mean_degC`, `global_mean_K`, `iceline_lat_nh`, `diag_lat_deg`,
`peak_diag_insol`.

---

## 5. Python API

The package re-exports the main entry points:

```python
from orbital_climate import (
    Config, load_config,
    EBM, run_equilibrium, run_perturbation,
    perturbed_config_from_scenario, save_perturbation_outputs,
    run_sweep,
)
```

- `run_equilibrium(config) -> EquilibriumResult` — spun-up profile + diagnostics.
- `run_perturbation(baseline, perturbed, n_years=60) -> PerturbationResult` — transient arrays.
- `perturbed_config_from_scenario(baseline, ecc, lon_perihelion_deg=None)` — copy a config with the orbit changed.
- `save_perturbation_outputs(result, out_dir)` — write CSV/PNG/summary.
- `run_sweep(base, ranges, out_dir, workers=1)` — parallel sweep → table path.

---

## 6. Worked examples

### Example 1 — Reproduce the 65 °N summer-insolation drop (CLI)

The headline Milankovitch result: injecting eccentricity (perihelion in 10 %,
aphelion out 10 % → e ≈ 0.117) with perihelion near Jan 3 cuts 65 °N summer
insolation by ~17 %.

```bash
# Present-day orbit
python -m orbital_climate.cli insolation --config input_climate.yaml
#   -> Peak daily-mean insolation at 65.0 deg ~ 478 W/m^2

# Perturbed orbit (pure eccentricity injection)
python -m orbital_climate.cli insolation --config input_climate.yaml --ecc 0.117
#   -> Peak daily-mean insolation at 65.0 deg ~ 400 W/m^2   (a 17% drop)
```

Add `--plot insol_perturbed.png` to save the latitude × season insolation map.

### Example 2 — Sudden orbital-change climate transient (CLI)

Spin up present-day Earth, then switch to e = 0.117 at t = 0 and track 80 years:

```bash
python -m orbital_climate.cli ebm \
    --config input_climate.yaml \
    --perturb-ecc 0.117 \
    --years 80 \
    --out-dir climate_runs
```

This prints the baseline equilibrium (global mean ≈ 288.15 K, NH ice edge, 65 N
peak insolation) and the transient endpoints, and writes `transient.csv`,
`transient.png`, and `summary.txt` under `climate_runs/<timestamp>/`.

> **Interpretation note.** With the default ocean-like uniform heat capacity
> (τ ≈ 1.6 yr) the seasonal temperature swing is strongly damped, so the
> transient shows a mild *annual-mean warming* (higher eccentricity raises the
> annual-mean insolation via Jensen's inequality) rather than the glacial
> inception described from the summer-insolation drop. Capturing inception would
> need a land/ocean heat-capacity contrast (small C over land → large summer
> response) and a snow/ablation treatment.

### Example 3 — Parameter sweep (YAML + CLI)

Add a `sweep:` block to a config file and run it across 5 workers:

```yaml
# my_sweep.yaml  (plus the usual physical parameters)
ecc: 0.0167
sweep:
  ecc:         [0.0167, 0.06, 0.117]
  diffusion_D: [0.4, 0.58, 0.8]
```

```bash
python -m orbital_climate.cli sweep --config my_sweep.yaml --workers 5 --out-dir climate_sweeps
```

Inspect the resulting table:

```python
import pandas as pd, glob
f = sorted(glob.glob("climate_sweeps/*/sweep_results.*"))[-1]
df = pd.read_parquet(f) if f.endswith(".parquet") else pd.read_csv(f)
print(df)
```

You should see the global mean rise with eccentricity and, at stronger
diffusion, the ice edge retreat poleward (more efficient heat transport).

### Example 4 — Feeding REBOUND-derived elements (Python API)

Take the post-flyby orbital elements from a REBOUND run and compute the climate
transient directly:

```python
from orbital_climate import Config, run_perturbation, save_perturbation_outputs

baseline  = Config()                       # present-day Earth
perturbed = Config(a_au=1.0017, ecc=0.117, # <- from the REBOUND black-hole flyby
                   lon_perihelion_deg=283.0)

result = run_perturbation(baseline, perturbed, n_years=100)
save_perturbation_outputs(result, "climate_runs/bh_flyby")

d = result.diagnostics
print(f"global mean: {d['global_mean'][0]:.2f} -> {d['global_mean'][-1]:.2f} degC")
print(f"NH ice edge: {d['iceline_lat_nh'][0]:.1f} -> {d['iceline_lat_nh'][-1]:.1f} deg")
```

---

## 7. Validation & tests

The physics is pinned by a `pytest` suite (run from the repo root):

```bash
python -m pytest orbital_climate/tests/ -v
```

Key validation targets exercised:
- **Analytic annual-mean insolation** `⟨S⟩ = S₀ / (4 a² √(1−e²))` — end-to-end Kepler + insolation check.
- **65 °N June-peak drop** ~480 → ~400 W/m² for e = 0.117, λ_p = 283°.
- **Present-day global mean = 288 K** (tuned coalbedo).
- **Local radiative equilibrium** with D = 0, no ice: `⟨T⟩ = (⟨Q⟩·a − A)/B` per latitude.
- **Global energy budget closes** at equilibrium; the diffusion operator conserves energy and annihilates a constant field.
- **Semi-implicit stability** from dt = 1 day to dt = 30 days.

---

## 8. Physics notes & assumptions

- **Grid.** Cell-centred in x = sin(latitude), so cells are equal-area and global
  means are simple cell means. The diffusion weight `(1 − x²)` vanishes at the
  poles, giving a natural no-flux boundary and exact global energy conservation.
- **Numerics.** Diffusion + linear OLR are treated implicitly (a constant
  tridiagonal solve, LU-factorised once); the nonlinear insolation×coalbedo
  source is explicit. The scheme is unconditionally stable for the stiff terms,
  so `dt_days` is an accuracy knob, not a stability limit.
- **"Sudden" perturbation.** The transient assumes the orbit changes
  instantaneously at t = 0. This is well justified when the orbital change is
  fast compared to the ocean mixed-layer time (τ ≈ 1.6 yr) — as it is for a
  black-hole flyby, which reconfigures the orbit over days to a few months.
- **Known simplifications.** Linear OLR rather than σT⁴ (calibrated near 288 K —
  treat results outside ~250–300 K as ordinal, not quantitative); hard
  ice-albedo step at `T_ice`; fixed obliquity within a run (no precession or
  obliquity cycling).

---

## 9. Two-surface land/ocean mode

By default the model carries **one temperature per latitude** with a single
ocean-like heat capacity, which damps the seasonal cycle everywhere. Setting
`two_surface: true` switches to the North & Coakley (1979) formulation:

```
C_l ∂T_l/∂t = Q·a(T_l) − (A + B·T_l) + D∇²T_l + ν(T_o − T_l)
C_o ∂T_o/∂t = Q·a(T_o) − (A + B·T_o) + D∇²T_o − ν(T_o − T_l)·f/(1−f)
```

Each latitude carries **separate land and ocean temperatures**, coupled by a
zonal exchange term. The `f/(1−f)` factor (with `f` the land fraction) makes the
exchange **energy-conserving**: once weighted by their areas, the flux leaving
the ocean column equals the flux entering the land column.

**Why a blended heat capacity does not work.** The obvious cheap alternative —
making `C` a latitude-dependent area-weighted average of land and ocean values —
changes the thermal time constant at 60 °N by **under 10%**, because the
arithmetic mean is dominated by the much larger ocean value (2.1e8 vs 1.2e6).
Averaging `C` is the wrong operation; the two surfaces need separate
temperatures for the 175× contrast in thermal inertia to survive.

**Land fraction** comes from Earth's real zonal profile (Antarctica in the far
south, near-pure ocean in the southern mid-latitudes, most continental in the
northern mid-latitudes), integrating to ~29% global land cover. Override it with
`land_fraction_override` for idealized experiments.

**What it changes.** At 65 °N, present-day Earth:

| | single surface | two surface |
|---|---|---|
| seasonal range (blended) | 11.5 K | **34.9 K** |
| seasonal range over land | — | **~40 K** |
| seasonal range over ocean | — | **~9 K** |
| global mean T | 288.14 K | 288.07 K |

The global mean is essentially unchanged — equilibrium is set by the radiation
balance, not by heat capacity — so **no re-tuning of `coalbedo_a0` is needed**.
Only the seasonal cycle changes, which is exactly the intent.

`ν = 3.5 W/m²/°C` was calibrated so the 65 °N seasonal range is ~40 K over land
and ~9 K over ocean, matching observed continental vs. maritime seasonality.

**Extra outputs.** In two-surface mode `EquilibriumResult` additionally carries
`diag_land_summer_max` / `diag_land_winter_min`, the ocean equivalents, and the
`diag_land_seasonal_range` / `diag_ocean_seasonal_range` properties. They are
`None` in single-surface mode.

```bash
python -m orbital_climate.cli ebm --config input_climate.yaml   # with two_surface: true
```
```
