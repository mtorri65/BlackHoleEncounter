# Orbital Perturbation Climate Model — Development Report

A complete record of the climate-modelling work: what was built, why each design
choice was made, how every claim was validated, and which of my initial
assertions turned out to be wrong.

This is the **narrative and scientific record**. For day-to-day usage see
[`orbital_climate/MANUAL.md`](orbital_climate/MANUAL.md).

**Scope.** A 1-D latitude-resolved seasonal energy-balance model (EBM) of
Earth's climate under an arbitrary orbit, plus the machinery to drive it from
the REBOUND black-hole-flyby sweeps in this repository.

---

## Table of contents

1. [Goal and overall architecture](#1-goal-and-overall-architecture)
2. [Part I — Orbital physics core](#2-part-i--orbital-physics-core)
3. [Part II — The energy-balance model](#3-part-ii--the-energy-balance-model)
4. [Part III — Two-surface land/ocean extension](#4-part-iii--two-surface-landocean-extension)
5. [Part IV — Bridging REBOUND to the climate model](#5-part-iv--bridging-rebound-to-the-climate-model)
6. [Part V — Results](#6-part-v--results)
7. [Corrections and negative results](#7-corrections-and-negative-results)
8. [Bugs found and fixed](#8-bugs-found-and-fixed)
9. [Known limitations](#9-known-limitations)
10. [File inventory and reproduction](#10-file-inventory-and-reproduction)

---

## 1. Goal and overall architecture

The starting point was `orbital-climate-model-context.md`, describing a scenario
where Earth's perihelion moves in 10% and aphelion out 10% — effectively a pure
eccentricity injection to `e ≈ 0.117`, roughly 7× today's value and about 2× the
Milankovitch maximum.

The eventual goal was broader: run a climate analysis on **every** simulation in
a REBOUND BH-flyby sweep. That yields a two-stage pipeline:

```
REBOUND sweep  ──►  Earth's post-flyby orbital elements  ──►  seasonal EBM
  (672 runs)         (a, e, obliquity, λ_p, year length)      (climate outcome)
```

The two halves are cleanly separable: REBOUND determines *what orbit Earth ends
up on*; the EBM determines *what climate that orbit produces*. Everything below
is either one of those halves or the bridge between them.

Work proceeded in four stages, each validated before the next was built:
physics core → EBM → REBOUND bridge → land/ocean extension.

---

## 2. Part I — Orbital physics core

Built first, deliberately, so that sign and geometry errors were caught before
anything sat on top of them.

### 2.1 Kepler solver — [`orbital_climate/kepler.py`](orbital_climate/kepler.py)

Newton iteration on Kepler's equation:

```
E − e·sin E = M
```

Design choices:

* **Wrapping.** `M` is wrapped to `[−π, π]` for a well-conditioned Newton start,
  and the removed full turns are added back so the returned `E` stays on the
  same branch as the input.
* **True anomaly via the stable half-angle form.** `ν = 2·atan2(β·sin(E/2),
  cos(E/2))` with `β = √((1+e)/(1−e))`, rather than the `arccos` form. This
  lands in the correct quadrant and stays continuous through perihelion.
* An inverse `mean_anomaly(ν, e)` exists purely so round-trip identities can be
  tested.

### 2.2 Insolation — [`orbital_climate/insolation.py`](orbital_climate/insolation.py)

Standard astronomical daily-mean formulation:

```
solar longitude     λ  = ν + λ_p
declination         δ  = asin(sin ε · sin λ)
sunset hour angle   H₀ = acos(−tan φ · tan δ)          [clamped to ±1]
daily-mean flux     Q  = (S₀/π)/r² · (H₀ sin φ sin δ + cos φ cos δ sin H₀)
```

The `arccos` argument is clamped, which handles **polar day** (`H₀ = π`) and
**polar night** (`H₀ = 0`) as the two saturated branches — no special-casing.

**Key numerical decision:** time-averaging is done by sampling the **mean
anomaly** uniformly, because `M` advances linearly in time. This is what makes
the analytic annual-mean identity an honest end-to-end test of the whole
Kepler→insolation chain rather than a tautology.

### 2.3 Validation

| Target | Expected | Achieved |
|---|---|---|
| Analytic annual mean `⟨S⟩ = S₀/(4a²√(1−e²))` | exact | rel. error **2.3 × 10⁻⁶** |
| 65 °N summer peak, present day | ~480 W/m² | **478.0** |
| 65 °N summer peak, `e = 0.117` | ~400 W/m² | **399.5** |
| Peak drop | −17% | **−16.4%** |

Supporting tests: Kepler residual `< 10⁻¹¹` across `e ∈ [0, 0.9]`; circular
limit `E = M`; round-trip `M → E → ν → M`; perihelion/aphelion radii;
declination bounded by obliquity; polar day/night; equator-at-equinox flux
`= S₀/π`; global mean scaling as `1/a²`.

The global mean uses an **equal-area grid in `sin φ`**, so no `cos φ` weighting
bug can creep in — the plain cell mean *is* the area-weighted mean.

---

## 3. Part II — The energy-balance model

### 3.1 Formulation — [`orbital_climate/ebm.py`](orbital_climate/ebm.py)

North / Budyko / Sellers lineage:

```
C ∂T/∂t = Q(x,t)·a(x,T) − (A + B·T) + D ∂/∂x[(1 − x²) ∂T/∂x]
```

with `x = sin(latitude)`, `T` in °C, coalbedo `a = a₀ + a₂·P₂(x)` where ice-free
and `a_ice` where `T < T_ice`.

### 3.2 Discretisation

* **Equal-area cell-centred grid in `x`.** Cells have equal area, so global
  means are simple cell means.
* **Flux-divergence diffusion operator** with weight `w = 1 − x²` evaluated at
  cell *interfaces*. Because `w = 0` at both poles, the no-flux boundary
  condition is satisfied automatically and the operator conserves energy
  **exactly** — its area-weighted sum vanishes for any field.

### 3.3 Numerics — semi-implicit (IMEX)

The linear diffusion and OLR terms are treated **implicitly**; the nonlinear
`Q·a(T)` source **explicitly** at the old temperature:

```
M_imp = diag(C/Δt + B) − D·L         (LU-factorised once, reused every step)
rhs   = (C/Δt)·T + Q·a(T) − A
```

This is unconditionally stable for the stiff terms, so `dt_days` is an
**accuracy** knob, not a stability limit — verified consistent from
`dt = 1 day` to `dt = 30 days`.

### 3.4 Calibration

`coalbedo_a0` was scanned to hit present-day Earth:

| `a₀` | global mean |
|---|---|
| 0.670 | 287.05 K |
| **0.676** | **288.15 K** ✓ |
| 0.680 | 288.88 K |

Default parameters: `A = 203.3`, `B = 2.09` W/m²/°C, `D = 0.58` W/m²/°C,
`a₂ = −0.200`, `a_ice = 0.38`, `T_ice = −10 °C`, `C = 1.05 × 10⁸` J/m²/°C
(τ = C/B ≈ 1.6 yr), `n_lat = 180`.

### 3.5 Validation

| Test | Result |
|---|---|
| Diffusion annihilates a constant field | ✓ |
| Diffusion conserves energy | area-sum ~10⁻¹² |
| **Local radiative equilibrium** (`D = 0`, no ice): `⟨T⟩ᵢ = (⟨Q⟩ᵢ·a₀ − A)/B` | matches to **0.05 °C** |
| Global energy budget closes: `⟨absorbed⟩ = A + B⟨T⟩` | ✓ to 0.05 W/m² |
| Present-day global mean | **288.15 K** |
| Ice-albedo feedback moves ice edge equatorward under dimmer Sun | ✓ |
| Semi-implicit stability `dt = 1 → 30 days` | ✓ |

The local-radiative-equilibrium test is the strongest of these: with diffusion
disabled and the ice feedback suppressed, every latitude must independently
reach an analytically known temperature. It exercises the entire stepper.

### 3.6 Configurability refactor

`diag_lat_deg` (the Milankovitch diagnostic latitude, 65 °N by convention) was
originally **hardcoded** in the transient diagnostics. It was promoted to a
`Config` field — YAML-settable, CLI-overridable (`--diag-lat`), and sweepable.
Associated keys were renamed `T65N_*` → `Tdiag_*` since "65N" was no longer
accurate.

---

## 4. Part III — Two-surface land/ocean extension

### 4.1 The problem

With a single ocean-like heat capacity everywhere, the seasonal cycle is
heavily damped: 65 °N seasonal range came out at ~10 K versus ~40 K observed
over continental interiors. Since seasonal extremes are precisely what drives
glacial inception, this was the most load-bearing simplification remaining.

### 4.2 Why the obvious fix does not work

The tempting cheap fix is to make `C` latitude-dependent as an area-weighted
blend, `C = f·C_land + (1−f)·C_ocean`. The solver already accepted a `C` array,
so this was a two-line change. **It barely does anything:**

> At 60 °N (55% land): τ goes from 1.59 yr → **1.44 yr — a 9% change.**

The reason is that the arithmetic mean is **dominated by the ocean term**
(2.1 × 10⁸ vs 1.2 × 10⁶ J/m²/°C — a factor of **175**). Averaging heat
capacities is simply the wrong operation.

### 4.3 The correct formulation (North & Coakley 1979)

Each latitude carries **separate land and ocean temperatures**, coupled by a
zonal exchange term:

```
C_l ∂T_l/∂t = Q·a(T_l) − (A + B·T_l) + D∇²T_l + ν(T_o − T_l)
C_o ∂T_o/∂t = Q·a(T_o) − (A + B·T_o) + D∇²T_o − ν(T_o − T_l)·f/(1−f)
```

The **`f/(1−f)` factor is required for energy conservation**: once weighted by
their respective areas, the flux leaving the ocean column must equal the flux
entering the land column. Without it the model leaks energy.

Implementation details:

* State becomes a `2N` vector `[T_land, T_ocean]`; the implicit operator becomes
  a `2N × 2N` block matrix, still LU-factorised once (180×180 at `n_lat = 90` —
  trivial).
* **Ice-albedo is evaluated per surface**, so land can freeze while the ocean at
  the same latitude does not — physically important and now represented.
* `blend()` is **idempotent** (a state already blended passes through), so
  callers can blend defensively without tracking whether an upstream step
  already did.
* Land fraction is clamped to `[10⁻³, 1−10⁻³]` so `f/(1−f)` stays finite at
  Antarctica (`f = 1`).
* The whole mode is **opt-in** via `two_surface: true`, so all pre-existing
  behaviour and tests are preserved unchanged.

### 4.4 Land fraction

Earth's real zonal profile, interpolated from 10° band centres — Antarctica in
the far south, near-pure ocean in the southern mid-latitudes, most continental
in the northern mid-latitudes.

**Validation:** integrates to **0.290** global land fraction (Earth: 0.29).

### 4.5 Calibration of the coupling ν

| ν [W/m²/°C] | 65 °N land range | 65 °N ocean range |
|---|---|---|
| 1.0 | 64.1 K | 7.3 K |
| 3.0 | 43.3 K | 8.8 K |
| **3.5** (chosen) | **40.1 K** | **9.1 K** |
| 4.0 | 37.5 K | 9.3 K |
| 8.0 | 25.1 K | 10.2 K |

Target: ~40 K over land, ~8–9 K over ocean at 65 °N. **ν = 3.5** matches both.

### 4.6 Effect

| At 65 °N, present-day Earth | single surface | two surface |
|---|---|---|
| seasonal range (blended) | 11.5 K | **34.9 K** |
| seasonal range over land | — | **40.1 K** |
| seasonal range over ocean | — | **9.1 K** |
| global mean T | 288.14 K | **288.07 K** |

The global mean is **essentially unchanged** — equilibrium is set by the
radiation balance, not by heat capacity. Only the seasonal cycle changes, which
is exactly the intent.

---

## 5. Part IV — Bridging REBOUND to the climate model

### 5.1 The gap

Each REBOUND run writes `*__planets_run_deltas.csv` containing Earth's `a`, `e`
and `q` before and after the flyby — the orbit's **size and shape**. Driving a
*seasonal* climate model also needs its **orientation**: the longitude of
perihelion `λ_p` (which hemisphere gets summer at perihelion) and the obliquity.
Neither is in that file.

Also: a YAML `sweep:` block computes a **Cartesian product**, but these
parameters are **paired** — run *N*'s `(a, e, ε, λ_p)` belong together. Sweeping
662 values of `a` against 662 values of `e` would give 438,244 physically
meaningless combinations. Hence a bridge script rather than configuration.

### 5.2 The key discovery: the simulation frame

Testing rather than assuming: at `t = 0`, Earth's orbit is inclined **23.457°**
to the z-axis — exactly the obliquity. Therefore **the REBOUND simulation runs
in the equatorial J2000 frame, and the frame's z-axis *is* Earth's spin axis.**

This unlocks everything, because a point-mass integration applies **no torque**
to that axis: it stays fixed in inertial space while the BH tilts the *orbital
plane*. So:

```
obliquity      ε    = angle(ẑ, ĥ)                    ĥ = orbit normal
equinox dir    ê    = ĥ × ŝ⊥                          ŝ⊥ = ẑ projected into orbit plane
λ_p                 = angle from ê to the eccentricity vector, about ĥ
```

The equinox convention follows the insolation module: `sin δ = sin ε · sin λ`
with `λ = 0` at the northern vernal equinox, which places `λ = 0` at 90° from
`ŝ⊥` along the direction of orbital motion.

### 5.3 Validation of the recovery

Applied at `t = 0`, the recovery must reproduce present-day Earth — and does:

| Quantity | Recovered | Expected |
|---|---|---|
| `a` | 0.99981 AU | 0.9998140 (CSV) ✓ |
| `e` | 0.01717 | 0.0171676 (CSV) ✓ |
| obliquity | 23.457° | 23.44° ✓ |
| **λ_p** | **282.29°** | ~283° today ✓ |

Three independent cross-checks:

1. **Perihelion date.** The simulation's minimum heliocentric distance falls on
   **Jan 1** — and `λ_p ≈ 282°` means exactly "perihelion in early January."
2. **Zero variance.** All 672 runs return *identical* pre-flyby values to every
   digit, as they must, since all start from the same epoch.
3. **Bound count.** 662/672 Earth-bound, matching the independent count from the
   deltas CSVs exactly.

### 5.4 On the 1873 epoch and precession

The runs start at 1873-09-01 (so BH perihelion falls in ~2027) and integrate
112,420 days (~308 years). Does the early epoch invalidate the J2000 claim?

**No — J2000 is a fixed frame by definition**, the mean equator and equinox *at*
2000.0. Skyfield/JPL return ICRF coordinates for *any* query date; asking for
1873 gives 1873 positions *expressed in the J2000 frame*.

But there is a real approximation: Earth's **actual** spin axis precesses on a
25,772-year cycle at ~20″/yr, so it was ~0.70° from `ẑ_J2000` in 1873 and
~1.01° the other way by 2181. Using `ẑ` is effectively a mid-run average.

**Consequence:** the true `λ_p` drifts ~**4.3°** over the run relative to the
fixed-frame value. This is negligible here because (a) the BH-induced `λ_p`
spread is the full 0–360°, so 4.3° is ~1% of the signal; (b) REBOUND treats
Earth as a point mass with no spin, so precession **cannot** be modelled
regardless of frame choice; and (c) before/after both use the same fixed `ẑ`, so
*differences* remain internally consistent.

### 5.5 Extraction performance

Archives were disabled for this sweep, so the ~220 MB `*__orbits__*.xlsx` files
were the only source. Only the **first and last rows** are needed.

| Method | Per run | Full sweep |
|---|---|---|
| `openpyxl` full scan | 26 s | 4.8 h |
| **Stream-decompress zip, keep head/tail bytes** | **0.83 s** | **~10 min** |

A **30× speedup**, achieved by reading the sheet XML straight out of the xlsx
zip and retaining only the boundary bytes.

### 5.6 Two corrections applied per run

1. **Year length.** Kepler's third law: `P = a^1.5` years. Seasonal damping
   depends on year length relative to the ocean time constant (τ ≈ 1.6 yr), so
   this materially changes the answer. In the usable band, year length spans
   **184–1711 days** — a 9× range.
2. **Timestep.** `dt_days = P / 180` rather than a fixed value, so every run
   resolves its seasonal cycle with the same fidelity regardless of orbital
   period.

### 5.7 Validity band

The linear OLR (`A + B·T`) is calibrated near 288 K. Runs far outside that
regime still return numbers, but fictitious ones. Excluded by default:
`a ≥ 3 AU` or `e ≥ 0.9`.

| | count |
|---|---|
| Total runs | 672 |
| Earth unbound (ejected) | 10 |
| Outside validity band | 18 |
| **Usable** | **644 (96%)** |

### 5.8 Distribution of recovered elements (662 bound runs)

| | min | median | max |
|---|---|---|---|
| `a` [AU] | 0.632 | 1.038 | 1097.6 |
| `e` | 0.005 | 0.126 | 0.9992 |
| obliquity [°] | 6.1 | 23.29 | 74.4 |
| **λ_p** [°] | 0.2 | 172.1 | 359.4 |
| year length [d] | 184 | 386 | 13.3 M |

**λ_p spans essentially the full circle** — recovering it was necessary, not
optional. Freezing it at 283° would have been a genuine error.

---

## 6. Part V — Results

### 6.1 Equilibrium climate across the sweep

| | |
|---|---|
| Pre-flyby baseline | **288.19 K** (single-surface) |
| Median outcome | 278.97 K (**−9.2 K**) |
| Range | 187.8 → 559.7 K |
| Recognisably Earth-like (±10 K) | **138 runs (21%)** |
| Snowball (fully glaciated) | **250 runs (39%)** |

### 6.2 The ice-albedo bifurcation

The temperature distribution is **bimodal with a forbidden gap at ~235–260 K**.
Temperatures jump from ~265 K straight to ~230 K around `a ≈ 1.07 AU` with
nothing in between — the classic Budyko–Sellers ice-albedo catastrophe, which
emerged from the model rather than being imposed.

There is an **overlap zone** (`a = 1.073–1.122 AU`) where both temperate and
snowball outcomes occur, depending on the other orbital parameters.

### 6.3 What controls the outcome

| Predictor | correlation with T |
|---|---|
| **Annual-mean insolation `S₀/(4a²√(1−e²))`** | **+0.966** |
| semi-major axis `a` | −0.790 |
| longitude of perihelion `λ_p` | +0.593 → **spurious** |
| eccentricity | +0.110 |
| obliquity | +0.067 |

The `λ_p` correlation looked substantial but is **collinearity**: `λ_p` is 42%
correlated with `a` in this sweep. Controlling for annual insolation collapses
it to **−0.22**, while `S_mean` explains **+0.998**. This matches the physics
exactly — `λ_p` shifts *seasonal phasing*, not the annual mean. The residual
−0.22 is the genuine (small) ice-albedo/seasonal effect.

### 6.4 Orbital disruption ≠ climate disruption

Correlation between the orbital-impact `Score` (from
[`rank_run_impact.py`](rank_run_impact.py)) and `|ΔT|` is only **+0.586**. A run
can eject Mercury and Neptune while leaving Earth's climate nearly untouched, so
the two rankings are genuinely complementary rather than redundant.

### 6.5 Transient adjustment

| | single surface | two surface |
|---|---|---|
| years to equilibrium (median) | 11 | 15 |
| peak hemispheric asymmetry | 3.74 K median, 18.29 K max | 4.69 K median, 22.30 K max |
| …the same at equilibrium | **0.00 K** median | **0.00 K** median |
| temperature overshoot | **0 everywhere** | **0 everywhere** |
| ice-edge migration (median) | −6.7° | −7.9° |
| peak migration rate | up to 44°/yr | up to 30°/yr |

**The headline transient result:** **198 runs (31%)** show a peak hemispheric
asymmetry above 5 K that decays to under 1 K, typically peaking in **year 1**.
Nearly a third of the sweep passes through a substantially lopsided climate that
leaves **no trace whatsoever** in the equilibrium state — information the
equilibrium sweep structurally cannot produce.

642 of 644 runs settle within 60 years; median 11, 90th percentile 18. The
climate reaches its new state within one to two decades.

### 6.6 Two-surface results

Within the defensible 250–300 K band (220 runs), 65 °N seasonal range:

| | median |
|---|---|
| Single-surface | 10.3 K |
| Two-surface, blended | 23.2 K |
| **Two-surface, land** | **41.2 K** |
| Two-surface, ocean | 8.5 K |

Global mean essentially unchanged (280.5 → 282.6 K). Snowball count 240 vs 250.

### 6.7 The Milankovitch signal — and the regime dependence

In a **controlled experiment** (`a` fixed at 1.0, `e: 0.0167 → 0.117`, the
scenario where 65 °N June insolation drops 480 → 400 W/m²):

| Model | 65 °N summer peak change |
|---|---|
| Single surface | **−0.04 K** — nothing |
| Two-surface, **land** | **−7.53 K** |
| Two-surface, ocean | −0.02 K |

The two-surface model recovers the glacial-inception signal the single-surface
model entirely missed. With a 7-day land time constant the surface tracks the
insolation drop almost directly, while the 3.2-year ocean damps it to nothing.
For reference, the context document cites "~4–6 K direct orbital summer cooling
at high northern latitudes" — the right magnitude and, crucially, the right sign.

**However, this does not generalise to the whole sweep.** Across all 220
habitable-band runs:

| 65 °N summer-peak change | single surface | two-surface land |
|---|---|---|
| median | **−6.31 K** | **−3.37 K** |
| runs cooling > 4 K | 138 (63%) | 104 (47%) |

The single-surface model shows *more* cooling. The reason: in this sweep Earth's
semi-major axis varies from 0.63 to 2.8 AU, so **annual-mean insolation changes
dominate**. A heavily damped single-surface world has summer ≈ annual mean and
inherits the full annual-mean cooling; land, sitting well above the annual mean
in summer, is partly buffered.

**Isolating the mechanism** — restricting to the 72 runs where `a` stayed within
2% of 1.0 AU (orbital *shape* change only, annual mean preserved at +0.45 K):

| | median |
|---|---|
| Single-surface summer | **+0.48 K** (warming) |
| Two-surface **land** summer | **−0.98 K** (cooling) |

**A sign flip.** With the annual mean held fixed, the single-surface model
reports summer *warming* while land correctly reports *cooling*.

**Accurate statement:** the two-surface model is right in both regimes, but the
correction it applies depends on what changed. For pure orbital-**shape**
changes it reveals a summer signal the old model gets backwards; for large
orbital-**size** changes it moderates an overestimate.

---

## 7. Corrections and negative results

Recorded deliberately — several initial claims did not survive testing.

### 7.1 Claims that were wrong

| Claim | Reality |
|---|---|
| "Obliquity change is a first-order effect" — based on one sampled run showing 23° → 74° | Full distribution: median change ~0; only **3%** of runs shift >10°, 0.3% >30°. The 74° case is an outlier. |
| "Adding two surfaces will require re-tuning `coalbedo_a0` and all results will shift" | **Not needed.** Global mean moved 288.14 → 288.07 K. Equilibrium is set by radiation balance, not heat capacity. |
| "The transient may land in a different attractor than `run_equilibrium` (which starts from a uniform 15 °C guess rather than Earth's real profile)" | **Tested 12 runs straddling the bifurcation: zero flips.** Equilibrium results are robust to initial conditions. |
| "The old model missed the glacial-inception signal" | **Regime-dependent.** True for pure eccentricity injection; across this sweep the single-surface model actually *over*-predicts summer cooling (63% vs 47%). |
| A hardcoded Kepler test constant of 1.1934205 | Recomputed independently: **1.1853242**. The original would have validated a broken solver. |
| Apparent `λ_p` → temperature correlation of +0.59 | **Spurious** — collinear with `a`. Controlling for insolation: −0.22. |

### 7.2 Genuine negative results

* **Temperature overshoot is exactly zero across all 644 runs**, in both single-
  and two-surface modes. The approach to equilibrium is monotonic. This was one
  of the stated motivations for running transients and it simply does not occur.
* **Initial-condition bistability does not occur** in the tested range, despite
  the model having a genuine bifurcation.
* **A blended heat capacity is not a useful approximation** (9% effect) — worth
  recording so the cheap approach is not re-attempted.

---

## 8. Bugs found and fixed

| Bug | Where | Consequence |
|---|---|---|
| Regex assumed XML attribute order (`Id` before `Target`); this writer emits the reverse | `extract_earth_elements.py` | Sheet lookup returned nothing; now parses XML properly |
| Double-blend: `ice_line_lat` re-blended an already-blended state | `ebm.py` | Shape error in two-surface mode; `blend()` made idempotent |
| `--limit` truncated the dataframe *before* statistics were computed | `climate_from_simulations.py` | Reported "654 outside validity band" when only 8 runs were requested |
| `ndarray.ptp()` removed in NumPy 2 | test suite | Test error; switched to `np.ptp()` |
| Hardcoded diagnostic-name list | `run_perturbation` | Land/ocean keys silently missing from transients |
| `rebound` 5.x renamed the `hash=` kwarg to `name=` | `solar_system_bh_rebound26.py` | All 1050 sweep combos failed instantly (pre-existing, fixed early) |
| Chained sweeps with `&&` | run scripts | A locked output file aborted the first pass and silently skipped the second |

---

## 9. Known limitations

1. **Linear OLR.** `A + B·T` is calibrated near 288 K. Only **33%** of runs land
   in the defensible 250–300 K band. Outside it, results are **ordinal, not
   quantitative** — at the 559.7 K extreme the linearisation implies ~800 W/m²
   outgoing where blackbody physics demands ~5500 W/m².
2. **Hard ice-albedo step** at `T_ice`. The *existence* of the bifurcation is
   robust EBM physics; its exact location is parameterisation-dependent.
3. **Fixed obliquity within a run** — no precession or obliquity cycling. The
   REBOUND runs cannot supply it (Earth is a point mass).
4. **Extreme-tail seasonal ranges** reach 311 K over land in small-`a`, high-`e`
   runs. Filter to the validity band before interpreting.
5. **Zonal-mean geometry.** Even with two surfaces, all land at a given latitude
   shares one temperature — no distinction between maritime and deep-continental
   climates at the same latitude.
6. **No carbonate-silicate thermostat, clouds, or ocean heat transport** beyond
   the diffusive parameterisation.
7. **`two_surface` cannot be meaningfully swept** by the sweep harness (it is a
   boolean).

---

## 10. File inventory and reproduction

### Package — [`orbital_climate/`](orbital_climate/)

| File | Role |
|---|---|
| `kepler.py` | Kepler solver, true anomaly, radius |
| `insolation.py` | Daily-mean insolation, analytic annual mean |
| `ebm.py` | EBM: grid, diffusion, coalbedo, semi-implicit stepper, spin-up, two-surface mode, Earth land-fraction profile |
| `experiment.py` | `run_equilibrium`, `run_perturbation`, output writers |
| `sweep.py` | Parallel Cartesian-product sweep → parquet/CSV |
| `config.py` | `Config` dataclass + YAML loading |
| `cli.py` | `insolation` / `ebm` / `sweep` subcommands |
| `tests/` | **54 tests** |
| `MANUAL.md` | Usage manual |

### Bridge scripts — repository root

| File | Role |
|---|---|
| `extract_earth_elements.py` | Recovers `a, e, ε, λ_p` per run from the orbits workbooks (30× faster than openpyxl) |
| `climate_from_simulations.py` | Runs the EBM on every simulation; equilibrium or transient; parallel |
| `rank_run_impact.py` | Ranks runs by orbital disruption (energy/eccentricity/perihelion metric) |
| `input_climate.yaml` | Shared configuration |

### Reproduction

```bash
# 0. Dependencies
python -m pip install numpy scipy pandas matplotlib pyyaml pyarrow pytest

# 1. Validate the physics
python -m pytest orbital_climate/tests/ -v            # 54 tests

# 2. Rank runs by orbital disruption
python rank_run_impact.py simulations/<STAMP> --plot impact.png

# 3. Recover Earth's post-flyby orbital elements (~10 min, 5 workers)
python extract_earth_elements.py simulations/<STAMP> --workers 5

# 4. Climate across the sweep — equilibrium
python climate_from_simulations.py simulations/<STAMP> \
    --config input_climate.yaml --two-surface --workers 5 --plot climate.png

# 5. Climate across the sweep — transient
python climate_from_simulations.py simulations/<STAMP> \
    --config input_climate.yaml --two-surface --transient --years 40 --workers 5
```

### Outputs

| File | Contents |
|---|---|
| `<STAMP>_impact_ranking.csv` | Orbital disruption score per run |
| `<STAMP>_earth_elements.csv` | Recovered `a, e, ε, λ_p`, year length, bound flag |
| `<STAMP>_climate.csv` | Single-surface equilibrium climate |
| `<STAMP>_climate_2surf.csv` | Two-surface equilibrium climate + land/ocean seasonal columns |
| `<STAMP>_climate_2surf_transient.csv` | Transient adjustment metrics |

---

## Appendix — validation summary

Every quantitative claim in this report, and where it was checked:

| Claim | Check | Result |
|---|---|---|
| Insolation pipeline correct | analytic `⟨S⟩ = S₀/(4a²√(1−e²))` | 2.3 × 10⁻⁶ rel. error |
| Milankovitch benchmark | 65 °N June peak, `e = 0.117` | 480 → 400 W/m² (−16.4%) |
| EBM stepper correct | local radiative equilibrium, `D = 0` | matches analytic to 0.05 °C |
| Diffusion conserves energy | area-weighted sum of `L·T` | ~10⁻¹² |
| Energy budget closes | `⟨absorbed⟩` vs `A + B⟨T⟩` | ✓ 0.05 W/m² |
| Present-day calibration | global mean | 288.15 K |
| λ_p recovery correct | applied at `t = 0` | 282.29° vs ~283° expected |
| λ_p semantics correct | perihelion date in simulation | Jan 1 ✓ |
| Extraction deterministic | pre-flyby values across 672 runs | zero variance |
| Bound/unbound classification | vs independent deltas CSVs | 662/672 both ✓ |
| Land fraction profile | area-weighted global mean | 0.290 (Earth 0.29) |
| Land/ocean calibration | 65 °N seasonal ranges | 40.1 K / 9.1 K vs ~40 / ~8–9 |
| Two surfaces preserve energetics | global mean shift | 288.14 → 288.07 K |
| Equilibrium robust to initial state | 12 runs at the bifurcation | zero flips |
