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

## Which sweep the numbers describe

**Read this before quoting any figure from this document.**

Results here come from **two different sweeps**, and the sections are marked
accordingly:

| | retired sweep | **current sweep** |
|---|---|---|
| tag | `20260724_230314` | **`20260811_184731`** |
| epoch | 1873-09-01 | **1885-09-01** |
| BH perihelion | ~2027 | **2047** |
| `v_inf` | 25 km/s | 25 km/s — **unchanged** |
| `omega` grid | 90° steps (4 values) | **15° steps (24 values)** |
| runs | 672 | **4032** |
| on disk? | **deleted 2026-08-15** | yes |

**The two sweeps differ in less than it looks.** `v_inf`, BH mass, and the `rp`,
inclination and `Omega` grids are identical. Only two things changed: the epoch
moved 12 years (carrying perihelion from 2027 to 2047, so the planets sit
differently at encounter), and the `omega` grid was refined 6-fold. Any
difference in results traces to those, **not** to encounter speed.

- **Parts I–IV (§2–§5)** describe how the model was built and validated. Physics,
  not sweep output — unaffected by the change, except §5.4 and §5.7–§5.8, which
  have been restated for the current sweep.
- **§6.1–§6.4, §7.6–§7.7** have been **regenerated against the current sweep**.
- **Transient, two-surface, Milankovitch and Sellers results have been removed**
  (they occupied §6.5–§6.8). They were measured on the retired sweep and were never
  re-run, since each needs a separate pass of the climate model with non-default
  options. Rather than leave stale numbers behind a warning label, the sections are
  gone and §6.5 now records what they covered, what survives of them, and how to
  regenerate. **This report currently has no results for transient behaviour, the
  land/ocean split, or nonlinear OLR.**
- **§8 (corrections)** deliberately preserves superseded numbers where they record
  a claim that turned out wrong. That is the point of the section — do not "fix"
  those.

The current sweep's derived products live in the repository root as
`simulations/20260811_184731_{climate,mars_climate,earth_elements,impact_ranking,bh_captures}.csv`.
Every regenerated figure below is computed from those files.

---

## Table of contents

1. [Goal and overall architecture](#1-goal-and-overall-architecture)
2. [Part I — Orbital physics core](#2-part-i--orbital-physics-core)
3. [Part II — The energy-balance model](#3-part-ii--the-energy-balance-model)
4. [Part III — Two-surface land/ocean extension](#4-part-iii--two-surface-landocean-extension)
5. [Part IV — Bridging REBOUND to the climate model](#5-part-iv--bridging-rebound-to-the-climate-model)
6. [Part V — Results](#6-part-v--results) — incl. [§6.9 the capture census](#69-where-the-planets-end-up--the-capture-census)
7. [Part VI — Mars: a condensing atmosphere](#7-part-vi--mars-a-condensing-atmosphere)
8. [Corrections and negative results](#8-corrections-and-negative-results)
9. [Bugs found and fixed](#9-bugs-found-and-fixed)
10. [Known limitations](#10-known-limitations)
11. [Data volume and archival](#11-data-volume-and-archival)
12. [File inventory and reproduction](#12-file-inventory-and-reproduction)
13. [References](#13-references)

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
 (4032 runs)         (a, e, obliquity, λ_p, year length)      (climate outcome)
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

The linear OLR `A + B·T` shown here is the default. A nonlinear alternative
(Sellers 1969) was added later once its validity range proved too narrow for
this sweep — see §10 for the failure analysis. The sweep-wide consequences were
measured only on the retired sweep and have been removed pending regeneration (§6.5).

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

**Extension for nonlinear OLR.** A nonlinear `OLR(T)` would ordinarily force the
matrix to be re-factorised every step. Instead the implicit operator *keeps* its
linear `B·T` relaxation term and the explicit source carries `B·T − OLR(T)`:

```
rhs = (C/Δt)·T + Q·a(T) − OLR(T) + B·T
```

The `B·T` contributions **cancel exactly at convergence**, so `olr_B` becomes a
pure numerical preconditioner while the converged state satisfies the true
balance `D∇²T + Q·a(T) − OLR(T) = 0`. With `olr_model = "linear"` the source
collapses to `Q·a(T) − A` and the scheme is identical to the original — which is
why all 54 pre-existing tests passed unchanged when the option was added.

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
| seasonal range (blended) | 11.5 K | **23.0 K** |
| seasonal range over land | — | **40.1 K** |
| seasonal range over ocean | — | **9.1 K** |
| global mean T | 288.14 K | **288.42 K** |

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
3966 values of `a` against 3966 values of `e` would give 15.7 million physically
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
2. **Near-zero variance.** All runs must return the same pre-flyby values, since
   they all start from the same epoch. In the retired sweep they were identical to
   every digit. In the current sweep they agree only to ~10⁻⁵ relative — 1697
   distinct values of `a` across 4032 runs, spanning 1.75×10⁻⁵ AU. This is **not
   a regression in the recovery**: it is float32 rounding introduced by the
   Parquet switch (§11), where the retired sweep's xlsx stored full-precision
   float64 as text. The check still passes, but as a tolerance rather than an
   equality — see §8 for why the tolerance is larger than §11 originally claimed.
3. **Bound count.** 3966/4032 Earth-bound in the current sweep (662/672 in the
   retired one), matching the independent count from the deltas CSVs exactly.

### 5.4 On the 19th-century epoch and precession

The current runs start at **1885-09-01** (so BH perihelion falls in 2047) and
integrate 115,342 days (~316 years). The retired sweep started 1873-09-01 for a
2027 perihelion. The argument below is unchanged by the 12-year shift — it is
about the size of the precession term, which moves by well under a degree — and
is stated for the current epoch.

Does the early epoch invalidate the J2000 claim?

**No — J2000 is a fixed frame by definition**, the mean equator and equinox *at*
2000.0. Skyfield/JPL return ICRF coordinates for *any* query date; asking for
1885 gives 1885 positions *expressed in the J2000 frame*.

But there is a real approximation: Earth's **actual** spin axis precesses on a
25,772-year cycle at ~20″/yr, so it was ~0.64° from `ẑ_J2000` in 1885 and
~1.12° the other way by 2201. Using `ẑ` is effectively a mid-run average.

**Consequence:** the true `λ_p` drifts ~**4.3°** over the run relative to the
fixed-frame value. This is negligible here because (a) the BH-induced `λ_p`
spread is the full 0–360°, so 4.3° is ~1% of the signal; (b) REBOUND treats
Earth as a point mass with no spin, so precession **cannot** be modelled
regardless of frame choice; and (c) before/after both use the same fixed `ẑ`, so
*differences* remain internally consistent.

### 5.5 Extraction performance

> **Historical.** This describes the retired sweep, where the only source was
> xlsx. The current sweep writes Parquet directly and the head/tail rows come
> out of it in milliseconds, so the optimisation below no longer runs — see §11.
> It is kept because the measurement is what justified changing the engine.

Archives were disabled for that sweep, so the ~220 MB `*__orbits__*.xlsx` files
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

*Current sweep (`20260811_184731`).*

| | count |
|---|---|
| Total runs | 4032 |
| Earth unbound (ejected) | 66 |
| Excluded (unbound or outside band) | 128 |
| **Usable** | **3904 (97%)** |

Earth survives bound in **98.4%** of runs, against 98.5% in the retired sweep —
refining the `omega` grid and moving the epoch did not measurably change Earth's
odds of staying.

### 5.8 Distribution of recovered elements (3966 bound runs)

*Current sweep.*

| | min | median | max |
|---|---|---|---|
| `a` [AU] | 0.458 | 0.985 | 82.8 |
| `e` | 0.003 | 0.127 | 0.990 |
| obliquity [°] | 3.6 | 23.66 | 131.2 |
| **λ_p** [°] | 0.8 | 222.4 | 359.4 |
| year length [d] | 113 | 357 | 275,000 |

Two differences from the retired sweep are worth noting. The median `a` is now
**0.985 AU against 1.038** — Earth ends up *inward* of where it started on
median, rather than outward, which propagates directly into §6.1. And obliquity
now reaches **131°**, past 90°, meaning some runs flip Earth's spin axis into
retrograde; the old sweep's maximum was 74°. The extreme `a` tail is also far
shorter, 82.8 AU against 1097.6.

**None of this is an encounter-speed effect** — both sweeps run at
`v_inf` = 25 km/s. The 6-fold refinement of the `omega` grid is the likely
driver: `omega` sets where periapsis falls relative to the planets, the old
sweep sampled it at only four values 90° apart, and §6.3 shows the sweep's
outcomes are dominated by where Earth ends up rather than by how hard it was
hit. Four coarse samples of that angle can easily land a median in a different
place than twenty-four. The epoch shift contributes too, since the planets sit
elsewhere at encounter. **This has not been separated into the two causes**, and
doing so would need a sweep at the old `omega` resolution on the new epoch.

**λ_p spans essentially the full circle** — recovering it was necessary, not
optional. Freezing it at 283° would have been a genuine error.

---

## 6. Part V — Results

*§6.1–§6.4 regenerated against the current sweep (`20260811_184731`), linear OLR,
single surface — the same configuration as the retired-sweep figures they replace,
so the comparison is like for like.*

### 6.1 Equilibrium climate across the sweep

| | current (4032) | retired (672) |
|---|---|---|
| Pre-flyby baseline | **288.19 K** | 288.19 K |
| Median outcome | **293.26 K (+5.1 K)** | 278.97 K (−9.2 K) |
| Range | **185.2 → 976.4 K** | 187.8 → 559.7 K |
| Recognisably Earth-like (±10 K) | **810 runs (20.7%)** | 138 (21%) |
| Snowball (fully glaciated) | **1182 runs (30.3%)** | 250 (39%) |
| Runaway flagged | **894 runs (22.9%)** | — |

**The median sign flipped.** The retired sweep cooled Earth by 9.2 K on median;
the current one *warms* it by 5.1 K. This is not a model change — it follows
directly from §5.8, where median `a` moved from 1.038 to 0.985 AU. Warming
follows from Earth sitting closer to the Sun; the climate model is doing nothing
surprising.

**What moved the median is a sampling change, not a physical one.** Both sweeps
use `v_inf` = 25 km/s, the same BH mass, and the same `rp`, inclination and
`Omega` grids. What differs is the `omega` resolution (4 values → 24) and the
12-year epoch shift. So this is a caution about reading sweep medians as physics:
a median can move 14 K because the grid under it was re-sampled.

The proportion of Earth-like outcomes is essentially unchanged (20.7% vs 21%),
which is reassuring: the *spread* of outcomes is a property of the sweep geometry,
while the *centre* tracks where Earth ends up. Snowball fraction falls (30.3% vs
39%) for the same inward-shift reason, and the hot tail correspondingly extends
to 976 K — far outside any defensible range, which is what the runaway flag is
for.

### 6.2 The ice-albedo bifurcation

The temperature distribution is **bimodal with a forbidden gap at 231–260 K** —
*zero* runs out of 3904 land in 235–260 K. Temperatures jump from ~265 K straight
to ~230 K around `a ≈ 1.07 AU` with nothing in between — the classic
Budyko–Sellers ice-albedo catastrophe, which emerged from the model rather than
being imposed.

There is an **overlap zone** (`a = 1.071–1.120 AU`) where both temperate and
snowball outcomes occur, depending on the other orbital parameters.

This is the most robust result in the document: the gap and the overlap zone
reproduced to within 0.003 AU across two sweeps differing in epoch, encounter
speed and run count by a factor of six. It is a property of the climate model,
not of the sampling.

### 6.3 What controls the outcome

| Predictor | current | retired |
|---|---|---|
| **Annual-mean insolation `S₀/(4a²√(1−e²))`** | **+0.978** | +0.966 |
| semi-major axis `a` | −0.751 | −0.790 |
| longitude of perihelion `λ_p` | +0.493 → **spurious** | +0.593 |
| eccentricity | +0.204 | +0.110 |
| obliquity | +0.190 | +0.067 |

The `λ_p` correlation looked substantial but is **collinearity**: `λ_p` is 42%
correlated with `a` in this sweep. Controlling for annual insolation collapses
it to **−0.22**, while `S_mean` explains **+0.998**. This matches the physics
exactly — `λ_p` shifts *seasonal phasing*, not the annual mean. The residual
−0.22 is the genuine (small) ice-albedo/seasonal effect.

### 6.4 Orbital disruption ≠ climate disruption

Correlation between the orbital-impact `Score` (from
[`rank_run_impact.py`](rank_run_impact.py)) and `|ΔT|` is only **+0.673**
(+0.586 in the retired sweep). A run can eject Mercury and Neptune while leaving
Earth's climate nearly untouched, so the two rankings are genuinely complementary
rather than redundant. The coupling is somewhat tighter in the current sweep but
the conclusion is unchanged: knowing how violently a run disturbed the system
still explains under half the variance in what it did to Earth's climate.

### 6.5 Transient, two-surface, Milankovitch and Sellers — *removed, pending regeneration*

Four results sections stood here: transient adjustment, the two-surface
land/ocean results, the Milankovitch signal and its regime dependence, and the
sweep re-run under the Sellers nonlinear OLR. Every count and percentage in them
was measured on the retired 672-run sweep.

**They have been removed rather than left in place behind a warning.** A flagged
number is still a number, and this report is quoted from. Nothing here now
describes a sweep that no longer exists.

*(They were §6.5–§6.8. The numbering gap before §6.9 is deliberate: it marks the
excision, and it keeps external references to §6.9 — the capture census, cited
from [`SCENARIO_post_flyby_system.md`](SCENARIO_post_flyby_system.md) — resolving.)*

What they established that does *not* depend on those statistics is preserved:

- Temperature overshoot is exactly zero and the approach to equilibrium is
  monotonic; equilibrium results are robust to initial conditions; a blended heat
  capacity is not a useful approximation (§8, where all three are recorded as
  negative results, with their retired-sweep provenance marked).
- The linear OLR's cold-end failure, the Sellers form that fixes it, and the
  Simpson–Nakajima runaway ceiling are described in §10.
- The two-surface formulation itself, its calibration and its validation are
  Part III (§4) and are unaffected — they are model construction, not sweep
  output.

**Regenerating them** means three further passes of
[`climate_from_simulations.py`](climate_from_simulations.py) over the current
sweep with non-default options (`--two-surface`, `--olr-model sellers`,
`--transient`); see §12 for the commands. Until then this report has **no current
results** for transient behaviour, the land/ocean split, or nonlinear OLR, and
should not be cited for any of the three.
### 6.9 Where the planets end up — the capture census

*Current sweep. New section: this analysis did not exist when the report was
first written.*

An ejection and a capture look identical in the per-run deltas — both show a body
that is no longer bound to the Sun. [`find_bh_captures.py`](find_bh_captures.py)
separates them with a two-body energy test (`ε = v²/2 − μ/r`) against the black
hole, evaluated on the final state vectors. Across 4032 runs × 9 bodies (36,288
body-outcomes):

| outcome | count | share |
|---|---|---|
| still bound to the Sun | 34,116 | 94.0% |
| unbound from both — **free** | 2,140 | 5.9% |
| **captured by the black hole** | **32** | **0.09%** |

**Capture is rare — roughly one body-outcome in a thousand, and 1.5% of all
losses.** This is the quantitative answer to a question the earlier sweep could
only pose anecdotally: §7.7 of the first edition of this report was built around
"the one case in the sweep where Mercury is captured," which invited the reading
that capture is a characteristic outcome. It is not. A body stripped from the Sun
overwhelmingly just leaves. Ending up bound to the perturber requires the
encounter to remove almost exactly the right amount of energy, and the sweep shows
how narrow that window is.

---

## 7. Part VI — Mars: a condensing atmosphere

The Earth model was extended to Mars, driven by the same sweep. The exercise
was a test of how much of the model was physics and how much was Earth.

### 7.1 What transferred unchanged

Most of it. The Kepler solver, the insolation chain, the equal-area grid, the
diffusion operator and the IMEX stepper are all planet-agnostic and were reused
without modification. What is Earth-specific is only the **parameterisations**.

### 7.2 Radiation got *simpler*

The expectation was that a second planet would need more approximation. The
opposite happened.

Mars's absorbed flux is 110.4 W/m² and `σ(210 K)⁴ = 110.3 W/m²`, so the required
emissivity is **1.00** against Earth's 0.60. Mars's greenhouse effect is ~5 K,
not 33 K, so `OLR = σT⁴` is not a fit at all — it is the physics.

| | Earth | Mars |
|---|---|---|
| OLR form | `A + B·T`, fitted | `σT⁴`, exact |
| Required emissivity | 0.60 | **1.00** |
| Validity window | 230–300 K | **none — it is a law** |

The single largest approximation in the Earth model, and the source of its
narrowest validity limit, simply does not arise. Radiative damping also falls
out rather than being fitted: `B = 4σT³ = 2.10 W/m²/K`, coincidentally almost
identical to Earth's fitted 2.09.

Implemented as `olr_model: "graybody"` with `olr_emissivity`.

### 7.3 The CO₂ cycle — the genuinely new physics

Mars's atmosphere condenses. Roughly a quarter of it freezes onto the winter
pole each year and sublimates back in spring, coupling three things an Earth EBM
keeps separate: surface temperature, surface albedo, and atmospheric mass.

Three couplings, all represented in [`orbital_climate/mars.py`](orbital_climate/mars.py):

* **Latent buffering.** Where the surface would cool below the frost point, CO₂
  condenses instead and the latent heat pins the temperature *at* that point.
  Polar winter temperature is therefore set by thermodynamics, not by the energy
  balance.
* **Mass exchange.** Condensed CO₂ leaves the atmosphere, lowering surface
  pressure, which lowers the frost point (Clausius–Clapeyron). Condensation is
  therefore **self-limiting** — a negative feedback, unlike ice-albedo.
* **Albedo.** Coalbedo follows the *presence of frost*, not a temperature
  threshold, so winter frost persists into spring while the surface warms.

`MarsEBM` carries a second state variable (frost mass per unit area) alongside
temperature, so it is a subclass rather than a flag — the Earth model is
untouched.

**Why it is not optional.** Without latent heat the modelled polar winter runs
away to **80 K** against an observed ~148 K — an error of ~70 K in exactly the
quantity that governs the seasonal cap cycle.

### 7.4 Calibration and validation

`D` was refitted to Mars's much weaker transport: **0.002**, roughly 300× smaller
than Earth's 0.58, which is what a 6 mbar atmosphere warrants.

| Quantity | Model | Observed |
|---|---|---|
| Polar winter minimum | **146.8 K** | ~148 K |
| Surface pressure (annual mean) | **596 Pa** | ~600 Pa |
| Seasonal atmospheric swing | **22%** | ~25% |
| Peak cap thickness | **0.57 m** | ~0.5–1 m |
| Global mean | 203.4 K | ~210 K |
| CO₂ conservation | **exact** | — |

Note the global mean sits *below* the 210 K zeroth-order estimate, and correctly
so: `OLR ∝ T⁴` is convex, so a planet with strong gradients radiates more
efficiently than a uniform one at the same mean, giving `⟨T⟩ < T_eff`.

### 7.5 Element recovery, generalised

[`extract_mars_elements.py`](extract_mars_elements.py) generalises §5's recovery
to any planet. Earth's case was special — the frame's z-axis *is* its spin axis —
so the general form takes the spin axis from an **IAU pole table** instead.
Validated by returning Mars's present-day `a = 1.5237`, `e = 0.0934`,
`obliquity = 25.18°` against the true 1.524 / 0.0934 / 25.19.

One cost worth recording: it reads the **float32** Parquet tree rather than the
original float64 workbooks, so pre-flyby values vary by ~10⁻¹¹ across runs rather
than Earth's ~10⁻²². Harmless at four significant figures, but it is the price
paid for the 5.7× compression, and belongs beside the saving.

### 7.6 Sweep results

*Regenerated against the current sweep (`20260811_184731`).*

3826 usable runs of 4032 (the remainder have Mars unbound or outside the band).

**Atmospheric collapse is flagged, not reported.** 598 runs (15.6%) freeze the
atmosphere out at some point in the year, 8 (0.2%) for the whole year. They are
systematically the cold, distant cases — median `a` **2.43 AU against 1.49** for
valid runs, and less than half the annual-mean insolation.

Restricted to the 3228 runs the model can speak about:

| | median | 5th | 95th |
|---|---|---|---|
| Global mean T | 204.9 K | 175.8 | 234.4 |
| Equatorial seasonal range | 44.6 K | 13.3 | 154.6 |
| Pressure swing | 0.28 | 0.12 | 0.83 |

Every one of these is within a few percent of the retired sweep's values
(206.2 K, 42.6 K, 0.30) despite six times the runs and a different encounter
speed — Mars's climate response is set by where it lands, and the distribution of
where it lands is similar.

**Eccentricity drives Martian seasons, not obliquity** — and this inverts the
Earth result:

| | current | retired |
|---|---|---|
| Eccentricity | **+0.957** | +0.771 |
| Obliquity | +0.028 | +0.071 |

The current sweep makes the point far more sharply: with six times the sampling
the eccentricity correlation tightens to +0.957 while obliquity falls to
essentially zero.

On Earth obliquity dominates (§4). The reason is thermal inertia: with
`τ ≈ 6.6 days` Mars tracks the instantaneous `1/r²` forcing almost perfectly, so
the perihelion–aphelion distance swing arrives undamped. Earth's ocean averages
that away and leaves only the tilt signal. **The same parameter can dominate on
one planet and be negligible on another, for reasons that have nothing to do
with the parameter itself.**

**890 runs (28%) see the equator exceed 273 K**, and among those it stays above
freezing for a median 19% of the year.

**Above freezing is not the same as liquid water.** That distinction was added
after this section was first written and is the single most consequential
correction in the document — see §8. Requiring the surface to be both above the
triple point *and* below the local boiling point:

| test | runs with any liquid water | share |
|---|---|---|
| above triple point only (the weak test) | 880 | 23.0% |
| **liquid water possible, at the equator** | **823** | **21.5%** |
| **liquid water possible, at the best latitude** | **1250** | **32.7%** |

Even the best case is a *fraction of a year at one latitude*: the maximum across
the whole sweep is **7.8% of the year**. No run in 3826 produces Mars with
year-round liquid water anywhere on its surface.

### 7.7 A worked scenario

*The adopted scenario is examined in full in
[`SCENARIO_mars_window.md`](SCENARIO_mars_window.md), narrated in
[`SCENARIO_timeline.md`](SCENARIO_timeline.md), and its aftermath worked out in
[`SCENARIO_post_flyby_system.md`](SCENARIO_post_flyby_system.md). The retired
sweep's worked example — a Mercury-capture run — had its own write-up, now
deleted: capture is no longer a distinguishing feature, since 32 runs of the
current sweep capture a body (26 Mercury, 6 Venus) including the adopted run
itself. See §6.9.*

The current sweep's adopted run is
`…rp0p75__vinf25__inc30__toff59132__Om0__om30`, selected for the Mars
habitability window rather than for a capture:

| | value | percentile in sweep |
|---|---|---|
| Post-flyby `a` | 1.336 AU | — |
| Eccentricity | 0.217 | — |
| Obliquity | 23.8° | — |
| Global mean T | 217.5 K | 77th |
| Equatorial seasonal range | 53.3 K | 57th |
| Peak temperature anywhere | 278.7 K | 69th |
| **Liquid water, best latitude (−21.5°)** | **7.2% of the year** | **99th** |

**The selection is honest about what it selected for.** This run is unremarkable
on every axis except the one it was chosen for — 57th percentile for seasonal
range, 69th for peak temperature — and sits at the 99th percentile for the
liquid-water fraction. That is the correct shape for a deliberately chosen
scenario: it is an outlier in the target variable and typical in everything else,
which is a much weaker claim than "the flyby makes Mars habitable."

It also remains true that the sweep, not the run, is what licenses any general
statement. A single-run study would support "the flyby warms Mars while freezing
Earth"; across 3826 runs the median Mars outcome is 204.9 K and no run anywhere
achieves year-round liquid water. This is the false generalisation §9 warns
about, and the reason the worked scenario appears *after* the sweep results
rather than instead of them.

---

## 8. Corrections and negative results

Recorded deliberately — several initial claims did not survive testing.

### 7.1 Claims that were wrong

| Claim | Reality |
|---|---|
| "Obliquity change is a first-order effect" — based on one sampled run showing 23° → 74° | Full distribution: median change ~0; only **3%** of runs shift >10°, 0.3% >30°. The 74° case is an outlier. |
| "Adding two surfaces will require re-tuning `coalbedo_a0` and all results will shift" | **Not needed.** Global mean moved 288.14 → 288.42 K. Equilibrium is set by radiation balance, not heat capacity. |
| "The transient may land in a different attractor than `run_equilibrium` (which starts from a uniform 15 °C guess rather than Earth's real profile)" | **Tested 12 runs straddling the bifurcation: zero flips.** Equilibrium results are robust to initial conditions. |
| "The old model missed the glacial-inception signal" | **Regime-dependent.** True for pure eccentricity injection; across this sweep the single-surface model actually *over*-predicts summer cooling (63% vs 47%). |
| A hardcoded Kepler test constant of 1.1934205 | Recomputed independently: **1.1853242**. The original would have validated a broken solver. |
| Apparent `λ_p` → temperature correlation of +0.59 | **Spurious** — collinear with `a`. Controlling for insolation: −0.22. |
| "A second planet will need *more* approximation than Earth" | **Backwards.** Mars's radiation is *simpler*: emissivity 1.00 makes `σT⁴` exact rather than fitted, removing the Earth model's narrowest validity limit entirely. |
| A Mars test asserting the year's minimum temperature against a single end-of-year frost point | **The test was wrong, not the model.** The frost point *moves* as pressure falls; the minimum matched the frost point at minimum pressure exactly. Rewritten to compare step by step, plus an assertion that the frost point moves at all. |
| "float32 costs ~10⁻⁶ AU and ~10⁻⁶ degrees" (§11) | **Optimistic by 1–4 orders of magnitude.** Measured across 4032 runs whose pre-flyby elements are physically identical — so any spread is pure numerical error — the actual scatter is 1.8×10⁻⁵ AU in `a` and 3.1×10⁻² degrees in `λ_p`. The conclusion survives (still ~140× under the accepted precession uncertainty) but the quoted figure did not. It also downgraded the §5.3 zero-variance check from an equality to a tolerance. |
| `liquid_water_possible` tested only that the surface was above the water triple point | **Incomplete — it ignored boiling.** On a thin-atmosphere Mars the surface can sit above 273 K while the pressure is far below water's saturation pressure, so the water sublimates rather than pooling. The test now requires `T > 273.16 K` **and** `p > p_sat(T)`. This changed the answer by a large factor and **changed which runs won**: the best run under the old test showed 33% of the year with "liquid water", against 1% under the correct one. The weaker test is retained as `above_water_triple_point` for comparison, and the two still disagree — 23.0% of runs against 21.5% at the equator (§7.6). A unit test asserting the old behaviour had to be renamed rather than deleted, because it was asserting the bug. |

### 7.2 Genuine negative results

* **Temperature overshoot is exactly zero**, in both single- and two-surface
  modes. The approach to equilibrium is monotonic. This was one of the stated
  motivations for running transients and it simply does not occur.
  *(Measured across the retired sweep's 644 usable runs; not re-tested since —
  see §6.5. A qualitative result, but the evidence for it is not current.)*
* **Initial-condition bistability does not occur** in the tested range, despite
  the model having a genuine bifurcation. *(12 runs straddling the bifurcation,
  retired sweep.)*
* **A blended heat capacity is not a useful approximation** (9% effect) — worth
  recording so the cheap approach is not re-attempted. *(Retired sweep.)*
* **Obliquity barely affects Martian seasonal amplitude** (r = +0.03 in the
  current sweep, against +0.96 for eccentricity) — the reverse of Earth, and a
  reminder that a parameter's importance is a property of the system, not of the
  parameter.
* **No run in 3826 gives Mars year-round liquid water.** The best case is 7.8% of
  the year at one latitude. The flyby can open a habitability *window*; it cannot
  produce a habitable Mars. This is a negative result about the scenario itself,
  not about the model.
* **Capture is not a characteristic outcome.** 32 of 36,288 body-outcomes end up
  bound to the black hole — 1.5% of bodies that leave the Sun, 0.09% overall
  (§6.9). The first edition's worked example was a capture case, which
  overstated how typical that is.
* **The median climate outcome is not a robust quantity.** It moved from −9.2 K
  to +5.1 K between two sweeps of the same scenario family (§6.1) purely because
  median post-flyby `a` shifted by 0.05 AU. The *distribution* of outcomes is
  stable; its centre is not. Do not quote the median as a property of "a black
  hole flyby."

---

## 9. Bugs found and fixed

| Bug | Where | Consequence |
|---|---|---|
| Regex assumed XML attribute order (`Id` before `Target`); this writer emits the reverse | `extract_earth_elements.py` | Sheet lookup returned nothing; now parses XML properly |
| Positional cell parsing assumed a uniform 12 numeric columns; the BH sheet has 10 (no `a_tidal`) | `convert_orbits_to_parquet.py` | Would have silently shifted every value; caught by the time-column integrity check, now addressed by column letter |
| Double-blend: `ice_line_lat` re-blended an already-blended state | `ebm.py` | Shape error in two-surface mode; `blend()` made idempotent |
| `--limit` truncated the dataframe *before* statistics were computed | `climate_from_simulations.py` | Reported "654 outside validity band" when only 8 runs were requested |
| `ndarray.ptp()` removed in NumPy 2 | test suite | Test error; switched to `np.ptp()` |
| Hardcoded diagnostic-name list | `run_perturbation` | Land/ocean keys silently missing from transients |
| `rebound` 5.x renamed the `hash=` kwarg to `name=` | `solar_system_bh_rebound26.py` | All 1050 sweep combos failed instantly (pre-existing, fixed early) |
| Chained sweeps with `&&` | run scripts | A locked output file aborted the first pass and silently skipped the second |

---

## 10. Known limitations

1. **OLR parameterisation.** The default linear form `A + B·T` (Budyko 1969) is
   calibrated near 288 K and fails differently at each end:

   * **Cold end (qualitative failure).** It reaches `OLR = 0` at **175.9 K** and
     goes negative below. Between ~176 K and ~230 K it behaves *backwards*: a
     freezing, drying atmosphere should emit *more* (toward blackbody), but the
     linear form drives emission toward zero. It therefore under-cools, and
     equilibria in this range are **too warm**.
   * **Hot end (missing phase transition).** Real moist atmospheres cannot
     exceed the **Simpson–Nakajima ceiling** (~300 W/m²; Nakajima, Hayashi & Abe
     1992) — beyond it no equilibrium exists and the planet runs away. The
     linear form has no ceiling and crosses that threshold at ~319 K, i.e.
     `a < 0.876 AU` on a circular orbit.

   Defensible range: roughly **230–300 K**. Across the sweep that is only
   **35%** of runs; 37% sit in the severe-cold regime and **16% (102 runs) are
   in the runaway regime**, where the model reports temperatures for worlds that
   would have lost their oceans entirely.

   **Both failures are now addressed** (see MANUAL §10): `olr_model: sellers`
   selects the Sellers (1969) nonlinear form, which matches the linear one at
   288 K to 0.3% but tends correctly to blackbody emission when frozen; and
   every equilibrium now reports `absorbed_mean` plus a `runaway` boolean
   against the Simpson–Nakajima ceiling. The runaway condition is **flagged, not
   capped** — capping would manufacture a stable state that does not exist.
   Sellers does not fix the hot end (above ~350 K it degenerates to a constant
   emissivity `0.5·σT⁴`), which is precisely why the flag is needed alongside it.

   Results quoted elsewhere in this report were computed with the **linear**
   form and retain the limitation as stated.
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

### Mars model (§7)

8. **Fixed emissivity and fixed CO₂ inventory are one assumption, not two.**
   This was originally recorded here as two separate limitations. That was
   wrong, and the coupling matters.

   `ε = 1.0` is exact at 600 Pa and cannot represent a thickened atmosphere
   developing a real greenhouse — so it appears to be a limit that bites in the
   *hot* direction. But the inventory is fixed at 200 kg/m², which caps surface
   pressure at `200 × 3.71 = 742 Pa` — **1.24× present-day Mars**, and the
   sweep confirms no run exceeds it. Over that range a grey-gas greenhouse
   varies from +3.8 K to +4.7 K, so a pressure-dependent emissivity would change
   nothing.

   **Neither limit is reachable while the other holds.** Relaxing emissivity
   alone buys under 1 K; relaxing the inventory alone lets pressure rise while
   the greenhouse stays deaf to it. Only the pair does anything — and then it
   does a great deal: with the inventory free, mean temperature runs 208.6 K at
   200 kg/m², 228.8 K at 1000, and 289.5 K at 5000, because warmer means less
   condensation means higher pressure means warmer still.

   Both are therefore **unflagged limits of the same physical assumption**:
   that Mars's exchangeable CO₂ is small and its greenhouse negligible. Both are
   true of present-day Mars and neither is checked by the model.

   Explored on the `mars-thick-atmosphere` branch (not merged), which adds
   `olr_model="graygas"` and finds a genuine but *narrow* tipping point — two
   stable states spanning ~1.4% of insolation. See §8 for the failed first
   attempt to measure it.
9. **No regolith or escape reservoir.** CO₂ moves only between atmosphere and
   seasonal caps. Mars's permanent cap, adsorbed regolith CO₂, and loss to space
   are all absent, which is what pins the inventory in the first place.
10. **Atmospheric collapse is flagged but not modelled.** 86 runs (13%) freeze
    out at some point; the honest output there is "the atmosphere freezes out",
    not a temperature.
11. **No dust.** Real Martian albedo and opacity swing with global dust storms,
    which dominate its interannual variability.

---

## 11. Data volume and archival

> **Superseded by an engine change.** Everything below describes the retired
> sweep, when the engine wrote xlsx and conversion was a *post-hoc* rescue. The
> engine now writes Parquet directly (`orbits_format: parquet`, the default) and
> supports two-rate logging — a dense window around perihelion, coarse elsewhere —
> so the current 4032-run sweep never creates the 150 GB in the first place. The
> xlsx path and `convert_orbits_to_parquet.py` remain only for reading older
> sweeps. The analysis of *why* the format was so wasteful is retained because it
> is what motivated the engine change.

The engine wrote one `*__orbits__*.xlsx` per run at ~220 MB, making the 672-run
sweep **150.7 GB**. Almost all of that is format overhead, not information:
every number is stored as XML text, five `*_str` columns duplicate numeric ones
as text, and four more columns are derivable from the state vector.

[`convert_orbits_to_parquet.py`](convert_orbits_to_parquet.py) reduced this to
**26.4 GB — a 5.71× reduction** — with all 672 runs verified and zero failures.

**What is kept:** `t_days, x, y, z, vx, vy, vz, disp_helio_au`. Dropped as
recoverable: `r_helio` (= |r − r_sun|), `a_tidal` (computable from the BH state),
the `*_m`/`*_m_s2` unit conversions, and the five `*_str` duplicates.
`disp_helio_au` is **deliberately kept** — it is the displacement against the
BH-free baseline integration and cannot be reconstructed from this file alone.

**Reading strategy.** `openpyxl` needs ~13 s per sheet (29 h for the sweep).
Instead the sheet XML is stream-decompressed straight from the xlsx zip and
parsed with one regex: ~0.6 s per sheet, **21× faster**.

**Precision.** float32 was measured, not assumed: re-deriving Earth's orbital
elements from float32 state vectors differs from float64 by ~10⁻⁶ AU and
~10⁻⁶ degrees — negligible beside the ~4.3° λ_p precession uncertainty already
accepted in §5.4.

> **That measurement was optimistic.** Now that the engine writes float32 Parquet
> for the whole sweep, the spread can be read straight off the pre-flyby elements,
> which are physically identical across all 4032 runs and so isolate the numerical
> error exactly. Actual scatter: **1.8×10⁻⁵ AU** in `a` (18× the stated figure)
> and **3.1×10⁻² degrees** in `λ_p` (four orders of magnitude above it). The
> original measurement evidently sampled a favourable case. The *conclusion* is
> unaffected — 0.03° is still ~140× smaller than the accepted 4.3° precession
> uncertainty, and float32 remains the right choice — but the error bar is not
> what this section claimed. See §5.3.

**Three layers of verification**, all required before a run is reported `ok`:

1. **Structural** — row numbers contiguous, and the time column must equal
   exactly `0..N−1`. Any cell misalignment breaks this immediately.
2. **Round-trip** — every file is read back and compared to what was written,
   including NaN-pattern agreement.
3. **Scientific** — orbital elements re-derived from the Parquet reproduce the
   originals to ~10⁻⁶ AU / ~10⁻³ degrees.

Layer 1 caught a real bug (see §8): the first implementation parsed cells
positionally, assuming a uniform 12 numeric columns per sheet. The **BH sheet
has only 10** — it omits `a_tidal`, since the black hole exerts no tidal
acceleration on itself. Positional parsing would have shifted every subsequent
value and written plausible-looking garbage. The reader now addresses cells by
spreadsheet column letter.

**Archival.** The scientific content actually drawn on — per-run `input.yaml`,
`planets_run_deltas.csv`, and the derived analysis CSVs — is only a few MB and
belongs in git. The trajectories are bulk intermediate data: compress, push to
object storage or a citable archive (Zenodo issues a DOI), and keep the local
copy cloud-only. They are also fully **regenerable** from the recipe
(`input.yaml` + engine git SHA + `de440s.bsp`), so storing them at all is a
compute-vs-storage tradeoff rather than a necessity.

**Prevention.** `output_interval_days: 1` is what produces 220 MB per run. A
future sweep will recreate the problem unless the cadence is coarsened away from
the encounter or the engine writes Parquet directly.

---

## 12. File inventory and reproduction

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
| `tests/` | **61 tests** |
| `mars.py` | Condensing-CO2 atmosphere (`MarsEBM`): frost point, latent buffering, mass exchange, collapse flag |
| `MANUAL.md` | Usage manual |

### Bridge scripts — repository root

| File | Role |
|---|---|
| `extract_earth_elements.py` | Recovers `a, e, ε, λ_p` per run from the orbit logs (reads Parquet directly; the 30×-faster xlsx path in §5.5 is retained for older sweeps) |
| `find_bh_captures.py` | Separates bodies captured by the BH from those merely ejected (§6.9) |
| `climate_from_simulations.py` | Runs the EBM on every simulation; equilibrium or transient; parallel |
| `rank_run_impact.py` | Ranks runs by orbital disruption (energy/eccentricity/perihelion metric) |
| `convert_orbits_to_parquet.py` | Converts the orbit workbooks to Parquet (5.7x smaller), with structural, round-trip and scientific verification |
| `extract_mars_elements.py` | Element recovery generalised to any planet via an IAU pole table |
| `climate_mars_from_simulations.py` | Runs the Mars model across a sweep |
| `input_climate.yaml` | Shared configuration (Earth) |
| `input_mars.yaml` | Mars configuration |

### Reproduction

The current sweep is `<STAMP>` = `20260811_184731`. Steps 2–5b below were run
against it and produced the `simulations/20260811_184731_*.csv` files this report
quotes — **with default options** (linear OLR, single surface, equilibrium). The
`--olr-model sellers --two-surface` and `--transient` variants shown have **not**
been run against the current sweep — running them is exactly what §6.5 asks for.

```bash
# 0. Dependencies
python -m pip install numpy scipy pandas matplotlib pyyaml pyarrow pytest

# 1. Validate the physics
python -m pytest orbital_climate/tests/ -v            # 61 tests

# 2. Rank runs by orbital disruption
python rank_run_impact.py simulations/<STAMP> --plot impact.png

# 2b. Capture census — which bodies leave *with* the black hole (§6.9)
python find_bh_captures.py simulations/<STAMP>

# 3. Recover Earth's post-flyby orbital elements (~10 min, 5 workers)
python extract_earth_elements.py simulations/<STAMP> --workers 5

# 4. Climate across the sweep — equilibrium
#    --olr-model sellers keeps frozen states physical; omit it for the
#    original linear (Budyko) form.
python climate_from_simulations.py simulations/<STAMP> \
    --config input_climate.yaml --olr-model sellers --two-surface \
    --workers 5 --plot climate.png

# 5. Climate across the sweep — transient
python climate_from_simulations.py simulations/<STAMP> \
    --config input_climate.yaml --olr-model sellers --two-surface \
    --transient --years 40 --workers 5

# 5b. Mars: recover its elements (any planet via --body), then run its climate
python extract_mars_elements.py simulations/<STAMP>_parquet --body Mars --workers 5
python climate_mars_from_simulations.py simulations/<STAMP>     --config input_mars.yaml --workers 5

# 6. (Storage) Compress the orbit workbooks ~5.7x, with verification.
#    Source files are never deleted; remove them yourself once backed up.
python convert_orbits_to_parquet.py simulations/<STAMP> --workers 5
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

## 13. References

**A caveat on this list.** These citations were compiled from working knowledge
rather than by consulting the papers, so **volume and page numbers in particular
should be verified against the primary sources** before any of this is quoted
elsewhere. Author, year and title are the parts to trust; the rest is a pointer.
This matters more than usual here, because §8 records a case where a plausible
number written down from memory would have certified a broken solver.

### Energy-balance climate models

* **Budyko, M. I.** (1969). *The effect of solar radiation variations on the
  climate of the Earth.* Tellus **21**, 611–619.
  — The linear OLR parameterisation `A + B·T` used throughout the Earth model,
  and the ice-albedo feedback.
* **Sellers, W. D.** (1969). *A global climatic model based on the energy balance
  of the earth–atmosphere system.* Journal of Applied Meteorology **8**, 392–400.
  — The nonlinear OLR form implemented as `olr_model="sellers"`.
* **North, G. R., & Coakley, J. A.** (1979). *Differences between seasonal and
  mean annual energy balance model calculations of climate and climate
  sensitivity.* Journal of the Atmospheric Sciences **36**, 1189–1204.
  — The two-surface land/ocean formulation of §4.
* **North, G. R., Cahalan, R. F., & Coakley, J. A.** (1981). *Energy balance
  climate models.* Reviews of Geophysics **19**, 91–121.
  — The benchmark review; source of the diffusion formulation and standard
  parameter values.

### Radiation and the limits of parameterised OLR

* **Koll, D. D. B., & Cronin, T. W.** (2018). *Earth's outgoing longwave
  radiation linear due to H₂O greenhouse effect.* PNAS **115**, 10293–10298.
  — Explains *why* Earth's OLR is quasi-linear at all, and where that breaks.
  The reference for §10.1's validity analysis.
* **Nakajima, S., Hayashi, Y.-Y., & Abe, Y.** (1992). *A study on the "runaway
  greenhouse effect" with a one-dimensional radiative–convective equilibrium
  model.* Journal of the Atmospheric Sciences **49**, 2256–2266.
  — The Simpson–Nakajima OLR ceiling used for the runaway flag.
* **Goody, R. M., & Yung, Y. L.** (1989). *Atmospheric Radiation: Theoretical
  Basis*, 2nd ed. Oxford University Press.
  — Grey-gas radiative transfer, the basis of `olr_model="graygas"`.
* **Pierrehumbert, R. T.** (2010). *Principles of Planetary Climate.* Cambridge
  University Press.
  — Source of the grey-atmosphere relation `T_surf⁴ = T_eff⁴(1 + 3τ/4)`, and of
  much of the framing in §7.

### Mars: the CO₂ cycle

* **Leighton, R. B., & Murray, B. C.** (1966). *Behavior of carbon dioxide and
  other volatiles on Mars.* Science **153**, 136–144.
  — The foundational paper. Predicted the seasonal CO₂ condensation cycle, and
  that polar temperatures would be buffered at the frost point, *before* it was
  observed. This is the physics implemented in `orbital_climate/mars.py`.
* **Gierasch, P. J., & Toon, O. B.** (1973). *Atmospheric pressure variation and
  the climate of Mars.* Journal of the Atmospheric Sciences **30**, 1502–1508.
  — Directly relevant to the `mars-thick-atmosphere` branch: identifies the
  pressure–greenhouse feedback and the possibility of multiple stable
  atmospheric states.
* **James, P. B., Kieffer, H. H., & Paige, D. A.** (1992). *The seasonal cycle of
  carbon dioxide on Mars.* In *Mars* (Kieffer et al., eds.), University of
  Arizona Press, 934–968.
  — Standard reference for the observed cycle and for the CO₂ vapour-pressure
  relation used to obtain the frost point.
* **Kieffer, H. H., et al.** (1976). *Infrared thermal mapping of the Martian
  surface and atmosphere: First results.* Science **193**, 780–786.
  — Viking observations establishing the ~148 K winter polar cap temperature
  against which the latent-buffering result is validated.
* **Soto, A., Mischna, M., Richardson, M., et al.** (2015). *Martian atmospheric
  collapse: Idealized GCM studies.* Icarus **250**, 553–569.
  — Atmospheric collapse in a full GCM; the phenomenon the collapse flag detects.

### Mars: thick atmospheres and early climate

* **Kasting, J. F.** (1991). *CO₂ condensation and the climate of early Mars.*
  Icarus **94**, 1–13.
* **Forget, F., & Pierrehumbert, R. T.** (1997). *Warming early Mars with carbon
  dioxide clouds that scatter infrared radiation.* Science **278**, 1273–1276.
* **Wordsworth, R. D.** (2016). *The climate of early Mars.* Annual Review of
  Earth and Planetary Sciences **44**, 381–408.
  — Review; the context for the thick-atmosphere results on the experimental
  branch, and for why a single grey optical depth is inadequate there.

### Orbital forcing and Milankovitch

* **Ward, W. R.** (1974). *Climatic variations on Mars: 1. Astronomical theory of
  insolation.* Journal of Geophysical Research **79**, 3375–3386.
  — Orbital forcing of Martian climate; counterpart to the Earth Milankovitch
  literature.
* **Laskar, J., Correia, A. C. M., Gastineau, M., et al.** (2004). *Long term
  evolution and chaotic diffusion of the insolation quantities of Mars.* Icarus
  **170**, 343–364.
  — Mars's chaotic obliquity, and why the 38° obliquity found in §7.7 is not
  exotic on long timescales.
* **Kasting, J. F., Whitmire, D. P., & Reynolds, R. T.** (1993). *Habitable zones
  around main sequence stars.* Icarus **101**, 108–128.
* **Kopparapu, R. K., et al.** (2013). *Habitable zones around main-sequence
  stars: New estimates.* The Astrophysical Journal **765**, 131.
  — The runaway-greenhouse and CO₂-condensation limits framing §10.1.

### Reference frames and planetary data

* **Archinal, B. A., et al.** (2018). *Report of the IAU Working Group on
  Cartographic Coordinates and Rotational Elements: 2015.* Celestial Mechanics
  and Dynamical Astronomy **130**, 22.
  — Source of the IAU pole orientations in `extract_mars_elements.py`; Mars's
  pole at α₀ = 317.68143°, δ₀ = 52.88650° is what makes the recovery of §7.5
  possible.
* **Putzig, N. E., & Mellon, M. T.** (2007). *Apparent thermal inertia and the
  surface heterogeneity of Mars.* Icarus **191**, 68–94.
  — Regolith thermal inertia, the basis for `heat_capacity = 1.2 × 10⁶`
  (τ ≈ 6.6 days).
* **NASA/GSFC Planetary Fact Sheets.**
  https://nssdc.gsfc.nasa.gov/planetary/factsheet/
  — Bulk orbital and physical parameters used as calibration targets.

### Numerical methods and the N-body simulation

* **Rein, H., & Liu, S.-F.** (2012). *REBOUND: An open-source multi-purpose
  N-body code for collisional dynamics.* Astronomy & Astrophysics **537**, A128.
* **Rein, H., & Spiegel, D. S.** (2015). *IAS15: A fast, adaptive, high-order
  integrator for gravitational dynamics, accurate to machine precision over a
  billion orbits.* MNRAS **446**, 1424–1437.
  — The integrator behind the flyby sweep.
* **Ascher, U. M., Ruuth, S. J., & Wetterton, B. T. R.** (1995).
  *Implicit–explicit methods for time-dependent partial differential equations.*
  SIAM Journal on Numerical Analysis **32**, 797–823.
  — The IMEX splitting of §3.3, including retaining a linear term implicitly as
  a preconditioner while carrying the nonlinear remainder explicitly.

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
| Extraction deterministic | pre-flyby values across 672 runs (xlsx/float64) | zero variance |
| …the same on the current sweep | 4032 runs (Parquet/float32) | agrees to 1.8×10⁻⁵ AU — a tolerance, not an equality (§5.3) |
| Bound/unbound classification | vs independent deltas CSVs | 662/672 and 3966/4032 both ✓ |
| Land fraction profile | area-weighted global mean | 0.290 (Earth 0.29) |
| Land/ocean calibration | 65 °N seasonal ranges | 40.1 K / 9.1 K vs ~40 / ~8–9 |
| Two surfaces preserve energetics | global mean shift | 288.14 → 288.42 K |
| Equilibrium robust to initial state | 12 runs at the bifurcation | zero flips |
| **Mars** — element recovery | applied at `t = 0` | a = 1.5237, e = 0.0934, ε = 25.18° |
| **Mars** — grey-body OLR justified | absorbed vs `σ(210 K)⁴` | 110.4 vs 110.3 W/m² → ε = 1.00 |
| **Mars** — latent buffering | polar winter minimum | 146.8 K vs observed ~148 K |
| **Mars** — CO₂ conservation | through a full year | exact (< 10⁻¹⁰) |
| **Mars** — frost point | inversion round-trip, and at 600 Pa | exact; 148 K vs observed ~148 K |
| **Mars** — seasonal cycle | pressure, swing, cap thickness | 596 Pa / 22% / 0.57 m vs ~600 / ~25% / ~0.5–1 m |
| **Mars** — reduces to base model | condensation suppressed | identical to `EBM` step for step |
| Linear OLR hard floor | solve `A + B·T = 0` | 175.9 K; negative below |
| Sellers matches linear at present day | point comparison at 288 K | 235.1 vs 234.3 W/m² (0.3%) |
| Sellers physical when frozen | Sellers/blackbody at 180–220 K | 0.89–0.97 (linear collapses to 0.14) |
| Sellers reproduces 288 K | `sellers_m` scan | 288.04 K at m = 0.51 |
| Runaway flag fires correctly | Earth moved inward | triggers between a = 0.90 and 0.85 AU (analytic 0.876) |
| Parquet preserves the science | elements re-derived from converted files | ~10⁻⁶ AU, ~10⁻³ deg — **optimistic, see §8/§11**; full-sweep measurement gives 1.8×10⁻⁵ AU, 3.1×10⁻² deg |
| Parquet conversion integrity | 672 runs, 3 verification layers | 672 ok, 0 failed |
| **Capture vs ejection census** | 36,288 body-outcomes, two-body energy test | 32 captured, 2140 free (§6.9) |
| **Mars liquid water** | triple point *and* boiling, 3826 runs | max 7.8% of year; 0 runs year-round (§7.6) |
