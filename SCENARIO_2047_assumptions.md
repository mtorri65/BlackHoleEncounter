# The 2047 scenario — assumptions and their provenance

Design record for the black-hole flyby sweep with **periapsis in July 2047**,
superseding the 2027 configuration used for `simulations/20260724_230314`.

> **State of this document.** Current for the configuration below
> (epoch 1885-09-01, `v_inf` = 25 km/s), run as
> `simulations/20260811_184731` — 4032 runs, complete. Sections 5 and 8 still
> quote measurements taken from superseded sweeps; the mechanisms they describe
> are unchanged and confirmed, but the specific tabulated values have not been
> re-measured against the current run. They are marked in place.

Every parameter below is either a free choice or is *forced* by a stated
constraint. This document records which is which, so that later work can tell
what may be varied freely and what cannot be changed without redoing an
analysis. Numbers are reproducible with
[`astrometric_detectability.py`](astrometric_detectability.py); the command that
produces each table is given with it.

Companion documents: [`SCENARIO_mars_window.md`](SCENARIO_mars_window.md) — the
chosen run from this sweep, studied in detail — and
[`SCENARIO_timeline.md`](SCENARIO_timeline.md), its narrative century, and
[`SCENARIO_post_flyby_system.md`](SCENARIO_post_flyby_system.md), which works out
what the surviving system looks like.

---

## 1. The configuration

```yaml
epoch: "1885-09-01T00:00:00Z"    # chosen: matches the constraint date below
bh_mass_msun: 0.1                # imposed
bh_vinf_kms: 25.0                # DERIVED -- see section 4
bh_tperi_offset_days: 59132      # imposed: 1885-09-01 + 59132 d = 2047-07-26
bh_rp_au: "0.25, 1.5, 0.25"      # imposed (swept), 6 values
bh_inc_deg: "-90, 90, 30"        # 7 values
bh_Omega_deg: "0, 270, 90"       # 4 values
bh_omega_deg: "0, 345, 15"       # 24 values -- DERIVED, see section 10
duration_days: 115342            # chosen: preserves the 2027 sweep's aftermath
```

The four imposed constraints were:

1. Uranus's orbital residuals **barely recognizable in 1885**
2. BH mass = 0.1 M☉
3. periapsis July 2047
4. periapsis distance swept over 0.25–1.50 AU

Constraint 1 is an *observability* constraint, not a dynamical one, and it is
the only reason `v_inf` is not free. Everything downstream of it inherits that
choice — see section 6 for what changes if it is relaxed.

`duration_days: 115342` gives 56,210 days after periapsis, identical to the 2027
sweep, so climate results remain directly comparable.

The epoch is set to 1885-09-01 so that t = 0 coincides with the date the
observability constraint is evaluated. Nothing physical depends on it — the
epoch only selects which real planetary configuration the integration starts
from — but it removes a standing source of confusion, since the BH's distance at
any date depends on the epoch through `bh_tperi_offset_days`.

---

## 2. Method: what "detectable" means here

`python astrometric_detectability.py horizon`

Two N-body integrations, identical but for the BH, differenced. Sun + 8 planets
from de440s at the epoch, IAS15.

Three deliberate choices, each of which materially affects the answer:

**The observable is angular, not radial.** 19th-century meridian-circle work
measures directions — right ascension and declination — and nothing else. There
is no ranging until radar. So the residual counted here is the angle between a
planet's heliocentric direction with and without the BH; the radial component of
the perturbation is invisible and is excluded.

**Residuals are detrended.** An astronomer facing unexplained residuals refits
the orbit. A constant offset is absorbed into the mean longitude at epoch and a
linear trend into the mean motion; neither is evidence of a perturber. Only what
survives removal of constant + linear counts. This is conservative in the
astronomer's favour: a full six-element refit would absorb *more*, so these
figures are an upper bound on what a real analysis would have retained.

**The baseline is 30 years.** Long enough for a secular signature to separate
from a refit, short enough to be a plausible single campaign.

Both detrended and raw values are reported by the tool; the raw ones ignore
refitting and should be read as an upper bound only.

---

## 3. Detection horizon: Uranus is the limiting planet

Peak detrended angular residual over a 30-year window, by the BH's distance at
the window's midpoint:

| BH distance | Mercury | Earth | Jupiter | Saturn | **Uranus** | Neptune |
|---:|---:|---:|---:|---:|---:|---:|
| 1519 AU | 0.0004″ | 0.0004″ | 0.011″ | 0.022″ | **0.067″** | 0.048″ |
| 941 AU | 0.002″ | 0.002″ | 0.046″ | 0.086″ | **0.277″** | 0.204″ |
| 741 AU | 0.004″ | 0.004″ | 0.094″ | 0.168″ | **0.561″** | 0.419″ |
| 579 AU | 0.009″ | 0.008″ | 0.200″ | 0.333″ | **1.16″** | 0.890″ |
| 434 AU | 0.021″ | 0.020″ | 0.490″ | 0.730″ | **2.72″** | 2.16″ |
| 333 AU | 0.051″ | 0.048″ | 1.16″ | 1.50″ | **5.98″** | 5.00″ |
| 246 AU | 0.144″ | 0.139″ | 3.33″ | 3.42″ | **14.8″** | 13.5″ |
| 173 AU | 0.565″ | 0.585″ | 13.2″ | 19.0″ | **44.3″** | 46.8″ |

Reading off Uranus: **1″ at ~610 AU, 3″ at ~420 AU, 30″ at ~196 AU.**

> **Read this table for its shape, not its absolute scale.** It was computed at
> `v_inf` = 25 with each window anchored at the epoch, and §4 shows the detrended
> metric moves by more than an order of magnitude with window position and
> length. The *relative* sensitivities across planets — which is what this
> section is about — are robust; the absolute arcsecond values are not.

**Why Uranus and not an inner planet.** A distant perturber acts as a tidal
quadrupole growing as $a^3$, while the orbital frequency that converts it into
an accumulated angle falls as $a^{-3/2}$; net sensitivity goes as $a^{3/2}$.
Mercury is ~100× less sensitive than Uranus — nowhere near enough to be rescued
by its far superior observations (transit timings).

**Why Uranus and not Neptune.** Neptune's raw signal is comparable, but a
30-year window covers only 18% of its orbit, so most of the perturbation is
degenerate with a mean-motion correction and is refitted away: its
detrended/raw ratio is 0.25 against Uranus's 0.94. Historically the case is
stronger still — in 1885 Neptune had been known for 39 years, and its orbit was
not determined well enough to argue about arcseconds. Uranus had been observed
since 1781 with prediscovery positions back to 1690, and its residuals were the
instrument that found Neptune.

**Threshold calibration.** Two anchors from the period. Le Verrier established
Mercury's perihelion anomaly at 38″/century in 1859, so a well-founded secular
anomaly of that scale was demonstrably within reach. The Uranus residuals that
drove the Neptune prediction reached ~30″ in the 1830s and ~120″ by 1845. A ~1″
detrended residual is therefore *heroic but not absurd*; ~30″ is unmistakable.

**The real 1885 limit was theory, not measurement.** To call a residual
anomalous, the gravitational model of the known planets must be better than the
residual. Le Verrier's tables (1858–1877) carried systematic errors of order an
arcsecond for the outer planets, with Neptune's own mass uncertain; Newcomb and
Hill's rigorous theory does not arrive until ~1895. A defensible practical
threshold is therefore a few arcseconds rather than one, which would move the
horizon inward to ~420 AU and lower the derived `v_inf` (section 4).

---

## 4. Deriving `v_inf` from the 1885 constraint

`python astrometric_detectability.py vinf --epoch 1885-09-01 --peri 2047-07-26 --when 1885`

With periapsis pinned to 2047-07-26, 1885 sits 162 years before arrival. That
fixes the travel time, so requiring a particular residual in 1885 fixes the
distance then, which fixes the approach speed.

**The answer depends strongly on how long an observing campaign is assumed**,
and that dependence is the single most important thing in this section. Uranus
detrended residual as of 1885, rp = 0.5 AU:

| `v_inf` | r(1885) | 20-yr | 30-yr | 60-yr | **100-yr** |
|---:|---:|---:|---:|---:|---:|
| 10 km/s | 375 AU | 0.41″ | 3.61″ | 8.15″ | **9.12″** |
| 15 km/s | 534 AU | 0.27″ | 1.35″ | 2.46″ | **3.22″** |
| 20 km/s | 698 AU | 0.16″ | 0.59″ | 0.84″ | **1.45″** |
| **25 km/s** | **866 AU** | 0.10″ | 0.28″ | 0.54″ | **0.76″** |
| 30 km/s | 1035 AU | 0.06″ | 0.14″ | 0.39″ | **0.45″** |

**The 100-year column is the defensible one.** Uranus was discovered in 1781 and
has prediscovery positions back to 1690, so an analyst in 1885 had roughly a
century of usable data, not thirty years. Read down it: **`v_inf` = 25 km/s puts
the residual at 0.76″**, just below the ~1″ that a Le Verrier-class analysis
might have claimed and well below what the era's planetary theory could have
defended. That is the constraint met.

At `v_inf` = 25 the BH sits at **863 AU at the 1885 epoch**.

### The metric is fragile, and this is a real limitation

The detrended residual is **not** a smooth function of distance. Holding
everything else fixed and varying only the window:

* **Window position**: a 30-year campaign slid across 1850–1935 gives anything
  from 0.42″ to 8.54″ — **a factor of 20** — while the raw (undetrended) signal
  climbs smoothly and monotonically as D⁻³ throughout.
* **Window length**: at a fixed 1885 end date, 20 → 100 years spans 0.41″ to
  9.12″ at `v_inf` = 10 — **a factor of 22**.
* **Epoch**: the same query at the retired 1873 epoch gives 0.03″ where the 1885
  epoch gives 0.28″.

The cause is the detrending, not the physics. Removing a constant and a linear
term is a proxy for an orbit refit, and Uranus's 84-year period means a campaign
of a few decades covers a fraction of an orbit over which the perturbation can
look nearly linear — in which case a refit absorbs almost all of it — or
strongly curved, in which case it survives. The nulls are refit degeneracies.

**Consequence for this section.** An earlier version of this document derived
`v_inf` = 10 km/s from a 30-year campaign ending in 1885, which happens to sit
near one of those nulls. That was not wrong arithmetic; it was an arbitrary
choice carrying the entire answer. The 100-year baseline is better justified on
historical grounds, and it gives 25 km/s.

**What would settle it properly** is replacing constant+linear detrending with a
genuine multi-element orbit refit — six elements plus planetary masses, fitted
to synthetic observations. That removes the degeneracy nulls and would give a
number that does not move when the window does. Until then, treat `v_inf` = 25
as a defensible choice rather than a derived value, and note that the whole
scenario is insensitive to it: §5 shows the encounter barely changes.

**Periapsis distance does not interact with this constraint.** Across
rp = 0.25 → 1.50 AU the required `v_inf` moves by well under 0.1 km/s, because
1.5 AU is nothing against 863 AU. rp may be swept freely without revisiting any
of the above.

---

## 5. Consequences for the encounter itself

> **PENDING REGENERATION.** The measurements below come from the superseded
> `v_inf` = 10 / epoch 1873 configuration. The *conclusion* — that gravitational
> focusing makes the encounter nearly independent of `v_inf` — is what justifies
> not worrying about §4's fragility, and it holds a fortiori at 25 km/s. The
> specific numbers need re-measuring against the current sweep.

`v_inf` = 10 km/s was expected to produce a far more violent flyby, since the
impulse approximation gives Δv ∝ 1/v. **It does not**, and the reason is worth
recording: at these periapsis distances gravitational focusing dominates, so the
BH's speed at periapsis is set by the Sun's potential rather than by `v_inf`.

| rp | v_peri at `v_inf`=10 | at `v_inf`=25 | ratio |
|---:|---:|---:|---:|
| 0.25 AU | 88.9 km/s | 91.8 km/s | 1.03 |
| 0.50 AU | 63.3 km/s | 67.3 km/s | 1.06 |
| 1.00 AU | 45.3 km/s | 50.8 km/s | 1.12 |
| 1.50 AU | 37.4 km/s | 43.9 km/s | 1.17 |

Confirmed against actual runs (rp = 0.5, inc = 0, Ω = ω = 0, periapsis 2047):

| body | `v_inf`=25: a, e | `v_inf`=10: a, e |
|---|---|---|
| Mercury | 0.319, 0.599 | 0.359, 0.274 |
| Venus | 0.848, 0.120 | 0.764, 0.074 |
| Earth | 1.292, 0.161 | 1.258, 0.132 |
| Mars | 1.197, 0.376 | 1.089, 0.527 |
| Jupiter | 3.951, 0.259 | 4.562, 0.141 |
| Saturn | **ejected** | **ejected** |
| Uranus | 17.88, 0.489 | 16.43, 0.350 |
| Neptune | 18.71, 0.705 | 17.51, 0.753 |

Damage is comparable in magnitude; the differences are geometric detail, not a
systematic strengthening. **Results from the 2027 sweep therefore remain broadly
representative** of the 2047 one at the same rp, which is a useful property.

---

## 6. Astrophysical plausibility — the weakest point in the scenario

`python astrometric_detectability.py rates`

A 0.1 M☉ black hole is sub-Chandrasekhar and cannot form by stellar collapse, so
it would have to be primordial — i.e. a dark-matter constituent, with halo
kinematics of 150–250 km/s relative to the Sun. The constraint above demands
10 km/s. **This is a real tension and it should be stated in any write-up.**

**No deceleration mechanism resolves it.** Time to change v by ~200 km/s:

| mechanism | timescale |
|---|---|
| dynamical friction, halo (ρ = 0.01 M☉/pc³) | 2.4×10⁸ Hubble times |
| dynamical friction, disk midplane (ρ = 0.1) | 2.4×10⁷ Hubble times |
| dynamical friction, inside a dense GMC (ρ = 100) | 2.4×10⁴ Hubble times |
| Bondi–Hoyle gas drag, ISM | 8.1×10⁸ Hubble times |

Chandrasekhar friction scales as 1/M; a Hubble-time brake needs ~10⁶ M☉ in the
disk. Gravitational-wave losses are irrelevant for a free-flying object. Capture
routes — into the Sun's birth cluster, or ejection from a nearby binary — do not
help either, because a cluster with a few km/s escape velocity cannot capture
something arriving at 250 km/s: those mechanisms *preserve* slowness, they do
not create it.

**The framing is wrong, though.** 10 km/s is relative to the *Sun*, not the
Galaxy. Nothing needs to decelerate; the BH's galactic orbit merely needs to
resemble the Sun's. That happens by **selection**, not by dynamics:

- **1 in 43,000** halo objects already has |v_rel| < 10 km/s (Standard Halo
  Model, σ = 156 km/s, Sun at 232 km/s).
- Slow encounters are *not* strongly disfavoured in rate. Low flux (×0.04) is
  nearly cancelled by gravitational focusing (×36 at rp = 0.5 AU): the rate
  ratio of 10 vs 270 km/s encounters is 1.29 at rp = 0.5 AU and 0.47 at 1.5 AU.

**Absolute rarity, stated plainly.** Even granting that 0.1 M☉ PBHs are 100% of
dark matter, encounters within 1.5 AU of the Sun run at **4.7×10⁻⁵ per 10 Gyr**
at any speed, and **5.1×10⁻¹⁰** once `v_inf` < 10 km/s is required. Microlensing
surveys already cap 0.1 M☉ PBHs at ~10–20% of dark matter, costing another
factor of ~10. The scenario was astronomically improbable before this constraint;
the velocity requirement makes it ~10⁵× more so.

### Three ways to respond

**(a) Accept it as a tail event.** Internally consistent, and the rate argument
above shows slow encounters are not additionally penalised once one occurs.

**(b) Drop the halo assumption** — the cleanest fix. A kinematically cold disk
population gives 10 km/s naturally: Gliese 710 will pass the Sun at ~14 km/s,
and the Sun's own peculiar velocity relative to the local standard of rest is
only ~18 km/s. The cost is giving up the identification with standard dark
matter; one would need a "dark disk", or simply decline to say what the object
is.

**(c) Relax the 1885 constraint.** This is the lever that created the problem,
and it is a narrative choice rather than physics. "Invisible in 1885" instead of
"barely recognizable" allows `v_inf` = 25 km/s (0.03″ in Uranus, comfortably
undetectable), restores 926 AU at epoch, and removes the kinematic tension
entirely.

**What will not work: a bound companion.** Free-fall from 400 AU takes ~707
years, but the 1885→2047 window allows only 162. At 400 AU the local escape
speed is 2.1 km/s while the BH moves at 10.2 — nearly five times escape. The
1885 and 2047 constraints *together force* a genuinely unbound flyby.

---

## 7. Engine convention discovered while checking this

**The BH's orbit is the point inversion of the orbit its labels describe:**
r → −r, v → −v. Verified exact to machine precision against a textbook
construction, at i = 0°, ±30°, −60° and 75°.

The cause is two lines in `solar_system_bh_rebound26.py`. The standard
hyperbolic radius is `r = a(1 − e cosh F)` with `a < 0`, which is positive. The
engine writes `r = a(e cosh F − 1)`, which has the right magnitude and the wrong
sign, so the position comes out antipodal. It then applies
`v_bh = -(R @ v_pqw)  # inbound branch` — the comment records that the BH was
observed leaving rather than arriving, and the fix was to flip the velocity.
Both vectors therefore end up negated.

**Nothing is broken.** The two-body equation `r̈ = −μ r/|r|³` is odd in `r`, so
if `r(t)` is a solution then so is `−r(t)`; point inversion maps orbits to
orbits. Under it:

- **h = r × v is unchanged** → the orbital plane, and hence `i` and `Ω`, are
  exactly as labelled
- **e_vec = (v × h)/μ − r̂ is negated** → periapsis lies on the opposite side
- `a`, `e`, `rp`, `v_inf`, `|r|`, `|v|` and the time of periapsis are all as
  requested

So the encounter has the right plane, shape, size and timing. What differs is
that the BH approaches from the opposite direction and rounds the far side of
the Sun.

In classical elements, for a labelled `i ≥ 0`, this reads **true ω =
`bh_omega_deg` + 180°** — measured directly from the eccentricity vector of the
engine's own initial state:

| `bh_omega_deg` | true ω (at i = 30°, Ω = 45°) |
|---:|---:|
| 0 | 180° |
| 90 | 270° |

For a *negative* labelled inclination the simple "+180°" no longer describes the
extracted number, because standard element extraction independently remaps the
plane: `i = −60°, Ω = 120°` comes back as `i = +60°, Ω = 300°`, and the `ω`
bookkeeping tangles with that remapping. The vector statement — **h preserved,
e_vec flipped** — is exact for every inclination and is the one to rely on.

**Benign for the sweeps.** The swept set {0, 90, 180, 270} is closed under
+180°, so the collection of orbits actually integrated is exactly the collection
intended; only the labels are permuted (the folder named `om0` holds the ω = 180°
orbit). Sweep statistics, the impact ranking and all climate results are
unaffected. Relative claims between two runs also survive, since every run
receives the same transformation.

### 7a. The sweep is exactly twofold redundant

Separately from the labelling issue above, the grid contains a genuine symmetry:

> **`(i, Ω, ω)` and `(−i, Ω+180°, ω+180°)` are the same physical orbit.**

The orbit normal is `(sin i sin Ω, −sin i cos Ω, cos i)`; negating `i` and adding
180° to `Ω` flips both signs in each of the first two components, leaving the
normal unchanged. Rotating the node and the argument of periapsis by 180°
together leaves periapsis pointing the same way in space. Same plane, same
direction, same apsidal orientation — the same encounter.

**Every one of the 4032 runs has such a partner**, so the sweep covers about
**2016 distinct configurations, each integrated twice.** Consequences:

- **Roughly half the compute is redundant** (~14 of ~28 hours). Constraining
  `Ω ∈ [0°, 180°)` or `i ≥ 0` would halve a future sweep at no cost in coverage.
- **Counts and percentages have an effective sample size of ~2016.** Means and
  medians are unbiased, but every tally quoted over "4032 runs" — 3904 usable,
  810 Earth-like, 32 captures — rests on half the independent information it
  appears to. The 32 captures are 16 configurations.

**The pairs are also an accidental convergence test.** Both members should be
identical, so their divergence measures the integration's numerical chaos
directly. Measured on the impact `Score`:

| | relative divergence |
|---|---|
| identical to 1e-12 | 4.1% of pairs |
| median | 7.5×10⁻¹⁰ |
| 90th percentile | 6.3×10⁻⁷ |
| max | 6.8×10⁻² |

The median run holds ~9 significant figures over 316 years, so individual
outcomes are sound. But divergence concentrates in the violent runs — the top
`Score` quartile diverges ~10× more than the quietest — and a few pairs disagree
at the percent level. **Treat any absolute claim about a single high-`Score` run
with that error bar in mind.**

*(A now-deleted document, `SCENARIO_mercury_capture.md`, noted this symmetry for
one run of the retired sweep and described the twins as "bit-identical". That is
true only for the 4.1% that agree to 1e-12: the mirrored initial state differs at
roundoff, and chaos amplifies it.)*

**Where it bites:** any absolute geometric statement about a *single named run*
— "in run `…om270` the BH approached from direction X". Note that plots reading
positions out of the simulation output (`plot_sky_tracks.py`,
`plot_local_sky.py`) show the true geometry and are correct; it is only the
mapping from a folder's ω label to that geometry that is off.

Left unchanged deliberately: correcting it would silently relabel every existing
run folder.

---

## 8. Nominal parameters are labels, not achieved geometry

> **PENDING REGENERATION** for the current configuration. The mechanism (below,
> under *Cause*) is exact and epoch-independent; the tabulated achieved values
> are from superseded sweeps. Three spot checks at `v_inf` = 25 with the 1885
> epoch, nominal periapsis 2047-07-26: i=0/Ω=0/ω=0 → 0.600 AU on 2047-08-27;
> i=30/Ω=0/ω=0 → 0.609 AU on 2047-08-29; i=60/Ω=270/ω=180 → **0.254 AU** on
> 2047-07-20. The spread is wider than at 10 km/s, as expected — the BH starts
> at 863 AU rather than 400, and the angular-momentum error scales as r·δv.

**`bh_rp_au` and `bh_tperi_offset_days` describe the BH's initial two-body orbit
at the epoch, not the encounter it actually has.** Both drift during the long
infall, outward in distance and later in time.

Measured across the existing 2027 sweep (`v_inf` = 25, epoch 1873, nominal
periapsis 2027-07-26):

| requested rp | achieved rp | achieved date | drift |
|---:|---:|---|---:|
| 0.25 AU | 0.385 AU | 2027-08-17 | +22 d |
| 0.50 AU | 0.637 AU | 2027-08-19 | +24 d |
| 0.75 AU | 0.878 AU | 2027-08-19 | +24 d |
| 1.00 AU | 1.116 AU | 2027-08-19 | +24 d |
| 1.25 AU | 1.355 AU | 2027-08-19 | +24 d |
| 1.50 AU | 1.594 AU | 2027-08-18 | +23 d |

and for the 2047 configuration (`v_inf` = 10, nominal 2047-07-26):

| requested rp | achieved rp | achieved date | drift |
|---:|---:|---|---:|
| 0.25 AU | 0.331 AU | 2047-08-29 | +34 d |
| 0.50 AU | 0.599 AU | 2047-09-03 | +39 d |
| 1.50 AU | 1.644 AU | 2047-09-19 | +55 d |

**The sign of the distance error depends on orbital orientation, not just on
infall length.** Every run in the two tables above shares Ω = 0, and they all
arrive wide. Sampling other orientations at rp = 0.5 AU shows the opposite is
equally possible:

| i | Ω | ω | achieved rp | achieved date |
|---:|---:|---:|---:|---|
| 0 | 0 | 0 | 0.599 AU | 2047-09-03 |
| 30 | 0 | 0 | 0.608 AU | 2047-09-06 |
| 60 | 0 | 0 | 0.586 AU | 2047-08-30 |
| 60 | 270 | 180 | **0.416 AU** | 2047-08-19 |
| 90 | 270 | 180 | **0.422 AU** | 2047-08-21 |

Ω = 0 arrives ~0.09 AU wide; Ω = 270 arrives ~0.08 AU narrow, and about two
weeks earlier. So the achieved periapsis spans 0.42–0.61 AU across a sweep that
requested 0.50 throughout — roughly ±20%. An earlier version of this section
described the drift as a uniform outward bias, which was an artefact of every
sampled run sharing one orientation.

Two consequences.

**The timing drift is why the epoch is not a free choice.** The lag grows with
the length of the infall, so moving the epoch changes the achieved periapsis date
even when the nominal date is held fixed: from an 1885 epoch the same nominal
2047-07-26 arrives on 2047-09-24 rather than 2047-09-03. That 21-day shift is
enough to move Earth's post-flyby semi-major axis from 1.259 to 1.025 AU — the
response to a timing shift is smooth, not chaotic, at about 0.011 AU per day.
**The epoch is therefore part of the scenario definition, not a bookkeeping
detail.**

**The drift is exploited deliberately in the 2047 configuration.** Section 1
asks for a periapsis "around 1 September 2047" and sets the nominal date to
26 July; the +34 to +55 day drift lands the achieved encounter between 29 August
and 19 September, median ~5 September. Setting the nominal date to 1 September
instead would put the real encounter in mid-October.

**Cause — measured, not inferred.** The BH is initialised on an exact two-body
orbit *relative to the Sun*, but from 400 AU it is not orbiting the Sun: it is
orbiting the whole solar system, whose mass sits at the barycentre. The Sun is
offset from that point and moving at ~13 m/s relative to it, mostly Jupiter's
doing, and that velocity is never removed.

Measuring the same orbit against both centres settles it. For the scenario run
(`rp0p5…Om270__om180`), osculating periapsis through the infall:

| year | vs the **Sun** | vs the **barycentre** |
|---:|---:|---:|
| 1874 | 0.4931 AU | **0.4176 AU** |
| 1880 | 0.3180 AU | **0.4176 AU** |
| 1950 | 0.3722 AU | **0.4176 AU** |
| 2000 | 0.4044 AU | **0.4176 AU** |
| 2046 | 0.4150 AU | **0.4174 AU** |

Against the barycentre the periapsis is **constant to 0.0003 AU over 170
years**; against the Sun it swings by 0.24 AU. The swing is the Sun wobbling on
Jupiter's 12-year period, not the orbit changing. The achieved closest approach,
0.4157 AU, is the barycentric value — the initial state encoded it from day one.

The leverage is that a nearly radial orbit carries almost no angular momentum:
h = 0.018 AU²/day here. At 400 AU a 13 m/s transverse discrepancy contributes
r·δv ≈ 3×10⁻³ AU²/day, up to 18% of the total, from a velocity error of one part
in a thousand. Since r_p = h²/[μ(1+e)], the measured −9% in h becomes −17% in
periapsis distance. It also explains the orientation dependence above: whether
the Sun's barycentric velocity adds to or subtracts from the BH's angular
momentum depends on the approach geometry, so Ω = 0 comes out wide and Ω = 270
narrow from the same requested value.

**The fix is known and small, and has not been applied.** Constructing the BH's
initial state relative to the Sun+planets barycentre, with
`mu = G(M_sun + M_planets + M_bh)`, would make `bh_rp_au` mean what it says: a
handful of lines in `build_sim`. It is not applied because it changes the
initial conditions of every run, so the entire sweep — and every derived
product — would have to be regenerated, to no benefit for any question asked so
far. The drift is systematic and now fully measured, so it is cheaper to read
the parameters as labels than to re-run 672 integrations.

---

## 9. What is assumed and not tested

- **The 30-year baseline** is a stand-in for a real observing campaign. Actual
  1885 knowledge of Uranus rested on ~195 years of positions of very uneven
  quality; the effective sensitivity is not a simple function of span.
- **Detrending by constant + linear** approximates an orbit refit. A real
  six-element refit absorbs more, so the derived `v_inf` is an upper bound.
- **Theory error is not modelled**, only measurement geometry. Section 3 argues
  it dominated in 1885, which would lower `v_inf` further.
- **Inner-planet figures are sampled at ~30-day intervals**, which undersamples
  Mercury and Venus. They are ~100× below Uranus and do not affect any
  conclusion, but the individual numbers should not be quoted.
- **The BH's own detectability** by other channels — astrometric deflection of
  background stars, microlensing — is a separate analysis thread in this
  repository (`GAIA_dr3_v1.py` and successors) and is not folded in here.


---

## 10. Angular grid resolution

`bh_omega_deg` is sampled at 15°, six times finer than the other angles, and that
asymmetry was measured rather than guessed.

**ω is the parameter the outcome is rough in.** Scanning it at 15° with
everything else fixed (rp = 0.5, i = 30°, Ω = 0), *inside a single 90° cell*
between ω = 45° and 75°, Earth goes from untouched to ejected and the number of
unbound planets from zero to four. The four values a 90° grid samples span Earth
a = 0.82–2.10 where the full scan spans −6.68 to 3.11, and every one of them
unbinds one or two planets while the finer scan finds cases with zero and with
four. **A 90° grid is not merely coarse, it is biased** — it samples the middle
of the outcome distribution and misses both tails, which affects every ejection
statistic derived from it.

**Ω and inclination are not.** The same scan in Ω shows the four sampled values
already recover nearly the full range (0.770–1.716 against 0.770–1.787).
Inclination at 30° recovers **99–100%** of what a 15° grid finds, at every ω
tested including the rough zone.

The reason is physical: **ω sets the direction from the Sun to periapsis**, so
sweeping it walks the encounter past one planet after another. Ω rotates the line
of nodes and i tilts the plane; neither aims the flyby at anything.

**A rejected intermediate.** A first scan suggested inclination was oversampled
and could be coarsened to pay for ω — but that scan was taken at ω = 0, which is
precisely where inclination does not matter. Repeating it across ω = 0, 30, 60,
90 showed a 45° inclination grid collapses to **20%** of the range at ω = 60°,
and a 60° grid loses a quarter of it at ω = 90° (the peak sits at i = −60°, which
neither samples). The saving would have cost more than it bought.

**An exact symmetry halves the grid.** Runs pair off under
**(i, Ω, ω) → (−i, Ω+180°, ω+180°)**, matching to six decimal places throughout
the impact ranking. Roughly half of the 4032 runs are therefore physically
redundant. Sampling i over 0–90° instead of ±90° would buy back a factor of two,
at the cost of losing a built-in consistency check that has already proved
useful for spotting bookkeeping errors.
