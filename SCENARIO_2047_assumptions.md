# The 2047 scenario — assumptions and their provenance

Design record for the black-hole flyby sweep with **periapsis in July 2047**,
superseding the 2027 configuration used for `simulations/20260724_230314`.

Every parameter below is either a free choice or is *forced* by a stated
constraint. This document records which is which, so that later work can tell
what may be varied freely and what cannot be changed without redoing an
analysis. Numbers are reproducible with
[`astrometric_detectability.py`](astrometric_detectability.py); the command that
produces each table is given with it.

Companion documents: [`SCENARIO_mars_window.md`](SCENARIO_mars_window.md) — the
chosen run from this sweep, studied in detail — and
[`SCENARIO_timeline.md`](SCENARIO_timeline.md), its narrative century.
[`SCENARIO_mercury_capture.md`](SCENARIO_mercury_capture.md) covers the retired
2027 sweep and is kept for the record.

---

## 1. The configuration

```yaml
bh_mass_msun: 0.1                # imposed
bh_vinf_kms: 10.0                # DERIVED -- see section 4
bh_tperi_offset_days: 63515      # imposed: 1873-09-01 + 63515 d = 2047-07-26
bh_rp_au: "0.25, 1.5, 0.25"      # imposed (swept)
duration_days: 119725            # chosen: preserves the 2027 sweep's aftermath
```

The four imposed constraints were:

1. Uranus's orbital residuals **barely recognizable in 1885**
2. BH mass = 0.1 M☉
3. periapsis July 2047
4. periapsis distance swept over 0.25–1.50 AU

Constraint 1 is an *observability* constraint, not a dynamical one, and it is
the only reason `v_inf` is not free. Everything downstream of it inherits that
choice — see section 6 for what changes if it is relaxed.

`duration_days: 119725` gives 56,210 days after periapsis, identical to the 2027
sweep, so climate results remain directly comparable. The alternative, 127030,
would place periapsis at the exact midpoint of the integration as the old config
did, at ~6% more cost.

---

## 2. Method: what "detectable" means here

`python astrometric_detectability.py horizon`

Two N-body integrations, identical but for the BH, differenced. Sun + 8 planets
from de440s at the 1873-09-01 epoch, IAS15.

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

`python astrometric_detectability.py vinf --peri 2047-07-26 --when 1885`

With periapsis pinned to 2047-07-26, 1885 sits 162.6 years before arrival. That
fixes the travel time, so requiring a particular residual in 1885 fixes the
distance then, which fixes the approach speed. Uranus residual over an
1855–1885 window, rp = 0.5 AU:

| `v_inf` | r(1855) | r(1885) | detrended | raw |
|---:|---:|---:|---:|---:|
| 8 km/s | 369 AU | 317 AU | 1.60″ | 3.07″ |
| **10 km/s** | **440 AU** | **376 AU** | **0.92″** | 1.94″ |
| 12 km/s | 514 AU | 437 AU | 0.55″ | 1.27″ |
| 14 km/s | 591 AU | 501 AU | 0.34″ | 0.86″ |
| 16 km/s | 668 AU | 566 AU | 0.22″ | 0.60″ |
| 19 km/s | 786 AU | 665 AU | 0.11″ | 0.37″ |
| 22 km/s | 905 AU | 765 AU | 0.06″ | 0.23″ |
| 25 km/s | 1025 AU | 866 AU | 0.03″ | 0.15″ |

**"Barely recognizable" at ~1″ gives `v_inf` = 10 km/s.** Read the phrase as 2″
and it becomes ~7.5; as 0.5″, ~12.5. **The defensible range is 8–13 km/s.**

At `v_inf` = 10 the BH sits at **400 AU at the 1873 epoch**, against 926 AU for
the old 25 km/s configuration.

**Periapsis distance does not interact with this constraint.** Across
rp = 0.25 → 1.50 AU the required `v_inf` moves by 0.04 km/s, because 1.5 AU is
nothing against 400 AU. rp may be swept freely without revisiting any of the
above.

---

## 5. Consequences for the encounter itself

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
receives the same transformation — including the mirror-twin identity asserted in
[`SCENARIO_mercury_capture.md`](SCENARIO_mercury_capture.md).

**Where it bites:** any absolute geometric statement about a *single named run*
— "in run `…om270` the BH approached from direction X". Note that plots reading
positions out of the simulation output (`plot_sky_tracks.py`,
`plot_local_sky.py`) show the true geometry and are correct; it is only the
mapping from a folder's ω label to that geometry that is off.

Left unchanged deliberately: correcting it would silently relabel every existing
run folder.

---

## 8. Nominal parameters are labels, not achieved geometry

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
