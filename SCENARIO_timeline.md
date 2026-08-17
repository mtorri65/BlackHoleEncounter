# A century of warning, 1885–2100

A narrative timeline of the scenario in
[`SCENARIO_mars_window.md`](SCENARIO_mars_window.md): a 0.1 M☉ black hole
reaching perihelion at **0.662 AU on 20 August 2047**, approaching at 25 km/s.

Every date and distance below is read from the simulation output. The framing is
narrative; the numbers are not.

The century is in three acts, and the first one is longer than it looks. The
object is gravitationally detectable from around 1900 and unambiguously so
within a generation — so the story does not begin in 2047. It begins with a
residual in Uranus's position, and the hundred and fifty years between that
residual and the encounter are the part with people in it.

*(Rebuilt for the current sweep — epoch 1885-09-01, `v_inf` = 25 km/s. An earlier
version of this document described a `v_inf` = 10 km/s configuration in which the
hole was less than half as far out at every epoch; all of Act I's distances and
signal strengths are correspondingly larger there. An earlier retired scenario
placed perihelion in 2027; its write-up has been deleted as superseded and
survives only in git history.)*

---

## Act I — The warning (1885–2046)

### The shape of the approach

| year | distance | | year | distance |
|---:|---:|---|---:|---:|
| 1885 | 863 AU | | 2010 | 205 AU |
| 1900 | 787 AU | | 2020 | 152 AU |
| 1930 | 628 AU | | 2030 | 98 AU |
| 1955 | 496 AU | | 2040 | 45 AU |
| 1980 | 364 AU | | 2044 | 22 AU |
| 2000 | 258 AU | | Jan 2047 | 4.9 AU |

Approach speed is **25 km/s ≈ 5.3 AU/yr**, essentially constant beyond 50 AU:
`v∞` dominates the solar potential across the whole outer approach.

The object emits nothing. It is 591 metres across, has a Hawking temperature of
6×10⁻⁷ K, and sits where there is nothing to accrete. It is optically dark at
every epoch, to every instrument, for the entire period. Its only signature is
gravitational — and gravitationally it is enormous.

### 1885 · 863 AU · at threshold

By construction, the scenario is tuned so the hole sits *just* at the edge of
nineteenth-century detectability. At 863 AU it induces a residual of about
**0.76 arcseconds** in Uranus's position over a century-long campaign — below
the ~1″ a Le Verrier-class analysis might have claimed, and well below what the
era's planetary theory could have defended. See §4 of
[`SCENARIO_2047_assumptions.md`](SCENARIO_2047_assumptions.md), which also
explains why that number is far less robust than it looks.

The limit here is theory, not measurement. That matters, because theory error is
the one term in the budget guaranteed to improve.

### 1900–1955 · the anomaly becomes real

Two independent things now happen, and they multiply. The signal rises because
the object approaches; the noise falls because Newcomb and Hill's rigorous
outer-planet theory arrives ~1895 and removes the systematic floor.

| | 1885 | 1900 | 1930 | 1955 |
|---|---:|---:|---:|---:|
| distance | 863 AU | 787 AU | 628 AU | 496 AU |
| signal (D⁻³) | 0.76″ | 1.00″ | 1.97″ | 4.00″ |
| theory error | ~1″ | ~0.3″ | ~0.2″ | ~0.15″ |
| **SNR** | **0.8** | **3.3** | **9.9** | **27** |

Of the 4.4× gain in significance from 1885 to 1900, about **3.3× is the
astronomers and only 1.3× is the object.** Newcomb–Hill dominates that first
decade and a half.

The 1.3× is worth stating explicitly, because the leverage is not obvious. The
object closes at ~5.3 AU/yr, so 1885→1900 brings it only **8.8% nearer**. Cubed,
that becomes a **32% stronger perturbation**: (863/787)³ = 1.32. The factor of
three in the exponent is the engine of this whole act.

That asymmetry is also why the two terms trade places. Theory error improves only
until it reaches the measurement floor, and then it stops. D⁻³ does not stop, and
accelerates as D shrinks:

    d/dt ln(M/D³) = 3v/D ≈ 2.5 %/yr at 630 AU

By the 1930s the approach contributes more per decade than any refinement of the
tables. **The perturbation doubles between 1900 and 1930 with no change of
instrument.** An accelerating, direction-fixed, non-Keplerian residual is not
absorbable by any refit of the known planets.

Historically this lands squarely on the Planet X searches — Lowell's from 1905,
Tombaugh's culminating in 1930. In the real solar system those searches chased
residuals that were spurious. Here they are real, and Pluto is the wrong object:
the anomaly survives its discovery and keeps growing.

### The 1930s · it is not a planet

A distant perturber's tidal quadrupole yields a **direction** and the combination
**M/D³**. Mass and distance are degenerate. Two measurements break it.

**The growth rate gives v/D.** Thirty years of residuals fix 3v/D, hence the
approach speed in units of the distance.

**The direction drifts, and that gives the impact parameter.** Angular momentum
h = r_p·v_p = 0.662 AU × 59.8 km/s, so the line of sight to the perturber sweeps
at h/D² ≈ **4.4 arcsec/yr at 630 AU** — more than two arcminutes over a
thirty-year campaign, and rising as D⁻².

Combining them changes the character of the problem entirely. At 630 AU the local
escape speed from the Sun is 1.7 km/s. The measured speed is **25 km/s: fifteen
times escape, 220 times in energy.**

**The object is not a planet. It is unbound, and it is arriving.**

An analysis of that era, carrying 15% errors, puts arrival at 628 AU ÷ 5.3 AU/yr
≈ 118 years hence — **the middle of the twenty-first century, ±18 years.** The
scenario acquires a deadline more than a century before the deadline.

### 1962–1990 · ranging ends the argument

Venus radar (1961) is *not* decisive here — at 459 AU the accumulated range
signal is a few kilometres against ~50 km precision, and it detects nothing. That
is a change from the closer configuration this document once described, and it
gives the argument twenty more years to run on astrometry alone.

Spacecraft ranging settles it. **Mariner 9 (1971)** makes it a 10²σ measurement;
**Viking's seven-metre lander ranging (1976–82)** makes it 10⁴σ.

By 1980, M, D, v and the impact parameter are pinned directly, with no degeneracy
left and no reliance on assumption. **The 2047 encounter is predictable in detail
from the Viking data set** — sixty-seven years out.

### 1990–2018 · three channels saturate

**VLBI.** The Sun's acceleration toward the hole produces a dipole proper-motion
pattern in the quasar frame of magnitude a/c — **159 µas/yr in 1995, 491 µas/yr
by 2018.** The analogous Galactic aberration drift, at 5.8 µas/yr, is already
measured and built into ICRF3. This is 27–85× larger and points straight at the
object.

**Ranging.** Cassini-class Saturn ranging (~25 m) against a ~1.4×10⁴ km signal is
a 6×10⁵σ detection. By 2030 it is 10⁸σ.

**Gaia.** The Einstein radius is θ_E ≈ **1.02 arcseconds**, which is enormous.
Weak deflection falls as θ_E²/θ, so a background star a full degree away is still
displaced by ~290 µas against Gaia's ~25 µas — a coherent astrometric distortion
field tens of degrees across, dragged over the sky by the object's parallax. This
is the moment the perturber stops being a direction in a residual field and
becomes **a point on the sky**.

### The detection table

| year | D | best instrument | observable | signal | precision | **SNR** |
|---:|---:|---|---|---:|---:|---:|
| 1885 | 863 AU | meridian circle, Le Verrier tables | Uranus angle | 0.76″ | ~1″ | **0.8** |
| 1900 | 787 AU | photographic plates, Newcomb–Hill | Uranus angle | 1.00″ | ~0.3″ | **3.3** |
| 1930 | 628 AU | Lowell-era plates, refined theory | Uranus angle | 1.97″ | ~0.2″ | **10** |
| 1955 | 496 AU | machine numerical ephemerides | Uranus angle | 4.00″ | ~0.15″ | **27** |
| 1962 | 459 AU | Venus radar ranging | Venus range | ~5 km | ~50 km | **0.1** |
| 1971 | 412 AU | Mariner 9 orbiter tracking | Mars range | 13 km | ~100 m | **1×10²** |
| 1980 | 364 AU | Viking lander ranging | Mars range | 84 km | 7 m | **1×10⁴** |
| 1995 | 285 AU | VLBI / ICRF | quasar aberration | 159 µas/yr | ~1 µas/yr | **2×10²** |
| 2010 | 205 AU | Cassini ranging | Saturn range | 1.4×10⁴ km | 25 m | **6×10⁵** |
| 2018 | 162 AU | Gaia DR2 | lensing deflection | 288 µas at 1° | 25 µas | **12** |
| 2030 | 98 AU | contemporary ranging | Saturn range | 1.3×10⁵ km | ~1 m | **1×10⁸** |

### April 2044 · Neptune at 19.5 AU · September 2046 · Uranus at 16.2 AU

The first two encounters, both distant. The hole crossed Neptune's orbital
distance in August 2042 and is now accelerating hard: it reaches 11 AU on
2 January 2046 and closes that last stretch in twenty months — 4.9 AU by the end
of 2046, inside Earth's orbit by 22 July 2047.

### 2046 · 11 AU · why the elements still read normal

Here is the thing that is easy to state wrongly. **In January 2047 — seven months
before perihelion, with the black hole 4.9 AU from the Sun — Earth's semi-major
axis still reads 1.000 AU and Mars's 1.523.** The orbital elements are, on paper,
entirely ordinary until the encounter itself.

That is true, and it is not in tension with anything above. The two describe
different quantities. A tidal field of this strength produces a large *positional*
perturbation while leaving the *secular elements* almost untouched: the Saturn
signal at 2010 is ~1.4×10⁴ km of accumulated range drift, which is 10⁻⁴ AU
against a semi-major axis of 9.573 AU. An element table printed to four decimals
shows nothing. A ranging residual shows a 10⁵σ anomaly. Both are correct.

So the century's real asymmetry is this: for a hundred and fifty years the solar
system is measurably falling toward something, and for a hundred and fifty years
the planets keep their orbits to four decimal places. **Everything is known and
nothing has happened yet.**

---

## Act II — The passage (May 2047 – November 2048)

### 10 May 2047 · Jupiter at 4.92 AU

Distant, and Jupiter comes through it almost unscathed — 5.206 AU to 5.011, with
eccentricity rising to 0.170. The largest planet is barely moved.

### 28 June 2047 · Mars at 0.810 AU

**The closest planetary approach of the whole encounter, and the one that makes
the scenario.** Mars is pushed *inward*: 1.524 AU to 1.336, eccentricity from
0.093 to 0.217, year from 687 days to 564.

**Earth and Mars have swapped order.** Mars is now the inner of the two.

### 20 August 2047 · perihelion, 0.662 AU

The black hole reaches its closest approach to the Sun — **0.662 AU, inside
Venus's orbit** — moving at **59.8 km/s**.

### 26 August 2047 · Mercury changes hands

Six days after perihelion, Mercury's orbital energy with respect to the black
hole goes negative and stays there. It is no longer a planet of the Sun.

Its new orbit is a = 0.372 AU, e = 0.923, **period 0.72 years**. On 3 July 2048
it makes its first periapsis passage around the hole, at **0.025 AU** — 3.7
million km from an object 591 metres across, and fourteen times the distance at
which it would be torn apart. **Mercury survives intact and leaves as a moon of
something dark.**

### 28 August – 11 September 2047 · Venus, Earth and the Moon

Venus at 1.210 AU on the 28th; Earth at 1.291 AU on 10 September; the Moon at
1.290 AU the day after.

**No near miss for Earth — and decisive anyway.** Within months its semi-major
axis has gone from 1.000 to 1.486 AU and its eccentricity from 0.017 to 0.236.
The mechanism is not a close pass but a **differential tug**: the hole pulls the
Sun and the Earth by different amounts in different directions, and what survives
is the difference.

**The Moon comes through it** — still bound, still ordinary, orbiting a world
about to freeze. Its orbit is essentially untouched: the semi-major axis ends at
about 386 000 km against today's 384 400, and the Earth–Moon separation stays
inside its normal 357 000–407 000 km perigee/apogee spread for the entire
integration. (The run ends with the Moon near apogee at 402 000 km, which is
where it happens to be, not a widened orbit.) The binding is verified on relative
state vectors at every output step in
[`SCENARIO_post_flyby_system.md`](SCENARIO_post_flyby_system.md) §4.

### 25 November 2048 · Saturn at 7.41 AU

The last encounter, and the strangest outcome. Saturn is not ejected; it is
flung from 9.573 AU onto an orbit whose semi-major axis reads 166 AU in 2048,
4100 AU by 2060, and 1650 AU a century later, with eccentricity 0.994. It
remains, formally, a planet of the Sun. It would not survive the next passing
star.

### The state by the end of 2048

| | Before | After |
|---|---|---|
| **Mercury** | 0.387 AU, e 0.206 | **captured by the black hole** |
| Venus | 0.723 AU, e 0.007 | 0.989 AU, e 0.200 |
| **Earth** | 1.000 AU, e 0.017 | **1.486 AU, e 0.236** |
| **Mars** | 1.524 AU, e 0.093 | **1.336 AU, e 0.217** |
| Jupiter | 5.206 AU, e 0.048 | 5.011 AU, e 0.170 |
| **Saturn** | 9.573 AU, e 0.055 | **~1650 AU, e 0.994** |
| Uranus | 19.299 AU, e 0.051 | 11.485 AU, e 0.679 |
| Neptune | 30.090 AU, e 0.008 | 29.927 AU, e 0.516 |

*The table above is read from the run's `__planets_run_deltas.csv`.
[`SCENARIO_post_flyby_system.md`](SCENARIO_post_flyby_system.md) tabulates the
same run from `bh_captures.csv`, and the two differ slightly — Earth 1.4874 against
1.4862 AU, Jupiter 5.011 against 5.004 — because the elements are extracted at
different points. For Saturn the gap looks larger, 1,650 against 1,627 AU, but
that is the same small difference amplified by e = 0.994, where the semi-major
axis is hypersensitive to orbital energy. Perihelia and eccentricities agree.*

*What these elements imply — Saturn's effective expulsion (aphelion 3,244 AU on a
65,600-year period), Venus inheriting Earth's orbit, the seven orbit-crossing
pairs that make the arrangement unstable beyond the integration, the survival of
the Earth–Moon system, the scattering of the asteroid belt onto Earth- and
Mars-crossing orbits, and the end of annular eclipses — is worked through in
[`SCENARIO_post_flyby_system.md`](SCENARIO_post_flyby_system.md).*

---

## Act III — Consequence (2049–2100)

### Earth freezes, in four years

Earth's year is now **662 days** — 1.81 calendar years — and it receives about
45% of the sunlight it does today.

| Earth orbit | ≈ calendar | Global mean | Northern ice edge |
|---:|---:|---|---|
| 1 | 2049 | 259.6 K | **21.9° latitude** |
| **2** | **2051** | **216.1 K** | **pole to pole** |
| 3 | 2052 | 197.6 K | frozen |
| 5 | 2056 | 185.8 K | frozen |
| 80 | 2192 | **182.1 K (−91.1 °C)** | equilibrium |

By the end of Earth's *first* orbit the permanent ice edge has reached the
subtropics. **On the second orbit — 2051, four calendar years after the
encounter — the ice closes over the equator.** Everything after that is a frozen
planet settling toward its **182 K** floor.

**It is permanent, and the reason is that freezing and thawing have different
thresholds.** Earth tips into the snowball once its annual-mean sunlight falls
below about **0.95×** today's, but climbing back out takes **1.28×** — some **28%
more sunlight than present-day Earth receives**. That gap is the ice-albedo
feedback: a white planet reflects most of what arrives, so restoring the original
sunlight is nowhere near enough to undo the original freeze. The new orbit
delivers well under half. There is no path back within the century, or any
century.

*(Caveat: a thermal result. Oceans reach freezing on this schedule;
kilometre-thick ice sheets take far longer. The model has no ice-sheet dynamics.
The two thresholds above are Sellers-OLR values, consistent with the temperatures
in this section; under the linear law they would be 0.91× and 1.17×.)*

> **Which outgoing-radiation law?** These figures use the **Sellers** nonlinear
> OLR, which is the defensible choice here: §10 of
> [`CLIMATE_MODEL_REPORT.md`](CLIMATE_MODEL_REPORT.md) documents the simpler
> linear (Budyko) law failing at exactly the cold end this run occupies, because
> its emission collapses toward zero and a frozen planet stops shedding heat
> properly. Under the linear law the same run reads 261.4 / 224.1 / 210.9 / 205.3
> / **204.7 K** — some 23 K warmer at the floor, and the figure earlier editions
> of these documents quoted. The effect is systematic, a median 15.97 K across all
> 1,182 snowball runs in the sweep (§6.8). **The freeze schedule is identical
> either way** — ice still closes over the equator on the second orbit; only the
> floor moves.

### Mars, meanwhile, has a season

Mars's year is now 564 days, and its distance from the Sun swings between 1.046
and 1.626 AU. Obliquity is essentially unchanged at 23.8°, and the longitude of
perihelion has barely moved.

The result is not a warm planet. The global mean is 217.5 K, the equator's annual
mean 235 K, the winter pole 171 K. A quarter of the atmosphere freezes onto the
winter cap and returns each year, cycling the surface pressure between 537 and
714 Pa. The atmosphere never collapses.

But at **21.5° south**, for **41 days of each Martian year**, the surface holds
between 273.2 and 275.0 K while the pressure sits near 700 Pa — and under those
conditions liquid water is thermodynamically permitted. It will not freeze, and
it will not boil.

It is the first time in the planet's history that has been true.

**And it happens because Mars stays cold.** The peak temperature anywhere on the
planet, at any point in its year, is 278.7 K. At ~700 Pa water boils at 275.3 K,
so the entire liquid range is barely two kelvin wide — and a planet that swings
to 340 K crosses that band twice a year at speed. Across the sweep the wettest
worlds are the ones whose maxima sit between 280 and 290 K; past 290 K a warmer
Mars is a *drier* one. This Mars creeps to the edge of freezing and lingers
there. The window is long because the planet is barely warm enough, not because
it is warm.

The margin is correspondingly thin. A modest error in the CO₂ inventory, the
albedo or the heat transport closes it entirely. And the window says only that
liquid water *could* exist — whether any is present to do so is a question the
model cannot answer.

### 2100

The black hole is 283 AU away and receding, with Mercury in a 0.72-year orbit
around it. Earth is a frozen world at 182 K and will remain one. Saturn is
somewhere past a thousand astronomical units on an orbit it will not keep. And
Mars — smaller, colder, and now the third planet rather than the fourth — has six
weeks each year in its southern subtropics during which water could, in
principle, be liquid.

Within a single century the solar system has lost one planet outright, put
another on an orbit no passing star will let it keep, frozen a third, and handed
a fourth a narrow, precarious opening. **None of it was a surprise.**

---

## What is and is not modelled

**Solid:** all dynamics (REBOUND, IAS15), every date and distance above, and the
post-flyby orbital elements, which are validated by reproducing present-day
values when applied at t = 0.

**Modelled with stated limits:** the climate outcomes. The Earth model omits
ice-sheet dynamics; the Mars model omits water entirely, along with the diurnal
cycle, topography and dust. See `CLIMATE_MODEL_REPORT.md`.

**Reconstructed, not simulated:** everything in Act I. The detection analysis is
analytic — tidal amplitudes, lensing geometry and aberration computed from the
trajectory, compared against representative period instrument precisions. It is
not the output of an ephemeris-fitting experiment, and the appendix states where
it is weakest.

**Not modelled at all:** tidal or relativistic effects, and everything about the
biosphere — including, pointedly, what a civilisation does with a century and a
half of warning. Act I establishes that the warning exists. What happens inside it
is not a question this repository can answer.

**The asteroid belt is absent from this run** (`n_belt = 0`), so nothing above
reflects it. A companion run at identical parameters
(`simulations/20260815_171552/`) adds 4,000 belt tracers and finds the belt thrown
onto Earth- and Mars-crossing orbits; see §6 of
[`SCENARIO_post_flyby_system.md`](SCENARIO_post_flyby_system.md). That opens an
era of impacts which none of the dates above account for.

**And this is one run of 4032.** It ranks 63rd on orbital disruption — not an
extreme event — but it is one of only four that freeze the Earth, hand a planet
to the black hole, and open a Martian water window at the same time, and the best
of those four: its 7.2% water window against 1.7% for the alternative. A typical
flyby in this sweep does none of it. The century described here is a possibility,
not a prediction.

*Read those counts with §7a of
[`SCENARIO_2047_assumptions.md`](SCENARIO_2047_assumptions.md) in hand: the sweep
is exactly twofold redundant, so 4032 runs are ~2016 distinct configurations, and
the four qualifying runs are two mirror pairs — **two configurations, not four.**
This run's twin ranks 64th with an identical score.*

---

## Appendix · detection method and assumptions

**Trajectory.** Radial infall, dr/dt = −√(v∞² + 2μ/r), with μ = G(M☉ + M_BH),
checked against the run's own logged distances.

**Tidal amplitude.** Δa = 2GM·Δr/D³, the differential acceleration between Earth
and the target body. The common-mode fall of the whole system toward the hole is
unobservable and is excluded throughout.

**Suppression factor.** The free-drift estimate ½·Δa·T² overstates the
*unabsorbed* residual, because an ephemeris refit takes up part of the signal.
Calibrated against the Planet Nine literature: 10 M⊕ at 700 AU over a 13-year
Cassini baseline gives 747 m free-drift against a published signal of ~100 m and
a threshold of ~20 m, so **a factor of 8** is applied to every ranging figure.
This is the least secure number here, and nearly irrelevant: at a factor of 100
the Cassini-era SNR is still 10⁴.

**Astrometric scaling.** Anchored to §4 of the assumptions document (0.76″ at
863 AU, 1885, over a century-long campaign) and extrapolated as D⁻³.

**The astrometric anchor is fragile, and this is the weakest link.** §4 shows the
detrended residual moves by a factor of ~20 with the *position* of the assumed
observing window and ~22 with its *length*, because removing a constant and a
linear trend is only a proxy for an orbit refit and Uranus's 84-year period means
a few decades of data can look nearly linear. Every astrometric figure in Act I
inherits that. The ranging and aberration channels do not — they are smooth,
monotonic in D⁻³ and D⁻², and colossal — so the *conclusion* that the object is
unmissable from the 1970s onward is robust even though the 1885 threshold is not.

**Instrument precisions** are representative period values, not results of
specific published analyses. Nineteenth-century meridian work is theory-limited
rather than measurement-limited; Venus radar improves from ~50 km (1961) to
~1 km (1970); Viking lander ranging reaches ~7 m; Cassini ~20–30 m; Gaia
end-of-mission ~25 µas.

**Not addressed:** the response of the Kuiper belt after the hole crosses
Neptune's distance around 2042, which is a further independent signal and is not
folded in.

**What cannot detect it:** imaging, at any epoch. Being a black hole makes the
object optically invisible while making it, at fixed orbit, far louder
gravitationally than any planet — 0.1 M☉ is 3,300× Planet Nine's mass, and even
at 863 AU its M/D³ exceeds Planet Nine's by ~2,000×. The same M/D³ that makes the
flyby consequential is what makes it observable, and both scale identically.
