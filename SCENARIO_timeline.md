# A century of warning, 1885–2097

A narrative timeline of the scenario in
[`SCENARIO_mars_window.md`](SCENARIO_mars_window.md): a 0.1 M☉ black hole
reaching perihelion at **0.416 AU on 19 August 2047**, approaching at 10 km/s.

Every date and distance below is read from the simulation output. The framing is
narrative; the numbers are not.

The century is in three acts, and the first one is longer than it looks. The
object is gravitationally detectable from the 1880s and unambiguously so from
around 1900 — so the story does not begin in 2047, or in 1997. It begins with a
residual in Uranus's position, and the hundred and sixty years between that
residual and the encounter are the part with people in it.

*(This replaces an earlier timeline written for the retired 2027 sweep; that
scenario is preserved in [`SCENARIO_mercury_capture.md`](SCENARIO_mercury_capture.md).
It also absorbs the former `SCENARIO_detection_timeline.md`, whose analysis is now
Act I and the appendix.)*

---

## Act I — The warning (1885–2046)

### The shape of the approach

A radial-infall integration (`v∞` = 10 km/s, `μ = G(M☉ + 0.1 M☉)`, periapsis
0.416 AU on 2047-08-19) reproduces the run's own output to better than 1%:

| year | distance | | year | distance |
|---:|---:|---|---:|---:|
| 1885 | 376 AU | | 2020 | 77 AU |
| 1900 | 344 AU | | 2035 | 41 AU |
| 1930 | 278 AU | | 2040 | 27 AU |
| 1955 | 224 AU | | 2044 | 16 AU |
| 1997 | 130 AU | | 2046 | 8.7 AU |
| 2000 | 123 AU | | Jan 2047 | 4.5 AU |

Approach speed is **10.3 km/s ≈ 2.2 AU/yr**, essentially constant beyond 50 AU:
`v∞` dominates the solar potential across the whole outer approach.

The object emits nothing. It is 295 metres across, has a Hawking temperature of
6×10⁻⁷ K, and sits inside the heliosphere where there is nothing to accrete. It
is optically dark at every epoch, to every instrument, for the entire century.
Its only signature is gravitational — and gravitationally it is enormous.

### 1885 · 376 AU · at threshold

By construction, the scenario is tuned so that the hole sits *just* at the edge
of nineteenth-century detectability. At 376 AU it induces a residual of about
**0.9 arcseconds** in Uranus's position over a thirty-year campaign — the level
at which a Le Verrier-class analysis might have claimed something, and below what
the era's planetary theory could have defended. See §4 of
[`SCENARIO_2047_assumptions.md`](SCENARIO_2047_assumptions.md).

The limit here is theory, not measurement. That matters, because theory error is
the one term in the budget guaranteed to improve.

### 1895–1930 · the anomaly becomes real

Two independent things now happen, and they multiply. The signal rises because
the object approaches; the noise falls because Newcomb and Hill's rigorous
outer-planet theory arrives ~1895 and removes the systematic floor.

| | 1885 | 1900 | 1930 |
|---|---:|---:|---:|
| distance | 376 AU | 344 AU | 278 AU |
| signal (D⁻³) | 0.92″ | 1.21″ | 2.27″ |
| theory error | ~1″ | ~0.3″ | ~0.2″ |
| **SNR** | **0.9** | **4.0** | **11** |

Of the 4.4× gain in significance from 1885 to 1900, about **3.3× is the
astronomers and only 1.3× is the object.** Newcomb–Hill dominates that first
decade and a half.

The 1.3× is worth stating explicitly, because the leverage is not obvious. The
object closes at ~2.2 AU/yr, so 1885→1900 brings it only **8.6% nearer**. Cubed,
that becomes a **31% stronger perturbation**: (376.1/343.6)³ = 1.31. The factor
of three in the exponent is the engine of this whole act.

That asymmetry is also why the two terms trade places. Theory error improves only
until it reaches the measurement floor, and then it stops. D⁻³ does not stop, and
accelerates as D shrinks:

    d/dt ln(M/D³) = 3v/D ≈ 2.3 %/yr at 280 AU

By the 1930s the approach contributes more per decade than any refinement of the
tables — 1930→1955 gains a factor of 1.9 from geometry alone against theory
improving only 0.2″ → 0.15″. **The perturbation nearly doubles between 1900 and
1930 with no change of instrument.** An accelerating, direction-fixed,
non-Keplerian residual is not absorbable by any refit of the known planets.

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
h = r_p·v_p = 0.416 AU × 69.2 km/s, so the line of sight to the perturber sweeps
at h/D² ≈ **16 arcsec/yr at 280 AU** — eight arcminutes over a thirty-year
campaign, and rising as D⁻².

Combining them changes the character of the problem entirely. At 280 AU the local
escape speed from the Sun is 2.5 km/s. The measured speed is **10.3 km/s: four
times escape, seventeen times in energy.**

**The object is not a planet. It is unbound, and it is arriving.**

An analysis of that era, carrying 15% errors, puts arrival at 278 AU ÷ 2.2 AU/yr
≈ 128 years hence — **the mid-twenty-first century, ±20 years.** The scenario
acquires a deadline roughly a century before the deadline.

### 1961–1990 · ranging ends the argument

Radar ranging to Venus (1961) is marginal on arrival and decisive within a decade
as precision falls from ~50 km to ~1 km. Mariner 9 (1971) makes it a 500σ
measurement. Viking's seven-metre lander ranging (1976–82) makes it 10⁴σ.

By 1980, M, D, v and the impact parameter are pinned directly, with no degeneracy
left and no reliance on assumption. **The 2047 encounter is predictable in detail
from the Viking data set** — sixty-seven years out. The open question is no
longer whether but how close, and §8 of the assumptions document shows why that
matters: Earth's post-flyby semi-major axis moves ~0.011 AU per day of timing
shift, so the entire century turns on refining r_p from ~0.5 AU to a few per cent.

### 1990–2013 · three channels saturate

**VLBI.** The Sun's acceleration toward the hole produces a dipole proper-motion
pattern in the quasar frame of magnitude a/c — **710 µas/yr in 1995, 2,200 µas/yr
by 2020.** The analogous Galactic aberration drift, at 5.8 µas/yr, is already
measured and built into ICRF3. This is 100–400× larger and points straight at the
object.

**Ranging.** Cassini-class Saturn ranging (~25 m) against a ~10⁵ km signal is a
10⁶σ detection.

**Gaia.** The Einstein radius is θ_E = √(4GM D)/c ÷ D ≈ **1.3–1.5 arcseconds**,
which is enormous. Weak deflection falls as θ_E²/θ, so a background star a full
degree away is still displaced by ~570 µas against Gaia's ~25 µas — a coherent
astrometric distortion field tens of degrees across, dragged over the sky by the
object's ~2,000″ parallax. This is the moment the perturber stops being a
direction in a residual field and becomes **a point on the sky**. Genuine
microlensing events run from ~1/yr at high Galactic latitude to ~10³/yr against
the bulge.

### The detection table

Signal is the accumulated observable over a plausible campaign at each epoch;
ranging figures carry the suppression factor described in the appendix.

| year | D | best instrument | observable | signal | precision | **SNR** |
|---:|---:|---|---|---:|---:|---:|
| 1885 | 376 AU | meridian circle, Le Verrier tables | Uranus angle | 0.92″ | ~1″ | **0.9** |
| 1900 | 344 AU | photographic plates, Newcomb–Hill | Uranus angle | 1.21″ | ~0.3″ | **4** |
| 1930 | 278 AU | Lowell-era plates, refined theory | Uranus angle | 2.27″ | ~0.2″ | **11** |
| 1955 | 224 AU | machine numerical ephemerides | Uranus angle | 4.37″ | ~0.15″ | **29** |
| 1962 | 208 AU | Venus radar ranging | Venus range | 57 km | ~50 km | **1** |
| 1971 | 188 AU | Mariner 9 orbiter tracking | Mars range | 52 km | ~100 m | **5×10²** |
| 1980 | 168 AU | Viking lander ranging | Mars range | 292 km | 7 m | **4×10⁴** |
| 1995 | 135 AU | VLBI / ICRF | quasar aberration | 710 µas/yr | ~1 µas/yr | **7×10²** |
| 2010 | 100 AU | Cassini ranging | Saturn range | 1.1×10⁵ km | 25 m | **4×10⁶** |
| 2018 | 82 AU | Gaia DR2 | lensing deflection | 570 µas at 1° | 25 µas | **23** |
| 2030 | 53 AU | contemporary ranging | Saturn range | 4.2×10⁵ km | ~1 m | **4×10⁸** |
| 2040 | 27 AU | contemporary ranging | Saturn range | 7.8×10⁵ km | ~1 m | **8×10⁸** |

### October 2044 · closest approach to Uranus, 14.96 AU

The first of the encounters, and a distant one. Uranus is the only body whose
closest approach falls clearly before perihelion. The hole crossed Neptune's
orbital distance around 2036 and is now accelerating hard: it will cover the
final 8 AU in a single year.

### 2046 · 8.7 AU · why the elements still read normal

Here is the thing that is easy to state wrongly. **In January 2047 — seven months
before perihelion, with the black hole 4.5 AU from the Sun — Earth's semi-major
axis reads 1.000 AU and Mars's 1.523.** The orbital elements are, on paper,
entirely ordinary until the encounter itself.

That is true, and it is not in tension with anything above. The two describe
different quantities. A tidal field of this strength produces a large *positional*
perturbation while leaving the *secular elements* almost untouched: the Saturn
signal at 2010 is ~10⁵ km of accumulated range drift, which is 7×10⁻⁴ AU against a
semi-major axis of 9.5290 AU. An element table printed to four decimals shows
nothing. A ranging residual shows a 10⁶σ anomaly. Both are correct.

So the century's real asymmetry is this: for a hundred and forty years the solar
system is measurably falling toward something, and for a hundred and forty years
the planets keep their orbits to four decimal places. **Everything is known and
nothing has happened yet.**

---

## Act II — The passage (August 2047 – August 2048)

The dynamics that matter take eight days. The consequences take a year to finish
arriving.

### 12 August 2047 · Saturn at 9.71 AU

Distant. Saturn ends the encounter barely moved — 9.53 AU to 10.69.

### 19 August 2047 · perihelion, 0.416 AU · Mercury at 0.541 AU

The black hole reaches its closest approach to the Sun — **0.416 AU, well inside
Mercury's orbit** — moving at **69 km/s**. Mercury passes within 0.541 AU of it
the same day and, remarkably, survives: it ends at 0.338 AU with `e = 0.217`,
scarred but bound.

### 21–26 August 2047 · Venus, then Earth and the Moon

Venus at 0.683 AU on the 21st; the Moon at 0.851 AU on the 25th; Earth at
0.850 AU on the 26th.

**Earth's closest approach is 0.85 AU — no near miss, and decisive anyway.**
Within months its semi-major axis has gone from 1.000 to 1.644 AU and its
eccentricity from 0.017 to 0.426. The mechanism is not a close pass but a
**differential tug**: the hole pulls the Sun and the Earth by different amounts in
different directions, and what survives is the difference.

**The Moon comes through it.** It ends 374 000 km from Earth against today's
384 000 — still bound, still ordinary, orbiting a world about to freeze.

### 12 September 2047 · Mars at 1.827 AU

The most distant of the inner-planet encounters, and the one that matters most.
Mars is pushed *inward*: 1.524 AU to 1.372, with eccentricity rising from 0.093
to 0.222 and obliquity from 25.2° to 30.3°.

**Earth and Mars have swapped order.** Mars is now the inner of the two.

### 6 August 2048 · Jupiter at 0.603 AU

A year after perihelion, the departing hole passes closer to Jupiter than it did
to Earth. Jupiter is thrown from 5.21 AU to 6.21 with `e = 0.477` — the largest
proportional change of any surviving planet.

### The state at the end of 2048

| | Before | After |
|---|---|---|
| Mercury | 0.387 AU, e 0.206 | 0.338 AU, e 0.217 |
| Venus | 0.723 AU, e 0.007 | 0.933 AU, e 0.244 |
| **Earth** | 1.000 AU, e 0.017 | **1.644 AU, e 0.426** |
| **Mars** | 1.524 AU, e 0.093 | **1.372 AU, e 0.222** |
| Jupiter | 5.207 AU, e 0.048 | 6.211 AU, e 0.477 |
| Neptune | 30.02 AU, e 0.008 | **206 AU, e 0.864** — for now |

---

## Act III — Consequence (2049–2097)

### Earth freezes, in four years

Earth's year is now **770 days** — 2.11 calendar years — and it receives
139 W/m², about 41% of what it does today.

| Earth orbit | ≈ calendar | Global mean | Northern ice edge |
|---:|---:|---|---|
| 1 | 2049 | 255.5 K | **1° latitude** |
| **2** | **2051** | **216.3 K** | **pole to pole** |
| 3 | 2053 | 205.2 K | frozen |
| 5 | 2058 | 201.5 K | frozen |
| 100 | 2258 | **201.2 K** | equilibrium |

By the end of Earth's *first* orbit the permanent ice edge has reached the
equator's doorstep. **On the second orbit — 2051, four calendar years after the
encounter — the ice closes over the equator.** Everything after that is a frozen
planet settling toward its 201 K floor.

It is permanent. Escaping this state needs about 29% more sunlight than
present-day Earth receives; this Earth receives 41% of it. There is no path back
within the century, or any century.

*(Caveat: a thermal result. Oceans reach freezing on this schedule;
kilometre-thick ice sheets take far longer. The model has no ice-sheet dynamics.)*

### April 2058 · Neptune is taken

For a decade Neptune has been travelling an absurd orbit — bound, but reaching
past 200 AU and taking centuries to complete. It is still, formally, a planet of
the Sun.

Then the outbound black hole, 25 AU from the Sun and receding, passes within
**3.200 AU of it**. Neptune's eccentricity goes to 1.79 and it leaves.

The solar system loses a planet eleven years after the black hole's closest
approach, to an encounter most of a decade in the making.

### Mars, meanwhile, has a season

Mars's year is now 587 days, and its distance from the Sun swings between 1.068
and 1.676 AU. Obliquity has risen to 30.3°, and the longitude of perihelion has
shifted to 216° — which places **southern summer at perihelion**.

The result is not a warm planet. The global mean is 212.9 K, the equator's annual
mean 230 K, the winter pole 173 K. Thirty per cent of the atmosphere freezes onto
the winter cap and returns each year, cycling the surface pressure between 501 and
716 Pa. The atmosphere never collapses.

But at **34.5° south**, for **33 days of each Martian year**, the surface holds
between 273.2 and 273.7 K while the pressure sits near 700 Pa — and under those
conditions liquid water is thermodynamically permitted. It will not freeze, and it
will not boil.

It is the first time in the planet's history that has been true.

**The margin is razor-thin.** At 700 Pa water boils at 274.5 K. The entire window
— freezing to boiling — is under one and a half kelvin, and the surface spends 33
days inside it because its temperature happens to level off there, not because
anything holds it. A small error in the CO₂ inventory, the albedo or the heat
transport closes it completely.

And the window says only that liquid water *could* exist. Whether any is present
to do so is a question the model cannot answer.

### 2097

The black hole is 124 AU away and receding, carrying Neptune. Earth is a frozen
world at 201 K and will remain one. Jupiter runs an eccentric orbit it will keep
for the rest of the Sun's life. And Mars — smaller, colder, and now the third
planet rather than the fourth — has a month each year in its southern subtropics
during which water could, in principle, be liquid.

Within a single century the solar system has lost one planet, frozen another, and
handed a third a narrow, precarious opening. None of it was a surprise.

---

## What is and is not modelled

**Solid:** all dynamics (REBOUND, IAS15), every date and distance above, and the
post-flyby orbital elements, which are validated by reproducing present-day values
when applied at t = 0.

**Modelled with stated limits:** the climate outcomes. The Earth model omits
ice-sheet dynamics; the Mars model omits water entirely, along with the diurnal
cycle, topography and dust. See `CLIMATE_MODEL_REPORT.md` §10.

**Reconstructed, not simulated:** everything in Act I. The detection analysis is
analytic — tidal amplitudes, lensing geometry and aberration computed from the
trajectory, compared against representative period instrument precisions. It is
not the output of an ephemeris-fitting experiment, and the appendix states where
it is weakest.

**Not modelled at all:** tidal or relativistic effects, the asteroid belt
(`n_belt = 0` in this run), and everything about the biosphere — including,
pointedly, what a civilisation does with a hundred years of warning. Act I
establishes that the warning exists. What happens inside it is not a question this
repository can answer.

**And this is one run of 672.** It ranks 87th on orbital disruption — not an
extreme event — but it is one of only 32 that freeze the Earth, eject a planet and
open a Martian water window at the same time. A typical flyby in this sweep does
none of those. The century described here is a possibility, not a prediction.

---

## Appendix · detection method and assumptions

**Trajectory.** Radial infall, dr/dt = −√(v∞² + 2μ/r), with μ = G(M☉ + M_BH).
Reproduces the run's own distances to <1% beyond ~10 AU. The transverse component
is negligible at these distances (h = 0.018 AU²/day, per assumptions §8).

**Tidal amplitude.** Δa = 2GM·Δr/D³, the differential acceleration between Earth
and the target body. The common-mode fall of the whole system toward the hole is
unobservable and is excluded throughout.

**Suppression factor.** The free-drift estimate ½·Δa·T² overstates the
*unabsorbed* residual, because an ephemeris refit takes up part of the signal.
Calibrated against the Planet Nine literature: 10 M⊕ at 700 AU over a 13-year
Cassini baseline gives 747 m free-drift against a published signal of ~100 m and a
threshold of ~20 m, so **a factor of 8** is applied to every ranging figure. This
is the least secure number here. It is also nearly irrelevant: at a factor of 100
the Cassini-era SNR is still 10⁵.

**Astrometric scaling.** Anchored to assumptions §4 (0.92″, 1885, 376 AU) and
extrapolated as D⁻³.

**Two open problems in the source documents.**

*The astrometric normalisation is inconsistent.* Assumptions §3 gives 2.72″ at
434 AU; §4 gives 0.92″ at a ~408 AU window midpoint. Under the D⁻³ scaling that
§3's own table exhibits, these differ by a factor of ~3.5 and cannot both describe
this trajectory. Act I uses §4, since it is the calibration point from which `v∞`
is derived. If §3 is correct instead, every astrometric figure rises by 3.5× and
the 1885 threshold disappears — the object would be comfortably detectable then
rather than marginally so.

*The 1885 SNR of 0.9 rests entirely on the ~1″ theory-error estimate for
Le Verrier's tables.* At 0.5″ it becomes a 2σ detection at the design epoch and
the "just at the edge" framing does not survive. Both problems point the same way:
the astrometric normalisation wants re-deriving rather than inheriting.

**Instrument precisions** are representative period values, not results of specific
published analyses. Nineteenth-century meridian work is theory-limited rather than
measurement-limited; Venus radar improves from ~50 km (1961) to ~1 km (1970);
Viking lander ranging reaches ~7 m; Cassini ~20–30 m; Gaia end-of-mission ~25 µas.

**Not addressed:** the response of the Kuiper belt after the hole crosses
Neptune's distance around 2036, which is a further independent signal and is not
folded in.

**What cannot detect it:** imaging, at any epoch. Being a black hole makes the
object optically invisible while making it, at fixed orbit, far louder
gravitationally than any planet — 0.1 M☉ is 3,300× Planet Nine's mass at a sixth
the distance, a factor of ~7×10⁵ in M/D³. This is also why no repositioning
rescues the original "unnoticed approach" framing: suppressing the Cassini-era
signal below threshold needs M/D³ smaller by ~10⁶, which means an asteroid, or a
distance incompatible with a 2047 arrival. The same M/D³ that makes the flyby
consequential is what makes it observable, and both scale identically.
