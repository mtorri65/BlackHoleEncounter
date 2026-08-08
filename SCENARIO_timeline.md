# A century, 1977–2077

A narrative timeline of the scenario in
[`SCENARIO_mercury_capture.md`](SCENARIO_mercury_capture.md): a 0.1 M☉ black hole
passing 0.637 AU from the Sun, reaching perihelion on **19 August 2027**.
(0.5 AU is the *requested* periapsis; §8 of `SCENARIO_2047_assumptions.md`
explains why the achieved value differs.)

Every date and distance below is read from the simulation output. The framing is
narrative; the numbers are not.

---

## Act I — The approach (1977–2026)

### 1977 · 271 AU

The black hole is nine times further away than Neptune. Nothing about the solar
system is measurably different, and nothing will be for decades. Every orbital
element in the system reads exactly as it does today: Earth at 1.000 AU with
`e = 0.017`, Mars at 1.524 with `e = 0.093`, Mercury at 0.387 with `e = 0.206`.

The object emits nothing. Its only observable signature would be gravitational —
the astrometric deflection of background stars as it drifts across the sky. That
detection problem is the subject of a separate analysis thread in this repository
and is not modelled here.

### 1997 · 165 AU · 2007 · 111 AU

Half a century of approach at roughly 5 AU per year. The planets' orbital
elements remain flat to four decimal places. Whatever is going to happen has not
started.

### 2017 · 58 AU · 2022 · 30 AU

The hole crosses Neptune's orbital distance around 2022. It is now inside the
solar system in the loosest sense. Planetary elements still show nothing —
through mid-2027 the numbers are indistinguishable from unperturbed.

### 2025 · 13.6 AU

Past Saturn's distance. Two years remain.

### September 2026 · closest approach to Saturn, 6.35 AU

The first of the encounters, and a distant one — and by a wide margin the
earliest. Saturn, Venus (31 July 2027) and Jupiter (5 August) all reach closest
approach *before* the hole reaches perihelion on 19 August; Earth, Mercury and
Mars all do so after. Which side of perihelion a planet falls on is simply a
matter of where it sits in its orbit as the hole sweeps through.

### Mid-2026 · 7.6 AU · Early 2027 · inside Jupiter's orbit

Twelve months out. Gravitational focusing is now accelerating the hole hard: it
covers the last 7 AU in a year, and will cover the final 1 AU in weeks.

Even now the orbital elements have barely moved. **In May 2027 — three months
before perihelion — Earth's semi-major axis reads 0.999 AU and its eccentricity
0.010.** The solar system is, on paper, still ordinary.

---

## Act II — The passage (July–October 2027)

Everything happens in about four months.

### 31 July 2027 · Venus at 0.059 AU

Venus passes within **0.059 AU** of the hole, about 15 times the Earth–Moon
distance. Venus does not survive the encounter as a planet of the Sun: within
weeks its heliocentric eccentricity exceeds 4, and it is gone.

### August 2027 · Jupiter at 4.62 AU (5th) · perihelion, 0.637 AU (19th)

Jupiter goes first, and distantly — though its mass means it perturbs the hole
nearly as much as the reverse.

On **19 August** the black hole reaches its own closest approach to the Sun:
**0.637 AU, well inside Earth's orbit**, moving at 61 km/s.

That 0.637 AU is not the 0.5 AU this run requested, and the difference is not an
error. `bh_rp_au` fixes the periapsis of the *initial two-body orbit* at the 1873
epoch; across the 154-year infall the achieved closest approach drifts outward
and late. Every run in this sweep does it, by a consistent amount — requested
0.25 → 1.50 AU arrive at 0.385, 0.637, 0.878, 1.116, 1.355 and 1.594 AU, all
22–24 days after the nominal date of 26 July. §8 of
[`SCENARIO_2047_assumptions.md`](SCENARIO_2047_assumptions.md) has the details.

### September 2027 · Earth at 1.127 AU

Earth's closest approach is more than an astronomical unit — no near miss. It is
nonetheless decisive. Within a month Earth's semi-major axis has gone from
0.999 to **1.400 AU** and its eccentricity from 0.010 to **0.253**.

The mechanism is not a close pass but a **differential tug**: the hole pulls the
Sun and the Earth by different amounts in different directions, and what survives
is the difference.

### October 2027 · Mercury at 0.092 AU (7th), Mars at 0.104 AU (11th)

The two closest planetary encounters of the entire event, four days apart.

**Mercury does not merely get ejected — it is captured.** It leaves bound to the
black hole on an orbit of `a = 0.128 AU`, `e = 0.251`, period 53 days. Of 386
planetary ejections across the full 672-run sweep, only eight end this way; this
is the cleanest of them.

**Mars is thrown inward.** Its semi-major axis first overshoots to 2.56 AU in
September, then settles by mid-2028 at **1.295 AU** with eccentricity 0.429.

### The state by the end of 2027

| | Before | After |
|---|---|---|
| Mercury | 0.387 AU, e 0.206 | **captured by the black hole** |
| Venus | 0.723 AU, e 0.007 | **ejected**, e > 4 |
| Earth | 0.999 AU, e 0.017 | **1.455 AU, e 0.220** |
| Mars | 1.524 AU, e 0.093 | **1.295 AU, e 0.429** |

**Earth and Mars have swapped places.** Earth is now the outer of the two.

---

## Act III — Departure and consequence (2028–2077)

The hole recedes on the mirror image of its approach: 6.8 AU by mid-2028,
30 AU by 2032, 111 AU by 2047, and **270 AU by 2077** — back where it was in 1977,
now with Mercury and Venus in tow.

The dynamics are finished within months. **The climate takes decades**, and this
is where the century's real story is.

### Earth freezes

Earth's year is now **639 days**. The transient below is in calendar time.

| | Global mean | Northern ice edge |
|---|---|---|
| **2029** | 263.0 K | **27° latitude** |
| **2031** | 231.8 K | **pole to pole** |
| 2033 | 212.0 K | frozen |
| 2038 | 190.5 K | frozen |
| 2045 | 185.1 K | frozen |
| **2062** | **184.1 K** | equilibrium |

By 2029 — two years after the passage — permanent ice reaches the latitude of
Cairo and northern India. **By 2031 the ice reaches the equator.** The remaining
three decades merely cool a frozen planet toward its 184 K floor.

The collapse is fast because both feedbacks push together: less sunlight cools
the surface, ice forms, the brighter surface reflects more, and it cools further.

**It is permanent.** Escaping a snowball in this model requires roughly 29% *more*
sunlight than present-day Earth receives. This Earth receives about half. There
is no path back within the century, or any century.

*(Caveat: this is a thermal result. The oceans reach freezing temperatures on
this schedule; kilometre-thick ice sheets take far longer. The model has no
ice-sheet dynamics.)*

### Mars becomes seasonally habitable

Mars, meanwhile, has moved the other way. Its year is now 538 days, and with
`e = 0.429` its distance from the Sun swings between 0.74 and 1.85 AU — a **6.4×
flux ratio**. With a thermal time constant of about a week, Mars tracks that
swing almost instantaneously.

The result is a planet of violent extremes: a mean temperature 12 K warmer than
before, an equatorial seasonal range of **119 K**, and peak surface temperatures
of **310 K at the equator and 338 K in the southern subtropics**.

**Sixty-one per cent of the atmosphere now freezes onto the winter pole and
returns each year**, against 22% today — so the surface pressure cycles between
288 and 737 Pa.

That cycle is what matters. Because perihelion warms the surface *and* sublimates
the caps, the atmosphere is thickest exactly when the planet is warmest — and
during roughly **79 days of each Martian year, at around 30° south, both the
temperature and the pressure exceed water's triple point simultaneously.**

For the first time in its history, liquid water is thermodynamically permitted on
Mars — seasonally, in one hemisphere, for about a seventh of the year.

*(Whether water would actually be present is a different question the model
cannot answer; see §6 of the scenario study. And the window depends on an assumed
CO₂ inventory with only a 21% margin — below it, the window closes entirely.)*

### 2077

The black hole is 270 AU away and receding, carrying Mercury. Venus is gone.
Earth is a frozen world at 184 K, and will remain one. Mars has an atmosphere
that collapses and rebuilds every year, and a southern summer in which water
could — thermodynamically — be liquid.

Within a single century the solar system has lost two planets, frozen a third,
and given a fourth its first opportunity for liquid water.

---

## What is and is not modelled

**Solid:** all dynamics (REBOUND, IAS15), every date and distance above, and the
post-flyby orbital elements, which are validated by reproducing present-day
values when applied at t = 0.

**Modelled with stated limits:** the climate outcomes. The Earth model omits
ice-sheet dynamics; the Mars model omits water, the diurnal cycle, topography and
dust. See `CLIMATE_MODEL_REPORT.md` §10.

**Not modelled at all:** the detectability of the hole during the approach, any
tidal or relativistic effects, the Moon's fate, the asteroid belt (`n_belt = 0`
in this run), and everything about the biosphere.

**And this is one run of 672.** It sits in the tail — 90th percentile for Martian
seasonal amplitude, 95th for peak temperature, and one of only eight capture
events. A typical flyby in this sweep does far less. The century described here is
a possibility, not a prediction.
