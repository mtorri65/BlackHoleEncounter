# The solar system after the flyby

What the surviving system actually *looks like* once the black hole has gone —
the consequences that follow from the post-flyby orbital elements rather than
from the encounter itself.

[`SCENARIO_timeline.md`](SCENARIO_timeline.md) narrates the passage and tabulates
the elements at the end of 2048; [`SCENARIO_mars_window.md`](SCENARIO_mars_window.md)
covers the Mars habitability result. This document picks up where those stop and
asks what the numbers imply: which bodies are effectively gone, which orbits
intersect, whether the Earth–Moon system survives, and what the sky looks like
from the ground.

**Run:** `20260811_184731__rp0p75__vinf25__inc30__toff59132__Om0__om30`
(sweep `20260811_184731`, BH 0.1 M☉, achieved perihelion 0.662 AU on
2047-08-20). Every figure below is computed from that run's output — the
`bh_captures.csv` elements and, where relative geometry matters, the orbit log
itself.

---

## 1. The surviving system

Sorted by perihelion. Mercury is absent: it is captured by the black hole and
leaves with it.

| body | perihelion | aphelion | e | period |
|---|---:|---:|---:|---:|
| Venus | 0.79 | 1.19 AU | 0.200 | 1.0 yr |
| Mars | 1.05 | 1.63 AU | 0.217 | 1.5 yr |
| Earth | 1.13 | 1.84 AU | 0.236 | 1.8 yr |
| Uranus | 3.69 | 19.28 AU | 0.679 | 38.9 yr |
| Jupiter | 4.16 | 5.85 AU | 0.169 | 11.2 yr |
| **Saturn** | **9.71** | **3243.65 AU** | **0.994** | **65,608 yr** |
| Neptune | 14.49 | 45.34 AU | 0.515 | 163.6 yr |

Two things are immediately visible. The ordering of the planets is no longer the
ordering of their semi-major axes — Mars now comes inside Earth, and Uranus
inside Jupiter. And Saturn is no longer on a planetary orbit at all.

## 2. Saturn is expelled in all but name

Saturn goes from a = 9.54 AU to **a = 1627 AU, e = 0.994**. Its perihelion barely
moves (9.71 AU, close to where it started), but its aphelion is thrown out to
**3,244 AU** on a **65,600-year** period.

Formally it is still bound: `e < 1`, so an energy test classifies it as `solar`
and [`find_bh_captures.py`](find_bh_captures.py) does not count it among the
losses. Functionally it is gone. It leaves the planetary system on a
near-parabolic trajectory and does not return on any timescale the scenario
covers.

> **This is a trap in reading the capture census.** §6.9 of
> [`CLIMATE_MODEL_REPORT.md`](CLIMATE_MODEL_REPORT.md) reports 5.9% of
> body-outcomes as "free" and 0.09% as captured, using `e < 1` as the test. That
> is the right test for the question it asks — *is this body still gravitationally
> bound to the Sun* — but it is a poor proxy for *is this body still part of the
> solar system*. Saturn here is counted as retained.

So this run loses **two** bodies by two different mechanisms: Mercury is captured
by the black hole and departs with it; Saturn is flung onto an orbit that removes
it from the system without formally unbinding it.

## 3. Venus inherits Earth's orbit

Venus moves *outward* from 0.723 AU to **0.989 AU** — within 1% of Earth's
original semi-major axis — while Earth is thrown out to 1.486 AU. Their year
lengths swap character too: Venus now takes 0.983 years to orbit, almost exactly
Earth's old period.

| | before | after |
|---|---:|---:|
| a | 0.723 AU | **0.989 AU** |
| e | 0.007 | **0.200** |
| perihelion / aphelion | 0.719 / 0.728 AU | **0.791 / 1.186 AU** |
| annual-mean insolation | 1.911 × Earth-today | **1.044 × Earth-today** |

Venus's annual-mean insolation falls **45%**, landing essentially at Earth's
present value. The two planets exchange places in the most literal sense the
scenario allows.

**What this does to Venus's climate is not established here, and the project
cannot answer it.** The energy-balance model was only ever driven for Earth and
Mars, and it is the wrong instrument regardless: Venus's 735 K surface is held up
by a 92-bar CO₂ atmosphere in a runaway greenhouse, precisely the regime §10 of
the climate report documents as outside the model's range. Cutting insolation by
45% does not obviously reverse a runaway — the greenhouse is largely
self-sustaining once the water is gone — but that is an open question, not a
result. For scale only, ignoring the atmosphere entirely, the bare equilibrium
temperature at Venus's albedo of 0.77 falls from 227 K to 195 K.

## 4. The Earth–Moon system survives intact

Earth keeps the Moon. This was checked directly rather than inferred, because the
obvious inference is wrong:

- Specific orbital energy of the Moon relative to Earth is **negative at all
  7,374 output steps** — never unbound, not even transiently during the
  encounter.
- Final lunar semi-major axis **386,011 km**, against 384,400 km today.
- Earth–Moon separation stays within **356,539–406,910 km** across the whole
  integration, a normal perigee/apogee spread.
- Earth's Hill radius on its new orbit is 1.70 million km at perihelion —
  **4.4× the lunar distance**, so the Moon retains ample margin even at 1.486 AU.

> **The osculating-elements trap.** `bh_captures.csv` reports *heliocentric*
> elements for every body, including satellites. Earth's row reads a = 1.486 AU
> and the Moon's reads a = 1.445 AU, a difference of 6 million km, which looks
> like a stripped Moon. It is not. Those are instantaneous heliocentric elements
> computed from each body's own state vector, and the Moon carries ~1.02 km/s of
> orbital motion around Earth on top of Earth's ~26.9 km/s heliocentric velocity.
> Because semi-major axis depends on v² through vis-viva, a ~4% velocity
> difference yields a ~3% difference in apparent heliocentric `a`. **Any satellite
> will appear to be on a distinct orbit from its primary.** Testing whether a moon
> is retained requires relative state vectors from the orbit log; the census file
> cannot answer it.

## 5. The system is left orbit-crossing

The post-flyby configuration is not a stable set of nested orbits. **Seven of the
21 body pairs have overlapping radial ranges:**

| pair | overlap |
|---|---|
| Venus × Mars | 1.05 – 1.19 AU |
| Venus × Earth | 1.13 – 1.19 AU |
| Mars × Earth | 1.13 – 1.63 AU |
| Uranus × Jupiter | 4.16 – 5.85 AU |
| Uranus × Saturn | 9.71 – 19.28 AU |
| Uranus × Neptune | 14.49 – 19.28 AU |
| Saturn × Neptune | 14.49 – 45.34 AU |

Every pair in the inner system crosses. **Uranus is the principal destabiliser
in the outer system** — driven inward from 19.19 AU to 11.49 AU and pumped to
e = 0.679, its perihelion of 3.69 AU falls inside Jupiter's orbit while its
aphelion of 19.28 AU reaches Neptune's. It crosses three of the four surviving
giants.

Radial overlap is a necessary but not sufficient condition for encounters —
mutual inclinations and the relative phasing of the orbits matter, and neither is
analysed here. But it establishes that the arrangement is not obviously stable,
and it is the standard precondition for a system that will keep evolving through
close encounters.

**This is the scenario's hard time limit.** The integration runs 316 years
(115,342 days), which resolves the encounter and its immediate aftermath and
nothing beyond. Secular and encounter-driven evolution of a crossing system plays
out over 10⁵–10⁶ years. So every downstream claim — including the Mars
habitability window — carries an unstated ceiling: **the configuration that
produces it is not shown to be durable, and this project cannot show how long it
lasts.**

## 6. Total solar eclipses become universal

A pleasing consequence, and a clean qualitative change. Today's eclipses are
marginal — the Moon and Sun are so nearly the same apparent size that whether an
eclipse is total or annular depends on where both bodies sit in their orbits. The
flyby destroys that coincidence.

The Moon's orbit is untouched (§4), so its apparent diameter is unchanged. Earth
moves out, so **the Sun shrinks**:

| | Moon | Sun |
|---|---|---|
| today | 29.4′ – 33.5′ | 31.4′ – 32.5′ |
| after | 29.4′ – 33.5′ | **17.4′ – 28.2′** |

**The ranges no longer overlap.** The smallest possible Moon (29.4′, at apogee)
exceeds the largest possible Sun (28.2′, at Earth's perihelion) by 4%. The worst
case ratio goes from 0.903 today — Moon too small, hence annular eclipses — to
**1.042**. Annular eclipses become impossible; every central solar eclipse is
total.

The shadow geometry confirms this independently. The lunar umbral cone is 374,532
km long at 1 AU, *shorter* than the Moon's mean distance, which is the geometric
reason the umbra so often fails to reach the ground. On the new orbit the cone
reaches 425,042 km at Earth's perihelion and 688,217 km at aphelion — always
overshooting Earth.

| configuration | umbra width |
|---|---:|
| today, best case | ~167 km |
| after, Earth perihelion + Moon apogee (worst) | 148 km |
| after, Earth aphelion + Moon perigee (best) | **1,675 km** |

The *worst* case on the new orbit is about as good as today's best; the best case
is a shadow wider than Europe.

Totality also lasts far longer. The Moon must clear a margin of
(θ_moon − θ_sun)/2, which is 1.03′ today but **8.05′** at Earth's aphelion,
giving ~30 minutes of geocentric totality against today's 4.1. Scaling by the
~1.8× that Earth's rotation contributes — the factor that turns 4.1 geocentric
minutes into the 7.5 observed today — suggests **roughly 50 minutes of totality**.
Treat that as an order-of-magnitude estimate: it ignores the shadow-path geometry
a proper calculation would include.

Frequency is broadly unchanged per unit time. Earth's year stretches to 662 days,
so there are 23.2 lunations per orbit instead of 12.4, and eclipse seasons still
occur twice per orbit.

The finest eclipses in the solar system's history would go unwatched: Earth is a
snowball at 205 K by this point (§Act III of the timeline). And Mars, which *does*
get its habitability window under exactly these conditions, has no moon remotely
large enough to eclipse anything.

---

## What this document does and does not establish

**Established, from the run's output:** the post-flyby elements and everything
directly derivable from them — periods, apsides, insolation, radial overlaps,
apparent angular sizes, umbral geometry. The Earth–Moon binding, tested on
relative state vectors at every output step.

**Estimated:** totality duration (order-of-magnitude, as flagged above).

**Not established:** Venus's climate; whether the crossing configuration actually
produces encounters, and on what timescale; anything at all beyond the 316-year
integration. Mutual inclinations are not analysed, so §5 shows only that the
orbits overlap in radius, not that the bodies come near each other.
