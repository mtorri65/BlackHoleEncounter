# A window on Mars

A single run from the 2047 sweep in which the black hole freezes the Earth,
leaves with Mercury, and gives Mars — for about six weeks a year — surface
conditions under which liquid water is thermodynamically permitted.

**Run:**

```
simulations/20260811_184731/20260811_184731__rp0p75__vinf25__inc30__toff59132__Om0__om30
```

Selected from 4032 runs against three criteria: Mars made habitable, Earth and
Moon rendered uninhabitable, another planet ejected — ideally as a satellite of
the black hole. **This run satisfies all three including the last**, which no run
in the previous, coarser sweep could.

Everything below is measured from the run's own output. The configuration and
why each parameter has the value it does are in
[`SCENARIO_2047_assumptions.md`](SCENARIO_2047_assumptions.md); the narrative
version of this century is in [`SCENARIO_timeline.md`](SCENARIO_timeline.md).

---

## 1. The encounter

| | |
|---|---|
| Black hole mass | 0.1 M☉ |
| v∞ | 25 km/s |
| Requested periapsis | 0.75 AU |
| **Achieved periapsis** | **0.662 AU on 20 August 2047**, at 59.8 km/s |
| Orientation | i = 30°, Ω = 0°, ω = 30° (labels — see §7) |
| Rank on orbital disruption | 63rd of 4032 (Score 6.976) |

Closest approach of each body to the hole:

| body | distance | date |
|---|---:|---|
| Neptune | 19.46 AU | 2044-04-05 |
| Uranus | 16.17 AU | 2046-09-18 |
| Jupiter | 4.92 AU | 2047-05-10 |
| **Mars** | **0.810 AU** | 2047-06-28 |
| Venus | 1.210 AU | 2047-08-28 |
| Earth | 1.291 AU | 2047-09-10 |
| Moon | 1.290 AU | 2047-09-11 |
| **Mercury** | **0.025 AU** | **2048-07-03** |
| Saturn | 7.41 AU | 2048-11-25 |

Mercury's entry is not a flyby. By July 2048 it is already in orbit **around the
black hole**, and 0.025 AU is the periapsis of that orbit (§3).

---

## 2. What survives

| body | a before | a after | e after | fate |
|---|---:|---:|---:|---|
| **Mercury** | 0.387 | — | — | **captured by the black hole** |
| Venus | 0.723 | 0.989 | 0.200 | bound to Sun |
| **Earth** | 1.000 | **1.486** | **0.236** | bound, frozen |
| Moon | — | — | — | **still bound to Earth**, 402 000 km |
| **Mars** | 1.524 | **1.336** | **0.217** | bound to Sun |
| Jupiter | 5.206 | 5.011 | 0.170 | bound to Sun |
| **Saturn** | 9.573 | **1650.6** | **0.994** | bound, barely |
| Uranus | 19.299 | 11.485 | 0.679 | bound to Sun |
| Neptune | 30.090 | 29.927 | 0.516 | bound to Sun |

**Earth and Mars swap order.** Earth moves out to 1.486 AU, Mars in to 1.336 —
Mars becomes the inner of the two, and the one with a temperate season.

The Moon survives, at 402 000 km against today's 384 000, and freezes with Earth.
Note that `__planets_run_deltas.csv` is misleading on this point: it reports the
Moon's osculating *solar* orbit, which differs from Earth's by the Moon's
~1 km/s orbital velocity and can look like a separate object when it is not.

**Saturn is the near-miss.** At a = 1650 AU with e = 0.994 it is still formally
bound, but its semi-major axis wanders between 165 and 4100 AU over the century
after the encounter — a body the next passing star would remove without effort.

---

## 3. Mercury leaves with the black hole

Mercury crosses to a bound orbit around the hole on **26 August 2047**, six days
after perihelion. It does not escape afterwards.

| | |
|---|---|
| semi-major axis about the BH | **0.372 AU** |
| eccentricity | **0.923** |
| period | **0.72 years** |
| periapsis about the BH | 0.029 AU |
| first periapsis passage | 2048-07-03, at 0.025 AU |

**Mercury survives the passage intact.** At 0.025 AU it is 3.74 million km from
a 0.1 M☉ black hole whose Schwarzschild radius is 295 metres, and **14 times the
rigid-body Roche limit** of 260 000 km. It is not disrupted; it simply becomes a
moon of something dark, on a 0.72-year orbit, and leaves the solar system on it.

Of 4032 runs, 32 produce a capture — 0.8% — and every one of them captures
Mercury.

---

## 4. Earth freezes, in four years

Earth's new orbit: **a = 1.486 AU, e = 0.236, year 662 days** (1.81 calendar
years), obliquity essentially unchanged at 23.9°.

| orbit | ≈ calendar | global mean | NH ice edge |
|---:|---:|---:|---|
| 1 | 2049 | 259.6 K | **21.9° latitude** |
| **2** | **2051** | **216.1 K** | **pole to pole** |
| 3 | 2052 | 197.6 K | frozen |
| 5 | 2056 | 185.8 K | frozen |
| 80 | 2192 | **182.1 K (−91.1 °C)** | equilibrium |

**Ice closes over the equator on Earth's second orbit — four calendar years after
the encounter.** Both feedbacks push together: less sunlight cools the surface,
ice forms, the brighter surface reflects more, and it cools further.

It is permanent, because freezing and thawing have different thresholds. Earth
tips into the snowball once its annual-mean sunlight falls below about 0.95×
today's, but escaping again needs **1.28× — roughly 28% *more* sunlight than
present-day Earth receives**. The gap is the ice-albedo feedback, and the new
orbit delivers well under half.

*(A thermal result only. The oceans reach freezing on this schedule;
kilometre-thick ice sheets take far longer. The model has no ice-sheet
dynamics.)*

*(These figures use the **Sellers** nonlinear outgoing-radiation law, which is the
defensible choice in the frozen regime — see §10 and §6.8 of
[`CLIMATE_MODEL_REPORT.md`](CLIMATE_MODEL_REPORT.md). Under the simpler linear law
the same run reads 261.4 / 224.1 / 210.9 / 205.3 / **204.7 K**, some 23 K warmer at
the floor; that is the figure earlier editions quoted. The freeze schedule is
identical either way. The two insolation thresholds above are Sellers values too,
consistent with these temperatures; under the linear law they would be 0.91× and
1.17×.)*

---

## 5. Mars gets a window

| | before | after |
|---|---:|---:|
| a | 1.524 AU | **1.336 AU** |
| e | 0.093 | **0.217** |
| obliquity | 25.18° | 23.81° |
| longitude of perihelion | 250.5° | 256.8° |
| year | 687 d | **564 d** |
| perihelion / aphelion | — | **1.046 / 1.626 AU** |

Climate on the new orbit:

| | |
|---|---:|
| global mean | 217.5 K |
| equator / pole (annual mean) | 234.6 K / 171.3 K |
| **peak anywhere, all year** | **278.7 K** |
| surface pressure | 537 – 714 Pa (24.7% cycles) |
| peak polar cap | 881 kg/m² |
| atmospheric collapse | **none** |

**The result: 41 days of each 564-day Martian year, at 21.5° S, both temperature
and pressure permit liquid water.** It arrives in two stretches — 25 days and
16 days — with the surface between 273.2 and 275.0 K while water boils at
274.9–275.3 K.

That is 7.2% of the year, the **99th percentile** of the whole sweep, against a
maximum anywhere of 7.8%.

### Why a colder Mars is a wetter one

The obvious reading — that this run works because Mars gets warm — is wrong, and
the sweep says so. Among all 1250 runs with a liquid-water window:

| peak T anywhere | runs | median wet fraction |
|---|---:|---:|
| 273–280 K | 154 | 3.89% |
| **280–290 K** | 218 | **4.44%** |
| 290–310 K | 296 | 1.67% |
| 310–400 K | 444 | 1.11% |
| > 400 K | 138 | 1.11% |

The relationship is an inverted U, and the correlation of peak temperature with
wet fraction is **−0.46**. Past about 290 K a warmer Mars is a *drier* one.

The reason is the width of the window. At ~700 Pa water freezes at 273.2 K and
boils at 275.3 K — a band **two kelvin wide**. A surface that swings to 340 K
crosses it twice a year at speed. This run's Mars never exceeds 278.7 K
anywhere; it creeps to the edge of freezing and lingers there. **The window is
long because the planet is barely warm enough, not because it is warm.**

This also means the run is invisible to a naive diagnostic. Its equator-only
liquid-water fraction is **0.0%**, and so is the triple-point-only figure — the
equator never crosses freezing at all. The entire result lives at 21.5° S and
requires both looking off-equator *and* applying the boiling condition.

---

## 6. How typical is this?

Of 4032 runs:

- **1238** give Mars a true liquid-water window with the atmosphere intact
- **2142** leave Earth uninhabitable (frozen, runaway, or ejected)
- **1642** eject some planet other than Earth or the Moon
- **254** satisfy all three
- **4** also capture a planet — two distinct scenarios, this being the better

The run has a mirror twin at `…incm30__toff59132__Om180__om210`, identical to
six decimal places. That reflects an exact symmetry of the sweep,
**(i, Ω, ω) → (−i, Ω+180°, ω+180°)**, visible throughout the ranking: roughly
half of the 4032 runs are physically redundant.

The runner-up among the four is
`rp0p25__vinf25__inc90__toff59132__Om0__om15` — Mercury captured *and* Saturn
ejected, Earth frozen at 211 K, but Mars gets only 7 days of water.

---

## 7. Two ways the parameters lie

**Achieved periapsis is not `bh_rp_au`.** This run requested 0.75 AU and
achieved 0.662. The engine builds the BH's initial state as a two-body orbit
relative to the *Sun*, but over the infall it orbits the whole solar system,
whose mass sits at the barycentre; the Sun's ~13 m/s barycentric motion is never
removed, and for a nearly radial orbit that is a large fractional error in
angular momentum. Since r_p ∝ h², the effect is amplified. §8 of the assumptions
document has the measurement.

**The argument of periapsis is offset by 180° from its label**, because the
engine builds the state with a negative radius and then negates the velocity.
The orbit is the point inversion of the one the labels describe — same plane,
same i, Ω, shape and timing, periapsis on the opposite side. So this run's true
ω is 210°, not 30°.

---

## 8. What is and is not modelled

**Solid.** All dynamics: REBOUND with IAS15, every date and distance in §1–3,
and the post-flyby elements, which reproduce present-day values when the same
extraction is applied at t = 0.

**Modelled, with stated limits.** The climate outcomes. Both are energy-balance
models on a latitude grid: no clouds, no ocean circulation, no weather. The
Earth model has no ice-sheet dynamics. The Mars model has no water cycle at all
— no reservoir, no evaporative cooling, no transport — and reports daily means,
while Mars's real diurnal range is 60–100 K.

**The liquid-water claim is necessary, not sufficient.** It says that for 41 days
a year the surface at 21.5° S is above freezing *and* the ambient pressure
exceeds water's saturation vapour pressure, so liquid would neither freeze nor
boil. It does **not** say water is present.

The margin is thin, and §5 explains why that is structural rather than bad luck:
the surface sits at 273.2–275.0 K and boils at 274.9–275.3 K, so **the whole
window is under two kelvin wide**. A modest error in the CO₂ inventory, the
albedo or the diffusion coefficient closes it. What makes this run the best of
4032 is that its Mars sits inside that band for longer than any other, not that
it has any margin.

**Not modelled at all.** The detectability of the hole during approach (a
separate thread, and Act I of the timeline), tidal or relativistic effects, the
asteroid belt (`n_belt = 0`), and everything biological.
