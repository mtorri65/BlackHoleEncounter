# A window on Mars

A single run from the 2047 sweep in which the black hole freezes the Earth,
takes Neptune, and leaves Mars — for about a month a year — with surface
conditions under which liquid water is thermodynamically permitted.

**Run:**

```
simulations/20260807_225002/20260807_225002__rp0p5__vinf10__inc60__toff63515__Om270__om180
```

Selected from 672 runs against three criteria: Mars made habitable, Earth and
Moon rendered uninhabitable, another planet ejected. Thirty-two runs satisfy all
three; this one gives Mars the longest liquid-water window of any of them.

Everything below is measured from the run's own output. The scenario
configuration and why each parameter has the value it does are in
[`SCENARIO_2047_assumptions.md`](SCENARIO_2047_assumptions.md); the narrative
version of this century is in [`SCENARIO_timeline.md`](SCENARIO_timeline.md).

---

## 1. The encounter

| | |
|---|---|
| Black hole mass | 0.1 M☉ |
| v∞ | 10 km/s |
| Requested periapsis | 0.5 AU |
| **Achieved periapsis** | **0.416 AU on 19 August 2047**, at 69.2 km/s |
| Orientation | i = 60°, Ω = 270°, ω = 180° (labels — see §7) |
| Rank on orbital disruption | **87th of 672** (Score 5.005) |

The achieved periapsis is *closer* than requested, not further. Every case in §8
of the assumptions document came out wider, because they all shared one orbital
orientation; this one does the opposite. See §7 below.

Closest approach of each body to the hole:

| body | distance | date |
|---|---:|---|
| Uranus | 14.96 AU | 2044-10-20 |
| Saturn | 9.71 AU | 2047-08-12 |
| **Mercury** | **0.541 AU** | 2047-08-19 |
| Venus | 0.683 AU | 2047-08-21 |
| Moon | 0.851 AU | 2047-08-25 |
| Earth | 0.850 AU | 2047-08-26 |
| Mars | 1.827 AU | 2047-09-12 |
| **Jupiter** | **0.603 AU** | **2048-08-06** |
| **Neptune** | **3.200 AU** | **2058-04-24** |

Two of those are worth noticing. Jupiter's closest approach comes **a year
after** perihelion, and it is closer than Earth's. Neptune's comes **eleven
years after**, on the hole's way out, and is what finally ejects it (§5).

---

## 2. What survives

| body | a before | a after | e after | fate |
|---|---:|---:|---:|---|
| Mercury | 0.387 | 0.338 | 0.217 | bound to Sun |
| Venus | 0.723 | 0.933 | 0.244 | bound to Sun |
| **Earth** | 1.000 | **1.644** | **0.426** | bound, frozen |
| Moon | — | — | — | **still bound to Earth**, 374 000 km |
| **Mars** | 1.524 | **1.372** | **0.222** | bound to Sun |
| Jupiter | 5.207 | 6.211 | 0.477 | bound to Sun |
| Saturn | 9.529 | 10.689 | 0.287 | bound to Sun |
| Uranus | 19.267 | 12.403 | 0.449 | bound to Sun |
| **Neptune** | 30.020 | **−17.568** | **2.462** | **ejected** |

**Earth and Mars swap order.** Earth moves out to 1.644 AU, Mars in to 1.372 —
Mars becomes the inner of the two, and the one with a temperate season.

The Moon survives. Its separation from Earth at the end of the integration is
374 000 km against today's 384 000; the pair is still bound, and goes into the
deep freeze together. Note that the heliocentric elements in
`__planets_run_deltas.csv` are misleading on this point: they give the Moon's
osculating *solar* orbit, which differs from Earth's by the Moon's ~1 km/s
orbital velocity, and can look like a separate object when it is not.

Nothing is captured by the black hole. Of the eight captures in the whole sweep,
every one is Mars — and a Mars orbiting a black hole receives no light at all,
so "captured" and "habitable" are mutually exclusive. That part of the original
brief cannot be satisfied.

---

## 3. Earth freezes, in four years

Earth's new orbit: **a = 1.644 AU, e = 0.426, year 770 days** (2.11 calendar
years), mean insolation **139 W/m²** against today's ~342.

Running the transient from the pre-flyby climate:

| orbit | ≈ calendar | global mean | NH ice edge |
|---:|---:|---:|---|
| 1 | 2049 | 255.5 K | **1° latitude** |
| 2 | 2051 | 216.3 K | **pole to pole** |
| 3 | 2053 | 205.2 K | frozen |
| 5 | 2058 | 201.5 K | frozen |
| 100 | 2258 | **201.2 K** | equilibrium |

**Ice reaches the equator on Earth's second orbit — within four calendar years
of the encounter.** The collapse is fast because both feedbacks push the same
way: less sunlight cools the surface, ice forms, the brighter surface reflects
more, and it cools further.

It is permanent. Escaping a snowball in this model needs roughly 29% *more*
sunlight than present-day Earth receives; this Earth receives 41%.

*(A thermal result only. The oceans reach freezing on this schedule;
kilometre-thick ice sheets take far longer. The model has no ice-sheet
dynamics.)*

---

## 4. Mars gets a window

| | before | after |
|---|---:|---:|
| a | 1.524 AU | **1.372 AU** |
| e | 0.093 | **0.222** |
| obliquity | 25.18° | **30.33°** |
| longitude of perihelion | 250.5° | 216.2° |
| year | 687 d | **587 d** |
| perihelion / aphelion | 1.38 / 1.67 AU | **1.068 / 1.676 AU** |

Climate at the new orbit:

| | |
|---|---:|
| global mean | 212.9 K |
| equator / pole (annual mean) | 230.0 K / 172.6 K |
| range anywhere in the year | 146.4 – 285.4 K |
| surface pressure | 501 – 716 Pa (30% of the atmosphere cycles) |
| peak polar cap | 841 kg/m² |
| atmospheric collapse | **none** — the atmosphere never freezes out |

**The result: 33 days of each 587-day Martian year, at 34.5° S, both the
temperature and the pressure permit liquid water.** A single contiguous window,
not a scatter of moments. Through it the surface holds at 273.2–273.7 K while
the pressure runs 672–716 Pa.

Note where it is. **At the equator the figure is zero** — the equator's annual
mean is 230 K and it never crosses freezing at all. The window exists only in
the southern subtropics, because obliquity of 30° combined with `λ_p = 216°`
puts southern summer near perihelion. An equator-only diagnostic reports nothing
happening here.

---

## 5. Neptune's delayed ejection

Neptune is not thrown out during the encounter. It is thrown out eleven years
later, and the intervening decade is the strangest part of the run.

| year | Neptune a | e |
|---:|---:|---:|
| 2040 | 29.31 AU | 0.028 |
| 2047 | 25.22 AU | 0.194 |
| 2049 | 205.7 AU | 0.864 |
| 2052 | 111.4 AU | 0.741 |
| 2055 | 123.7 AU | 0.756 |
| **2058** | **−40.9 AU** | **1.793** |
| 2150 | −17.60 AU | 2.490 |

The passage flings Neptune onto a wildly eccentric but still *bound* orbit
reaching past 200 AU. It stays there for a decade. Then, in **April 2058**, the
receding black hole — by now 25 AU from the Sun and heading out — passes within
**3.200 AU of Neptune** and removes it for good.

---

## 6. How typical is this?

It is not. Of 672 runs:

- **218 (35%)** give Mars a true liquid-water window at some latitude
- **208 of 654 (32%)** freeze the Earth into a snowball
- **112 (17%)** eject some planet other than Earth or the Moon
- **32** satisfy all three criteria together
- **0** also produce a capture

This run has the longest Martian window among those 32 (5.6% of the year); the
best in the entire sweep is 8.9%, in runs that fail the other two criteria. On
raw orbital disruption it ranks only **87th of 672** — the interest here is in
the *combination*, not the violence.

The run has a mirror twin, `…incm60__toff63515__Om90__om0`, which is the same
encounter reflected through the ecliptic and produces identical results.

---

## 7. Two ways the parameters lie

**Achieved periapsis is not `bh_rp_au`, and the sign of the error depends on
orientation.** Measured across rp = 0.5 runs in this sweep:

| i | Ω | ω | achieved rp | achieved date |
|---:|---:|---:|---:|---|
| 0 | 0 | 0 | 0.599 AU | 2047-09-03 |
| 30 | 0 | 0 | 0.608 AU | 2047-09-06 |
| 60 | 0 | 0 | 0.586 AU | 2047-08-30 |
| **60** | **270** | **180** | **0.416 AU** | **2047-08-19** |
| 90 | 270 | 180 | 0.422 AU | 2047-08-21 |

Ω = 0 orientations arrive ~0.09 AU wide; Ω = 270 orientations arrive ~0.08 AU
*narrow*. Section 8 of the assumptions document originally described this as a
uniform outward drift, which was an artefact of every sampled run sharing one
orientation. It has been corrected.

**The argument of periapsis is offset by 180° from its label.** The engine
builds the black hole's state with a negative radius and then negates the
velocity, so the orbit is the point inversion of the one its labels describe —
same plane, same `i`, same `Ω`, same shape and timing, but periapsis on the
opposite side. Harmless for sweep statistics, since the swept set is closed
under +180°, but the `om180` in this run's folder name means a true ω of 0°.

---

## 8. What is and is not modelled

**Solid.** All dynamics: REBOUND with IAS15, every date and distance in §1, and
the post-flyby elements, which reproduce present-day values when the same
extraction is applied at t = 0.

**Modelled, with stated limits.** The climate outcomes. Both models are
energy-balance models on a latitude grid: no clouds, no ocean circulation, no
weather. The Earth model has no ice-sheet dynamics. The Mars model has no water
cycle whatsoever — no reservoir, no evaporative cooling, no transport — and
reports daily means, while Mars's real diurnal range is 60–100 K.

**The liquid-water claim is necessary, not sufficient.** It says that for 33
days a year the surface at 34.5° S is above freezing *and* the pressure exceeds
water's saturation vapour pressure, so liquid would neither freeze nor boil. It
does **not** say water is present. Whether any survived four billion years and
the flyby is a question this model cannot address.

The margin is thin in a way worth stating plainly: during the window the surface
sits at 273.2–273.7 K and water boils at 274.5 K. **The whole window is under
1.5 K from both ends.** A modest error in the CO₂ inventory, the albedo, or the
diffusion coefficient closes it.

**Not modelled at all.** The detectability of the hole during approach (a
separate thread in this repository), tidal or relativistic effects, the asteroid
belt (`n_belt = 0`), and everything biological.
