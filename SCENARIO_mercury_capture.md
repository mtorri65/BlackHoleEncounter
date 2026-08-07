# Scenario study — the run where Mercury leaves with the black hole

A detailed case study of one simulation from the sweep in
[`CLIMATE_MODEL_REPORT.md`](CLIMATE_MODEL_REPORT.md). Everything here is a claim
about **this run**, not about flybys in general; §7 quantifies how atypical it is.

**Run:**

```
simulations/20260724_230314/20260724_230314__rp0p5__vinf25__inc30__toff56210__Om0__om0
```

A mirror-image twin exists at `…incm30__toff56210__Om180__om180` — the same
encounter geometry reflected through the ecliptic, producing bit-identical
results. Physically there is one scenario, recorded twice.

---

## 1. The encounter

| Parameter | Value |
|---|---|
| Black hole mass | 0.1 M☉ |
| Periapsis distance | **0.5 AU** |
| v∞ | 25 km/s |
| Inclination | 30° |
| Ω, ω | 0°, 0° |
| Time of periapsis | 56 210 days after epoch (≈ 2027) |
| Integration | 1873-09-01, 112 420 days (≈ 308 yr) |

Gravitational focusing matters: at 0.5 AU periapsis the hole is moving at
~67 km/s, so the encounter with the inner system lasts weeks, not years. From
the climate model's point of view — where the ocean mixed layer responds over
~1.6 years — the orbital change is effectively **instantaneous**.

The run ranks **3rd of 672** on the orbital-disruption metric (Score 8.36), with
two planets unbound: **Mercury and Venus**.

---

## 2. Mercury: captured

Of 386 planetary ejections across the whole sweep, **8 end up gravitationally
bound to the departing black hole** — about 2% of ejections, ~1% of runs. Two of
those are this scenario and its mirror.

Mercury's state at the end of the integration, relative to the black hole:

| | |
|---|---|
| Semi-major axis | **0.128 AU** |
| Eccentricity | **0.251** |
| Orbital period | **53 days** |
| Separation at t_end | 0.099 AU |
| Distance from the Sun | **820 AU** and receding |

This is a **properly circularised bound orbit**, not a marginal one — `e = 0.25`
is comfortably elliptical. Mercury is a satellite of a black hole, on a two-month
orbit, leaving the solar system.

The Venus captures elsewhere in the sweep are far more tenuous (`e` = 0.89 and
0.99). This is the one clean case.

**Caveat.** The capture test is a two-body energy check at the final timestep:
unbound from the Sun (`½v² − GM☉/r > 0`) *and* bound to the hole
(`½v²_rel − GM_BH/r_rel < 0`). It has not been verified by integrating further,
so long-term stability is assumed rather than demonstrated. The tight, modestly
eccentric orbit makes it plausible; it is not proof.

---

## 3. Earth: a snowball within two years

| | Before | After |
|---|---|---|
| Semi-major axis | 0.99981 AU | **1.4514 AU** |
| Eccentricity | 0.01717 | 0.2205 |
| Obliquity | 23.46° | 23.67° |
| Longitude of perihelion | 282.3° | 99.4° |
| Year length | 365 d | **639 d** |

Earth is thrown **outward by 45%**, receiving roughly half its former sunlight.
Note also that `λ_p` swings by 183° — the hemispheres exchange which one gets
summer at perihelion — though this is academic given what follows.

The transient, starting from Earth's pre-flyby equilibrium:

| Year | Global mean | NH ice edge |
|---|---|---|
| 1 | 263.0 K | 27.1° |
| **2** | **231.8 K** | **ice to the equator** |
| 5 | 194.4 K | — |
| 10 | 185.1 K | — |
| 20+ | **184.1 K** | equilibrium |

**Year 1 already has ice at 27° latitude** — the latitude of Florida, Cairo,
northern India. By **year 2 the planet is frozen pole to pole.**

The collapse is fast because both feedbacks push the same way: the insolation
drop cools the surface, ice forms, the brighter surface reflects more, and it
cools further. This is the Budyko–Sellers ice-albedo runaway, but not teetering
near a threshold — driven hard past it.

**There is no way back.** Escaping a snowball in this model requires ~29% *more*
sunlight than present-day Earth receives; this Earth gets about half. The state
is permanent, absent the carbonate–silicate thermostat the model omits.

**Timescale caveat.** Two years is a *thermal* result — the surface reaches
freezing temperatures that fast, which follows from the ocean mixed layer's time
constant and is credible. It does **not** mean kilometre-thick ice sheets exist
in year 2. The model has no ice-sheet dynamics, no snow accumulation, and no sea
ice thickness — only a temperature threshold flipping the albedo. The oceans
begin freezing within ~1–2 years; the geology takes far longer to catch up.

---

## 4. Mars: pulled inward, and made violently seasonal

| | Before | After |
|---|---|---|
| Semi-major axis | 1.5237 AU | **1.2948 AU** |
| Eccentricity | 0.0934 | **0.4346** |
| Obliquity | 25.18° | **38.34°** |
| Longitude of perihelion | 250.5° | 247.4° |
| Year length | 687 d | 538 d |

**Earth and Mars swap places.** Earth is flung out to 1.45 AU; Mars is pulled in
to 1.29 AU. After the encounter Mars is the inner of the two.

Climate consequences:

| | Today | After |
|---|---|---|
| Global mean | 203.4 K | **215.9 K** (+12.5) |
| Equatorial seasonal range | 23 K | **119 K** |
| Equatorial summer peak | 230.7 K | **310.1 K** |
| Peak temperature anywhere | 243.9 K | **392.7 K** |
| Surface pressure | 596 Pa | 288–737 Pa |
| Seasonal atmospheric swing | 22% | **61%** |
| Peak seasonal cap | 902 kg/m² | 978 kg/m² |

The mean warming (+12.5 K) is not the story. **The seasonal violence is.** With
`e = 0.435` the Sun–Mars distance swings from 0.73 to 1.86 AU — a **6.4× flux
ratio** — and Mars, with a thermal time constant of ~6.6 days, tracks that
forcing almost instantaneously rather than averaging it away as Earth's ocean
does.

**61% of the atmosphere now condenses and re-sublimates every year**, against
22% today. At aphelion Mars has less than half the air it has at perihelion.

---

## 5. Habitability of post-flyby Mars

The question this scenario was ultimately posed to answer. Computed at
`n_lat = 90` with 360 steps per year.

### The condition

Liquid water requires **both** conditions simultaneously:

```
T > 273.16 K     AND     p > 611.657 Pa      (water's triple point)
```

The pressure condition is easy to forget and is often the binding one. Mars sits
essentially *on* the triple point — its mean surface pressure is ~610 Pa — so
"above freezing" alone says nothing: below 611.657 Pa, ice sublimates directly
to vapour however warm the surface gets.

### Result

| | |
|---|---|
| Pressure above the triple point | **47% of the year** |
| Most favourable latitude | **−30°** (southern subtropics) |
| Window there | **79 of 538 days (14.7%)** |
| Peak daily-mean T in window | **338 K** |
| Pressure during window | 615–736 Pa |

By latitude:

| Latitude | Fraction of year | Peak daily-mean T |
|---|---|---|
| **−30°** | **14.7%** | 338.1 K |
| 0° | 13.3% | 310.0 K |
| +15° | 4.7% | 294.4 K |
| +30° | 0% | 280.1 K |
| +45° | 0% | 268.2 K |

The southern preference comes from `λ_p = 247°`, which places **perihelion in
southern summer**. Combined with 38° obliquity, that hemisphere receives the
close approach and the strong tilt together.

### Why it works: the conditions peak together

This is the interesting physical result. On a planet with a **condensing**
atmosphere, the two requirements are positively coupled:

> perihelion → warmer → polar caps sublimate → atmosphere thickens →
> pressure rises above the triple point

The atmosphere is thickest (736 Pa) exactly when the surface is warmest. A planet
with a non-condensing atmosphere gets no such help. **The 61% seasonal pressure
swing — which reads as a liability — is what makes the window possible at all.**

### Sensitivities

**The CO₂ inventory is load-bearing and the margin is thin.** The pressure
condition requires an inventory above `611.657 / 3.71 = 164.9 kg/m²`. The model
assumes **200 kg/m²** — a margin of only **21%**.

If Mars's exchangeable CO₂ is smaller than assumed, **the window does not shrink,
it closes**, because the pressure condition is then never satisfied anywhere at
any time. One unverified parameter separates "seasonally habitable" from "never".

**Factors that would widen it:**

| | Window (somewhere on the planet) |
|---|---|
| Daily means, pure water (as computed) | 14.7% |
| With a +30 K diurnal peak | 21.4% |
| NaCl brines (252 K) | 18.9% |
| **Perchlorate brines (199–206 K)** | **47.2%** |

Perchlorates are known to be present on Mars and are what the literature actually
discusses. With them, **pressure alone becomes the limiting condition**, satisfied
47% of the year.

---

## 6. What can and cannot be claimed

**Defensible:**

> This orbit moves Mars from a state where liquid water is essentially never
> thermodynamically permitted to one where it is permitted seasonally in the
> southern subtropics, because perihelion warming and CO₂ cap sublimation
> coincide. The conclusion depends critically on the assumed CO₂ inventory.

**Not defensible:** that there *would be* liquid water. The condition is
**necessary, not sufficient**. Specifically absent from the model:

1. **Any water at all.** No inventory, no transport, no vapour pressure. Whether
   water is present at −30° is not addressed.
2. **Evaporative cooling**, which at these pressures is violent and self-limiting.
3. **The diurnal cycle.** All temperatures are daily means; Mars's real diurnal
   range is 60–100 K, so daily maxima exceed these substantially.
4. **Topography.** The model is zonal-mean. Real Mars spans 20+ km of relief:
   Hellas at ~1200 Pa would be dramatically more favourable, the highlands at
   ~100 Pa impossible. **That variation is probably larger than the seasonal
   signal computed here, and is entirely invisible to this model.**
5. **Dust**, which dominates real Martian interannual variability.

The model is a **screening tool for comparing orbits**, not an instrument for
settling habitability. Its 654 runs share every simplification, so *differences*
between them are meaningful where absolute values are not.

---

## 7. How representative is this scenario?

**Not very** — which is precisely why the sweep exists.

| Quantity | This run | Percentile in sweep |
|---|---|---|
| Mars global mean T | 215.9 K | 72nd |
| Mars pressure swing | 61% | 76th |
| **Mars equatorial seasonal range** | **119 K** | **90th** |
| **Mars peak temperature** | **392.7 K** | **95th** |
| Mercury captured | yes | 1 of 8 across 672 runs |

Across the full sweep, **Mars warming is a coin flip**: 320 runs warmer, 338
cooler, median change **−1.7 K**. Reporting *"the flyby warms Mars while freezing
Earth"* as a general result would be a false generalisation drawn from a
single tail case.

There is also a systematic reason the two findings coincide here and rarely
elsewhere: **capturing Mercury requires the hole to pass deep into the inner
system**, which is the same geometry that flings Earth outward. The captures are
dynamically dramatic and climatically dull — in all 8, Earth ends as a snowball.
This run is the exception only in that Mars happens to land somewhere
interesting.

---

## 8. Reproduction

```bash
# Orbital elements for both planets (Earth uses the workbooks; Mars the Parquet tree)
python extract_earth_elements.py simulations/20260724_230314 --workers 5
python extract_mars_elements.py simulations/20260724_230314_parquet --body Mars --workers 5

# Earth climate (two-surface, Sellers OLR)
python climate_from_simulations.py simulations/20260724_230314 \
    --config input_climate.yaml --olr-model sellers --two-surface --workers 5

# Mars climate (condensing CO2 atmosphere)
python climate_mars_from_simulations.py simulations/20260724_230314 \
    --config input_mars.yaml --workers 5
```

The per-latitude habitability analysis of §5 is not part of the standard
pipeline; it was computed directly from `MarsEBM.run_year_co2` output using
`orbital_climate.mars.liquid_water_possible`.

---

## 9. Provenance

Every number here was computed from the simulation outputs, not estimated.
Validation of the underlying models is in `CLIMATE_MODEL_REPORT.md` §10 and its
appendix; the known limitations are §10, of which the ones most relevant to this
scenario are the fixed CO₂ inventory (§5 above), the daily-mean temperatures, and
the zonal-mean geometry.

The Mars model reproduces present-day Mars to: polar winter minimum 146.8 K
(observed ~148 K), surface pressure 596 Pa (~600 Pa), seasonal swing 22% (~25%),
and cap thickness 0.57 m (~0.5–1 m). The Earth model reproduces a 288.15 K global
mean and the 65 °N Milankovitch benchmark to 478.0 → 399.5 W/m².
