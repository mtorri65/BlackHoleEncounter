# A window on Mars

A single run from the 2047 sweep in which the black hole freezes the Earth,
leaves with Mercury, and gives Mars — for about six weeks a year — surface
conditions under which liquid water is thermodynamically permitted.

**Run:**

```
simulations/20260811_184731/20260811_184731__rp0p75__vinf25__inc30__toff59132__Om0__om30
```

*(That sweep's raw runs are no longer on disk. Figures from §1 on come from a
re-run of the same configuration with every day logged,
`simulations/20261107_205226/20261107_205226__rp0p75__vinf25__inc30__toff59132__Om0__om30`,
made on 2026-11-07 under REBOUND 5.2.2 with the engine's Earth and Moon masses
corrected to the de440 values (see §9). It reproduces §2's elements — Earth–Moon
barycentre 1.486, Mars 1.336, Saturn 1650.6 AU. The sky positions, launch
windows and noon times were first computed on an earlier daily re-run
(`20261030_232613`, old masses); planetary positions in the two agree to within
40 000 km through 2055, so those figures stand, and only the Moon's position
about Earth and the hours of Mercury's close passes changed. The hours in §1 and
§4 come from a further re-integration of September 2047 at 10-minute steps.)*

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
| Orientation | i = 30°, Ω = 0°, ω = 30° (labels — see §8) |
| Rank on orbital disruption | 63rd of 4032 (Score 6.976) |

Closest approach of each body to the hole:

| body | distance | date |
|---|---:|---|
| Neptune | 19.46 AU | 2044-04-05 |
| Uranus | 16.17 AU | 2046-09-18 |
| Jupiter | 4.92 AU | 2047-05-10 |
| **Mars** | **0.810 AU** | 2047-06-28 |
| Venus | 1.210 AU | 2047-08-28 |
| Earth | 1.291 AU | 2047-09-10, 05:43 UTC |
| Moon | 1.294 AU | 2047-09-10, 12:25 UTC |
| **Mercury** | **0.025 AU** | **2048-07-03** |
| Saturn | 7.41 AU | 2048-11-25 |

Mercury's entry is not a flyby. By July 2048 it is already in orbit **around the
black hole**, and 0.025 AU is the periapsis of that orbit (§3).

Earth's closest approach is the nearest of the whole encounter: 1.2913 AU
(193 million km) at a relative speed of 68.6 km/s. The Moon, on the far side of
Earth from the hole that morning, passes 0.0024 AU (360 000 km) farther out,
seven hours later. *(Earlier editions gave Earth 05:34 UTC and the Moon 1.290 AU
at 19:10 UTC; both came from a run whose Moon was 22 days out of phase — §9.)*

### Seen from Earth

**Mars passes behind the Sun on the day it passes the hole.** Its elongation —
its angle from the Sun as seen from Earth — bottoms out at **0.87° on 27–28 June
2047**, within a day of its 0.810 AU closest approach to the hole. It is on the
far side of the Sun, 2.59 AU from Earth, at magnitude +1.5, and it clears the
Sun's disc (0.26° in radius) rather than passing behind it.

| date | elongation | side | constellation |
|---|---:|---|---|
| 14 Jun | 3.66° | east (evening) | Taurus |
| 18 Jun | 2.64° | east | Taurus |
| 21 Jun | 1.91° | east | Gemini (enters) |
| 24 Jun | 1.25° | east | Gemini |
| 26 Jun | 0.94° | east | Gemini |
| **27 Jun** | **0.87°** | east | Gemini |
| **28 Jun** | **0.87°** | west (morning) | Gemini |
| 30 Jun | 1.09° | west | Gemini |
| 2 Jul | 1.47° | west | Gemini |
| 4 Jul | 1.92° | west | Gemini |
| 7 Jul | 2.64° | west | Gemini |
| 10 Jul | 3.39° | west | Gemini |

Mars is within 2° of the Sun from 21 June to 4 July and within 5° from 9 June
to 16 July, so it is lost in the glare for about five weeks. Throughout, the
hole is 18–20° from Mars on the sky.

**Mercury goes where Mercury cannot go.** Already bound to the hole by
September (§3), it climbs the evening sky with it:

| date | elongation (east) | RA / Dec (J2000) | constellation | distance from Sun | mag |
|---|---:|---|---|---:|---:|
| 1 Sep | 21.6° | 11h57m / −0.3° | Virgo | 0.50 AU | +0.2 |
| 5 Sep | 24.4° | 12h20m / −3.5° | Virgo | 0.54 | +0.4 |
| 10 Sep | 27.8° | 12h48m / −7.6° | Virgo | 0.59 | +0.7 |
| 15 Sep | 31.2° | 13h16m / −11.5° | Virgo | 0.65 | +0.9 |
| 20 Sep | 34.5° | 13h44m / −15.2° | Virgo | 0.72 | +1.2 |
| 25 Sep | 37.7° | 14h13m / −18.5° | Virgo | 0.79 | +1.4 |
| 30 Sep | 40.7° | 14h41m / −21.4° | Libra (from 27 Sep) | 0.86 | +1.6 |

Mercury's orbit normally keeps it within about 28° of the Sun. It passes that
limit around **11 September** and keeps going, 12–15° from the hole all month.
Any observer who knows the sky can see this; nothing else in it says so as
plainly.

*(Positions are J2000. Coordinates of date in 2047 differ by about 0.65° of
precession. Magnitudes use the standard phase law for each planet and are
approximate.)*

---

## 2. What survives

| body | a before | a after | e after | fate |
|---|---:|---:|---:|---|
| **Mercury** | 0.387 | — | — | **captured by the black hole** |
| Venus | 0.723 | 0.989 | 0.200 | bound to Sun |
| **Earth** | 1.000 | **1.486** | **0.236** | bound, frozen |
| Moon | — | — | — | **still bound to Earth**, a ≈ 383 400 km |
| **Mars** | 1.524 | **1.336** | **0.217** | bound to Sun |
| Jupiter | 5.206 | 5.011 | 0.170 | bound to Sun |
| **Saturn** | 9.573 | **1650.6** | **0.994** | bound, barely |
| Uranus | 19.299 | 11.485 | 0.679 | bound to Sun |
| Neptune | 30.090 | 29.927 | 0.516 | bound to Sun |

**Earth and Mars swap order.** Earth moves out to 1.486 AU, Mars in to 1.336 —
Mars becomes the inner of the two, and the one with a temperate season.

The Moon survives, its orbit essentially untouched, and freezes with Earth. Its
semi-major axis about Earth averages 383 400 km over 1990–2000 and 383 410 km over
2100–2200. *(Earlier editions gave "402 000 km against today's 384 000", which
read as a widened orbit. That figure was just the Moon's distance at the run's
final step, a point on an orbit that ranges over 356 000–407 000 km.)*
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
| apoapsis about the BH | 0.716 AU |
| first periapsis passage | 2048-07-03, at 0.025 AU |

### Distance from the hole

In the first year Mercury drifts out to 0.74 AU, as far as it ever gets, and
falls back for its first close pass:

| date | distance | million km | range during that month |
|---|---:|---:|---|
| 20 Aug 2047 (hole's perihelion) | 0.325 AU | 48.6 | 0.283–0.325 |
| 1 Sep 2047 | 0.283 | 42.3 | 0.283–0.383 |
| 1 Oct | 0.388 | 58.0 | 0.388–0.521 |
| 1 Nov | 0.525 | 78.5 | 0.525–0.627 |
| 1 Dec | 0.629 | 94.1 | 0.629–0.699 |
| 1 Jan 2048 | 0.701 | 104.9 | 0.701–0.734 |
| 1 Feb | 0.735 | 110.0 | 0.732–0.738 (farthest, 12 Feb: 0.738) |
| 1 Mar | 0.731 | 109.4 | 0.687–0.731 |
| 1 Apr | 0.684 | 102.3 | 0.593–0.684 |
| 1 May | 0.589 | 88.1 | 0.423–0.589 |
| 1 Jun | 0.415 | 62.1 | 0.088–0.415 |
| 1 Jul | 0.067 | 9.9 | 0.025–0.382 (**closest, 3 Jul: 0.025 AU**) |
| 1 Aug 2048 | 0.390 | 58.3 | 0.390–0.569 |

After that the orbit settles and does not change. From 2050 to the end of the
run in 2201, Mercury swings from **0.0285 AU (4.3 million km) to 0.716 AU
(107 million km)** and back every **262 days**. Its average distance is 0.53 AU
(79 million km). Every year from 2049 on spans that full range at least once, and
the distance on any given calendar date is just a matter of where Mercury is in
its orbit.

The first close passes, for dating scenes:

| closest to the hole (UTC) | distance | farthest from the hole | distance |
|---|---:|---|---:|
| 3 Jul 2048, ~04:00 | 0.0250 AU | 12 Feb 2048 | 0.738 AU |
| 24 Mar 2049, ~04:00 | 0.0272 | 12 Nov 2048 | 0.719 |
| 12 Dec 2049, ~00:00 | 0.0278 | 3 Aug 2049 | 0.717 |
| 31 Aug 2050, ~14:00 | 0.0281 | 22 Apr 2050 | 0.717 |
| 21 May 2051, ~02:00 | 0.0282 | 10 Jan 2051 | 0.716 |
| 7 Feb 2052, ~13:00 | 0.0283 | 29 Sep 2051 | 0.716 |
| 26 Oct 2052, ~23:00 | 0.0284 | 18 Jun 2052 | 0.716 |

Seen from Earth, the hole is 95° west of the Sun (morning sky) at the 24 March 2049
pass and 64° east (evening sky) at the December 2049 pass, both times in Aquarius,
10.6 and 15.6 AU away. At those distances Mercury's whole orbit is tiny on the sky:
at a close pass it is at most 9′ (March) or 6′ (December) from the hole.

Later close passes follow every 262.6 days at 0.0285 AU, 214 of them by 2201. The
closest distance creeps out from 0.025 to 0.0285 AU over the first few orbits as
the Sun's pull fades. Meanwhile the hole, with Mercury, moves away from the Sun
at about 5.3 AU a year: 3 AU on 1 January 2048, 32 AU by 2053, 283 AU by 2100 and
817 AU by 2201.

*(Distances are from the daily re-run noted under the run path. Close passes are
taken from the orbit's shape, since at periapsis Mercury moves about 0.05 AU a
day, too fast for daily samples to catch the minimum: each pass is timed from
Mercury's orbit about the hole a day and two days beforehand. The two estimates
agree to the minute, but the model itself does not pin the hour that finely:
correcting Earth's and the Moon's masses by 0.017% (§9) moved the passes 1–6
hours later, the shift growing with each orbit, because Mercury's capture
amplifies small differences. Treat the hours as good to several hours, and a
pass near midnight (12 Dec 2049) as able to fall on either date. Distances did not
change. For late 2047 that
orbit-shape calculation is not reliable: the Sun, still close to the hole,
distorts it so much that it predicts close passes of 0.008–0.016 AU that never
happen. The measured distance rises steadily through those months.)*

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

### The year lengthens

Earth's year grows during the encounter rather than at a single instant. These are
osculating elements — the ellipse Earth would follow if the hole vanished at that
moment — for the Earth–Moon barycentre. Earth's centre alone reads about
0.0017 AU higher, because of the Moon's monthly wobble.

| | a | year | e |
|---|---:|---:|---:|
| 1 Jun 2047 | 1.009 AU | 370 d | |
| 20 Aug (hole at perihelion) | 1.191 | 475 d | |
| 10 Sep 00:00 UTC | 1.3756 | 589.3 d | 0.237 |
| 10 Sep 12:00 | 1.3784 | 591.1 d | 0.238 |
| 11 Sep 00:00 | 1.3811 | 592.8 d | 0.239 |
| 11 Sep 12:00 | 1.3837 | 594.5 d | 0.240 |
| 12 Sep 00:00 | 1.3862 | 596.1 d | 0.241 |
| 1 Oct | 1.435 | 628 d | |
| mid-2048 onward | 1.487 | 662 d | 0.236 |

Over 10–11 September the semi-major axis grows by about 0.0055 AU a day, and
**the year lengthens by about 0.15 days every hour** — 3.4 days per day. The
growth is already slowing: it peaked on 25 August, at 0.0102 AU a day, and the
closest approach on 10 September comes after the fastest change, not with it.

Earth's distance from the Sun on the new orbit runs from **1.134 AU at perihelion
to 1.839 AU at aphelion**. The first aphelion is **15 July 2048** (about 08:00 UTC,
1.8385 AU), the first perihelion **12 June 2049** (1.134 AU); after that they
alternate every 662 days (aphelion 9 May 2050, perihelion 4 April 2051, …).

### Noon leaves the clock

Earth's spin is untouched, so a 24-hour clock still keeps time with the stars. But
the Sun now takes 662 days, not 365, to work its way round the sky against them,
and the clock assumes the old rate. **Each day the Sun reaches its highest point a
little earlier by the clock: 1.76 minutes a day on average, so noon goes round the
whole clock face once every 2.2 calendar years.** The rate swings with Earth's
distance from the Sun, from about 2½ minutes a day near aphelion to about half a
minute near perihelion.

Solar noon at Bristol (51.45° N, 2.59° W), by the clock:

| date | solar noon (UTC) | noon height | Sun's diameter | sunlight vs today |
|---|---:|---:|---:|---:|
| 12 Jun 2047 | 12:09 | 61.7° | 31.6′ | 98% |
| 12 Sep 2047 | 12:04 | 42.9° | 28.9′ | 81% |
| 12 Dec 2047 | 10:25 | 18.5° | 22.8′ | 51% |
| 12 Mar 2048 | 07:38 | 15.1° | 19.0′ | 35% |
| 12 Jun 2048 | 03:59 | 22.1° | 17.5′ | 30% |
| 12 Sep 2048 | 23:56 | 33.9° | 17.8′ | 31% |
| 12 Dec 2048 | 20:11 | 48.2° | 19.8′ | 38% |
| 12 Mar 2049 | 17:34 | 61.5° | 24.2′ | 57% |
| **12 Jun 2049** | **16:46** (17:46 BST) | **50.9°** | **28.2′** | **78%** |
| 22 Jun 2049 | 16:41 | 47.6° | 28.1′ | 77% |

By September 2048 the midday Sun stands over Bristol at midnight by the clock. In
June 2049, at perihelion, the Sun is as large as it ever gets on the new orbit,
28.2′ against a normal June's 31.6′, and it stands where it normally would in late
April or mid-August, sinking a degree and a half a week: the calendar no longer
tracks the seasons either. Whether the story's society keeps the old clock or
resets it to the Sun is a choice for the narrative; these times assume the old
24-hour clock (UTC). Through the encounter itself the drift is still only a few
minutes, so sunrise and sunset times in 2047 are close to normal.

*(Sun positions from the daily re-run; heights are for the Sun's centre with
standard refraction, which adds only 0.02° at these heights.)*

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

## 6. Getting to Mars

How often a ship can go from Earth to Mars, before the encounter and after.
Every figure here comes from [`launch_windows.py`](launch_windows.py), run on the
daily re-run noted under the run path. For each departure day and flight time it
solves for the Sun-centred transfer arc between Earth on the departure day and
Mars on the arrival day. Each arc is costed as the burn out of a 300 km circular
Earth orbit plus the burn into a 300 km circular Mars orbit, all propulsive.
**On that measure today's windows cost 5.7–6.2 km/s**, and the search reproduces
the real ones from the pre-encounter years of the same run.

### Before the encounter

| best departure | window (within 0.5 km/s) | flight | arrives | cost | Mars leads Earth by |
|---|---|---:|---|---:|---:|
| 30 Oct 2026 | 24 Sep – 27 Nov | 312 d | Sep 2027 | 5.67 km/s | 67° |
| 24 Nov 2028 | 19 Oct – 20 Dec | 300 d | Sep 2029 | 5.87 | 74° |
| 24 Dec 2030 | 16 Nov – 29 Jan | 284 d | Oct 2031 | 6.22 | 73° |
| 18 Apr 2033 | 23 Mar – 16 May | 200 d | Nov 2033 | 6.05 | 29° |
| 27 Jun 2035 | 6 May – 29 Jul | 202 d | Jan 2036 | 5.74 | 29° |
| 15 Aug 2037 | 30 Jun – 24 Sep | 346 d | Jul 2038 | 6.09 | 36° |
| 18 Sep 2039 | 15 Aug – 22 Oct | 340 d | Aug 2040 | 5.78 | 49° |
| 19 Oct 2041 | 13 Sep – 20 Nov | 320 d | Sep 2042 | 5.65 | 62° |
| 16 Nov 2043 | 9 Oct – 10 Dec | 304 d | Sep 2044 | 5.77 | 72° |
| 11 Dec 2045 | 3 Nov – 10 Jan | 290 d | Sep 2046 | 6.07 | 75° |

Windows come every 25–28 months and stay open about two months, on the real
calendar, because the hole does not move the planets measurably until the
mid-2040s. **The December 2045 window is the last one.** The hole is 11 AU out
when it opens, and that ship would arrive in September 2046 with the hole
inbound. This calculation counts only the Sun's pull on the spacecraft, so that
flight's cost is indicative. After it, no flight of a year or less is affordable
until 2052.

### After the encounter

The new orbits do not drift. From 2048 to 2200, Earth's year holds at
660–662 days and Mars's at 563–564; the two orbits are tilted 2.9° to each other.
Because the years are so close, the planets return to the same relative
position only every **10.4–10.6 years** (the synodic period, against 2.1 today).

Windows below are for flights of a year or less. A window is "open" while it
costs **6.0 km/s or less**, i.e. no more than a typical window today. Dates are
from a 2-day grid.

| best departure | open | flight | arrives | cost | Mars leads Earth by | years since previous |
|---|---|---:|---|---:|---:|---:|
| **14 Jul 2053** | Mar 2053 – Nov 2054 | 313 d | 23 May 2054 | 4.86 km/s | −18° | — |
| 7 Oct 2063 | Mar 2063 – Jul 2064 | 362 d | Oct 2064 | 4.94 | −30° | 10.2 |
| 27 Jul 2074 | May 2073 – Jan 2075 | 306 d | May 2075 | 4.73 | −11° | 10.8 |
| 27 May 2084 | Feb 2084 – Oct 2085 | 362 d | May 2085 | 4.78 | −32° | 9.8 |
| 9 Mar 2095 | Nov 2094 – Jul 2096 | 326 d | Jan 2096 | 4.84 | −21° | 10.8 |
| 5 Jun 2105 | Oct 2104 – Mar 2106 | 362 d | Jun 2106 | 4.95 | −30° | 10.2 |
| 25 Mar 2116 | Dec 2114 – Aug 2116 | 274 d | Dec 2116 | 4.72 | −9° | 10.8 |
| 16 Jan 2126 | Sep 2125 – Jun 2127 | 362 d | Jan 2127 | 4.77 | −30° | 9.8 |
| 28 Oct 2136 | Jul 2136 – Mar 2138 | 314 d | Sep 2137 | 4.87 | −19° | 10.8 |
| 24 Jan 2147 | Jun 2146 – Oct 2147 | 362 d | Jan 2148 | 4.89 | −27° | 10.2 |
| 15 Nov 2157 | Aug 2156 – Apr 2158 | 266 d | Aug 2158 | 4.72 | −7° | 10.8 |
| 10 Sep 2167 | May 2167 – Jan 2169 | 362 d | Sep 2168 | 4.77 | −31° | 9.8 |
| 24 Jun 2178 | Mar 2178 – Nov 2179 | 322 d | May 2179 | 4.84 | −20° | 10.8 |
| 25 Sep 2188 | Feb 2188 – Jul 2189 | 362 d | Sep 2189 | 4.96 | −31° | 10.3 |
| 14 Jul 2199 | Apr 2198 – Dec 2199 | 314 d | May 2200 | 4.74 | −11° | 10.8 |

A 362-day flight means the one-year cap is binding: those cycles would choose a
slightly longer, slightly cheaper flight if allowed. Some cycles also have a short shoulder window a few months either side (for
example Sep 2052 – Jan 2053, at best 5.66 km/s).

What changes:

- **Windows come once a decade.** The gaps run 10.2, 10.8, 9.8, 10.8 years and
  then repeat, a 41.6-year pattern. The very cheapest windows recur every fourth
  cycle: 2074, 2116, 2157, 2199.
- **Each window stays open much longer:** 16–20 months at today's cost, against
  about two months now.
- **Each trip is cheaper, by about 1 km/s.** The orbits are now neighbours, so a
  ship leaves Earth at 0.8–2.1 km/s and meets Mars at 0.4–1.2 km/s, against about
  3 km/s at each end today. The cost is 4.7–5.0 km/s.
- **Mars must be behind Earth at launch, by 7–32°**, as seen from the Sun. It is
  now the inner and faster planet, so it starts behind and catches up during the
  flight. Today it must lead by 29–75°. The near-Hohmann transfers, at 266–326 days,
  depart with Mars 7–21° behind.
- **Flights take 9–12 months**, about the same as today.

**Between windows, the door is closed only to short flights.** With the flight
capped at a year, off-window departures cost up to 29 km/s. If a flight of
1½–2½ years is acceptable, **46% of all departure dates cost 6.0 km/s or less,
and the worst date of the century costs 9.0 km/s.** A desperate mission can go
off-cycle by flying longer, without needing much more fuel.

**Right after the encounter**, from 2048 to mid-2052, every flight of a year or
less costs 6.0–25 km/s. A slower flight could leave sooner, from September 2051 at
5.09 km/s on a 538-day transfer.

**July 2053 in detail.** The best day is **14 July**: 313 days, arriving
**23 May 2054**, 4.86 km/s in total. The ship leaves Earth's pull at 1.73 km/s and
meets Mars at 1.05 km/s (a launch energy, C3, of 3.0 km²/s², against 9–17 for
today's windows). Little rides on the exact date: every day from 1 to 31 July
costs 4.861–4.883 km/s, so a launch slip of weeks barely matters.

**The planets also come much closer.** Earth and Mars pass nearest once a cycle,
at 0.11–0.25 AU. The closest is 0.107 AU (16 million km) in 2085; today's best
opposition is about 0.37 AU.

---

## 7. How typical is this?

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

## 8. Two ways the parameters lie

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

## 9. What is and is not modelled

**Solid.** All dynamics: REBOUND with IAS15, every date and distance in §1–3,
and the post-flyby elements, which reproduce present-day values when the same
extraction is applied at t = 0.

**Modelled, with stated limits.** The climate outcomes. Both are energy-balance
models on a latitude grid: no clouds, no ocean circulation, no weather. The
Earth model has no ice-sheet dynamics. The Mars model has no water cycle at all
— no reservoir, no evaporative cooling, no transport — and reports daily means,
while Mars's real diurnal range is 60–100 K.

**The Moon's phase was wrong until 2026-11-07.** The engine's Earth and Moon
masses were together 0.017% short of the ephemeris they start from, so the
simulated month ran about 13 minutes long and the Moon fell steadily behind: by
2047 its phases were 22 days late. Earth's orbit about the Sun was unaffected,
but anything that depends on where the Moon is around Earth was not — the Moon's
closest approach to the hole in §1 among them. The masses now match de440; the
re-run's lunar phases stay within an hour of the real Moon's through 2047, and
the figures here come from it.

**The launch windows (§6) are mission-design estimates, not integrated
trajectories.** The planets' positions come from the run, but each transfer arc
feels only the Sun, which is accurate once the hole is far away (14 AU by late
2049). The costs assume 300 km parking orbits at both ends and no aerobraking,
and ignore launch-site, plane-change and planetary-departure details. They are
right for comparing one date with another and the new era with today; they are
not a flight plan.

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
