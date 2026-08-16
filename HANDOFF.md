# Handoff — state as of 2026-08-15

Written at the end of a long working session, for whoever (or whatever) picks
this up next. Everything below was checked against the repository at the time of
writing rather than recalled.

---

## Read this first: verify, don't remember

The single most expensive mistake of the last session was **writing facts into a
document from memory of earlier conversation instead of checking the repo.** One
table in `CLIMATE_MODEL_REPORT.md` claimed the retired sweep used
`v_inf` = 10 km/s and had been deleted. Both were false — it used 25 km/s and was
still on disk — and those two errors then seeded two more: a wrong causal
explanation for a result, and a wrong replacement for it after the first was
caught. Four corrections, one root cause, and a single `grep` of an `input.yaml`
would have prevented all of them.

Concretely, before asserting any of the following, go and look:

- what a sweep's parameters were → `simulations/<STAMP>/<run>/*__input.yaml`
- what is on disk → `ls simulations/`, not what a document says
- what a number is → recompute it from the derived CSVs
- what a script does → its docstring, which is current; **not**
  `Script_Summaries.pdf`, which is not

---

## Current state

**Git:** `main` at `c9a2728`, in sync with `origin/main`
(github.com/mtorri65/BlackHoleEncounter). Nothing uncommitted except
`impact_2047.png`, which was modified before the session and was never touched.

**The scenario.** A 0.1 M☉ black hole, `v_inf` = 25 km/s, epoch 1885-09-01,
perihelion **0.662 AU on 2047-08-20**. The adopted run is

```
simulations/20260811_184731/20260811_184731__rp0p75__vinf25__inc30__toff59132__Om0__om30
```

**Sweep data.** `simulations/20260811_184731/` — 4032 runs, 16.9 GB, current and
the basis for everything. Its derived products are 6 CSVs in `simulations/`
(`_climate`, `_mars_climate`, `_earth_elements`, `_mars_elements`,
`_impact_ranking`, `_bh_captures`). **`simulations/` is gitignored** — none of
this is in version control.

**Retired sweep.** `20260724_230314` (epoch 1873, perihelion 2027, 672 runs) had
its raw runs deleted on 2026-08-15, reclaiming 172.7 GB. **Its 10 derived CSVs
were deliberately kept** — they are the only surviving record of the two-surface,
Sellers and transient climate analyses, and they make the epoch-vs-grid
comparison in §5.8 of the climate report possible without re-running anything.
Do not delete them.

**Belt run.** `simulations/20260815_171552/` — one run, adopted-scenario
parameters plus a 4000-tracer belt, produced with all engine fixes in place. Its
results are **not yet written into any document** (see below).

---

## What is settled

- **The engine's belt path had four bugs**, all fixed and committed in `06e3630`,
  none of which affected any published result because `n_belt` was 0. Three were
  in the hand-rolled element→state conversion (velocity perpendicular to radius
  so every tracer started at an apsis; inclination applied without shrinking the
  in-plane components; mean anomaly used as true anomaly), the fourth in
  reporting (`mu = G(Msun + M_bh)` where heliocentric elements need `G*Msun`).
  `sim.N_active` was also never set, costing 52× in speed.
- **The sweep is exactly twofold redundant.** `(i, Ω, ω)` and
  `(−i, Ω+180°, ω+180°)` are the same physical orbit, so 4032 runs cover ~2016
  configurations. Every count has an effective sample size of ~2016. Documented
  in `SCENARIO_2047_assumptions.md` §7a, which also uses the mirror pairs as a
  free convergence test.
- **The median climate outcome is epoch-dependent, not robust.** It moved from
  −9.2 K to +5.1 K between sweeps; §5.8 of the climate report separates this into
  +12.64 K from the 12-year epoch shift and +1.65 K from the ω grid refinement.
  Quote the distribution, not its centre.
- **The Earth–Moon system survives** the flyby, verified on relative state
  vectors at all 7,374 output steps. Note that `bh_captures.csv` reports
  *heliocentric* elements, so any satellite appears to be on a distinct orbit
  from its primary — that file structurally cannot answer whether a moon is
  retained.

---

## What is open

1. **The belt run's results are undocumented.** `simulations/20260815_171552/`
   has correct numbers — 0 planet-crossers before, 1334 Mars-crossing and 92
   Earth-crossing after, 0.90% unbound, mean impact waiting time 5,458 years —
   and no document mentions them. A short section in
   `SCENARIO_post_flyby_system.md` would be the natural home.
2. **`SCENARIO_timeline.md` Act I** predates the §6 rate-figure regeneration in
   `SCENARIO_2047_assumptions.md`. Worth a consistency pass; nothing is known to
   be wrong.
3. **`CLIMATE_MODEL_REPORT.md` §6.5** records four removed sections (transient,
   two-surface, Milankovitch, Sellers) awaiting regeneration against the current
   sweep. Needs three extra passes of `climate_from_simulations.py` with
   non-default options. The §6.5→§6.9 numbering gap is deliberate.
4. **`Script_Summaries.pdf`** is dated 2026-07-23 and 16 of 42 scripts have
   changed or appeared since. No markdown source exists.
5. **Book assets are untracked**: `bh_sky_track_full.png`,
   `bh_encounter_storyboard.png`, `bh_sky_track_2047.png`,
   `bh_sky_track_1885_1900.png`, `bh_observing_window.png`, and three MP4s.
   Decide whether they belong in the repo.
6. **`sky_stars.csv` / `sky_constellation_lines.csv`** are gitignored by the
   blanket `*.csv` rule, so a fresh clone has no star data and the sky charts
   silently fall back to plain labels. Re-run `fetch_constellation_data.py`.

---

## Context

The user is **writing a science-fiction novel** set in this scenario, and intends
to feed documents to a separate chat for that purpose. `SCENARIO_story_bible.md`
was written for that: constraints, dated beats, what an observer sees, and an
explicit boundary between what the simulation establishes and what is free
invention. The recommended set to share is the bible plus
`SCENARIO_timeline.md`, plus `bh_sky_track_full.png` and
`bh_encounter_storyboard.png`.

A **brown-dwarf variant** (0.075 M☉ instead of a black hole) was explored and
**rejected** — the user asked that it not be raised again. It fails because such
an object is luminous: it would have been a naked-eye star since roughly 1252 AD,
destroying the detection narrative and mooting the Gaia thread.

Two documents were deleted this session as superseded:
`SCENARIO_mercury_capture.md` (capture is no longer distinctive — 32 runs capture
a body) and `SCENARIO_timeline_old.md` / `SCENARIO_detection_timeline.md` (both
described the retired 10 km/s configuration). All survive in git history.
