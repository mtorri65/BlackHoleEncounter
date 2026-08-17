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

1. ~~**The belt run's results are undocumented.**~~ **Done** — written up as §6
   of `SCENARIO_post_flyby_system.md`, with all figures recomputed from
   `simulations/20260815_171552/`'s before/after element files rather than
   transcribed.

   Two corrections came out of that recomputation, both worth carrying forward:

   - The crossing counts quoted in the earlier draft of this note (1334
     Mars-crossing, 92 Earth-crossing) use **today's** planetary orbits as the
     thresholds. Against the *post-flyby* orbits the figures are **1,247
     Mars-crossing and 1,617 Earth-crossing** — Earth moved out to 1.135–1.840 AU,
     straight into the perturbed belt, so the old threshold undercounts it ~18×.
   - `postprocess_belt_sizes_and_hazard.py` hardcodes `q_after < 1.0` as
     "Earth-crossing" (lines 41 and 170). The shipped 5,458-year waiting time is
     therefore computed against a planet that is no longer at 1 AU; recomputed
     against Earth's actual post-flyby orbit it is **206 years**. Worth
     parameterizing that threshold if the script is touched again — deliberately
     left alone here so the committed `__hazard_summary.txt` still matches what
     the script produces.
2. ~~**`SCENARIO_timeline.md` Act I** predates the §6 rate-figure regeneration.~~
   **Done.** The stated worry was unfounded — Act I never quotes §6's rate
   figures, and its astrometric anchor (0.76″ at 863 AU, 100-year baseline)
   already matches §4. But the pass found seven other things, all fixed:

   - Detection-table distances for **1962 (474→459 AU)** and **1971 (426→412 AU)**
     were ~3% high; every other distance in Act I checks out against the run.
   - The BH was called **"295 metres across"** twice — that is the Schwarzschild
     *radius*; the diameter is 591 m, as the story bible already said.
   - "the final 11 AU **in a single year**" — it takes twenty months.
   - The Moon was given as "402 000 km against today's 384 000", implying a
     widened orbit. 402 000 km is just where the Moon sits at the run's end; the
     semi-major axis is ~386 000 km, essentially unchanged.
   - "4.2× gain" was inconsistent with its own two factors (3.3 × 1.32 = 4.4).
   - 2044 distance 23 → 22 AU.
   - The four runs that freeze Earth, capture a planet and open a Mars window are
     **two mirror pairs = two configurations**, per §7a. Verified independently;
     the count of four is right, the independence is not.

   Also noted, not an error in either document: `__planets_run_deltas.csv` and
   `bh_captures.csv` disagree slightly for the same run (Saturn 1,650 vs 1,627 AU,
   amplified by e = 0.994). The timeline now carries a footnote explaining it.

   **Not regenerated:** the detection table's signal/precision/SNR columns are
   derived quantities, not simulation output. A 3% distance shift moves them ~10%,
   inside the quoted precision, so they were left alone — but they were not
   recomputed against the current configuration.
3. **`CLIMATE_MODEL_REPORT.md` §6.5** — **data regenerated 2026-08-15; the
   write-up is still outstanding.** The three passes have been run against the
   current sweep (2h19m total) and their CSVs are on disk:

   - `20260811_184731_climate_sellers.csv` — Sellers OLR, single surface
   - `20260811_184731_climate_2surf.csv` — two-surface, linear OLR
   - `20260811_184731_climate_2surf_transient.csv` — Sellers + two-surface, 40 yr

   Column schemas match the retired sweep's equivalents exactly, so the retired
   analyses can be reproduced directly. **§6.5 still reads "removed, pending
   regeneration" and must be rewritten into its four sections** — the numbers now
   exist, the prose does not. Headline findings, all on 3,904 usable runs:

   - **Sellers acts only on the cold end.** Median shift on already-snowball runs
     −15.97 K; on non-snowball runs +0.20 K. 182 runs flip into snowball, none
     out (non-runaway snowball 39.3% → 45.3%, n = 3,010).
   - **Two-surface is a seasonal-amplitude story, not a mean-temperature one.**
     Global means barely move; land seasonal range is 4.43× ocean (36.1 vs 8.2 K),
     and the diagnostic-latitude range doubles, 10.5 → 20.8 K. That doubling is
     what the Milankovitch section needs.
   - **Overshoot is ~zero, confirmed on the current sweep.** Median 0.0000 K, max
     0.0293 K, 8 runs above 0.01 K. §8 currently says "exactly zero" with retired
     provenance — tighten to "at most 0.03 K" and re-provenance it.
   - **`--years 40` is adequate**: median 14 years to equilibrium, 90th percentile
     23, only 1.0% hit the cap.

   §12's command block was corrected at the same time — see the `--out` warning
   there; the old block would have overwritten `_climate.csv`.

   **Write-up completed 2026-08-16.** §6.5–§6.8 are now four real sections again
   and the numbering gap before §6.9 is closed. Also updated: the report header,
   §3.1's forward reference, the derived-products list, and §8's overshoot and
   glacial-inception entries. Outcomes worth knowing:

   - **§6.6 and the §6.7 controlled experiment reproduce the retired sweep almost
     exactly** (land seasonal range 41.3 vs 41.2 K; controlled land cooling
     −7.51 vs −7.53 K). Both are model properties, so they are epoch-independent.
   - **§6.7's sweep-wide medians flipped sign** with the epoch, like §6.1. The
     ordering — single-surface over-predicts cooling relative to land — survives.
   - **One retired claim was withdrawn.** The retired §6.7 "sign flip" rested on
     an `a`-within-2% subset that does *not* hold annual-mean insolation fixed
     (+3.05 W/m² off baseline). Selecting on insolation instead shows no sign
     flip. Recorded in §8's wrong-claims table.
   - **Overshoot is not identically zero** — max 0.0293 K over 3,904 runs. §8
     reworded from "exactly zero" to "at most 0.03 K".

   **Both gaps closed 2026-08-16** with two further passes (43 min and 61 min),
   giving the full 2×2 of OLR form × surface treatment plus both transients. §6.5
   now carries the single-vs-two-surface comparison and §6.8 the glacial-inception
   result. Two further findings from those passes:

   - **OLR and surface treatment barely interact.** Snowball counts across the
     2×2: linear 1,182 / 1,174 and Sellers 1,364 / 1,346 (single / two surface).
     The OLR moves the count by ~+177, the surface split by ~−13. Varying one
     factor at a time therefore loses almost nothing.
   - **The retired glacial-inception comparison had a moving denominator.**
     Computing it within the 250–300 K band is unfair across OLR forms, because
     Sellers moves runs *out* of the band (1,135 → 911). On the band definition
     the population appears to fall, 28% → 22%; on a fixed non-runaway denominator
     it rises, 49% → 51%. §6.8 uses the fixed denominator and says why.

   Measured pass costs, for anyone re-running: 25 / 30 / 43 min (equilibrium:
   sellers, 2surf, sellers+2surf) and 61 / 83 min (transient: 1surf, 2surf).
   §12 now carries these. **Estimating from `--limit` samples underpredicts by
   1.4–2.2×** — the first runs in the elements file are cheaper than the sweep
   average, so extrapolate from full passes, not samples.

   ### Earth's frozen temperature was restated: 205 K → 182 K

   **Resolved 2026-08-16 — the documents now lead with the Sellers value.** For
   the adopted run:

   | | linear (Budyko) | **Sellers — now used** |
   |---|---:|---:|
   | Earth equilibrium | 204.70 K (−68.4 °C) | **182.09 K (−91.1 °C)** |

   The run is a snowball, squarely in the regime §10 of the climate report says
   the linear law gets wrong — its emission collapses toward zero, so a frozen
   planet stays artificially warm. Systematic, not a one-run quirk: median
   −15.97 K across all 1,182 snowball runs.

   The **full decline trajectory** was regenerated under both forms, not just the
   endpoint, so the tables are internally consistent (orbits 1/2/3/5/80):

   | OLR | orbit 1 | 2 | 3 | 5 | 80 |
   |---|---:|---:|---:|---:|---:|
   | linear | 261.4 | 224.1 | 210.9 | 205.3 | 204.7 K |
   | **Sellers** | 259.6 | 216.1 | 197.6 | 185.8 | **182.1 K** |

   The linear column reproduces the previously published table **exactly**, which
   is what validates the Sellers column beside it. NH ice edge at orbit 1 moves
   24.6° → 21.9°. **The freeze schedule is unchanged** — ice closes over the
   equator on the second orbit under both.

   Updated: `SCENARIO_timeline.md`, `SCENARIO_story_bible.md`,
   `SCENARIO_mars_window.md`, `SCENARIO_post_flyby_system.md`, and §6.8 of the
   climate report. Each keeps the linear figures as a documented comparison, since
   older notes quote them.

   **The 29% escape threshold was checked and survives — but it was a *Sellers*
   number all along.** No derivation for it existed anywhere: not in the climate
   report, not in the code, only the bare claim in two scenario documents. It was
   reconstructed by holding the post-flyby orbit fixed, scaling `S0`, and spinning
   the EBM up **from a frozen initial state** (`EBM.spin_up(T0=...)`), then
   comparing the resulting annual-mean insolation to present-day Earth's
   340.3 W/m². The transition is sharp — 247 → 321 K across one 0.4% step.

   | | freezes below | escapes above |
   |---|---:|---:|
   | linear | 0.910 × today | 1.171 × today (+17.1%) |
   | **Sellers** | **0.945 × today** | **1.284 × today (+28.4%)** |

   +28.4% matches the documented 29%; the linear value (+17.1%) does not. So the
   figure was already consistent with the Sellers configuration the documents now
   use, and **an earlier note in this handoff calling it "a linear-OLR number"
   was wrong** — both scenario documents have been corrected accordingly.

   Two things worth keeping from this:

   - **The model has a wide hysteresis loop** — freeze at 0.945×, thaw at 1.284×.
     That is *why* the snowball is permanent, and both documents now say so
     instead of quoting the escape threshold alone. Note this does **not**
     contradict §8's "initial-condition bistability does not occur": that test
     compared a uniform warm start against Earth's real profile, i.e. two points
     on the same warm branch, not a warm start against a frozen one.
   - **The method is undocumented in the report.** If the escape threshold matters
     again, write it into `CLIMATE_MODEL_REPORT.md` rather than re-deriving it a
     third time.

   *Note, pre-existing and untouched:* §8's subsections are numbered 7.1/7.2,
   duplicating Part VI's. Refer to that material as "§8" until it is renumbered.
4. **`Script_Summaries.pdf`** is dated 2026-07-23 and 16 of 42 scripts have
   changed or appeared since. No markdown source exists.
5. ~~**Book assets are untracked.**~~ **Done** — the five PNGs and three MP4s are
   committed, along with `.claude/settings.json` and **all 21 top-level
   `simulations/*.csv`**, which had been excluded by the `simulations/` and
   `*.csv` ignore rules and existed only on local disk.

   The retired sweep's 10 CSVs were the reason to do it: its raw runs were deleted
   on 2026-08-15 to reclaim 172.7 GB, so unlike the current sweep's products they
   **cannot be regenerated at any price**. The current sweep's 11 are regenerable
   from `simulations/20260811_184731/`, but only at ~4 hours of compute for the
   five climate passes alone.

   > **`.gitignore` was deliberately left unchanged; the files were force-added.**
   > Broadening the rule would sweep in the per-run `.csv` files under
   > `simulations/<stamp>/<run>/`, which are large and genuinely disposable.
   > **Any future derived CSV therefore needs `git add -f`** or it will be silently
   > skipped — this is the easiest way to lose new sweep products.

6. **`sky_stars.csv` / `sky_constellation_lines.csv` remain excluded, on purpose.**
   Unlike the derived sweep products they are a reproducible ~250 KB download:
   re-run `fetch_constellation_data.py`. The cost of a fresh clone missing them is
   that `sky_backdrop.py` falls back to plain constellation labels and the star
   field silently disappears, which reads like a plotting bug rather than missing
   data. Worth committing them too if that trap bites anyone again.

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
