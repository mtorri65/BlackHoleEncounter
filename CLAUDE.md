# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project overview

~42 standalone Python scripts in the repo root, simulating a close flyby of a passing
black hole through the solar system (REBOUND N-body integration), plus post-processing,
visualization, a climate thread that asks what the perturbed orbits do to Earth and Mars,
and a separate analysis thread that checks whether the BH's astrometric signature would be
detectable in real Gaia DR3 data. Every script is run directly from the repo root.

The one exception to "no packages" is `orbital_climate/`, an actual importable package
(`ebm.py`, `mars.py`, `insolation.py`, `kepler.py`, `sweep.py`, `experiment.py`,
`config.py`, `cli.py`) holding the energy-balance climate model. Root-level scripts import
from it; it is not run directly except through `orbital_climate/cli.py`. Its physics
background is in `orbital-climate-model-context.md`.

`Script_Summaries.pdf` in the repo root has a full purpose/inputs/outputs/techniques
write-up for every script — consult it before re-deriving what a script does from scratch.
**It is out of date**: it dates from 2026-07-23, and 16 of the 42 scripts have been added
or substantially changed since — the whole scenario-selection and climate thread, all the
sky-track/export tooling, and edits to the core engine itself. Treat a missing entry as
"not yet written up", not "no such script", and check a script's own docstring before
trusting the PDF's account of it.

## Environment / commands

- Python 3.13, no virtualenv or requirements file — dependencies are installed straight
  into the system Python. Key third-party packages: `rebound`, `astropy`, `astroquery`,
  `skyfield`, `pandas`, `numpy`, `matplotlib`, `openpyxl`, `PIL`, `PyYAML`.
- No test suite, no linter/formatter config, no build system. To sanity-check a script
  after editing, use `python -m py_compile <script>.py` (catches syntax errors without
  running the (often slow / data-dependent) script itself).
- Run the core simulation: `python solar_system_bh_rebound26.py input.yaml`
  (add `--resume-dir simulations/<timestamp>` to continue an interrupted parameter sweep).
- Most other scripts take a run folder, parent (batch) folder, or archive `.bin` file as a
  positional CLI arg — check each script's `argparse` block. A few have no CLI at all and
  instead hardcode input paths at the top of the file under a "CONFIGURE THIS"-style
  comment (`replay_archive.py`, `run_animation.py`, `convert_npz_to_csv.py`,
  `overlay_plot_v2.py`) — edit those constants before running.
- MP4-producing scripts (`animate_snapshots3D_1.py`, `make_inner_3au_videos_from_archive.py`,
  and the core engine when `write_mp4: true`) require `ffmpeg` on PATH. `animate_sky_track.py`
  is the exception — it prefers the binary bundled with the `imageio-ffmpeg` package and only
  falls back to PATH, so it works with no system-wide install.
- **One-time setup for the sky charts**: `python fetch_constellation_data.py` downloads the
  Hipparcos catalogue and Stellarium's constellation figures and distils them to
  `sky_stars.csv` + `sky_constellation_lines.csv` (~250 KB total; the 53 MB source is
  discarded). `.gitignore` excludes `*.csv`, **so these two files are never committed and a
  fresh clone will not have them** — re-run the fetch. Without them `sky_backdrop.py` prints
  a hint and falls back to plain constellation name labels, so charts still render; the star
  field just silently disappears, which is easy to misread as a plotting bug.

## Architecture

### Core engine drives everything

`solar_system_bh_rebound26.py` is the center of the project. Nearly every other
script either consumes its output files or duplicates a fragment of its logic. It:

- Reads a YAML config (`input.yaml`) describing the BH's hyperbolic orbit — periapsis,
  v-infinity, inclination, longitude of ascending node, argument of periapsis, time-of-
  periapsis offset — each optionally given as a `min,max,step` sweep range — plus
  integration and output-writing settings. `epoch` is t = 0 for every `t_days` in the
  output and the reference for `bh_tperi_offset_days`, so changing it changes the
  scenario, not just labels. The current config uses 1885-09-01; the retired 2027 sweep
  used 1873-09-01. Scripts that turn `t_days` back into dates read the epoch from each
  run's own `__input.yaml` rather than assuming one.
- Initializes the Sun, planets and Moon from a Skyfield JPL ephemeris (de440s) at a UTC
  epoch. The Moon used to be a synthetic fixed offset from Earth with an orbital velocity
  but no compensating change to Earth's — which injected ~12.4 m/s of spurious momentum
  and stretched Earth's year to 365.570 days. Runs made before that fix (all of
  `simulations/20260724_230314`) carry a seasonal drift of +0.33 d/yr; see the docstring
  of `plot_local_sky.py`.
- Integrates with REBOUND's IAS15 adaptive integrator.
- Logs the trajectory at `output_interval_days`, optionally at two rates: setting
  `output_dense_window_days` logs every `output_dense_interval_days` within that many days
  of `bh_tperi_offset_days` and at the coarse rate elsewhere. The dynamics are confined to
  a few years around periapsis, so 1-day/30-day cuts output ~20x with no loss.
- Sweeps the full Cartesian product of any BH parameter ranges (`itertools.product`),
  writing one output subfolder per parameter combination — wide ranges produce many runs.
- Can resume an interrupted sweep via `--resume-dir`.

"Tidal acceleration" throughout the whole pipeline is a simple point-mass approximation
(`2*G*M_BH/R^3 * r_helio`), not a full tidal tensor or a post-Newtonian/GR correction.
Perturbation from the BH is measured by diffing against a parallel BH-free baseline
integration, not by adding relativistic terms to the equations of motion.

### Run-folder output convention

Each run lands in
`simulations/<timestamp>/<timestamp>__rp<..>__vinf<..>__inc<..>__toff<..>__Om<..>__om<../`,
with every file prefixed by the subfolder name:

- `__input.yaml` — copy of the config used for that run
- `__archive.bin` — REBOUND SimulationArchive (only if `archive_enable: true`)
- `snapshot_<t>.npz` — per-body Cartesian state per dump (only if `write_npz: true`)
- `__orbits__<rp>__<mass>.parquet` — per-body position/velocity time series, long
  format (`body, t_days, x_au, y_au, z_au, vx, vy, vz, disp_helio_au`, float32).
  The default since the `orbits_format` option was added; the same schema
  `convert_orbits_to_parquet.py` produces from older sweeps (as `orbits.parquet`).
- `__orbits__<rp>__<mass>.xlsx` — the legacy form of the same log, with derived
  columns and `_str` duplicates. Only written for `orbits_format: xlsx` or `both`.
  Roughly 5x larger; a full 672-run sweep at 1-day cadence came to 150 GB.
- `__bh_radec__<rp>__<mass>.xlsx` — BH geocentric RA/Dec track, feeds the Gaia analysis scripts
- `__planets_run_deltas.csv`, `__belt_run_before/after.csv` — orbital elements before/after
- `uncertainties/uncertainty_input__<body>__*.csv` — feeds `estimate_BH_parameters_uncertainties*.py`
- plots/animations: `__heliocentric_plot.png`, `__heliocentric_inner_3AU.png`, belt
  scatter/histogram PNGs, `__hazard_summary.txt`, MP4s

Downstream scripts locate their inputs either via a run/parent-folder CLI arg, or (several
Gaia-analysis scripts) by auto-discovering the most recently modified matching file under
`simulations/` with `glob`/`rglob`. Because of this, most post-processing scripts are not
independently runnable — they implicitly depend on a prior core-engine run having already
produced the expected files in the expected place.

### Five tiers of scripts

1. **Core engine** — `solar_system_bh_rebound26.py`. Run this first.
2. **Post-processing / visualization** (consume core-engine output): orbit comparison
   (`post_orbit_compare_geocentric_fixed4.py`), BH-parameter-recovery uncertainty via a
   Fisher-matrix method (`estimate_BH_parameters_uncertainties.py` and the vectorized,
   preferred `_fast.py` sibling), belt impact-hazard estimation
   (`postprocess_belt_sizes_and_hazard.py`), Earth-temperature sweep summary
   (`summarize_earth_temperatures1.py`), an NPZ-based repair utility for older runs
   (`rebuild_planets_run_deltas_from_npz.py`), and archive/NPZ-based visualization/animation
   (`plot_bh_perihelion.py`, `make_inner_3au_videos_from_archive.py`,
   `animate_snapshots3D_1.py`, `run_animation.py`, `convert_npz_to_csv.py`, the two
   `gallery_heliocentric_*` interactive viewers, `inspect_archive_hashes.py`,
   `replay_archive.py`, `reconstruct_trajectories_from_archive.py`,
   `count_planets_run_rows.py`).
3. **Gaia DR3 detectability analysis** (independent thread — cross-references the core
   engine's BH RA/Dec track against real Gaia astrometry to assess whether the BH's
   point-lens microlensing signature would be observable under a synthetic Vera Rubin/LSST
   observing cadence): `GAIA_dr3_v1.py` is the entry point (live Gaia TAP query; produces
   the pass-catalog CSVs everything else in this tier reads), then
   `astrometric_shift.py` / `astrometric_shift_option1.py` (per-star shift plots — time
   series vs. sky-plane track), `detectable_shift_on_date.py`,
   `distance_from_bh_on_first_detect.py`, `overlay_bh_on_vr_image.py`, `overlay_plot_v2.py`.
   `centroid_motion.py` is not independently runnable — it's a code fragment that assumes
   variables already defined by the `astrometric_shift*.py` scripts.
4. **Scenario selection & climate** (the "which run is interesting, and what would living
   there be like" thread; runs in this order): `rank_run_impact.py` ranks a whole sweep by
   how much each run disturbed the system; `find_bh_captures.py` separates planets that
   merely leave from those that leave *with* the BH, which the per-run deltas cannot
   distinguish on their own; `extract_earth_elements.py` / `extract_mars_elements.py`
   recover full post-flyby elements, since `__planets_run_deltas.csv` carries only a, e, q;
   then `climate_from_simulations.py` / `climate_mars_from_simulations.py` drive the
   `orbital_climate/` energy-balance model over every run. Mars's habitability verdict comes
   from `orbital_climate/mars.py::liquid_water_possible`, which requires the surface to be
   above the triple point **and** below boiling at the local pressure — a weaker
   triple-point-only test is kept alongside as `above_water_triple_point` for comparison,
   and the two disagree by a large factor, so do not substitute one for the other.
5. **Sky-track charts** (where the BH appears on the sky, for the narrative documents):
   `animate_sky_track.py` is the entry point, producing MP4/GIF animations paced either by
   constant apparent speed (`--pace arc`, the default) or linearly in time; `--from-year` /
   `--to-year` restrict the window. `sky_backdrop.py` draws the star field and constellation
   figures behind it, and `fetch_constellation_data.py` builds the data those need (see
   Environment above — this is a required one-time step). `export_bh_track.py` and
   `export_bh_track_cdc.py` emit the track for external planetarium software instead, as
   generic CSV and Cartes du Ciel user-object format. `plot_sky_tracks.py` is older and
   unrelated in output: static RA/Dec tracks of the BH *and the planets* during the flyby.

### Sky-chart conventions

Three things about the sky charts are easy to get wrong silently:

- **Frame.** The simulation runs in equatorial J2000, and Hipparcos positions are ICRS,
  so stars and BH track share a frame directly and **no precession is applied**. Charts
  drawn for a historical epoch in "apparent" coordinates would need it — J2000 and
  equinox-of-date differ by ~1.6° at 1885, far above plotting tolerance. `export_bh_track.py`
  therefore writes *both*, and the chart's own setting decides which column to use.
- **Orientation.** `animate_sky_track.py` sets it with `RA_LEFT_EDGE` / `RA_RIGHT`. The
  default is the star-atlas convention — RA increasing to the *left*, axis reading 24ʰ→0ʰ —
  so the chart can be compared against a printed atlas. Flipping `RA_RIGHT = True` reads
  more naturally as a graph but is a mirror image of the sky, which is not obvious by
  inspection unless you check an asymmetric asterism.
- **Seam.** Independent of orientation: the RA wrap point is placed where the track does not
  cross it. The current scenario spans RA 5.05ʰ–23.17ʰ, so seaming at 0ʰ draws it unbroken;
  an 8ʰ seam cut it into two disconnected pieces. Constellation figures that straddle the
  seam are drawn twice, shifted, and name positions use a mean of unit vectors so they do
  not land halfway across the sky.

`time_dilation_tides.py` is a standalone physics illustration (Schwarzschild tidal/
time-dilation effects near the BH) with no dependency on the rest of the pipeline.

`astrometric_detectability.py` is standalone too — it runs its own paired N-body
integrations rather than reading sweep output. It asks whether 19th-century astrometry
could have seen the BH approaching, and inverts that into the `bh_vinf_kms` a scenario
needs. Its BH-orbit construction deliberately mirrors the engine's, sign conventions
included; the derived assumptions live in `SCENARIO_2047_assumptions.md`.

Note when reading any single run's geometry: the engine builds the BH state with a
negative `r` and then negates the velocity (`# inbound branch`), so the orbit is the
point inversion of what its labels describe — same plane, same `i`, `Ω`, `rp`, `e` and
arrival time, but periapsis on the opposite side (true ω = `bh_omega_deg` + 180° for
labelled `i ≥ 0`). Harmless for sweeps, since {0, 90, 180, 270} is closed under +180°,
and harmless for relative claims between runs; it matters only for absolute geometric
statements about one named run. See §7 of `SCENARIO_2047_assumptions.md`.

### Particle lookup pattern

Archive- and snapshot-based scripts use a hash-first, index-fallback pattern to look up
bodies (`get_particle`-style helpers), so they work against both newer hash-named REBOUND
archives and older anonymous ones. Follow this pattern rather than assuming particle hashes
always resolve.

### Versioning

The core engine is one of a versioned series (`solar_system_bh_reboundNN[_patched].py`,
seen up to v26) — `input.yaml` notes in a header comment which engine version it's
compatible with. Check that comment before assuming a config field exists across engine
versions.
