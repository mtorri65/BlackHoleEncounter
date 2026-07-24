# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project overview

A flat collection of ~26 standalone Python scripts simulating a close flyby of a passing
black hole through the solar system (REBOUND N-body integration), plus post-processing,
visualization, and a separate analysis thread that checks whether the BH's astrometric
signature would be detectable in real Gaia DR3 data. There is no package structure — every
script is run directly from the repo root.

`Script_Summaries.pdf` in the repo root has a full purpose/inputs/outputs/techniques
write-up for every script — consult it before re-deriving what a script does from scratch.

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
  and the core engine when `write_mp4: true`) require `ffmpeg` on PATH.

## Architecture

### Core engine drives everything

`solar_system_bh_rebound26.py` is the center of the project. Nearly every other
script either consumes its output files or duplicates a fragment of its logic. It:

- Reads a YAML config (`input.yaml`) describing the BH's hyperbolic orbit — periapsis,
  v-infinity, inclination, longitude of ascending node, argument of periapsis, time-of-
  periapsis offset — each optionally given as a `min,max,step` sweep range — plus
  integration and output-writing settings.
- Initializes the Sun and planets from a Skyfield JPL ephemeris (de440s) at a UTC epoch;
  the Moon is added as a simplified fixed offset from Earth, not a real lunar ephemeris.
- Integrates with REBOUND's IAS15 adaptive integrator.
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
- `__orbits__<rp>__<mass>.xlsx` — per-body position/velocity/tidal-accel time series
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

### Three tiers of scripts

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
   `replay_archive.py`, `count_planets_run_rows.py`).
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

`time_dilation_tides.py` is a standalone physics illustration (Schwarzschild tidal/
time-dilation effects near the BH) with no dependency on the rest of the pipeline.

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
