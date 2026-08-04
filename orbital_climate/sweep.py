"""Parallel sweep harness over EBM configurations.

Given a base :class:`~orbital_climate.config.Config` and a mapping of parameter
names to value lists, this runs the Cartesian product of overrides, spins each
configuration up to equilibrium, and collects scalar diagnostics into one table
(written as parquet if available, else CSV).

Each combo is fully independent, so combos are distributed across worker
processes with a :class:`concurrent.futures.ProcessPoolExecutor`, mirroring the
pattern used by the REBOUND driver in this repository.

Example
-------
    from orbital_climate.config import Config
    from orbital_climate.sweep import run_sweep
    run_sweep(Config(), {"ecc": [0.0167, 0.06, 0.117],
                         "diffusion_D": [0.4, 0.58, 0.8]},
              out_dir="climate_sweeps", workers=5)
"""

from __future__ import annotations

import dataclasses
import itertools
from datetime import datetime
from pathlib import Path
from concurrent.futures import ProcessPoolExecutor, as_completed

from .config import Config
from .experiment import run_equilibrium


def _combo_label(combo: dict) -> str:
    """Filesystem-safe label for one parameter combination."""
    parts = []
    for k, v in combo.items():
        vs = f"{v:g}" if isinstance(v, (int, float)) else str(v)
        vs = vs.replace("-", "m").replace(".", "p")
        parts.append(f"{k}{vs}")
    return "__".join(parts) if parts else "base"


def _sweep_worker(base_dict: dict, combo: dict) -> dict:
    """Run one equilibrium and return a flat row of parameters + diagnostics.

    Runs in a worker process; takes/returns only plain dicts so everything is
    trivially picklable.
    """
    cfg = dataclasses.replace(Config(**base_dict), **combo)
    res = run_equilibrium(cfg)
    row = dict(combo)
    row.update({
        "spinup_years": res.spinup_years,
        "global_mean_degC": res.global_mean,
        "global_mean_K": res.global_mean + 273.15,
        "iceline_lat_nh": (float("nan") if res.iceline_lat_nh is None
                           else res.iceline_lat_nh),
        "diag_lat_deg": res.diag_lat_deg,
        "peak_diag_insol": res.peak_diag_insol,
        "nh_mean_degC": res.nh_mean,
        "sh_mean_degC": res.sh_mean,
        "diag_summer_max_degC": res.diag_summer_max,
        "diag_winter_min_degC": res.diag_winter_min,
        "diag_seasonal_range_K": res.diag_seasonal_range,
        "absorbed_mean_Wm2": res.absorbed_mean,
        "runaway": res.runaway,
    })
    return row


def run_sweep(
    base: Config,
    ranges: dict[str, list],
    out_dir: str | Path = "climate_sweeps",
    workers: int = 1,
    stamp: str | None = None,
) -> Path:
    """Run the Cartesian product of ``ranges`` and write a diagnostics table.

    Parameters
    ----------
    base : Config
        Baseline configuration; each combo overrides a subset of its fields.
    ranges : dict[str, list]
        Mapping ``field_name -> list of values`` to sweep. Field names must be
        valid :class:`Config` fields.
    out_dir : path
        Parent directory; a timestamped subdirectory is created inside it.
    workers : int
        Number of worker processes. ``<= 1`` runs sequentially in-process.
    stamp : str, optional
        Override the timestamp label (mainly for tests / reproducibility).

    Returns
    -------
    Path
        The run directory containing ``sweep_results.parquet`` (or ``.csv``).
    """
    valid = {f.name for f in dataclasses.fields(Config)}
    bad = set(ranges) - valid
    if bad:
        raise ValueError(f"Unknown Config field(s) in sweep ranges: {sorted(bad)}")

    stamp = stamp or datetime.now().strftime("%Y%m%d_%H%M%S")
    run_dir = Path(out_dir) / stamp
    run_dir.mkdir(parents=True, exist_ok=True)

    keys = list(ranges.keys())
    combos = [dict(zip(keys, vals)) for vals in itertools.product(*(ranges[k] for k in keys))]
    base_dict = base.to_dict()

    rows: list[dict] = []
    if workers <= 1 or len(combos) <= 1:
        for combo in combos:
            rows.append(_sweep_worker(base_dict, combo))
    else:
        with ProcessPoolExecutor(max_workers=workers) as ex:
            futures = {ex.submit(_sweep_worker, base_dict, c): c for c in combos}
            for fut in as_completed(futures):
                rows.append(fut.result())

    _write_table(rows, keys, run_dir)
    return run_dir


def _write_table(rows: list[dict], sweep_keys: list[str], run_dir: Path) -> Path:
    """Write rows as parquet (preferred) or CSV; keep sweep columns leftmost."""
    try:
        import pandas as pd
    except Exception:  # pragma: no cover - pandas is a project dependency
        pd = None

    # Stable column order: swept parameters first, then diagnostics.
    diag_keys = [k for k in (rows[0].keys() if rows else []) if k not in sweep_keys]
    columns = sweep_keys + diag_keys

    if pd is not None:
        df = pd.DataFrame(rows, columns=columns)
        # Sort by the swept parameters for a tidy, reproducible table.
        df = df.sort_values(sweep_keys).reset_index(drop=True)
        try:
            path = run_dir / "sweep_results.parquet"
            df.to_parquet(path, index=False)
        except Exception:
            path = run_dir / "sweep_results.csv"
            df.to_csv(path, index=False)
        return path

    # Fallback: hand-rolled CSV.
    import csv
    path = run_dir / "sweep_results.csv"
    with open(path, "w", newline="", encoding="utf-8") as fh:
        writer = csv.DictWriter(fh, fieldnames=columns)
        writer.writeheader()
        for row in rows:
            writer.writerow(row)
    return path
