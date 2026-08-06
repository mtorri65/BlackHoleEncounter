"""Run the Mars climate model on every simulation in a BH-flyby sweep.

Mars counterpart to ``climate_from_simulations.py``. Takes the orbital elements
recovered by ``extract_mars_elements.py`` and, for each run, spins the
condensing-CO2 EBM (``orbital_climate.mars.MarsEBM``) up to its repeating
seasonal cycle.

Differences from the Earth driver, all physical rather than structural:

* **Grey-body OLR.** Mars's greenhouse is ~5 K, so ``sigma T^4`` is the physics
  rather than a fit, and the 230-300 K validity window that constrains the
  Earth model does not apply. The validity limits here are the CO2 cycle's,
  not the radiation's.
* **No land/ocean split.** Mars has no ocean; a single regolith surface with
  tau ~ 6.6 days is correct.
* **CO2 condensation.** The headline diagnostics are the seasonal pressure
  swing and cap mass, which have no Earth analogue.

Usage
-----
    python climate_mars_from_simulations.py simulations/<STAMP> --workers 5
"""

from __future__ import annotations

import argparse
import dataclasses
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path

import numpy as np
import pandas as pd

from orbital_climate.config import Config, load_config
from orbital_climate.mars import MarsEBM

STEPS_PER_YEAR = 180
KELVIN = 273.15


def _evaluate(base_dict: dict, row: dict) -> dict:
    """Spin one run's Mars climate up to equilibrium. Runs in a worker process."""
    base = Config(**base_dict)
    yr = float(row["days_per_year_after"])
    cfg = dataclasses.replace(
        base,
        a_au=float(row["a_au_after"]),
        ecc=float(row["ecc_after"]),
        obliquity_deg=float(row["obliquity_deg_after"]),
        lon_perihelion_deg=float(row["lon_perihelion_deg_after"]),
        days_per_year=yr,
        dt_days=yr / STEPS_PER_YEAR,
    )
    out = {"run": row["run"], "a_au": cfg.a_au, "ecc": cfg.ecc,
           "obliquity_deg": cfg.obliquity_deg,
           "lon_perihelion_deg": cfg.lon_perihelion_deg, "days_per_year": yr}
    try:
        m = MarsEBM(cfg)
        T, mf, M, info = m.spin_up_co2()
        T, mf, M, rec = m.run_year_co2(T, mf, M0=M)
    except Exception as exc:                       # noqa: BLE001 - record, don't abort
        out["error"] = f"{type(exc).__name__}: {exc}"
        return out

    Ta = rec["T"].mean(axis=0)
    i_eq = int(np.argmin(np.abs(m.lat_deg)))
    p = rec["pressure_Pa"]
    p_max = float(p.max())
    out.update({
        "T_global_K": float(Ta.mean()) + KELVIN,
        "T_equator_K": float(Ta[i_eq]) + KELVIN,
        "T_pole_K": float(Ta[-1]) + KELVIN,
        "T_max_anywhere_K": float(rec["T"].max()) + KELVIN,
        "T_min_anywhere_K": float(rec["T"].min()) + KELVIN,
        "equator_seasonal_range_K": float(rec["T"][:, i_eq].max() - rec["T"][:, i_eq].min()),
        "pressure_mean_Pa": float(p.mean()),
        "pressure_min_Pa": float(p.min()),
        "pressure_max_Pa": float(p.max()),
        # Fraction of the atmosphere cycling in and out of the caps each year.
        "pressure_swing_frac": float(1.0 - p.min() / p_max) if p_max > 0 else np.nan,
        "peak_cap_kg_m2": float(rec["m_frost"].max()),
        # Fraction of the year the equator is above the melting point of water.
        "frac_year_equator_above_273K": float(
            np.mean(rec["T"][:, i_eq] + KELVIN > 273.15)),
        "spinup_years": int(info["years"]),
        "S_mean_Wm2": cfg.S0 / (4.0 * cfg.a_au ** 2 * np.sqrt(1.0 - cfg.ecc ** 2)),
    })
    # A permanently collapsed atmosphere is a qualitatively different world.
    out["atmosphere_collapsed"] = bool(out["pressure_mean_Pa"] < 0.05 * base.co2_inventory_kg_m2 * base.surface_gravity)
    return out


def main() -> int:
    p = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    p.add_argument("parent_dir", type=Path, help="Sweep folder, e.g. simulations/<STAMP>")
    p.add_argument("--elements", type=Path, default=None)
    p.add_argument("--config", type=Path, default=Path("input_mars.yaml"))
    p.add_argument("--out", type=Path, default=None)
    p.add_argument("--workers", type=int, default=1)
    p.add_argument("--max-a", type=float, default=5.0, help="Exclude a >= this [AU].")
    p.add_argument("--max-e", type=float, default=0.9, help="Exclude e >= this.")
    p.add_argument("--diffusion-d", type=float, default=0.002,
                   help="Meridional transport, calibrated for Mars's thin atmosphere.")
    p.add_argument("--limit", type=int, default=None)
    args = p.parse_args()

    elements = args.elements or args.parent_dir.parent / f"{args.parent_dir.name}_mars_elements.csv"
    if not elements.exists():
        raise SystemExit(
            f"Elements CSV not found: {elements}\n"
            f"Run:  python extract_mars_elements.py {args.parent_dir}_parquet --workers 5")

    el = pd.read_csv(elements)
    n_total = len(el)
    bound = el[el["bound_after"].astype(bool)]
    usable = bound[(bound["a_au_after"] < args.max_a)
                   & (bound["ecc_after"] < args.max_e)
                   & bound["lon_perihelion_deg_after"].notna()]
    stats = dict(n_total=n_total, n_unbound=n_total - len(bound),
                 n_excluded=len(bound) - len(usable), n_usable=len(usable))
    if args.limit:
        usable = usable.head(args.limit)

    base = dataclasses.replace(load_config(args.config), co2_cycle=True,
                               diffusion_D=args.diffusion_d)
    base_dict = base.to_dict()
    records = usable.to_dict("records")

    rows = []
    if args.workers <= 1:
        for i, r in enumerate(records, 1):
            rows.append(_evaluate(base_dict, r))
            if i % 25 == 0 or i == len(records):
                print(f"  {i}/{len(records)} runs", flush=True)
    else:
        with ProcessPoolExecutor(max_workers=args.workers) as ex:
            futs = [ex.submit(_evaluate, base_dict, r) for r in records]
            for i, f in enumerate(as_completed(futs), 1):
                rows.append(f.result())
                if i % 25 == 0 or i == len(records):
                    print(f"  {i}/{len(records)} runs", flush=True)

    df = pd.DataFrame(rows).sort_values("run").reset_index(drop=True)
    out = args.out or args.parent_dir.parent / f"{args.parent_dir.name}_mars_climate.csv"
    df.to_csv(out, index=False)

    print(f"\nRuns: {stats['n_total']} total | {stats['n_unbound']} Mars unbound | "
          f"{stats['n_excluded']} outside band | {stats['n_usable']} usable")
    ok = df[df["T_global_K"].notna()] if "T_global_K" in df.columns else df.iloc[:0]
    if len(ok):
        print(f"\nMars today (this model): 203.4 K, 22% pressure swing\n")
        for lab, col, unit in [("global mean T", "T_global_K", "K"),
                               ("equatorial seasonal range", "equator_seasonal_range_K", "K"),
                               ("pressure swing", "pressure_swing_frac", "frac"),
                               ("peak cap mass", "peak_cap_kg_m2", "kg/m2")]:
            s = ok[col].dropna()
            print(f"  {lab:28s} median {s.median():9.2f}  "
                  f"5th {s.quantile(.05):9.2f}  95th {s.quantile(.95):9.2f}  [{unit}]")
        warm = int((ok["frac_year_equator_above_273K"] > 0).sum())
        print(f"\n  runs where the equator exceeds 273 K at some point: {warm} "
              f"({100 * warm / len(ok):.0f}%)")
        if "atmosphere_collapsed" in ok.columns:
            print(f"  runs with a collapsed atmosphere: {int(ok['atmosphere_collapsed'].sum())}")
    n_err = int(df["error"].notna().sum()) if "error" in df.columns else 0
    if n_err:
        print(f"  errors: {n_err}")
    print(f"\nWrote {out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
