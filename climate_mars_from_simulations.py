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

All temperatures reported here are **daily means**. Mars's real diurnal range is
60-100 K -- far larger than Earth's, because its thermal skin depth is
centimetres rather than an ocean mixed layer -- so daily maxima can exceed these
values substantially. The diurnal cycle is outside this model's scope, which
matters most for any threshold diagnostic sitting near a phase boundary.

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
from orbital_climate.mars import (MarsEBM, above_water_triple_point,
                                  liquid_water_possible)

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
    # Liquid-water fraction at every latitude, so the best band can be reported
    # alongside the equatorial value. rec["T"] is (time, lat); pressure is
    # global, so it broadcasts down the latitude axis.
    liquid_best = liquid_water_possible(rec["T"] + KELVIN, p[:, None]).mean(axis=0)
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
        # Fraction of the year the equator is above 273 K. On its own this is
        # NOT a statement about liquid water: Mars sits essentially on water's
        # triple point (611.657 Pa), so below that pressure ice sublimates
        # directly to vapour however warm the surface gets. Kept as a bare
        # temperature diagnostic; use the columns below for habitability.
        "frac_year_equator_above_273K": float(
            np.mean(rec["T"][:, i_eq] + KELVIN > 273.15)),
        # Three habitability diagnostics of increasing strictness, all reported
        # so the gap between them is visible rather than a matter of trust.
        #
        # (a) Above the triple point in both T and p. This is what most
        #     habitability summaries quote, and it is too generous: it counts
        #     surfaces hot enough to boil.
        "frac_year_above_triple_point": float(np.mean(
            above_water_triple_point(rec["T"][:, i_eq] + KELVIN, rec["pressure_Pa"]))),
        # (b) The same, but also requiring the ambient pressure to exceed
        #     water's saturation vapour pressure -- i.e. liquid that does not
        #     immediately boil. At ~700 Pa that window is only ~2 K wide, so
        #     this is typically several times smaller than (a).
        "frac_year_liquid_water_possible": float(np.mean(
            liquid_water_possible(rec["T"][:, i_eq] + KELVIN, rec["pressure_Pa"]))),
        # (c) The same again, but at whichever latitude does best rather than at
        #     the equator. On a tilted, eccentric orbit the warm band migrates
        #     off the equator -- in the strongest run of the 2047 sweep it sits
        #     near 35 deg S -- so an equator-only diagnostic understates the
        #     planet's best case, sometimes by a factor of several.
        "frac_year_liquid_water_best_lat": float(liquid_best.max()),
        "best_water_lat_deg": float(m.lat_deg[int(np.argmax(liquid_best))]),
        "spinup_years": int(info["years"]),
        "S_mean_Wm2": cfg.S0 / (4.0 * cfg.a_au ** 2 * np.sqrt(1.0 - cfg.ecc ** 2)),
        # Atmospheric collapse: the fraction of the year spent with the CO2
        # essentially all condensed. Above zero, the reported temperatures and
        # pressures are outside the model's domain for part of the year -- the
        # frost point degenerates as p -> 0, so those values describe nothing.
        "collapsed_fraction": float(rec["collapsed_fraction"]),
        "min_airborne_fraction": float(
            1.0 - rec["m_frost"].mean(axis=1).max() / cfg.co2_inventory_kg_m2),
    })
    out["atmosphere_collapsed"] = bool(out["collapsed_fraction"] > 0.0)
    out["permanently_collapsed"] = bool(out["collapsed_fraction"] > 0.99)
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
            f"Run:  python extract_mars_elements.py {args.parent_dir} --workers 5")

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
        # Four nested tests, each stricter than the last. Printed together
        # because the drop between them is the interesting part: a habitability
        # claim built on the first line is off by a large factor.
        warm = int((ok["frac_year_equator_above_273K"] > 0).sum())
        trip = int((ok["frac_year_above_triple_point"] > 0).sum())
        wet = int((ok["frac_year_liquid_water_possible"] > 0).sum())
        best = int((ok["frac_year_liquid_water_best_lat"] > 0).sum())
        print(f"\n  equator above 273 K at some point            : {warm} "
              f"({100 * warm / len(ok):.0f}%)")
        print(f"  ...AND above the triple-point pressure      : {trip} "
              f"({100 * trip / len(ok):.0f}%)")
        print(f"  ...AND not boiling (true liquid, equator)   : {wet} "
              f"({100 * wet / len(ok):.0f}%)")
        print(f"  ...same, at the best latitude rather than 0 : {best} "
              f"({100 * best / len(ok):.0f}%)")
        if warm > trip:
            print(f"    -> {warm - trip} are above freezing below the triple-point "
                  "pressure: ice")
            print("       sublimates straight to vapour, no liquid at any temperature.")
        if trip > wet:
            print(f"    -> {trip - wet} more are above the triple point but hot enough "
                  "to boil at")
            print("       ~700 Pa, where water boils at 275 K. The liquid window is ~2 K.")
        if best:
            q = ok["frac_year_liquid_water_best_lat"]
            print(f"    best-latitude liquid fraction: median {q[q > 0].median():.3f} "
                  f"of the year, max {q.max():.3f}")
        if "atmosphere_collapsed" in ok.columns:
            coll = int(ok["atmosphere_collapsed"].sum())
            perm = int(ok["permanently_collapsed"].sum())
            print(f"\n  ATMOSPHERE FREEZES OUT (model outside its domain):")
            print(f"    at some point in the year : {coll} ({100 * coll / len(ok):.0f}%)")
            print(f"    for the whole year        : {perm} ({100 * perm / len(ok):.0f}%)")
            if coll:
                print("    -> temperatures and pressures for these runs are not "
                      "meaningful; the honest result is 'the atmosphere freezes out'.")
    n_err = int(df["error"].notna().sum()) if "error" in df.columns else 0
    if n_err:
        print(f"  errors: {n_err}")
    print(f"\nWrote {out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
