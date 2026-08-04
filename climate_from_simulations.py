"""Run the seasonal climate model on every BH-flyby simulation in a sweep.

Bridges the two halves of this repository:

    REBOUND sweep  ->  Earth's post-flyby orbital elements  ->  seasonal EBM

``extract_earth_elements.py`` recovers Earth's full post-flyby orbit from each
run's ``*__orbits__*.xlsx`` (a, e, obliquity and longitude of perihelion -- the
latter two are *not* in ``planets_run_deltas.csv``). This script turns each of
those into an ``orbital_climate.Config``, spins the EBM up to its repeating
seasonal cycle, and collects the resulting climate per run.

Why a script and not a ``sweep:`` block
---------------------------------------
The EBM's own sweep harness takes the Cartesian product of parameter lists. Here
the parameters are *paired* -- run N's (a, e, obliquity, lambda_p) belong
together -- so the runs must be enumerated explicitly.

Two corrections applied per run
-------------------------------
* **Year length.** Kepler's third law: P = a^1.5 years. Seasonal damping depends
  on the year length relative to the ocean mixed-layer time constant
  (tau = C/B ~ 1.6 yr), so this materially changes the answer.
* **Timestep.** ``dt_days`` is set to P / ``steps_per_year`` rather than left
  fixed, so every run resolves its seasonal cycle with the same fidelity
  regardless of whether its year is 184 or 1711 days long.

Validity band
-------------
The EBM linearises outgoing longwave radiation as A + B*T, calibrated near
288 K. Runs whose post-flyby orbit falls far outside that regime (very large a,
near-parabolic e) would still return a number, but a fictitious one, so they are
excluded by default (``--max-a`` / ``--max-e``) and reported separately.

Usage
-----
    python climate_from_simulations.py simulations/20260724_230314 --workers 5
    python climate_from_simulations.py simulations/20260724_230314 \
        --config input_climate.yaml --plot climate.png
"""

from __future__ import annotations

import argparse
import dataclasses
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path

import numpy as np
import pandas as pd

from orbital_climate.config import Config, load_config
from orbital_climate.experiment import run_equilibrium, run_perturbation

DEFAULT_STEPS_PER_YEAR = 180


# ---------------------------------------------------------------------------
# Per-run climate evaluation
# ---------------------------------------------------------------------------
def _config_from_elements(base: Config, a_au: float, ecc: float, obliquity_deg: float,
                          lon_perihelion_deg: float, days_per_year: float,
                          steps_per_year: int) -> Config:
    """Build a Config from a set of orbital elements.

    ``dt_days`` is derived from the year length so every run resolves its
    seasonal cycle with the same number of steps, regardless of orbital period.
    """
    return dataclasses.replace(
        base,
        a_au=a_au,
        ecc=ecc,
        obliquity_deg=obliquity_deg,
        lon_perihelion_deg=lon_perihelion_deg,
        days_per_year=days_per_year,
        dt_days=days_per_year / steps_per_year,
    )


def config_for_run(base: Config, row: dict, steps_per_year: int, when: str = "after") -> Config:
    """Build a Config for one run from its recovered orbital elements.

    ``when`` selects the pre-flyby (``"before"``) or post-flyby (``"after"``)
    element set recorded by ``extract_earth_elements.py``.
    """
    a = float(row[f"a_au_{when}"])
    days_per_year = (float(row["days_per_year_after"]) if when == "after"
                     else 365.256363 * a ** 1.5)
    return _config_from_elements(
        base, a, float(row[f"ecc_{when}"]), float(row[f"obliquity_deg_{when}"]),
        float(row[f"lon_perihelion_deg_{when}"]), days_per_year, steps_per_year)


def _evaluate(base_dict: dict, row: dict, steps_per_year: int) -> dict:
    """Spin one run's climate up to equilibrium. Runs in a worker process."""
    base = Config(**base_dict)
    cfg = config_for_run(base, row, steps_per_year)
    out = {
        "run": row["run"],
        "a_au": cfg.a_au,
        "ecc": cfg.ecc,
        "obliquity_deg": cfg.obliquity_deg,
        "lon_perihelion_deg": cfg.lon_perihelion_deg,
        "days_per_year": cfg.days_per_year,
    }
    try:
        res = run_equilibrium(cfg)
    except Exception as exc:                       # noqa: BLE001 - record, don't abort
        out["error"] = f"{type(exc).__name__}: {exc}"
        return out

    snowball = res.iceline_lat_nh is None
    out.update({
        "T_global_degC": res.global_mean,
        "T_global_K": res.global_mean + 273.15,
        "iceline_lat_nh": np.nan if snowball else res.iceline_lat_nh,
        "snowball": snowball,
        "peak_diag_insol": res.peak_diag_insol,
        "spinup_years": res.spinup_years,
        # Analytic annual-mean insolation: depends only on a and e.
        "S_mean_Wm2": cfg.S0 / (4.0 * cfg.a_au ** 2 * np.sqrt(1.0 - cfg.ecc ** 2)),
        # Seasonal structure -- invisible in the annual mean above.
        "nh_mean_degC": res.nh_mean,
        "sh_mean_degC": res.sh_mean,
        "diag_summer_max_degC": res.diag_summer_max,
        "diag_winter_min_degC": res.diag_winter_min,
        "diag_seasonal_range_K": res.diag_seasonal_range,
        # Simpson-Nakajima: absorbed flux above the ceiling means no
        # equilibrium exists at all, so T is meaningless for these runs.
        "absorbed_mean_Wm2": res.absorbed_mean,
        "runaway": res.runaway,
    })
    if res.diag_land_summer_max is not None:
        # Two-surface mode: land carries the Milankovitch summer signal that the
        # ocean-damped blend averages away.
        out.update({
            "land_summer_max_degC": res.diag_land_summer_max,
            "land_winter_min_degC": res.diag_land_winter_min,
            "land_seasonal_range_K": res.diag_land_seasonal_range,
            "ocean_summer_max_degC": res.diag_ocean_summer_max,
            "ocean_seasonal_range_K": res.diag_ocean_seasonal_range,
        })
    return out


def _reduce_transient(g: dict, years: np.ndarray, tol: float = 0.1) -> dict:
    """Reduce a transient's per-year arrays to scalar 'how it got there' metrics."""
    gm = np.asarray(g["global_mean"], dtype=float)
    final = gm[-1]

    # Years to equilibrium: first year after which every value stays within tol.
    within = np.abs(gm - final) < tol
    idx = len(gm)
    for i in range(len(gm) - 1, -1, -1):
        if within[i]:
            idx = i
        else:
            break
    years_to_eq = int(years[idx]) if idx < len(years) else int(years[-1])

    # Overshoot: how far past the final state it travelled, in the direction of change.
    warming = final >= gm[0]
    overshoot = float(max(0.0, (gm.max() - final) if warming else (final - gm.min())))

    # Transient hemispheric asymmetry -- often vanishes at equilibrium.
    asym = np.abs(np.asarray(g["nh_mean"], float) - np.asarray(g["sh_mean"], float))
    i_asym = int(np.argmax(asym))

    ice = np.asarray(g["iceline_lat_nh"], dtype=float)
    finite = np.isfinite(ice)
    if finite.sum() >= 2:
        ice_start, ice_end = float(ice[finite][0]), float(ice[finite][-1])
        migration = ice_end - ice_start
        step = np.diff(ice[finite])
        peak_rate = float(np.max(np.abs(step))) if len(step) else 0.0
    else:
        ice_start = ice_end = migration = peak_rate = np.nan

    return {
        "years_to_equilibrium": years_to_eq,
        "T_overshoot_K": overshoot,
        "max_hemis_asymmetry_K": float(asym.max()),
        "max_hemis_asymmetry_year": int(years[i_asym]),
        "final_hemis_asymmetry_K": float(asym[-1]),
        "iceline_start_deg": ice_start,
        "iceline_end_deg": ice_end,
        "iceline_migration_deg": migration,
        "iceline_peak_rate_deg_per_yr": peak_rate,
        "T_year1_degC": float(gm[0]),
        "T_final_degC": float(final),
    }


def _evaluate_transient(base_dict: dict, row: dict, steps_per_year: int,
                        n_years: int) -> dict:
    """Spin up the pre-flyby orbit, switch to the post-flyby orbit, track the path."""
    base = Config(**base_dict)
    pre = config_for_run(base, row, steps_per_year, when="before")
    post = config_for_run(base, row, steps_per_year, when="after")
    out = {"run": row["run"], "a_au": post.a_au, "ecc": post.ecc,
           "obliquity_deg": post.obliquity_deg,
           "lon_perihelion_deg": post.lon_perihelion_deg,
           "days_per_year": post.days_per_year}
    try:
        res = run_perturbation(pre, post, n_years=n_years)
    except Exception as exc:                       # noqa: BLE001 - record, don't abort
        out["error"] = f"{type(exc).__name__}: {exc}"
        return out
    out.update(_reduce_transient(res.diagnostics, res.years))
    out["diag_seasonal_range_K"] = (
        res.diagnostics["Tdiag_summer_max"][-1] - res.diagnostics["Tdiag_winter_min"][-1])
    return out


# ---------------------------------------------------------------------------
# Driver
# ---------------------------------------------------------------------------
def run_all(
    elements_csv: Path,
    base: Config,
    workers: int = 1,
    max_a: float = 3.0,
    max_e: float = 0.9,
    steps_per_year: int = DEFAULT_STEPS_PER_YEAR,
    limit: int | None = None,
    transient: bool = False,
    n_years: int = 60,
) -> tuple[pd.DataFrame, dict]:
    el = pd.read_csv(elements_csv)

    n_total = len(el)
    bound = el[el["bound_after"].astype(bool)].copy()
    usable = bound[
        (bound["a_au_after"] < max_a)
        & (bound["ecc_after"] < max_e)
        & bound["lon_perihelion_deg_after"].notna()
        & bound["obliquity_deg_after"].notna()
    ].copy()
    # Stats describe the whole sweep; --limit only truncates what is evaluated.
    stats = {
        "n_total": n_total,
        "n_bound": len(bound),
        "n_unbound": n_total - len(bound),
        "n_usable": len(usable),
        "n_excluded_band": len(bound) - len(usable),
    }
    if limit:
        usable = usable.head(limit)
    stats["n_evaluated"] = len(usable)

    # Baseline: every run starts from the same pre-flyby orbit, so evaluate once.
    first = usable.iloc[0]
    baseline_row = {
        "run": "__baseline__",
        "a_au_after": first["a_au_before"],
        "ecc_after": first["ecc_before"],
        "obliquity_deg_after": first["obliquity_deg_before"],
        "lon_perihelion_deg_after": first["lon_perihelion_deg_before"],
        "days_per_year_after": 365.256363 * float(first["a_au_before"]) ** 1.5,
    }
    base_dict = base.to_dict()
    baseline = _evaluate(base_dict, baseline_row, steps_per_year)

    records = usable.to_dict("records")
    args = (steps_per_year, n_years) if transient else (steps_per_year,)
    fn = _evaluate_transient if transient else _evaluate

    rows: list[dict] = []
    if workers <= 1:
        for i, r in enumerate(records, 1):
            rows.append(fn(base_dict, r, *args))
            if i % 50 == 0 or i == len(records):
                print(f"  {i}/{len(records)} runs", flush=True)
    else:
        with ProcessPoolExecutor(max_workers=workers) as ex:
            futs = [ex.submit(fn, base_dict, r, *args) for r in records]
            for i, fut in enumerate(as_completed(futs), 1):
                rows.append(fut.result())
                if i % 50 == 0 or i == len(records):
                    print(f"  {i}/{len(records)} runs", flush=True)

    df = pd.DataFrame(rows).sort_values("run").reset_index(drop=True)
    if "T_global_degC" in df.columns:
        df["dT_vs_baseline_K"] = df["T_global_degC"] - baseline["T_global_degC"]
    if "T_final_degC" in df.columns:
        df["dT_vs_baseline_K"] = df["T_final_degC"] - baseline["T_global_degC"]
    # Change in high-latitude summer warmth over land: the Milankovitch
    # glacial-inception signal. Strongly negative => snow can survive the summer.
    for col, base_key in (("land_summer_max_degC", "land_summer_max_degC"),
                          ("diag_summer_max_degC", "diag_summer_max_degC")):
        if col in df.columns and base_key in baseline:
            df[f"d{col.replace('_degC', '')}_vs_baseline_K"] = df[col] - baseline[base_key]
    return df, {**stats, "baseline": baseline}


def join_impact_ranking(df: pd.DataFrame, ranking_csv: Path) -> pd.DataFrame:
    """Merge the orbital-disruption Score from rank_run_impact.py, if present."""
    if not ranking_csv.exists():
        return df
    rank = pd.read_csv(ranking_csv)
    keep = [c for c in ("run", "Score", "n_ejected", "bh_rp_au", "bh_inc_deg",
                        "bh_Omega_deg", "bh_omega_deg") if c in rank.columns]
    return df.merge(rank[keep], on="run", how="left")


# ---------------------------------------------------------------------------
# Plot
# ---------------------------------------------------------------------------
_SEQ_BLUE = ["#cde2fb", "#b7d3f6", "#9ec5f4", "#86b6ef", "#6da7ec", "#5598e7",
             "#3987e5", "#2a78d6", "#256abf", "#1c5cab", "#184f95", "#104281", "#0d366b"]
_INK_PRIMARY, _INK_SECONDARY, _INK_MUTED = "#0b0b0b", "#52514e", "#898781"
_GRIDLINE, _AXIS, _SURFACE = "#e1e0d9", "#c3c2b7", "#fcfcfb"
_CRITICAL = "#d03b3b"


def _style(ax):
    ax.set_facecolor(_SURFACE)
    for name, sp in ax.spines.items():
        sp.set_visible(name not in ("top", "right"))
        sp.set_color(_AXIS)
    ax.tick_params(colors=_INK_MUTED, labelsize=9)
    ax.xaxis.label.set_color(_INK_SECONDARY)
    ax.yaxis.label.set_color(_INK_SECONDARY)
    ax.grid(True, color=_GRIDLINE, linewidth=0.8)
    ax.set_axisbelow(True)


def make_plot(df: pd.DataFrame, baseline: dict, out_path: Path) -> None:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.colors import LinearSegmentedColormap

    ok = df[df["T_global_degC"].notna()]
    cmap = LinearSegmentedColormap.from_list("seq_blue", _SEQ_BLUE, N=256)

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(13, 5.2))
    fig.patch.set_facecolor(_SURFACE)

    # Panel A: equilibrium temperature vs semi-major axis, shaded by eccentricity.
    _style(ax1)
    sc = ax1.scatter(ok["a_au"], ok["T_global_K"], c=ok["ecc"], cmap=cmap,
                     s=40, alpha=0.9, linewidths=0.5, edgecolors=_SURFACE, zorder=3)
    ax1.axhline(baseline["T_global_K"], color=_CRITICAL, lw=1.5, ls="--", zorder=4)
    ax1.annotate(f"pre-flyby  {baseline['T_global_K']:.1f} K",
                 xy=(ok["a_au"].max(), baseline["T_global_K"]),
                 xytext=(-6, 6), textcoords="offset points",
                 ha="right", fontsize=8.5, color=_CRITICAL)
    ax1.set_xlabel("Earth semi-major axis after flyby $a$ [AU]")
    ax1.set_ylabel("Equilibrium global mean $T$ [K]")
    ax1.set_title("Climate outcome vs. post-flyby orbit", color=_INK_PRIMARY,
                  fontsize=11, loc="left")
    cb = fig.colorbar(sc, ax=ax1, fraction=0.046, pad=0.04)
    cb.set_label("eccentricity $e$", color=_INK_SECONDARY, fontsize=9)
    cb.ax.tick_params(colors=_INK_MUTED, labelsize=8)
    cb.outline.set_edgecolor(_AXIS)

    # Panel B: distribution of outcomes -- the ice-albedo bifurcation splits the
    # runs into a glaciated and a temperate branch with a gap between them.
    _style(ax2)
    lo, hi = 175, 420          # clip the handful of extreme-insolation outliers
    clipped = ok[(ok["T_global_K"] >= lo) & (ok["T_global_K"] <= hi)]
    n_out = len(ok) - len(clipped)

    # Shade where the linear OLR (A + B*T, calibrated near 288 K) is defensible.
    ax2.axvspan(250, 300, color="#cde2fb", alpha=0.55, zorder=1,
                label="linear-OLR validity band")
    ax2.hist(clipped["T_global_K"], bins=48, color="#2a78d6",
             edgecolor=_SURFACE, linewidth=0.5, zorder=3)
    ax2.axvline(baseline["T_global_K"], color=_CRITICAL, lw=1.5, ls="--", zorder=4)
    ax2.annotate(f"pre-flyby\n{baseline['T_global_K']:.1f} K",
                 xy=(baseline["T_global_K"], ax2.get_ylim()[1]),
                 xytext=(6, -12), textcoords="offset points",
                 fontsize=8.5, color=_CRITICAL, va="top")
    ax2.set_xlabel("Equilibrium global mean $T$ [K]")
    ax2.set_ylabel("Number of runs")
    title = "Distribution of climate outcomes"
    if n_out:
        title += f"  ({n_out} outside {lo}-{hi} K not shown)"
    ax2.set_title(title, color=_INK_PRIMARY, fontsize=11, loc="left")
    leg = ax2.legend(loc="upper left", frameon=False, fontsize=8.5)
    for t in leg.get_texts():
        t.set_color(_INK_SECONDARY)

    fig.tight_layout()
    fig.savefig(out_path, dpi=150, facecolor=_SURFACE)
    plt.close(fig)


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------
def main() -> int:
    p = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    p.add_argument("parent_dir", type=Path,
                   help="Sweep folder, e.g. simulations/20260724_230314")
    p.add_argument("--elements", type=Path, default=None,
                   help="Earth-elements CSV. Default: <parent>_earth_elements.csv")
    p.add_argument("--config", type=Path, default=None,
                   help="YAML with the shared EBM physics (orbital fields are overridden).")
    p.add_argument("--out", type=Path, default=None, help="Output CSV path.")
    p.add_argument("--plot", type=Path, default=None, help="Write a 2-panel plot here.")
    p.add_argument("--workers", type=int, default=1)
    p.add_argument("--max-a", type=float, default=3.0, help="Exclude a >= this [AU].")
    p.add_argument("--max-e", type=float, default=0.9, help="Exclude e >= this.")
    p.add_argument("--steps-per-year", type=int, default=DEFAULT_STEPS_PER_YEAR)
    p.add_argument("--n-lat", type=int, default=90, help="Latitude cells (speed knob).")
    p.add_argument("--limit", type=int, default=None, help="Only the first N usable runs.")
    p.add_argument("--transient", action="store_true",
                   help="Track the year-by-year adjustment from the pre-flyby climate "
                        "instead of only reporting the final equilibrium.")
    p.add_argument("--years", type=int, default=60,
                   help="Transient length in years (only with --transient).")
    p.add_argument("--olr-model", choices=("linear", "sellers"), default=None,
                   help="OLR parameterisation. 'sellers' is nonlinear and stays "
                        "physical when frozen; 'linear' is the Budyko form.")
    p.add_argument("--two-surface", action="store_true",
                   help="Enable the land/ocean two-surface model, so seasonal extremes "
                        "over land are resolved instead of being ocean-damped.")
    args = p.parse_args()

    elements = args.elements or args.parent_dir.parent / f"{args.parent_dir.name}_earth_elements.csv"
    if not elements.exists():
        raise SystemExit(
            f"Elements CSV not found: {elements}\n"
            f"Run:  python extract_earth_elements.py {args.parent_dir} --workers 5")

    base = load_config(args.config)
    base = dataclasses.replace(base, n_lat=args.n_lat)
    if args.two_surface:
        base = dataclasses.replace(base, two_surface=True)
    if args.olr_model:
        base = dataclasses.replace(base, olr_model=args.olr_model)

    df, info = run_all(elements, base, workers=args.workers, max_a=args.max_a,
                       max_e=args.max_e, steps_per_year=args.steps_per_year,
                       limit=args.limit, transient=args.transient, n_years=args.years)

    ranking = args.parent_dir.parent / f"{args.parent_dir.name}_impact_ranking.csv"
    df = join_impact_ranking(df, ranking)

    suffix = "_climate_transient" if args.transient else "_climate"
    out = args.out or args.parent_dir.parent / f"{args.parent_dir.name}{suffix}.csv"
    df.to_csv(out, index=False)

    b = info["baseline"]
    print(f"\nRuns: {info['n_total']} total | {info['n_unbound']} Earth unbound | "
          f"{info['n_excluded_band']} outside EBM validity band | "
          f"{info['n_usable']} usable | {info['n_evaluated']} evaluated")
    print(f"\nPre-flyby baseline: T = {b['T_global_K']:.2f} K, "
          f"ice edge {b['iceline_lat_nh']:.1f} deg")

    n_err = int(df["error"].notna().sum()) if "error" in df.columns else 0

    if args.transient:
        ok = df[df["T_final_degC"].notna()] if "T_final_degC" in df.columns else df.iloc[:0]
        if len(ok):
            print(f"\nTransient adjustment over {len(ok)} runs "
                  f"({args.years}-year integrations):")
            for label, col, unit in [
                ("years to equilibrium", "years_to_equilibrium", "yr"),
                ("peak hemispheric asymmetry", "max_hemis_asymmetry_K", "K"),
                ("  ...at equilibrium", "final_hemis_asymmetry_K", "K"),
                ("temperature overshoot", "T_overshoot_K", "K"),
                ("ice-edge migration", "iceline_migration_deg", "deg"),
                ("peak migration rate", "iceline_peak_rate_deg_per_yr", "deg/yr"),
                ("seasonal range at diag lat", "diag_seasonal_range_K", "K"),
            ]:
                if col in ok.columns:
                    s = ok[col].dropna()
                    if len(s):
                        print(f"  {label:28s} median {s.median():8.2f}  "
                              f"max {s.max():8.2f}  [{unit}]")
    else:
        ok = df[df["T_global_degC"].notna()] if "T_global_degC" in df.columns else df.iloc[:0]
        if len(ok):
            print(f"\nPost-flyby equilibrium temperature over {len(ok)} runs:")
            for q in (0.0, 0.25, 0.5, 0.75, 1.0):
                print(f"  {int(q * 100):3d}th pct : {ok['T_global_K'].quantile(q):8.2f} K "
                      f"({ok['dT_vs_baseline_K'].quantile(q):+7.2f} K vs baseline)")
            if "snowball" in ok.columns:
                print(f"\nFully glaciated (snowball): {int(ok['snowball'].sum())} runs")
            if "diag_seasonal_range_K" in ok.columns:
                s = ok["diag_seasonal_range_K"].dropna()
                print(f"Seasonal range at diagnostic latitude: "
                      f"median {s.median():.1f} K, max {s.max():.1f} K")
            if "runaway" in ok.columns:
                nr = int(ok["runaway"].sum())
                print(f"Runaway greenhouse (no equilibrium exists): {nr} runs "
                      f"({100 * nr / len(ok):.0f}%) -- their temperatures are not meaningful")
            if "land_seasonal_range_K" in ok.columns:
                sl = ok["land_seasonal_range_K"].dropna()
                so = ok["ocean_seasonal_range_K"].dropna()
                print(f"  over land : median {sl.median():.1f} K, max {sl.max():.1f} K")
                print(f"  over ocean: median {so.median():.1f} K, max {so.max():.1f} K")
                dcol = "dland_summer_max_vs_baseline_K"
                if dcol in ok.columns:
                    d = ok[dcol].dropna()
                    cooler = int((d < -4.0).sum())
                    print(f"\nHigh-latitude summer warmth over land vs pre-flyby:")
                    print(f"  median {d.median():+.2f} K, 5th pct {d.quantile(0.05):+.2f} K")
                    print(f"  runs cooling >4 K (glacial-inception regime): {cooler}"
                          f" ({100 * cooler / len(d):.0f}%)")
    if n_err:
        print(f"Errors: {n_err} runs")
    print(f"\nWrote {out}")

    if args.plot is not None and len(ok) and not args.transient:
        make_plot(df, b, args.plot)
        print(f"Wrote plot to {args.plot}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
