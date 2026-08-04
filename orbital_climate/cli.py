"""Command-line entry point for the orbital climate model.

Subcommands
-----------
    insolation   Insolation diagnostics (annual mean vs analytic; peak-latitude).
    ebm          Spin the EBM up to equilibrium; optionally run a sudden
                 orbital-change transient and write outputs.
    sweep        Run a parallel parameter sweep and write a diagnostics table.

Examples
--------
    python -m orbital_climate.cli insolation --config input_climate.yaml
    python -m orbital_climate.cli insolation --ecc 0.117 --plot insol.png
    python -m orbital_climate.cli ebm --config input_climate.yaml --perturb-ecc 0.117 \
        --years 60 --out-dir climate_runs
    python -m orbital_climate.cli sweep --config input_climate.yaml --workers 5
"""

from __future__ import annotations

import argparse
from datetime import datetime
from pathlib import Path

import numpy as np

from .config import load_config, Config
from .insolation import (
    daily_mean_insolation,
    global_annual_mean_insolation,
    analytic_global_annual_mean,
)

TWO_PI = 2.0 * np.pi


# ---------------------------------------------------------------------------
# Shared config-override handling
# ---------------------------------------------------------------------------
def _add_config_args(p: argparse.ArgumentParser) -> None:
    p.add_argument("--config", default=None, help="Path to a YAML config file.")
    p.add_argument("--S0", type=float, default=None, help="Solar constant at 1 AU [W/m^2].")
    p.add_argument("--a-au", type=float, default=None, help="Semi-major axis [AU].")
    p.add_argument("--ecc", type=float, default=None, help="Eccentricity [-].")
    p.add_argument("--obliquity", type=float, default=None, help="Obliquity [deg].")
    p.add_argument("--lon-perihelion", type=float, default=None,
                   help="Longitude of perihelion [deg].")
    p.add_argument("--diag-lat", type=float, default=None,
                   help="Diagnostic latitude for the seasonal/insolation peak [deg].")


def _overrides(args) -> dict:
    mapping = {
        "S0": getattr(args, "S0", None),
        "a_au": getattr(args, "a_au", None),
        "ecc": getattr(args, "ecc", None),
        "obliquity_deg": getattr(args, "obliquity", None),
        "lon_perihelion_deg": getattr(args, "lon_perihelion", None),
        "diag_lat_deg": getattr(args, "diag_lat", None),
    }
    return {k: v for k, v in mapping.items() if v is not None}


# ---------------------------------------------------------------------------
# insolation subcommand
# ---------------------------------------------------------------------------
def _peak_insolation(lat_deg, config, n_time):
    phi = np.radians(lat_deg)
    M = np.linspace(0.0, TWO_PI, n_time, endpoint=False)
    Q = daily_mean_insolation(phi, M, config)
    i = int(np.argmax(Q))
    return float(Q[i]), float(M[i])


def cmd_insolation(args) -> int:
    config = load_config(args.config, **_overrides(args))
    # Peak-insolation latitude: explicit --peak-lat wins, else the config's
    # diagnostic latitude (which --diag-lat / YAML can set).
    peak_lat = args.peak_lat if args.peak_lat is not None else config.diag_lat_deg
    gmean = global_annual_mean_insolation(config, n_time=args.n_time)
    analytic = analytic_global_annual_mean(config)
    peak_val, peak_M = _peak_insolation(peak_lat, config, args.n_time)

    print("Orbital / stellar configuration")
    print(f"  S0              = {config.S0:.3f} W/m^2")
    print(f"  a               = {config.a_au:.6f} AU")
    print(f"  eccentricity    = {config.ecc:.6f}")
    print(f"  obliquity       = {config.obliquity_deg:.4f} deg")
    print(f"  lon. perihelion = {config.lon_perihelion_deg:.4f} deg")
    print()
    print("Global-and-annual mean insolation")
    print(f"  numerical       = {gmean:.4f} W/m^2")
    print(f"  analytic        = {analytic:.4f} W/m^2  [S0 / (4 a^2 sqrt(1-e^2))]")
    print(f"  rel. error      = {abs(gmean - analytic) / analytic:.2e}")
    print()
    print(f"Peak daily-mean insolation at {peak_lat:.1f} deg latitude")
    print(f"  value           = {peak_val:.2f} W/m^2")
    print(f"  at day-of-orbit = {peak_M / TWO_PI * config.days_per_year:.1f}"
          f"  (perihelion at day 0)")

    if args.plot:
        _make_insolation_plot(config, args.plot)
        print(f"\nWrote insolation plot to {args.plot}")
    return 0


def _make_insolation_plot(config, path, n_time=361, n_lat=181):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    lat_deg = np.linspace(-90.0, 90.0, n_lat)
    phi = np.radians(lat_deg)
    M = np.linspace(0.0, TWO_PI, n_time, endpoint=False)
    Q = daily_mean_insolation(phi, M, config)

    fig, ax = plt.subplots(figsize=(8, 5))
    day = M / TWO_PI * config.days_per_year
    mesh = ax.pcolormesh(day, lat_deg, Q.T, shading="auto", cmap="inferno")
    ax.set_xlabel("Day of orbit (perihelion at 0)")
    ax.set_ylabel("Latitude [deg]")
    ax.set_title(f"Daily-mean insolation (e={config.ecc:.3f}, "
                 f"lon_peri={config.lon_perihelion_deg:.0f} deg)")
    fig.colorbar(mesh, ax=ax, label="W m$^{-2}$")
    fig.tight_layout()
    fig.savefig(path, dpi=130)
    plt.close(fig)


# ---------------------------------------------------------------------------
# ebm subcommand
# ---------------------------------------------------------------------------
def cmd_ebm(args) -> int:
    from .experiment import (
        run_equilibrium, run_perturbation,
        perturbed_config_from_scenario, save_perturbation_outputs,
    )

    config = load_config(args.config, **_overrides(args))

    eq = run_equilibrium(config)
    print("Baseline equilibrium")
    print(f"  spin-up years   = {eq.spinup_years}")
    print(f"  global mean T   = {eq.global_mean:.3f} degC = {eq.global_mean + 273.15:.2f} K")
    print(f"  NH ice edge     = "
          f"{'none' if eq.iceline_lat_nh is None else f'{eq.iceline_lat_nh:.2f} deg'}")
    print(f"  {eq.diag_lat_deg:.0f}N peak insol = {eq.peak_diag_insol:.2f} W/m^2")

    if args.perturb_ecc is None:
        return 0

    perturbed = perturbed_config_from_scenario(
        config, ecc=args.perturb_ecc, lon_perihelion_deg=args.perturb_lon_perihelion)
    result = run_perturbation(config, perturbed, n_years=args.years)

    stamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    out_dir = Path(args.out_dir) / stamp
    save_perturbation_outputs(result, out_dir)

    d = result.diagnostics
    print("\nSudden-perturbation transient")
    print(f"  perturbed e     = {perturbed.ecc:.4f}")
    print(f"  years           = {args.years}")
    print(f"  global mean T   = {d['global_mean'][0]:.3f} -> {d['global_mean'][-1]:.3f} degC")
    print(f"  {result.diag_lat_deg:.0f}N summer max  = {d['Tdiag_summer_max'][0]:.3f} -> "
          f"{d['Tdiag_summer_max'][-1]:.3f} degC")
    print(f"  NH ice edge     = {d['iceline_lat_nh'][0]:.2f} -> {d['iceline_lat_nh'][-1]:.2f} deg")
    print(f"\nWrote outputs to {out_dir}")
    return 0


# ---------------------------------------------------------------------------
# sweep subcommand
# ---------------------------------------------------------------------------
def cmd_sweep(args) -> int:
    import yaml
    from .sweep import run_sweep

    config = load_config(args.config, **_overrides(args))

    ranges = {}
    if args.config is not None:
        with open(args.config, "r", encoding="utf-8") as fh:
            raw = yaml.safe_load(fh) or {}
        ranges = raw.get("sweep", {}) or {}
    if not ranges:
        raise SystemExit(
            "No sweep ranges found. Add a top-level 'sweep:' mapping "
            "(field -> list of values) to the config file.")

    print(f"Sweeping {list(ranges.keys())} with up to {args.workers} worker(s)...")
    run_dir = run_sweep(config, ranges, out_dir=args.out_dir, workers=args.workers)
    print(f"Wrote sweep table to {run_dir}")
    return 0


# ---------------------------------------------------------------------------
# Parser assembly
# ---------------------------------------------------------------------------
def _build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(
        prog="orbital_climate.cli",
        description="Orbital perturbation climate model.",
    )
    sub = p.add_subparsers(dest="command", required=True)

    p_ins = sub.add_parser("insolation", help="Insolation diagnostics.")
    _add_config_args(p_ins)
    p_ins.add_argument("--peak-lat", type=float, default=None,
                       help="Latitude for the peak-insolation diagnostic [deg]. "
                            "Defaults to the config's diag_lat_deg.")
    p_ins.add_argument("--n-time", type=int, default=2000,
                       help="Orbital samples for time averages / peak search.")
    p_ins.add_argument("--plot", default=None, help="Write a latitude x season plot here.")
    p_ins.set_defaults(func=cmd_insolation)

    p_ebm = sub.add_parser("ebm", help="EBM equilibrium + optional perturbation.")
    _add_config_args(p_ebm)
    p_ebm.add_argument("--perturb-ecc", type=float, default=None,
                       help="If set, run a transient with this perturbed eccentricity.")
    p_ebm.add_argument("--perturb-lon-perihelion", type=float, default=None,
                       help="Optional perturbed longitude of perihelion [deg].")
    p_ebm.add_argument("--years", type=int, default=60, help="Transient length [years].")
    p_ebm.add_argument("--out-dir", default="climate_runs",
                       help="Parent directory for transient outputs.")
    p_ebm.set_defaults(func=cmd_ebm)

    p_sw = sub.add_parser("sweep", help="Parallel parameter sweep.")
    _add_config_args(p_sw)
    p_sw.add_argument("--workers", type=int, default=1, help="Worker processes.")
    p_sw.add_argument("--out-dir", default="climate_sweeps",
                      help="Parent directory for sweep outputs.")
    p_sw.set_defaults(func=cmd_sweep)

    return p


def main(argv=None) -> int:
    args = _build_parser().parse_args(argv)
    return args.func(args)


if __name__ == "__main__":
    raise SystemExit(main())
