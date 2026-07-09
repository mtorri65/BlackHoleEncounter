#!/usr/bin/env python3
"""
Postprocess belt_run_before.csv and belt_run_after.csv to:
- assign synthetic diameters to tracers (power-law size distribution)
- compute orbital deltas
- estimate Earth-crossing frequency
- estimate aggregate impact rate scaled to real asteroid population
- generate histogram of Earth-crossers by size bin (sample + scaled)
"""

import numpy as np
import pandas as pd
import argparse
import textwrap
import math
from pathlib import Path
import matplotlib.pyplot as plt


# ---------------------------
# Helper: assign synthetic sizes
# ---------------------------
def assign_sizes(n, dmin_km=0.05, dmax_km=300.0, slope=2.5, seed=1234):
    """Draw diameters (km) from a cumulative power law N(>D) ~ D^-slope."""
    rng = np.random.default_rng(seed)
    u = rng.random(n)
    emin, emax = dmin_km ** (-slope + 1), dmax_km ** (-slope + 1)
    diam_km = (emin + (emax - emin) * u) ** (1.0 / (-slope + 1))
    return diam_km


# ---------------------------
# Helper: make histogram + plots
# ---------------------------
def make_histogram(df, folder, f_cross_scale=1.0, bins_km=None, save_svg=False):
    """Generate histogram of Earth-crossers by diameter and save CSV + PNG (+SVG)."""
    if bins_km is None:
        # logarithmic binning from 0.1 km to 300 km
        bins_km = np.logspace(-1, 2.5, 25)

    crossers = df[(df["q_after"] < 1.0)]
    diam = crossers["diameter_km"]

    hist_counts, edges = np.histogram(diam, bins=bins_km)
    mids = np.sqrt(edges[:-1] * edges[1:])
    scaled_counts = hist_counts * f_cross_scale

    hist_df = pd.DataFrame({
        "diameter_bin_mid_km": mids,
        "count_sample": hist_counts,
        "count_scaled": scaled_counts
    })

    hist_csv = folder / "belt_crossers_size_histogram.csv"
    hist_png = folder / "belt_crossers_size_histogram.png"
    hist_svg = folder / "belt_crossers_size_histogram.svg"

    hist_df.to_csv(hist_csv, index=False)

    # --- Plot histogram ---
    fig, ax = plt.subplots(figsize=(7, 5))
    ax.loglog(mids, np.maximum(hist_counts, 1), "o-", label="Sample", color="tab:blue")
    ax.loglog(mids, np.maximum(scaled_counts, 1), "s--", label="Scaled (real)", color="tab:orange")
    ax.set_xlabel("Diameter [km]")
    ax.set_ylabel("Earth-crossers (count)")
    ax.set_title("Post-Encounter Earth-Crossers by Size")
    ax.legend()
    ax.grid(True, which="both", ls=":")
    plt.tight_layout()

    # Save in multiple formats
    plt.savefig(hist_png, dpi=200)
    if save_svg:
        plt.savefig(hist_svg)
    plt.close(fig)

    print(f"Histogram written to {hist_csv}")
    print(f"Plot saved to {hist_png}")
    if save_svg:
        print(f"SVG vector plot saved to {hist_svg}")


# ---------------------------
# Helper: locate CSVs (plain or prefixed)
# ---------------------------
def _pick_belt_paths(folder: Path):
    """
    Return (before_path, after_path) trying both plain names and prefixed names
    like '<folder.name>__belt_run_before.csv'.
    """
    prefix = folder.name + "__"
    # Find candidates
    before_candidates = list(folder.glob(f"{prefix}*belt_run_before.csv"))
    after_candidates  = list(folder.glob(f"{prefix}*belt_run_after.csv"))

    before_path = (before_candidates[0] if before_candidates
                   else folder / "belt_run_before.csv")
    after_path  = (after_candidates[0] if after_candidates
                   else folder / "belt_run_after.csv")

    if not before_path.exists() or not after_path.exists():
        raise FileNotFoundError(
            "belt_run_before.csv and belt_run_after.csv must both exist "
            "(plain or prefixed) in the specified folder"
        )
    print(f"Using BEFORE CSV : {before_path.name}")
    print(f"Using AFTER  CSV : {after_path.name}")
    return before_path, after_path


# ---------------------------
# Main
# ---------------------------
def main():
    ap = argparse.ArgumentParser(
        description="Analyze asteroid belt perturbations and estimate impact hazard."
    )
    ap.add_argument("--folder", required=True,
                    help="Simulation folder containing belt_run_before.csv and belt_run_after.csv (plain or prefixed)")
    ap.add_argument("--impact_prob", type=float, default=1e-8,
                    help="Annual impact probability per Earth-crossing ≥1 km asteroid (default=1e-8)")
    ap.add_argument("--real_gt1km", type=float, default=1.2e6,
                    help="Assumed real count of ≥1 km main-belt asteroids (default=1.2e6)")
    ap.add_argument("--horizon_years", type=float, default=10.0,
                    help="Horizon for cumulative probability P(≤T) (default=10 years)")
    ap.add_argument("--seed", type=int, default=1234,
                    help="Random seed for size assignment")
    ap.add_argument("--save_hist_svg", action="store_true",
                    help="Also save histogram as SVG vector graphic")
    args = ap.parse_args()

    folder = Path(args.folder)
    summary_path = folder / "hazard_summary.txt"
    delta_path = folder / "belt_run_deltas.csv"

    # Locate CSVs (handles plain and prefixed names)
    before_path, after_path = _pick_belt_paths(folder)

    print(f"Loading belt runs from: {folder}")
    df_before = pd.read_csv(before_path)
    df_after = pd.read_csv(after_path)

    if len(df_before) != len(df_after):
        raise ValueError("Mismatch in number of belt entries between before and after CSVs")

    n = len(df_before)
    diam_km = assign_sizes(n, seed=args.seed)

    df_merged = pd.DataFrame({
        "a_before": df_before["a_before"],
        "a_after": df_after["a_after"],
        "e_before": df_before["e_before"],
        "e_after": df_after["e_after"],
        "i_before": df_before["i_before"],
        "i_after": df_after["i_after"],
        "q_before": df_before["q_before"],
        "q_after": df_after["q_after"],
        "diameter_km": diam_km,
    })

    df_merged["da"] = df_merged["a_after"] - df_merged["a_before"]
    df_merged["de"] = df_merged["e_after"] - df_merged["e_before"]
    df_merged["di_deg"] = np.degrees(df_merged["i_after"] - df_merged["i_before"])
    df_merged.to_csv(delta_path, index=False)
    print(f"Wrote deltas and synthetic sizes to {delta_path}")

    # --- counts and scaling ---
    P_impact = args.impact_prob
    n_gt1km = int((df_merged["diameter_km"] > 1.0).sum())
    n_cross_gt1km = int(((df_merged["diameter_km"] > 1.0) & (df_merged["q_after"] < 1.0)).sum())
    f_cross = 0.0 if n_gt1km == 0 else (n_cross_gt1km / n_gt1km)

    N_real_gt1km = float(args.real_gt1km)
    N_real_cross = f_cross * N_real_gt1km
    R_real = N_real_cross * P_impact
    T_mean = float("inf") if R_real <= 0 else 1.0 / R_real
    P_T = 1.0 - math.exp(-R_real * args.horizon_years)

    txt = textwrap.dedent(f"""
        === Belt Impact Hazard Summary ===
        Folder: {folder}
        Tracers analyzed: {n}
        Assigned size slope: 2.5
        Seed: {args.seed}

        --- Sample stats ---
        ≥1 km tracers                : {n_gt1km}
        ≥1 km Earth-crossers (q<1 AU): {n_cross_gt1km}
        Fraction f_cross             : {f_cross:.6f}

        --- Real-belt scaling ---
        Assumed real ≥1 km count     : {N_real_gt1km:,.0f}
        Estimated real ≥1 km crossers: {N_real_cross:,.0f}

        --- Hazard metrics ---
        Per-object impact prob       : {P_impact:.1e} yr^-1
        Aggregate annual rate R_real : {R_real:.3e} yr^-1
        Mean waiting time            : {"∞" if not math.isfinite(T_mean) else f"{T_mean:,.0f} years"}
        P(≤{args.horizon_years:.0f} years)             : {P_T*100:.4f} %

        (Assumes simple Poisson process with constant R_real.)
    """).strip("\n")

    with open(summary_path, "w", encoding="utf-8") as f:
        f.write(txt)
    print(txt)
    print(f"\nSummary written to {summary_path}")

    # --- histogram ---
    f_cross_scale = (N_real_gt1km / n) if n > 0 else 1.0
    make_histogram(df_merged, folder, f_cross_scale=f_cross_scale, save_svg=args.save_hist_svg)


if __name__ == "__main__":
    main()
