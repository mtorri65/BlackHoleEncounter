#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Plot the astrometric microlensing-induced centroid motion of a selected Gaia star
in the sky plane (ΔRA vs ΔDec) over time, using your Rubin-like cadence (Option 2).

What you get:
- Figure 1: Centroid track in mas (ΔRA vs ΔDec), arrows show time direction.
- Figure 2 (optional): Same track colored by time (if you want it).

Inputs (command line):
  --gaia-id   Gaia source_id (required)
Optional:
  --sched     schedule CSV path
  --cat       Gaia pass catalog CSV path
  --bh-mass   BH mass in Msun (default 0.1)
  --k-sigma   detection threshold multiplier (default 3)
  --pair-max-min   max minutes between pair visits (default 90)
  --stride-days    keep one observing night every N days (default 3)
  --min-days-between-nights minimal spacing between kept nights (default 2)

Assumptions:
- BH track Excel file (*bh_radec*.xlsx) exists somewhere under the current directory.
- Schedule CSV contains columns: timestamp_utc, ra_deg, dec_deg
- Pass catalog contains: source_id, ra, dec, pmra, pmdec, sigma_ast_mas_per_visit, phot_g_mean_mag, max_snr
"""

import argparse
import math
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt


# -----------------------------
# Helpers
# -----------------------------
def find_bh_track() -> Path:
    for p in Path(".").rglob("*bh_radec*.xlsx"):
        return p
    raise FileNotFoundError("BH RA/Dec file not found (expected *bh_radec*.xlsx under cwd).")


def enforce_rubin_cadence(
    sched: pd.DataFrame,
    pair_visits_per_night: int = 2,
    pair_max_separation_min: int = 90,
    night_stride_days: int = 3,
    min_days_between_nights: int = 2,
    obs_night_utc_shift_hours: int = 12
) -> pd.DataFrame:
    """Downsample schedule to Rubin-like cadence: paired visits and night spacing."""
    s = sched.copy().sort_values("timestamp_utc").reset_index(drop=True)

    s["obs_night"] = (s["timestamp_utc"] - pd.Timedelta(hours=obs_night_utc_shift_hours)).dt.floor("D")

    keep_idx = []
    for night, g in s.groupby("obs_night", sort=True):
        g = g.sort_values("timestamp_utc")
        if g.empty:
            continue

        first = g.iloc[0]
        keep = [first.name]

        if pair_visits_per_night >= 2 and len(g) > 1:
            t0 = first["timestamp_utc"]
            dt_min = (g["timestamp_utc"] - t0).dt.total_seconds() / 60.0
            cand = g[(dt_min > 0) & (dt_min <= pair_max_separation_min)]
            if not cand.empty:
                keep.append(cand.iloc[0].name)

        keep_idx.extend(keep[:pair_visits_per_night])

    s_pair = s.loc[keep_idx].sort_values("timestamp_utc").reset_index(drop=True)

    reps = s_pair.groupby("obs_night", as_index=False).first().sort_values("obs_night").reset_index(drop=True)

    selected_nights = []
    last_kept = None
    for _, row in reps.iterrows():
        night = row["obs_night"]
        if last_kept is None:
            selected_nights.append(night)
            last_kept = night
            continue

        days_since = (night - last_kept).days
        if days_since < min_days_between_nights:
            continue
        if night_stride_days > 1 and days_since < night_stride_days:
            continue

        selected_nights.append(night)
        last_kept = night

    out = s_pair[s_pair["obs_night"].isin(selected_nights)].copy()
    out = out.sort_values("timestamp_utc").reset_index(drop=True)
    return out


def compute_centroid_components_mas(
    sched: pd.DataFrame,
    star: pd.Series,
    bh_tsec: np.ndarray,
    bh_range_au: np.ndarray,
    bh_mass_msun: float
):
    """
    Compute centroid shift components (ΔRA, ΔDec) in mas for each schedule row.
    Uses:
      delta = (u/(u^2+2))*theta_E, direction toward lens
    """
    # constants
    G = 6.67430e-11
    c = 299792458.0
    M_sun = 1.98847e30
    AU = 1.495978707e11

    # interpolate lens distance
    tsec = sched["timestamp_utc"].astype("int64").to_numpy() / 1e9
    Dl_au = np.interp(tsec, bh_tsec, bh_range_au)
    Dl_m = Dl_au * AU

    M = bh_mass_msun * M_sun
    thetaE_rad = np.sqrt((4 * G * M / c**2) / Dl_m)
    thetaE_arcsec = thetaE_rad * (180/math.pi) * 3600.0

    ra_l = sched["ra_deg"].to_numpy()
    dec_l = sched["dec_deg"].to_numpy()

    # Tangent plane centered on this schedule subset (stable)
    ra0 = float(np.median(ra_l))
    dec0 = float(np.median(dec_l))
    cosd0 = np.cos(np.deg2rad(dec0))

    def tan_xy(ra_deg, dec_deg):
        return (ra_deg - ra0) * cosd0, (dec_deg - dec0)

    lx, ly = tan_xy(ra_l, dec_l)

    # proper motion propagate star from Gaia epoch ~2016 to each visit time
    t = sched["timestamp_utc"]
    dt_yr = (t.dt.year + t.dt.dayofyear/365.25 - 2016.0).to_numpy()

    ra_s = star["ra"] + (star["pmra"] * dt_yr) / (np.cos(np.deg2rad(star["dec"])) * 3.6e6)
    dec_s = star["dec"] + (star["pmdec"] * dt_yr) / 3.6e6

    sx, sy = tan_xy(ra_s, dec_s)

    dx_arcsec = (sx - lx) * 3600.0
    dy_arcsec = (sy - ly) * 3600.0
    sep_arcsec = np.sqrt(dx_arcsec**2 + dy_arcsec**2)

    u = sep_arcsec / thetaE_arcsec
    shift_arcsec = (u / (u*u + 2.0)) * thetaE_arcsec
    shift_mas = shift_arcsec * 1000.0

    # components: along separation unit vector
    ux = np.divide(dx_arcsec, sep_arcsec, out=np.zeros_like(dx_arcsec), where=sep_arcsec > 0)
    uy = np.divide(dy_arcsec, sep_arcsec, out=np.zeros_like(dy_arcsec), where=sep_arcsec > 0)

    dRA_mas = ux * shift_mas
    dDec_mas = uy * shift_mas

    return dRA_mas, dDec_mas, shift_mas


# -----------------------------
# Main
# -----------------------------
def main():
    ap = argparse.ArgumentParser(description="Plot ΔDec vs ΔRA centroid track over time for a Gaia star.")
    ap.add_argument("--gaia-id", required=True, help="Gaia source_id (e.g. 6876643760779497216)")
    ap.add_argument("--sched", default="synthetic_rubin_schedule_REGENERATED_2026-05-04_to_2026-11-01_baseline3d_pairs30m.csv",
                    help="Schedule CSV with timestamp_utc, ra_deg, dec_deg")
    ap.add_argument("--cat", default="gaia_real_detection_catalog_pass_ge5nights.csv",
                    help="Gaia pass catalog CSV")
    ap.add_argument("--bh-mass", type=float, default=0.1, help="BH mass in Msun (default 0.1)")
    ap.add_argument("--pair-max-min", type=int, default=90, help="Max minutes between paired visits")
    ap.add_argument("--stride-days", type=int, default=3, help="Keep one observing night every N days")
    ap.add_argument("--min-days-between-nights", type=int, default=2, help="Minimum spacing between kept nights")
    ap.add_argument("--obs-night-shift-hours", type=int, default=12, help="UTC shift to define observing night")
    ap.add_argument("--arrow-step", type=int, default=1, help="Arrow every N points on the sky-track plot")
    args = ap.parse_args()

    target_id = str(args.gaia_id)

    # Load schedule
    sched = pd.read_csv(args.sched)
    sched["timestamp_utc"] = pd.to_datetime(sched["timestamp_utc"], utc=True)
    sched = sched.sort_values("timestamp_utc").reset_index(drop=True)

    # Load star catalog
    cat = pd.read_csv(args.cat, dtype={"source_id": "string"})
    star_df = cat[cat["source_id"] == target_id]
    if star_df.empty:
        raise ValueError(f"Target source_id {target_id} not found in {args.cat}")
    star = star_df.iloc[0]

    # Find BH track and load range vs time
    bh = pd.read_excel(find_bh_track())
    bh["t"] = pd.to_datetime(bh["utc_iso"], utc=True)
    bh = bh.sort_values("t").reset_index(drop=True)
    bh_tsec = bh["t"].astype("int64").to_numpy() / 1e9
    bh_range_au = bh["range_AU"].astype(float).to_numpy()

    # Enforce Rubin cadence (Option 2)
    sched_cad = enforce_rubin_cadence(
        sched,
        pair_visits_per_night=2,
        pair_max_separation_min=args.pair_max_min,
        night_stride_days=args.stride_days,
        min_days_between_nights=args.min_days_between_nights,
        obs_night_utc_shift_hours=args.obs_night_shift_hours,
    )

    # Compute centroid shift components for cadence schedule
    dRA_mas, dDec_mas, shift_mas = compute_centroid_components_mas(
        sched_cad, star, bh_tsec, bh_range_au, args.bh_mass
    )

    # Label
    label_star = (
        f"Gaia {target_id} | G={star['phot_g_mean_mag']:.2f}, "
        f"max SNR={star['max_snr']:.1f}"
    )

    # -----------------------------
    # Figure 1: ΔDec vs ΔRA with arrows
    # -----------------------------
    plt.figure(figsize=(7.5, 7.5))
    plt.plot(dRA_mas, dDec_mas, marker="o", ms=4, lw=1.5, label=label_star)

    # arrows indicating time direction
    step = max(1, int(args.arrow_step))
    for i in range(0, len(dRA_mas) - step, step):
        plt.annotate(
            "",
            xy=(dRA_mas[i + step], dDec_mas[i + step]),
            xytext=(dRA_mas[i], dDec_mas[i]),
            arrowprops=dict(arrowstyle="->", lw=1.0, alpha=0.6),
        )

    # mark maximum shift point
    imax = int(np.nanargmax(shift_mas))
    plt.scatter([dRA_mas[imax]], [dDec_mas[imax]], s=140, marker="*", zorder=5)
    plt.text(dRA_mas[imax], dDec_mas[imax], "  max", va="bottom")

    plt.axhline(0, lw=0.6)
    plt.axvline(0, lw=0.6)
    plt.gca().set_aspect("equal", adjustable="box")
    plt.xlabel("ΔRA (mas)")
    plt.ylabel("ΔDec (mas)")
    plt.title("Centroid motion in sky plane (ΔDec vs ΔRA)")
    plt.grid(True, alpha=0.4)
    plt.legend()
    plt.tight_layout()
    plt.show()

    # -----------------------------
    # Figure 2: optional — RA/Dec components vs time (helpful context)
    # -----------------------------
    plt.figure(figsize=(11, 5))
    plt.plot(sched_cad["timestamp_utc"], dRA_mas, marker="o", ms=4, lw=1.5, label="ΔRA (mas)")
    plt.plot(sched_cad["timestamp_utc"], dDec_mas, marker="o", ms=4, lw=1.5, label="ΔDec (mas)")
    plt.axhline(0, lw=0.6)
    plt.xlabel("Time (UTC)")
    plt.ylabel("Centroid shift component (mas)")
    plt.title(f"Centroid components vs time — {label_star}")
    plt.grid(True, alpha=0.4)
    plt.legend()
    plt.tight_layout()
    plt.show()


if __name__ == "__main__":
    main()
