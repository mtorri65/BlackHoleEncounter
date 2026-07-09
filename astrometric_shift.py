#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Astrometric microlensing time series for one Gaia star (Option 2: enforce Rubin cadence explicitly)

What this script does:
1) Loads your synthetic schedule CSV (dense timestamps, RA/Dec of BH track).
2) Loads your Gaia pass catalog and selects one star by source_id.
3) Computes the astrometric centroid shift per timestamp (|shift|, ΔRA, ΔDec) using the same
   centroid-shift formula you've been using.
4) Enforces a Rubin-like cadence by downsampling:
   - keep at most TWO visits per observing night (pair)
   - enforce a minimum spacing between selected nights (e.g. every 3 nights)
   - optionally keep only a seasonal window (already implied by your schedule file)
5) Plots the dense vs cadence-enforced time series for comparison.

Notes:
- Uses an "observing night" definition (UTC-12h) so paired visits around midnight UTC are grouped.
- Works with your existing files:
    synthetic_rubin_schedule_REGENERATED_2026-05-04_to_2026-11-01_baseline3d_pairs30m.csv
    gaia_real_detection_catalog_pass_ge5nights.csv
    and a BH track Excel file matching *bh_radec*.xlsx somewhere under the current directory.
"""

import math
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import argparse


# -----------------------------
# USER SETTINGS
# -----------------------------
SCHED_CSV = "synthetic_rubin_schedule_REGENERATED_2026-05-04_to_2026-11-01_baseline3d_pairs30m.csv"
PASS_CSV  = "gaia_real_detection_catalog_pass_ge5nights.csv"

#TARGET_ID = "6876643760779497216"  # Gaia source_id as string

BH_MASS_Msun = 0.1     # assumed BH mass used in theta_E
K_SIGMA = 3.0          # for plotting detection threshold line

# Rubin-like cadence enforcement knobs (Option 2)
PAIR_VISITS_PER_NIGHT = 2          # keep at most 2 visits per observing night
PAIR_MAX_SEPARATION_MIN = 90       # minutes; second visit must be within this window (e.g., 30–60 typical)
NIGHT_STRIDE_DAYS = 3              # keep at most one observing night every N days (e.g. 3–5 is WFD-ish)
MIN_DAYS_BETWEEN_NIGHTS = 2        # additional guardrail: if nights are too close, skip

# Define "observing night" as UTC shifted by 12h (groups local-night visits)
OBS_NIGHT_UTC_SHIFT_HOURS = 12

# -----------------------------
# Command-line arguments
# -----------------------------
parser = argparse.ArgumentParser(
    description="Astrometric microlensing time series for a Gaia star (Rubin-like cadence)"
)

parser.add_argument(
    "--gaia-id",
    required=True,
    help="Gaia source_id of the star to analyze (e.g. 6876643760779497216)"
)

args = parser.parse_args()

TARGET_ID = str(args.gaia_id)

# -----------------------------
# Helpers
# -----------------------------
def find_bh_track() -> Path:
    """Find the BH RA/Dec track Excel file under current directory."""
    for p in Path(".").rglob("*bh_radec*.xlsx"):
        return p
    raise FileNotFoundError("BH RA/Dec file not found (expected *bh_radec*.xlsx somewhere under cwd).")


def enforce_rubin_cadence(sched: pd.DataFrame) -> pd.DataFrame:
    """
    Enforce a Rubin-like cadence on the schedule:
    - group visits by observing night (UTC-12h)
    - select up to PAIR_VISITS_PER_NIGHT visits per night:
        * pick the earliest visit
        * pick the next visit within PAIR_MAX_SEPARATION_MIN minutes (if available)
    - then downselect nights to mimic WFD-ish cadence:
        * keep one night every NIGHT_STRIDE_DAYS
        * ensure MIN_DAYS_BETWEEN_NIGHTS spacing

    Returns a filtered schedule dataframe (subset of rows).
    """
    s = sched.copy().sort_values("timestamp_utc").reset_index(drop=True)

    # Observing night label
    s["obs_night"] = (s["timestamp_utc"] - pd.Timedelta(hours=OBS_NIGHT_UTC_SHIFT_HOURS)).dt.floor("D")

    # First: pick paired visits per observing night
    keep_idx = []

    for night, g in s.groupby("obs_night", sort=True):
        g = g.sort_values("timestamp_utc")
        if len(g) == 0:
            continue

        # Always keep earliest visit
        first = g.iloc[0]
        keep = [first.name]

        # Optionally keep a second visit within PAIR_MAX_SEPARATION_MIN minutes
        if PAIR_VISITS_PER_NIGHT >= 2 and len(g) > 1:
            t0 = first["timestamp_utc"]
            # candidates within the time window after t0
            dt_min = (g["timestamp_utc"] - t0).dt.total_seconds() / 60.0
            cand = g[(dt_min > 0) & (dt_min <= PAIR_MAX_SEPARATION_MIN)]
            if not cand.empty:
                keep.append(cand.iloc[0].name)

        # If you ever want >2 per night, you can extend this logic.
        keep_idx.extend(keep[:PAIR_VISITS_PER_NIGHT])

    s_pair = s.loc[keep_idx].sort_values("timestamp_utc").reset_index(drop=True)

    # Second: downselect nights (one night every NIGHT_STRIDE_DAYS + min separation)
    # Build one representative row per night (earliest kept)
    reps = s_pair.sort_values("timestamp_utc").groupby("obs_night", as_index=False).first()
    reps = reps.sort_values("obs_night").reset_index(drop=True)

    selected_nights = []
    last_kept = None

    for _, row in reps.iterrows():
        night = row["obs_night"]
        if last_kept is None:
            selected_nights.append(night)
            last_kept = night
            continue

        days_since = (night - last_kept).days
        if days_since < MIN_DAYS_BETWEEN_NIGHTS:
            continue
        if NIGHT_STRIDE_DAYS > 1 and (days_since < NIGHT_STRIDE_DAYS):
            continue

        selected_nights.append(night)
        last_kept = night

    s_final = s_pair[s_pair["obs_night"].isin(selected_nights)].copy()
    s_final = s_final.sort_values("timestamp_utc").reset_index(drop=True)

    # Cleanup helper columns for output
    return s_final


# -----------------------------
# Load data
# -----------------------------
sched = pd.read_csv(SCHED_CSV)
sched["timestamp_utc"] = pd.to_datetime(sched["timestamp_utc"], utc=True)
sched = sched.sort_values("timestamp_utc").reset_index(drop=True)

cat = pd.read_csv(PASS_CSV, dtype={"source_id": "string"})
star_df = cat[cat["source_id"] == TARGET_ID]
if star_df.empty:
    raise ValueError(f"Target source_id {TARGET_ID} not found in {PASS_CSV}")
star = star_df.iloc[0]

sigma_ast = float(star["sigma_ast_mas_per_visit"])
det_line = K_SIGMA * sigma_ast

label_star = (
    f"Gaia {TARGET_ID}  |  "
    f"G={star['phot_g_mean_mag']:.2f}, "
    f"max SNR={star['max_snr']:.1f}"
)

# BH track / range interpolation
bh = pd.read_excel(find_bh_track())
bh["t"] = pd.to_datetime(bh["utc_iso"], utc=True)
bh = bh.sort_values("t").reset_index(drop=True)

tsec_bh = bh["t"].astype("int64").to_numpy() / 1e9
range_au = bh["range_AU"].astype(float).to_numpy()

# Einstein radius helper
G = 6.67430e-11
c = 299792458.0
M_sun = 1.98847e30
AU = 1.495978707e11
M = BH_MASS_Msun * M_sun


def compute_shifts(s: pd.DataFrame):
    """Compute |shift|, ΔRA, ΔDec (mas) for each row in schedule dataframe s."""
    tsec = s["timestamp_utc"].astype("int64").to_numpy() / 1e9
    range_i = np.interp(tsec, tsec_bh, range_au)

    Dl_m = range_i * AU
    thetaE_rad = np.sqrt((4 * G * M / c**2) / Dl_m)
    thetaE_arcsec = thetaE_rad * (180/math.pi) * 3600.0

    ra_l = s["ra_deg"].to_numpy()
    dec_l = s["dec_deg"].to_numpy()

    # Tangent plane centered on THIS schedule subset for stability
    ra0 = float(np.median(ra_l))
    dec0 = float(np.median(dec_l))
    cosd0 = np.cos(np.deg2rad(dec0))

    def tan_xy(ra_deg, dec_deg):
        return (ra_deg - ra0) * cosd0, (dec_deg - dec0)

    lx, ly = tan_xy(ra_l, dec_l)

    # Proper motion propagation (Gaia ref ~2016)
    t = s["timestamp_utc"]
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

    # Components toward lens
    ux = np.divide(dx_arcsec, sep_arcsec, out=np.zeros_like(dx_arcsec), where=sep_arcsec > 0)
    uy = np.divide(dy_arcsec, sep_arcsec, out=np.zeros_like(dy_arcsec), where=sep_arcsec > 0)

    dRA_mas = ux * shift_mas
    dDec_mas = uy * shift_mas

    return shift_mas, dRA_mas, dDec_mas


# -----------------------------
# Compute dense series
# -----------------------------
shift_dense, dRA_dense, dDec_dense = compute_shifts(sched)

# -----------------------------
# Enforce Rubin cadence explicitly (Option 2)
# -----------------------------
sched_cad = enforce_rubin_cadence(sched)
shift_cad, dRA_cad, dDec_cad = compute_shifts(sched_cad)

print("Dense schedule rows:", len(sched))
print("Cadence-enforced rows:", len(sched_cad))
print("Unique observing nights (cadence):", sched_cad["obs_night"].nunique())

# -----------------------------
# Plot: |shift| comparison
# -----------------------------
plt.figure(figsize=(11, 5))

plt.plot(
    sched_cad["timestamp_utc"],
    shift_cad,
    marker="o", ms=5, lw=2,
    label=f"Rubin-like cadence — {label_star}"
)

# Detection line
plt.axhline(det_line, lw=1, linestyle="--", alpha=0.7, label=f"{K_SIGMA:.0f}σ (per-visit)")

# Mark maxima

imax_c = int(np.nanargmax(shift_cad))
plt.scatter([sched_cad.loc[imax_c, "timestamp_utc"]], [shift_cad[imax_c]], s=110, marker="*", zorder=5)
plt.text(sched_cad.loc[imax_c, "timestamp_utc"], shift_cad[imax_c], "  max (cadence)", va="bottom")

plt.ylabel("|Centroid shift| (mas)")
plt.xlabel("Time (UTC)")
plt.title("Astrometric microlensing |shift| vs time (Rubin-like cadence)")
plt.grid(True, alpha=0.4)
plt.legend()
plt.tight_layout()
plt.show()

# -----------------------------
# Plot: components (cadence-enforced only)
# -----------------------------
plt.figure(figsize=(11, 5))

plt.plot(
    sched_cad["timestamp_utc"],
    dRA_cad,
    marker="o", ms=5, lw=2,
    label="ΔRA (mas) — cadence-enforced"
)
plt.plot(
    sched_cad["timestamp_utc"],
    dDec_cad,
    marker="o", ms=5, lw=2,
    label="ΔDec (mas) — cadence-enforced"
)

plt.axhline(0, lw=0.6)
plt.ylabel("Centroid component shift (mas)")
plt.xlabel("Time (UTC)")
plt.title(f"Astrometric shift components vs time (cadence-enforced) — {label_star}")
plt.grid(True, alpha=0.4)
plt.legend()
plt.tight_layout()
plt.show()
