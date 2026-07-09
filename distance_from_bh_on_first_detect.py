#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import argparse
from pathlib import Path
import numpy as np
import pandas as pd


# -----------------------------
# Helpers
# -----------------------------
def find_bh_track_from_name(filename: str, base_dir: Path) -> Path:
    """
    Locate BH RA/Dec xlsx given only its filename.
    Expected layout:
      base_dir/simulations/<timestamp>/**/<filename>
    """
    timestamp = filename.split("__")[0]
    sim_root = base_dir / "simulations" / timestamp
    if not sim_root.exists():
        raise FileNotFoundError(f"Simulation folder not found: {sim_root}")

    matches = list(sim_root.rglob(filename))
    if not matches:
        raise FileNotFoundError(f"Could not find {filename} under {sim_root}")

    if len(matches) > 1:
        print("⚠️ Multiple matches found; using first:")
        for m in matches:
            print("   ", m)

    return matches[0]


def auto_find_any_bh_track(base_dir: Path) -> Path:
    """Fallback: find any *bh_radec*.xlsx under simulations."""
    sims = base_dir / "simulations"
    if not sims.exists():
        raise FileNotFoundError(f"No simulations folder found at: {sims}")
    matches = list(sims.rglob("*bh_radec*.xlsx"))
    if not matches:
        raise FileNotFoundError("Could not find any *bh_radec*.xlsx under ./simulations/")
    # pick the most recently modified
    matches.sort(key=lambda p: p.stat().st_mtime, reverse=True)
    return matches[0]


def angular_sep_deg(ra1_deg, dec1_deg, ra2_deg, dec2_deg):
    """
    Great-circle separation (degrees) using a numerically stable haversine formulation.
    Inputs can be numpy arrays.
    """
    ra1 = np.deg2rad(ra1_deg)
    dec1 = np.deg2rad(dec1_deg)
    ra2 = np.deg2rad(ra2_deg)
    dec2 = np.deg2rad(dec2_deg)

    d_ra = ra2 - ra1
    d_dec = dec2 - dec1

    sin_ddec = np.sin(d_dec / 2.0)
    sin_dra = np.sin(d_ra / 2.0)

    a = sin_ddec**2 + np.cos(dec1) * np.cos(dec2) * sin_dra**2
    c = 2.0 * np.arcsin(np.minimum(1.0, np.sqrt(a)))
    return np.rad2deg(c)


def parse_first_detect_date(series: pd.Series) -> pd.Series:
    """
    Parse first_detect_utc column which may look like '5/4/2026' or '2026-05-04'
    Returns pandas datetime64[ns, UTC] normalized to date.
    """
    dt = pd.to_datetime(series.astype(str).str.strip(), errors="coerce", utc=True, infer_datetime_format=True)
    # Some values might parse as NaT; try explicit month/day/year if needed
    if dt.isna().any():
        dt2 = pd.to_datetime(series.astype(str).str.strip(), errors="coerce", utc=True, format="%m/%d/%Y")
        dt = dt.fillna(dt2)
    return dt.dt.normalize()


# -----------------------------
# Main
# -----------------------------
def main():
    ap = argparse.ArgumentParser(
        description="For stars with first_detect_utc=2026-05-04, compute angular distance to BH on that date."
    )
    ap.add_argument("--pass-csv", default="gaia_real_detection_catalog_pass_ge5nights_EXCELSAFE.csv",
                    help="Input pass catalog CSV (EXCELSAFE is fine).")
    ap.add_argument("--bh-radec", default=None,
                    help="BH track xlsx filename (e.g. 20260121_114641__...__bh_radec__1p5__0p1.xlsx). "
                         "If omitted, script auto-finds the most recent *bh_radec*.xlsx under ./simulations.")
    ap.add_argument("--date", default="2026-05-04",
                    help="Date to evaluate BH position (UTC). Default 2026-05-04.")
    ap.add_argument("--out", default="gaia_pass_ge5nights_firstdetect_2026-05-04_sep_to_bh.csv",
                    help="Output CSV filename.")
    args = ap.parse_args()

    base_dir = Path.cwd()

    # ---- Load pass catalog (read everything as string first to avoid Excel-safe quirks)
    df = pd.read_csv(args.pass_csv, dtype=str)

    # Normalize source_id if EXCELSAFE uses a leading apostrophe
    if "source_id" in df.columns:
        df["source_id"] = df["source_id"].astype(str).str.strip().str.lstrip("'")
    else:
        raise KeyError(f"'source_id' column not found in {args.pass_csv}. Columns: {df.columns.tolist()}")

    # Ensure ra/dec numeric
    for col in ["ra", "dec"]:
        if col not in df.columns:
            raise KeyError(f"Missing required column '{col}' in {args.pass_csv}. Columns: {df.columns.tolist()}")
        df[col] = pd.to_numeric(df[col], errors="coerce")

    if "first_detect_utc" not in df.columns:
        raise KeyError(f"'first_detect_utc' not found in {args.pass_csv}. Columns: {df.columns.tolist()}")

    df["first_detect_date"] = parse_first_detect_date(df["first_detect_utc"])

    target_date = pd.to_datetime(args.date, utc=True).normalize()
    sub = df[df["first_detect_date"] == target_date].copy()

    print(f"Rows in pass catalog: {len(df)}")
    print(f"Rows with first_detect_utc date == {target_date.date()}: {len(sub)}")

    if sub.empty:
        print("No matching rows; nothing to do.")
        return

    # ---- Load BH track
    if args.bh_radec:
        bh_path = find_bh_track_from_name(args.bh_radec, base_dir)
    else:
        bh_path = auto_find_any_bh_track(base_dir)

    print(f"Using BH track: {bh_path}")
    bh = pd.read_excel(bh_path)

    # Your BH track columns:
    # time_days utc_iso RA_hours Dec_deg range_AU ...
    if "utc_iso" not in bh.columns:
        raise KeyError(f"BH track missing 'utc_iso'. Columns: {bh.columns.tolist()}")

    bh["t"] = pd.to_datetime(bh["utc_iso"], utc=True, errors="coerce")
    bh = bh.dropna(subset=["t"]).sort_values("t").reset_index(drop=True)

    if "Dec_deg" not in bh.columns:
        raise KeyError(f"BH track missing 'Dec_deg'. Columns: {bh.columns.tolist()}")

    if "RA_hours" in bh.columns:
        bh_ra_deg = pd.to_numeric(bh["RA_hours"], errors="coerce") * 15.0
    else:
        raise KeyError(f"BH track missing 'RA_hours'. Columns: {bh.columns.tolist()}")

    bh_dec_deg = pd.to_numeric(bh["Dec_deg"], errors="coerce")

    bh_tsec = bh["t"].astype("int64").to_numpy() / 1e9
    ra_arr = bh_ra_deg.to_numpy()
    dec_arr = bh_dec_deg.to_numpy()

    # Evaluate BH position at the target date (UTC midnight) via interpolation
    t_eval = target_date.to_datetime64().astype("datetime64[s]").astype(np.int64)
    t_eval = float(t_eval)  # seconds since epoch

    bh_ra_eval = np.interp(t_eval, bh_tsec, ra_arr)
    bh_dec_eval = np.interp(t_eval, bh_tsec, dec_arr)

    print(f"BH @ {target_date.date()} 00:00Z: RA={bh_ra_eval:.6f} deg, Dec={bh_dec_eval:.6f} deg")

    # ---- Compute separations
    sub["bh_ra_deg"] = bh_ra_eval
    sub["bh_dec_deg"] = bh_dec_eval

    sub["bh_sep_deg"] = angular_sep_deg(sub["ra"].to_numpy(), sub["dec"].to_numpy(),
                                        bh_ra_eval, bh_dec_eval)
    sub["bh_sep_arcsec"] = sub["bh_sep_deg"] * 3600.0

    # Write a tidy output
    keep_cols = [
        "source_id", "ra", "dec", "phot_g_mean_mag",
        "first_detect_utc", "epochs_detected",
        "dist_to_track_deg", "max_shift_mas", "max_snr",
        "bh_ra_deg", "bh_dec_deg", "bh_sep_deg", "bh_sep_arcsec"
    ]
    keep_cols = [c for c in keep_cols if c in sub.columns]  # keep only those present

    out_df = sub[keep_cols].sort_values("bh_sep_deg").reset_index(drop=True)
    out_df.to_csv(args.out, index=False)
    print(f"Wrote: {args.out}")


if __name__ == "__main__":
    main()
