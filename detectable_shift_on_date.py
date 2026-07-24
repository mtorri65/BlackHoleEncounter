#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import argparse
import math
from pathlib import Path

import numpy as np
import pandas as pd


# -----------------------------
# Constants
# -----------------------------
G = 6.67430e-11          # m^3 kg^-1 s^-2
c = 299792458.0          # m/s
M_sun = 1.98847e30       # kg
AU_m = 1.495978707e11    # m
RAD2ARCSEC = (180.0 / math.pi) * 3600.0


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
    """Fallback: find a BH track under ./simulations."""
    sims = base_dir / "simulations"
    if not sims.exists():
        raise FileNotFoundError(f"No simulations folder found at: {sims}")

    matches = list(sims.rglob("*bh_radec*.xlsx"))
    if not matches:
        raise FileNotFoundError("Could not find any *bh_radec*.xlsx under ./simulations/")

    # pick most recently modified
    matches.sort(key=lambda p: p.stat().st_mtime, reverse=True)
    return matches[0]


def angular_sep_deg(ra1_deg, dec1_deg, ra2_deg, dec2_deg):
    """
    Great-circle separation (deg) via haversine.
    Inputs may be numpy arrays.
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
    cang = 2.0 * np.arcsin(np.minimum(1.0, np.sqrt(a)))
    return np.rad2deg(cang)


def schwarzschild_radius_m(m_msun: float) -> float:
    """r_s = 2GM/c^2 in meters."""
    M = m_msun * M_sun
    return 2.0 * G * M / (c**2)


def einstein_radius_arcsec_simplified(m_msun: float, Dl_AU: float) -> float:
    """
    Simplified Einstein angle for Ds >> Dl:
      thetaE_rad = sqrt( (4GM/c^2) / Dl )
    Dl is lens distance (meters).
    Returns arcseconds.
    """
    M = m_msun * M_sun
    Dl_m = Dl_AU * AU_m
    thetaE_rad = math.sqrt((4.0 * G * M / (c**2)) / Dl_m)
    return thetaE_rad * RAD2ARCSEC


# -----------------------------
# Main
# -----------------------------
def main():
    ap = argparse.ArgumentParser(
        description="Compute BH-star separation and astrometric centroid shift on a single date."
    )
    ap.add_argument(
        "--pass-csv",
        default="gaia_real_detection_catalog_pass_ge5nights_EXCELSAFE.csv",
        help="Input pass catalog (EXCELSAFE ok). Must contain source_id, ra, dec."
    )
    ap.add_argument(
        "--date",
        default="2026-05-04",
        help="Date to evaluate (YYYY-MM-DD). Default 2026-05-04."
    )
    ap.add_argument(
        "--bh-radec",
        default=None,
        help="BH track xlsx filename (e.g. 20260121_114641__...__bh_radec__...xlsx). "
             "If omitted, auto-picks most recent *bh_radec*.xlsx under ./simulations."
    )
    ap.add_argument(
        "--bh-mass-msun",
        type=float,
        default=0.1,
        help="BH mass in solar masses. Default 0.1."
    )
    ap.add_argument(
        "--Dl-au",
        type=float,
        default=27.5,
        help="Lens distance Dl in AU for Einstein radius (simplified Ds>>Dl). Default 27.5 AU."
    )
    ap.add_argument(
        "--out",
        default="bh_shifts_on_2026-05-04_from_pass_catalog.csv",
        help="Output CSV filename."
    )
    args = ap.parse_args()

    base_dir = Path.cwd()
    eval_dt = pd.to_datetime(args.date, utc=True).normalize()
    t_eval_sec = float(eval_dt.to_datetime64().astype("datetime64[s]").astype(np.int64))

    # ---- (1) Read Gaia pass catalog and extract stars by source_id
    df = pd.read_csv(args.pass_csv, dtype=str)

    if "source_id" not in df.columns:
        raise KeyError(f"'source_id' not found in {args.pass_csv}. Columns: {df.columns.tolist()}")

    # EXCELSAFE may prefix apostrophe
    df["source_id"] = df["source_id"].astype(str).str.strip().str.lstrip("'")

    # We use RA/Dec already present in this file
    for col in ["ra", "dec"]:
        if col not in df.columns:
            raise KeyError(f"Missing column '{col}' in {args.pass_csv}. Columns: {df.columns.tolist()}")
        df[col] = pd.to_numeric(df[col], errors="coerce")

    df = df.dropna(subset=["ra", "dec"]).reset_index(drop=True)

    print(f"Loaded {len(df)} stars from {args.pass_csv}")

    # ---- (2) Determine BH RA/Dec on eval date from BH track
    if args.bh_radec:
        bh_path = find_bh_track_from_name(args.bh_radec, base_dir)
    else:
        bh_path = auto_find_any_bh_track(base_dir)

    print(f"Using BH track: {bh_path}")
    bh = pd.read_excel(bh_path)

    # Your BH columns: utc_iso, RA_hours, Dec_deg, range_AU, ...
    for col in ["utc_iso", "RA_hours", "Dec_deg"]:
        if col not in bh.columns:
            raise KeyError(f"BH track missing '{col}'. Columns: {bh.columns.tolist()}")

    bh["t"] = pd.to_datetime(bh["utc_iso"], utc=True, errors="coerce")
    bh = bh.dropna(subset=["t"]).sort_values("t").reset_index(drop=True)

    bh_tsec = bh["t"].astype("int64").to_numpy() / 1e9
    bh_ra_deg_arr = pd.to_numeric(bh["RA_hours"], errors="coerce").to_numpy() * 15.0
    bh_dec_deg_arr = pd.to_numeric(bh["Dec_deg"], errors="coerce").to_numpy()

    ra_bh = float(np.interp(t_eval_sec, bh_tsec, bh_ra_deg_arr))
    dec_bh = float(np.interp(t_eval_sec, bh_tsec, bh_dec_deg_arr))

    print(f"BH @ {eval_dt.date()} 00:00Z: RA={ra_bh:.6f} deg, Dec={dec_bh:.6f} deg")

    # ---- (3) Angular separation β between BH and each star on that date
    beta_deg = angular_sep_deg(df["ra"].to_numpy(), df["dec"].to_numpy(), ra_bh, dec_bh)
    beta_arcsec = beta_deg * 3600.0

    # ---- (4) Schwarzschild radius r_s for BH
    rs_m = schwarzschild_radius_m(args.bh_mass_msun)
    rs_km = rs_m / 1000.0
    rs_AU = rs_m / AU_m

    # ---- (5) Einstein radius theta_E with Ds >> Dl (Dl fixed to 27.5 AU by default)
    thetaE_arcsec = einstein_radius_arcsec_simplified(args.bh_mass_msun, args.Dl_au)
    thetaE_mas = thetaE_arcsec * 1000.0

    # ---- (6) Centroid shift δθ on that date using δθ = u/(u^2+2) θE, u = β/θE
    u = beta_arcsec / thetaE_arcsec
    delta_arcsec = (u / (u*u + 2.0)) * thetaE_arcsec
    delta_mas = delta_arcsec * 1000.0

    out = df.copy()
    out["eval_date_utc"] = str(eval_dt.date())
    out["bh_ra_deg"] = ra_bh
    out["bh_dec_deg"] = dec_bh
    out["beta_deg"] = beta_deg
    out["beta_arcsec"] = beta_arcsec
    out["bh_mass_msun"] = args.bh_mass_msun
    out["Dl_AU_for_thetaE"] = args.Dl_au
    out["schwarzschild_radius_km"] = rs_km
    out["schwarzschild_radius_AU"] = rs_AU
    out["thetaE_arcsec"] = thetaE_arcsec
    out["thetaE_mas"] = thetaE_mas
    out["centroid_shift_mas_on_date"] = delta_mas
    print(out.columns.tolist())
    row = out[out["source_id"] == "6876597890528117376"]
    print(f"beta_arcsec= {row['beta_arcsec']}")


    # Helpful ordering
    out = out.sort_values("centroid_shift_mas_on_date", ascending=False).reset_index(drop=True)

    # Ensure source_id is Excel-safe (force text)
    out["source_id"] = "'" + out["source_id"].astype(str)
    out.to_csv(args.out, index=False)
    print(f"Wrote: {args.out}")
    print(f"Schwarzschild radius: {rs_km:.6f} km  ({rs_AU:.6e} AU)")
    print(f"Einstein radius (Ds>>Dl): {thetaE_arcsec:.6f} arcsec  ({thetaE_mas:.3f} mas) at Dl={args.Dl_au} AU")


if __name__ == "__main__":
    main()
