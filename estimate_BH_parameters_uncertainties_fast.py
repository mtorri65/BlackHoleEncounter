#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import argparse
import math
import os
from dataclasses import dataclass
from typing import Optional

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

# ============================================================
# Directory handling
# Script runs from SolarSystem/
# Data are stored in SolarSystem/uncertainties/
# ============================================================

import os

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
BASE_DIR = os.path.join(SCRIPT_DIR, "uncertainties")
os.makedirs(BASE_DIR, exist_ok=True)

def fullpath(filename: str) -> str:
    return os.path.join(BASE_DIR, filename)

# ============================================================
# Data containers
# ============================================================

@dataclass
class Arrays:
    time_days: np.ndarray
    delta_m: np.ndarray
    bh_sun_d_au: np.ndarray


@dataclass
class RunSet:
    fid: Arrays
    lon_p: Arrays
    lon_m: Arrays
    lat_p: Arrays
    lat_m: Arrays
    mass_p: Optional[Arrays] = None
    mass_m: Optional[Arrays] = None


# ============================================================
# Utilities
# ============================================================

def safe_inv(a: np.ndarray) -> np.ndarray:
    try:
        return np.linalg.inv(a)
    except np.linalg.LinAlgError:
        return np.linalg.pinv(a)


def first_crossing(x: np.ndarray, threshold: float) -> Optional[int]:
    idx = np.where(np.isfinite(x) & (x >= threshold))[0]
    return int(idx[0]) if len(idx) else None


def load_arrays(filename: str) -> Arrays:
    path = fullpath(filename)
    df = pd.read_csv(path)

    required = ["time_days", "Delta_tidal (m)", "BH_Sun_d (AU)"]
    missing = [c for c in required if c not in df.columns]
    if missing:
        raise ValueError(f"{path} missing columns: {missing}")

    return Arrays(
        time_days=df["time_days"].to_numpy(dtype=float),
        delta_m=df["Delta_tidal (m)"].to_numpy(dtype=float),
        bh_sun_d_au=df["BH_Sun_d (AU)"].to_numpy(dtype=float),
    )


def check_same_timebase(ref: Arrays, arr: Arrays, name: str) -> None:
    if len(ref.time_days) != len(arr.time_days):
        raise ValueError(f"Time length mismatch for {name}")
    if not np.allclose(ref.time_days, arr.time_days, atol=1e-12, rtol=0.0):
        raise ValueError(f"Time grid mismatch for {name}")


def central_diff(y_plus: np.ndarray, y_minus: np.ndarray, step: float) -> np.ndarray:
    return (y_plus - y_minus) / (2.0 * step)


# ============================================================
# Core calculations
# ============================================================

def build_design_matrix(
    runs: RunSet,
    mass_fid_msun: float,
    dlon_deg: float,
    dlat_deg: float,
    dmass_msun: Optional[float] = None,
    use_mass_runs: bool = False,
):
    y0 = runs.fid.delta_m
    R_au = runs.fid.bh_sun_d_au

    # Derivatives wrt longitude and latitude
    d_y_d_lon = central_diff(runs.lon_p.delta_m, runs.lon_m.delta_m, dlon_deg)
    d_y_d_lat = central_diff(runs.lat_p.delta_m, runs.lat_m.delta_m, dlat_deg)

    # Derivative wrt mass
    if use_mass_runs:
        if runs.mass_p is None or runs.mass_m is None:
            raise ValueError("Mass derivative requested but mass± runs not provided.")
        if dmass_msun is None:
            raise ValueError("dmass_msun must be provided when using mass± runs.")
        d_y_d_mass = central_diff(runs.mass_p.delta_m, runs.mass_m.delta_m, dmass_msun)
    else:
        # Linear regime approximation
        d_y_d_mass = y0 / mass_fid_msun

    # Design matrix G: columns correspond to [mass, lon, lat]
    G = np.column_stack([d_y_d_mass, d_y_d_lon, d_y_d_lat])
    param_names = ["mass_msun", "lon_deg", "lat_deg"]

    return G, param_names, y0, R_au


def cumulative_fisher(
    G: np.ndarray,
    sigma_obs_m: float,
    min_points: int = 30,
    cadence_days: float = 1.0,
    corr_length_days: float = 1.0,
):
    """
    Fast cumulative covariance using prefix sums.

    For each time index i:
        F_i = sum_{k<=i} (g_k g_k^T) / sigma_eff^2
    """
    n, p = G.shape

    corr_factor = max(1.0, math.sqrt(corr_length_days / cadence_days))
    sigma_eff = sigma_obs_m * corr_factor
    inv_var = 1.0 / (sigma_eff * sigma_eff)

    # Prefix sums of outer products
    outer_terms = np.einsum("ni,nj->nij", G, G) * inv_var
    fisher_prefix = np.cumsum(outer_terms, axis=0)

    covs = np.full((n, p, p), np.nan, dtype=float)
    sigmas = np.full((n, p), np.nan, dtype=float)

    for i in range(min_points - 1, n):
        F = fisher_prefix[i]
        C = safe_inv(F)
        covs[i] = C
        sigmas[i] = np.sqrt(np.diag(C))

    return covs, sigmas


def cumulative_snr(
    y_true: np.ndarray,
    sigma_obs_m: float,
    min_points: int = 30,
    cadence_days: float = 1.0,
    corr_length_days: float = 1.0,
):
    corr_factor = max(1.0, math.sqrt(corr_length_days / cadence_days))
    sigma_eff = sigma_obs_m * corr_factor

    contrib = (y_true / sigma_eff) ** 2
    snr2 = np.cumsum(contrib)
    snr = np.sqrt(snr2)
    snr[:max(0, min_points - 1)] = np.nan
    return snr


def build_local_chi2_map(
    cov_lon_lat: np.ndarray,
    nsig: float = 4.0,
    ngrid: int = 201,
):
    sigma_lon = math.sqrt(cov_lon_lat[0, 0])
    sigma_lat = math.sqrt(cov_lon_lat[1, 1])

    lon_half = max(5.0, nsig * sigma_lon)
    lat_half = max(3.0, nsig * sigma_lat)

    lon = np.linspace(-lon_half, lon_half, ngrid)
    lat = np.linspace(-lat_half, lat_half, ngrid)
    LLON, LLAT = np.meshgrid(lon, lat)

    Cinv = safe_inv(cov_lon_lat)

    d = np.stack([LLON, LLAT], axis=-1)
    chi2 = np.einsum("...i,ij,...j->...", d, Cinv, d)

    return lon, lat, chi2


# ============================================================
# Plotting
# ============================================================

def plot_summary(
    time_years: np.ndarray,
    R_au: np.ndarray,
    sigmas: np.ndarray,
    snr: np.ndarray,
    out_prefix: str,
):
    fig, axes = plt.subplots(4, 1, figsize=(10, 12), sharex=True)

    axes[0].plot(time_years, R_au, lw=1.5)
    axes[0].set_ylabel("BH distance (AU)")
    axes[0].set_title("BH detection / localization uncertainties")
    axes[0].grid(True, alpha=0.3)

    sigma_mass_mj = sigmas[:, 0] / 9.545e-4
    axes[1].plot(time_years, sigma_mass_mj, lw=1.5)
    axes[1].set_ylabel(r"$\sigma(M)$ (M$_J$)")
    axes[1].set_yscale("log")
    axes[1].grid(True, alpha=0.3)

    axes[2].plot(time_years, sigmas[:, 1], lw=1.5, label=r"$\sigma(\lambda)$")
    axes[2].plot(time_years, sigmas[:, 2], lw=1.5, label=r"$\sigma(\beta)$")
    axes[2].set_ylabel("Angle uncertainty (deg)")
    axes[2].set_yscale("log")
    axes[2].legend()
    axes[2].grid(True, alpha=0.3)

    axes[3].plot(time_years, snr, lw=1.5)
    for thr in [3, 5, 10]:
        axes[3].axhline(thr, ls="--", lw=1)
    axes[3].set_ylabel("Cumulative SNR")
    axes[3].set_xlabel("Time since start of simulation (years)")
    axes[3].grid(True, alpha=0.3)

    plt.tight_layout()
    fig.savefig(fullpath(f"{out_prefix}_uncertainties.png"), dpi=160)
    plt.close(fig)


def plot_sky_map(
    lon: np.ndarray,
    lat: np.ndarray,
    chi2: np.ndarray,
    epoch_year: float,
    out_prefix: str,
):
    fig, ax = plt.subplots(figsize=(8, 6))
    levels = [2.30, 6.17, 11.8]  # 1σ, 2σ, 3σ for 2 parameters
    cf = ax.contourf(lon, lat, chi2, levels=60, cmap="viridis")
    cc = ax.contour(lon, lat, chi2, levels=levels, colors="white", linewidths=1.2)
    ax.clabel(cc, fmt={2.30: "1σ", 6.17: "2σ", 11.8: "3σ"})
    ax.set_xlabel(r"$\Delta \lambda$ (deg)")
    ax.set_ylabel(r"$\Delta \beta$ (deg)")
    ax.set_title(f"Local sky uncertainty at t = {epoch_year:.2f} yr")
    fig.colorbar(cf, ax=ax, label=r"$\Delta \chi^2$")
    plt.tight_layout()
    fig.savefig(fullpath(f"{out_prefix}_sky_map.png"), dpi=160)
    plt.close(fig)


# ============================================================
# Main
# ============================================================

def main():
    parser = argparse.ArgumentParser(description="Fast BH uncertainty estimator.")
    parser.add_argument("--fiducial", required=True)
    parser.add_argument("--lon-plus", required=True)
    parser.add_argument("--lon-minus", required=True)
    parser.add_argument("--lat-plus", required=True)
    parser.add_argument("--lat-minus", required=True)
    parser.add_argument("--mass-plus", default=None)
    parser.add_argument("--mass-minus", default=None)

    parser.add_argument("--dlon-deg", type=float, required=True)
    parser.add_argument("--dlat-deg", type=float, required=True)
    parser.add_argument("--mass-fid-msun", type=float, required=True)
    parser.add_argument("--dmass-msun", type=float, default=None)

    parser.add_argument("--sigma-obs-m", type=float, default=1.0)
    parser.add_argument("--cadence-days", type=float, default=1.0)
    parser.add_argument("--corr-length-days", type=float, default=7.0)
    parser.add_argument("--min-points", type=int, default=60)
    parser.add_argument("--epoch-year", type=float, default=None)
    parser.add_argument("--out-prefix", default="bh_case")

    args = parser.parse_args()

    # Load arrays once
    fid = load_arrays(args.fiducial)
    lon_p = load_arrays(args.lon_plus)
    lon_m = load_arrays(args.lon_minus)
    lat_p = load_arrays(args.lat_plus)
    lat_m = load_arrays(args.lat_minus)

    check_same_timebase(fid, lon_p, "lon_plus")
    check_same_timebase(fid, lon_m, "lon_minus")
    check_same_timebase(fid, lat_p, "lat_plus")
    check_same_timebase(fid, lat_m, "lat_minus")

    mass_p = mass_m = None
    use_mass_runs = args.mass_plus is not None and args.mass_minus is not None
    if use_mass_runs:
        mass_p = load_arrays(args.mass_plus)
        mass_m = load_arrays(args.mass_minus)
        check_same_timebase(fid, mass_p, "mass_plus")
        check_same_timebase(fid, mass_m, "mass_minus")

    runs = RunSet(
        fid=fid,
        lon_p=lon_p,
        lon_m=lon_m,
        lat_p=lat_p,
        lat_m=lat_m,
        mass_p=mass_p,
        mass_m=mass_m,
    )

    G, param_names, y_true, R_au = build_design_matrix(
        runs=runs,
        mass_fid_msun=args.mass_fid_msun,
        dlon_deg=args.dlon_deg,
        dlat_deg=args.dlat_deg,
        dmass_msun=args.dmass_msun,
        use_mass_runs=use_mass_runs,
    )

    covs, sigmas = cumulative_fisher(
        G=G,
        sigma_obs_m=args.sigma_obs_m,
        min_points=args.min_points,
        cadence_days=args.cadence_days,
        corr_length_days=args.corr_length_days,
    )

    snr = cumulative_snr(
        y_true=y_true,
        sigma_obs_m=args.sigma_obs_m,
        min_points=args.min_points,
        cadence_days=args.cadence_days,
        corr_length_days=args.corr_length_days,
    )

    time_days = fid.time_days
    time_years = time_days / 365.25

    # Report thresholds
    for thr in [3, 5, 10]:
        idx = first_crossing(snr, thr)
        if idx is None:
            print(f"SNR >= {thr}: not reached")
        else:
            print(
                f"SNR >= {thr}: "
                f"t = {time_years[idx]:.2f} yr, "
                f"BH distance = {R_au[idx]:.1f} AU"
            )

    # Final uncertainties
    finite = np.where(np.isfinite(sigmas[:, 0]))[0]
    if len(finite):
        i = finite[-1]
        print(f"Final sigma(M) = {sigmas[i, 0] / 9.545e-4:.4f} Mj")
        print(f"Final sigma(lon) = {sigmas[i, 1]:.4f} deg")
        print(f"Final sigma(lat) = {sigmas[i, 2]:.4f} deg")

    # Save CSV
    out_df = pd.DataFrame({
        "time_days": time_days,
        "time_years": time_years,
        "BH_Sun_d_AU": R_au,
        "snr": snr,
        "sigma_mass_msun": sigmas[:, 0],
        "sigma_mass_mj": sigmas[:, 0] / 9.545e-4,
        "sigma_lon_deg": sigmas[:, 1],
        "sigma_lat_deg": sigmas[:, 2],
    })
    out_csv = fullpath(f"{args.out_prefix}_sigmas.csv")
    out_df.to_csv(out_csv, index=False)

    # Plots
    plot_summary(
        time_years=time_years,
        R_au=R_au,
        sigmas=sigmas,
        snr=snr,
        out_prefix=args.out_prefix,
    )

    if args.epoch_year is not None:
        idx = int(np.argmin(np.abs(time_years - args.epoch_year)))
        if np.all(np.isfinite(covs[idx])):
            cov_lon_lat = covs[idx][1:3, 1:3]
            lon, lat, chi2 = build_local_chi2_map(cov_lon_lat)
            plot_sky_map(
                lon=lon,
                lat=lat,
                chi2=chi2,
                epoch_year=time_years[idx],
                out_prefix=args.out_prefix,
            )

    print(f"Saved: {out_csv}")
    print(f"Saved: {fullpath(args.out_prefix + '_uncertainties.png')}")
    if args.epoch_year is not None:
        print(f"Saved: {fullpath(args.out_prefix + '_sky_map.png')}")


if __name__ == "__main__":
    main()