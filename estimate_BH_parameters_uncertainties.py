#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Estimate BH parameter uncertainties from perturbed simulation runs.

This script computes:
- sigma_A(t)     : uncertainty on tidal amplitude A ~ M/R^3
- sigma_lon(t)   : uncertainty on heliocentric longitude (deg)
- sigma_lat(t)   : uncertainty on heliocentric latitude (deg)

using a linearized least-squares / Fisher-matrix approach.

Required input files:
    fiducial.csv
    lon_plus.csv
    lon_minus.csv
    lat_plus.csv
    lat_minus.csv

Optional:
    mass_plus.csv
    mass_minus.csv

Expected columns in each CSV:
    time_days
    Delta_tidal (m)
    BH_Sun_d (AU)

Author: OpenAI
"""

import argparse
import math
import os
from dataclasses import dataclass
from typing import Dict, List, Optional, Tuple

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

# ------------------------------------------------------------
# Working directory for all inputs and outputs
# ------------------------------------------------------------
import os

BASE_DIR = os.path.join("SolarSystem", "uncertainties")

# Ensure folder exists
os.makedirs(BASE_DIR, exist_ok=True)

def fullpath(filename):
    return os.path.join(BASE_DIR, filename)

# ============================================================
# Utilities
# ============================================================

def deg2rad(x: float) -> float:
    return np.deg2rad(x)


def rad2deg(x: float) -> float:
    return np.rad2deg(x)


def safe_inv(mat: np.ndarray) -> np.ndarray:
    """Robust inverse using pseudo-inverse if needed."""
    try:
        return np.linalg.inv(mat)
    except np.linalg.LinAlgError:
        return np.linalg.pinv(mat)


def first_crossing(x: np.ndarray, threshold: float) -> Optional[int]:
    idx = np.where(np.isfinite(x) & (x >= threshold))[0]
    return None if len(idx) == 0 else int(idx[0])


def ensure_same_timebase(reference_df: pd.DataFrame, other_df: pd.DataFrame, name: str) -> None:
    a = reference_df["time_days"].to_numpy(dtype=float)
    b = other_df["time_days"].to_numpy(dtype=float)
    if len(a) != len(b) or not np.allclose(a, b, atol=1e-9, rtol=0):
        raise ValueError(f"Time base mismatch in file: {name}")


def load_csv(path: str) -> pd.DataFrame:
    df = pd.read_csv(path)
    required = ["time_days", "Delta_tidal (m)", "BH_Sun_d (AU)"]
    for col in required:
        if col not in df.columns:
            raise ValueError(f"{path} is missing required column: {col}")
    return df


@dataclass
class RunSet:
    fiducial: pd.DataFrame
    lon_plus: pd.DataFrame
    lon_minus: pd.DataFrame
    lat_plus: pd.DataFrame
    lat_minus: pd.DataFrame
    mass_plus: Optional[pd.DataFrame] = None
    mass_minus: Optional[pd.DataFrame] = None


# ============================================================
# Core math
# ============================================================

def central_difference(y_plus: np.ndarray, y_minus: np.ndarray, step: float) -> np.ndarray:
    return (y_plus - y_minus) / (2.0 * step)


def build_design_matrix(
    runs: RunSet,
    dlon_deg: float,
    dlat_deg: float,
    mass_fid_msun: float,
    dmass_msun: Optional[float] = None,
    use_mass_derivative_from_runs: bool = False,
) -> Tuple[np.ndarray, List[str], np.ndarray, np.ndarray]:
    """
    Build design matrix G(t, p) where columns are partial derivatives of the observable
    Delta_tidal(m) with respect to the chosen parameters.

    Returns
    -------
    G : (N, P) array
        Partial derivatives wrt parameters.
    param_names : list[str]
        Names of fitted parameters.
    y_fid : (N,) array
        Fiducial signal (treated as the "true" signal for uncertainty work).
    R_au : (N,) array
        BH distance from Sun.
    """

    y0 = runs.fiducial["Delta_tidal (m)"].to_numpy(dtype=float)
    R_au = runs.fiducial["BH_Sun_d (AU)"].to_numpy(dtype=float)

    y_lon_p = runs.lon_plus["Delta_tidal (m)"].to_numpy(dtype=float)
    y_lon_m = runs.lon_minus["Delta_tidal (m)"].to_numpy(dtype=float)
    y_lat_p = runs.lat_plus["Delta_tidal (m)"].to_numpy(dtype=float)
    y_lat_m = runs.lat_minus["Delta_tidal (m)"].to_numpy(dtype=float)

    # Derivatives wrt sky angles (units: meters per degree)
    d_y_d_lon = central_difference(y_lon_p, y_lon_m, dlon_deg)
    d_y_d_lat = central_difference(y_lat_p, y_lat_m, dlat_deg)

    param_names = []
    cols = []

    # Amplitude parameter A: use mass derivative if available, else use linear scaling
    if use_mass_derivative_from_runs and runs.mass_plus is not None and runs.mass_minus is not None:
        if dmass_msun is None:
            raise ValueError("dmass_msun must be provided if using mass derivative from runs.")
        y_mass_p = runs.mass_plus["Delta_tidal (m)"].to_numpy(dtype=float)
        y_mass_m = runs.mass_minus["Delta_tidal (m)"].to_numpy(dtype=float)
        d_y_d_mass = central_difference(y_mass_p, y_mass_m, dmass_msun)  # meters per solar mass
        cols.append(d_y_d_mass)
        param_names.append("mass_msun")
    else:
        # In the linear regime, signal ~ mass, so derivative wrt mass is y0 / mass
        d_y_d_mass = y0 / mass_fid_msun
        cols.append(d_y_d_mass)
        param_names.append("mass_msun")

    cols.append(d_y_d_lon)
    param_names.append("lon_deg")

    cols.append(d_y_d_lat)
    param_names.append("lat_deg")

    G = np.column_stack(cols)
    return G, param_names, y0, R_au


def fisher_covariance_over_time(
    G: np.ndarray,
    sigma_obs_m: float,
    min_points: int = 30,
    cadence_days: float = 1.0,
    corr_length_days: float = 1.0,
) -> Tuple[np.ndarray, np.ndarray]:
    """
    Compute covariance C(t) = (G^T W G)^-1 over increasing time windows.

    A simple way to account for time-correlated noise is to inflate sigma by
    sqrt(corr_length_days / cadence_days) if corr_length_days > cadence_days.
    """

    n, p = G.shape
    covariances = np.full((n, p, p), np.nan, dtype=float)
    sigmas = np.full((n, p), np.nan, dtype=float)

    corr_factor = max(1.0, math.sqrt(corr_length_days / cadence_days))
    sigma_eff = sigma_obs_m * corr_factor
    w = 1.0 / (sigma_eff ** 2)

    for i in range(min_points - 1, n):
        Gi = G[: i + 1, :]
        F = w * (Gi.T @ Gi)   # Fisher information matrix
        C = safe_inv(F)
        covariances[i] = C
        sigmas[i] = np.sqrt(np.diag(C))

    return covariances, sigmas


def snr_over_time(
    y_true: np.ndarray,
    sigma_obs_m: float,
    min_points: int = 30,
    cadence_days: float = 1.0,
    corr_length_days: float = 1.0,
) -> np.ndarray:
    """
    Compute matched-filter-like cumulative SNR for the fiducial signal itself.
    """
    n = len(y_true)
    out = np.full(n, np.nan, dtype=float)

    corr_factor = max(1.0, math.sqrt(corr_length_days / cadence_days))
    sigma_eff = sigma_obs_m * corr_factor

    for i in range(min_points - 1, n):
        yi = y_true[: i + 1]
        out[i] = np.sqrt(np.sum((yi / sigma_eff) ** 2))

    return out


def chi2_sky_map_quadratic(
    sigma_lon_deg: float,
    sigma_lat_deg: float,
    rho_lon_lat: float,
    lon_halfwidth_deg: float = 20.0,
    lat_halfwidth_deg: float = 10.0,
    nlon: int = 201,
    nlat: int = 201,
) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
    """
    Build an approximate chi2 map near the best fit from the 2x2 covariance
    of (lon, lat), using the local quadratic approximation:
        Δχ² = dθ^T C^-1 dθ
    """

    C = np.array([
        [sigma_lon_deg ** 2, rho_lon_lat * sigma_lon_deg * sigma_lat_deg],
        [rho_lon_lat * sigma_lon_deg * sigma_lat_deg, sigma_lat_deg ** 2],
    ])
    Cinv = safe_inv(C)

    lon = np.linspace(-lon_halfwidth_deg, lon_halfwidth_deg, nlon)
    lat = np.linspace(-lat_halfwidth_deg, lat_halfwidth_deg, nlat)
    LLON, LLAT = np.meshgrid(lon, lat)

    d = np.stack([LLON, LLAT], axis=-1)
    chi2 = np.einsum("...i,ij,...j->...", d, Cinv, d)

    return lon, lat, chi2


# ============================================================
# Plotting
# ============================================================

def plot_uncertainties(
    time_years: np.ndarray,
    R_au: np.ndarray,
    sigmas: np.ndarray,
    param_names: List[str],
    snr: np.ndarray,
    out_prefix: str,
) -> None:
    fig, axes = plt.subplots(4, 1, figsize=(10, 12), sharex=True)

    # 1) BH distance
    axes[0].plot(time_years, R_au, lw=1.5)
    axes[0].set_ylabel("BH distance (AU)")
    axes[0].set_title("BH distance and parameter uncertainties vs time")
    axes[0].grid(True, alpha=0.3)

    # 2) mass uncertainty
    if "mass_msun" in param_names:
        i = param_names.index("mass_msun")
        mjup_per_msun = 1.0 / 9.545e-4
        axes[1].plot(time_years, sigmas[:, i] * mjup_per_msun, lw=1.5)
        axes[1].set_ylabel(r"$\sigma(M)$ (M$_J$)")
        axes[1].set_yscale("log")
        axes[1].grid(True, alpha=0.3)

    # 3) sky uncertainties
    if "lon_deg" in param_names and "lat_deg" in param_names:
        ilon = param_names.index("lon_deg")
        ilat = param_names.index("lat_deg")
        axes[2].plot(time_years, sigmas[:, ilon], label=r"$\sigma(\lambda)$", lw=1.5)
        axes[2].plot(time_years, sigmas[:, ilat], label=r"$\sigma(\beta)$", lw=1.5)
        axes[2].set_ylabel("Angle uncertainty (deg)")
        axes[2].set_yscale("log")
        axes[2].legend()
        axes[2].grid(True, alpha=0.3)

    # 4) cumulative detection significance
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
    lon_grid: np.ndarray,
    lat_grid: np.ndarray,
    chi2: np.ndarray,
    out_prefix: str,
    epoch_year: float,
) -> None:
    fig, ax = plt.subplots(figsize=(8, 6))
    levels = [2.30, 6.17, 11.8]  # ~1σ, 2σ, 3σ for 2 params
    cf = ax.contourf(lon_grid, lat_grid, chi2, levels=50, cmap="viridis")
    c = ax.contour(lon_grid, lat_grid, chi2, levels=levels, colors="white", linewidths=1.5)
    ax.clabel(c, fmt={2.30: "1σ", 6.17: "2σ", 11.8: "3σ"})
    ax.set_xlabel(r"$\Delta \lambda$ (deg)")
    ax.set_ylabel(r"$\Delta \beta$ (deg)")
    ax.set_title(f"Approximate local sky uncertainty at t = {epoch_year:.2f} yr")
    fig.colorbar(cf, ax=ax, label=r"$\Delta \chi^2$")
    plt.tight_layout()
    fig.savefig(fullpath(f"{out_prefix}_sky_map.png"), dpi=160)
    plt.close(fig)


# ============================================================
# Main
# ============================================================

def main() -> None:
    parser = argparse.ArgumentParser(description="Estimate BH direction and mass uncertainties from perturbed runs.")
    parser.add_argument("--fiducial", required=True)
    parser.add_argument("--lon-plus", required=True)
    parser.add_argument("--lon-minus", required=True)
    parser.add_argument("--lat-plus", required=True)
    parser.add_argument("--lat-minus", required=True)
    parser.add_argument("--mass-plus", default=None)
    parser.add_argument("--mass-minus", default=None)

    parser.add_argument("--dlon-deg", type=float, required=True, help="Longitude perturbation in degrees used in the runs.")
    parser.add_argument("--dlat-deg", type=float, required=True, help="Latitude perturbation in degrees used in the runs.")
    parser.add_argument("--mass-fid-msun", type=float, required=True, help="Fiducial BH mass in solar masses.")
    parser.add_argument("--dmass-msun", type=float, default=None, help="Mass perturbation in solar masses, if mass± runs exist.")

    parser.add_argument("--sigma-obs-m", type=float, default=1.0, help="Assumed 1σ measurement noise in meters.")
    parser.add_argument("--cadence-days", type=float, default=1.0, help="Sampling cadence in days.")
    parser.add_argument("--corr-length-days", type=float, default=1.0, help="Effective correlation length in days.")
    parser.add_argument("--min-points", type=int, default=30, help="Minimum number of time samples before solving.")
    parser.add_argument("--epoch-year", type=float, default=None, help="Epoch (years since start) at which to draw sky map.")
    parser.add_argument("--out-prefix", default="bh_fit")
    args = parser.parse_args()

    # Load data
    fid = load_csv(fullpath(args.fiducial))
    lonp = load_csv(fullpath(args.lon_plus))
    lonm = load_csv(fullpath(args.lon_minus))
    latp = load_csv(fullpath(args.lat_plus))
    latm = load_csv(fullpath(args.lat_minus))

    for name, df in [
        ("lon_plus", lonp),
        ("lon_minus", lonm),
        ("lat_plus", latp),
        ("lat_minus", latm),
    ]:
        ensure_same_timebase(fid, df, name)

    massp = massm = None
    if args.mass_plus is not None and args.mass_minus is not None:
        massp = load_csv(fullpath(args.mass_plus))
        massm = load_csv(fullpath(args.mass_minus))        
        ensure_same_timebase(fid, massp, "mass_plus")
        ensure_same_timebase(fid, massm, "mass_minus")

    runs = RunSet(
        fiducial=fid,
        lon_plus=lonp,
        lon_minus=lonm,
        lat_plus=latp,
        lat_minus=latm,
        mass_plus=massp,
        mass_minus=massm,
    )

    # Build design matrix
    use_mass_runs = (massp is not None and massm is not None)
    G, param_names, y_fid, R_au = build_design_matrix(
        runs=runs,
        dlon_deg=args.dlon_deg,
        dlat_deg=args.dlat_deg,
        mass_fid_msun=args.mass_fid_msun,
        dmass_msun=args.dmass_msun,
        use_mass_derivative_from_runs=use_mass_runs,
    )

    # Covariance vs time
    covs, sigmas = fisher_covariance_over_time(
        G=G,
        sigma_obs_m=args.sigma_obs_m,
        min_points=args.min_points,
        cadence_days=args.cadence_days,
        corr_length_days=args.corr_length_days,
    )

    # Cumulative detection significance
    snr = snr_over_time(
        y_true=y_fid,
        sigma_obs_m=args.sigma_obs_m,
        min_points=args.min_points,
        cadence_days=args.cadence_days,
        corr_length_days=args.corr_length_days,
    )

    time_days = fid["time_days"].to_numpy(dtype=float)
    time_years = time_days / 365.25

    # Report useful thresholds
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

    if "mass_msun" in param_names:
        i = param_names.index("mass_msun")
        mjup_per_msun = 1.0 / 9.545e-4
        finite = np.where(np.isfinite(sigmas[:, i]))[0]
        if len(finite) > 0:
            j = finite[-1]
            print(f"Final sigma(M): {sigmas[j, i] * mjup_per_msun:.4f} Mj")

    if "lon_deg" in param_names and "lat_deg" in param_names:
        ilon = param_names.index("lon_deg")
        ilat = param_names.index("lat_deg")
        finite = np.where(np.isfinite(sigmas[:, ilon]) & np.isfinite(sigmas[:, ilat]))[0]
        if len(finite) > 0:
            j = finite[-1]
            print(f"Final sigma(lon): {sigmas[j, ilon]:.4f} deg")
            print(f"Final sigma(lat): {sigmas[j, ilat]:.4f} deg")

    # Save uncertainty plot
    plot_uncertainties(
        time_years=time_years,
        R_au=R_au,
        sigmas=sigmas,
        param_names=param_names,
        snr=snr,
        out_prefix=args.out_prefix,
    )

    # Sky map at chosen epoch
    if args.epoch_year is not None and "lon_deg" in param_names and "lat_deg" in param_names:
        idx = int(np.argmin(np.abs(time_years - args.epoch_year)))
        C = covs[idx]
        ilon = param_names.index("lon_deg")
        ilat = param_names.index("lat_deg")

        if np.all(np.isfinite(C)):
            Csky = C[np.ix_([ilon, ilat], [ilon, ilat])]
            sigma_lon = math.sqrt(Csky[0, 0])
            sigma_lat = math.sqrt(Csky[1, 1])
            rho = Csky[0, 1] / (sigma_lon * sigma_lat) if sigma_lon > 0 and sigma_lat > 0 else 0.0

            lon_grid, lat_grid, chi2 = chi2_sky_map_quadratic(
                sigma_lon_deg=sigma_lon,
                sigma_lat_deg=sigma_lat,
                rho_lon_lat=rho,
                lon_halfwidth_deg=max(5.0, 4.0 * sigma_lon),
                lat_halfwidth_deg=max(3.0, 4.0 * sigma_lat),
            )
            plot_sky_map(
                lon_grid=lon_grid,
                lat_grid=lat_grid,
                chi2=chi2,
                out_prefix=args.out_prefix,
                epoch_year=time_years[idx],
            )

    # Save sigma table
    out_df = pd.DataFrame({
        "time_days": time_days,
        "time_years": time_years,
        "BH_Sun_d_AU": R_au,
        "snr": snr,
    })

    for k, name in enumerate(param_names):
        out_df[f"sigma_{name}"] = sigmas[:, k]

    out_csv = fullpath(f"{args.out_prefix}_sigmas.csv")
    out_df.to_csv(out_csv, index=False)
    print(f"Saved: {out_csv}")
    print(f"Saved: {args.out_prefix}_uncertainties.png")
    if args.epoch_year is not None:
        print(f"Saved: {args.out_prefix}_sky_map.png")


if __name__ == "__main__":
    main()