#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
post_orbit_compare_geocentric_fixed2.py

Creates orbit-compare plots from a single run folder that contains a prefixed
__orbits__*.xlsx workbook produced by solar_system_bh_rebound21.py.

For each planet (+ Moon):
  - Unperturbed baseline (no BH) heliocentric orbit for 1 period at epoch
  - Perturbed "last orbit" near --tend (default 36500 days), using period
    inferred from the end-state (requires vx/vy/vz columns)

Extra combined plots:
  - orbitcompare__Earth_Moon.png  (Earth heliocentric + Moon heliocentric-with-wiggles)
  - orbitcompare__ALL_epoch.png   (all bodies baseline @ epoch, one orbit each)
  - orbitcompare__ALL_final.png   (all bodies perturbed last-orbit window near --tend)

Extra geocentric plot:
  - orbitcompare__Moon_geocentric.png:
        baseline Moon orbit around Earth (1 lunar period) vs
        perturbed last lunar orbit around Earth near --tend.
    Period is computed from Moon-Earth *relative* state using mu=G*(M_earth+M_moon).
    If that fails, falls back to 27.321661 days (sidereal month).

Notes:
  - Requires the run folder to contain <prefix>__input.yaml so we can read epoch/kernel.
  - Requires vx_AU_per_d, vy_AU_per_d, vz_AU_per_d columns in the workbook sheets to
    infer end-state orbital periods robustly. If missing, we can add a fallback table.
"""

import argparse
import math
from pathlib import Path

import numpy as np
import pandas as pd
import yaml
import matplotlib.pyplot as plt

import rebound
from skyfield.api import load
from skyfield.constants import AU_KM

# -----------------------
# Constants / body order
# -----------------------
G_AU3_Msun_day2 = 0.0002959122082855911
MSUN = 1.0

PLANET_MASSES_MSUN = {
    "mercury": 1.651e-7, "venus": 2.447e-6, "earth": 3.003e-6, "mars": 3.227e-7,
    "jupiter": 9.5479e-4, "saturn": 2.8585e-4, "uranus": 4.3662e-5, "neptune": 5.1513e-5,
}
MOON_MASS_MSUN = 3.694e-8

# de440s keys
SKYFIELD_KEYS = {
    "sun": "sun",
    "mercury": "mercury",
    "venus": "venus",
    "earth": "earth",
    "mars": "mars barycenter",
    "jupiter": "jupiter barycenter",
    "saturn": "saturn barycenter",
    "uranus": "uranus barycenter",
    "neptune": "neptune barycenter",
}

PLANET_NAMES = ["Mercury","Venus","Earth","Mars","Jupiter","Saturn","Uranus","Neptune","Moon"]

# Keep consistent colors across plots
COLOR_MAP = {
    "Mercury": "tab:gray",
    "Venus": "tab:orange",
    "Earth": "tab:blue",
    "Mars": "tab:red",
    "Jupiter": "tab:brown",
    "Saturn": "tab:olive",
    "Uranus": "tab:cyan",
    "Neptune": "tab:purple",
    "Moon": "tab:green",
}

# -----------------------
# Small helpers
# -----------------------
def safe_num_for_filename(x: float) -> str:
    s = f"{x:g}"
    return s.replace("-", "m").replace(".", "p")

def _prefix_from_run_dir(run_dir: Path) -> str:
    # run_dir name looks like "<STAMP>__rp...__vinf...__inc...__toff...__Om...__om..."
    return run_dir.name

def _prefix_path(run_dir: Path, name: str) -> Path:
    return run_dir / f"{_prefix_from_run_dir(run_dir)}__{name}"

def _find_orbits_xlsx(run_dir: Path) -> Path:
    cands = sorted(run_dir.glob(f"{run_dir.name}__orbits__*.xlsx"))
    if not cands:
        # fallback: any orbits__*.xlsx
        cands = sorted(run_dir.glob("*__orbits__*.xlsx"))
    if not cands:
        raise FileNotFoundError(f"No __orbits__*.xlsx found in {run_dir}")
    return cands[0]

def _find_prefixed_input_yaml(run_dir: Path) -> Path | None:
    cand = run_dir / f"{run_dir.name}__input.yaml"
    if cand.exists():
        return cand
    # fallback: any *input.yaml
    cands = list(run_dir.glob("*input.yaml"))
    return cands[0] if cands else None

def _read_sheet(xls: pd.ExcelFile, name: str) -> pd.DataFrame:
    df = pd.read_excel(xls, sheet_name=name)
    # normalize column names if needed
    df.columns = [str(c).strip() for c in df.columns]
    return df

def _add_params_box(ax, params: dict, loc=(0.02, 0.98), max_lines=12):
    """
    Draw a small parameter text box inside the axes.

    Parameters
    ----------
    ax : matplotlib Axes
    params : dict
        Key-value pairs to show.
    loc : (x, y)
        Upper-left corner in axes fraction coordinates.
    max_lines : int
        Maximum number of lines to display.
    """
    lines = []
    for k, v in params.items():
        if v is None or (isinstance(v, float) and not np.isfinite(v)):
            continue
        lines.append(f"{k}: {v}")
        if len(lines) >= max_lines:
            break

    if not lines:
        return

    text = "\n".join(lines)

    ax.text(
        loc[0], loc[1],
        text,
        transform=ax.transAxes,
        va="top",
        ha="left",
        fontsize=9,
        bbox=dict(
            boxstyle="round",
            facecolor="white",
            edgecolor="gray",
            alpha=0.85
        )
    )

# -----------------------
# Orbit math helpers
# -----------------------
def _elem_from_state(x, y, z, vx, vy, vz, mu):
    r = math.sqrt(x*x + y*y + z*z)
    v2 = vx*vx + vy*vy + vz*vz

    hx = y*vz - z*vy
    hy = z*vx - x*vz
    hz = x*vy - y*vx
    h2 = hx*hx + hy*hy + hz*hz

    inv_a = (2.0 / max(r, 1e-30)) - (v2 / max(mu, 1e-30))
    a = (1.0 / inv_a) if abs(inv_a) > 1e-30 else math.copysign(1e30, inv_a)

    # eccentricity vector
    cx = vy*hz - vz*hy
    cy = vz*hx - vx*hz
    cz = vx*hy - vy*hx
    ex = cx/max(mu,1e-30) - x/max(r,1e-30)
    ey = cy/max(mu,1e-30) - y/max(r,1e-30)
    ez = cz/max(mu,1e-30) - z/max(r,1e-30)
    e = math.sqrt(ex*ex + ey*ey + ez*ez)

    return a, e, h2

def _kepler_period_days(a, mu):
    if not np.isfinite(a) or a <= 0:
        return float("nan")
    return 2.0*math.pi*math.sqrt(a**3 / mu)

def _infer_period_from_state_heliocentric(x, y, z, vx, vy, vz, mu_sun):
    a, e, _ = _elem_from_state(x,y,z,vx,vy,vz,mu_sun)
    return _kepler_period_days(a, mu_sun), a, e

def _infer_period_from_state_relative(dx, dy, dz, dvx, dvy, dvz, mu_primary):
    a, e, _ = _elem_from_state(dx,dy,dz,dvx,dvy,dvz,mu_primary)
    return _kepler_period_days(a, mu_primary), a, e

def _add_params_box(ax, params: dict, loc=(0.02, 0.98), max_lines=12):
    """
    Draw a small parameter text box inside the axes.

    Parameters
    ----------
    ax : matplotlib Axes
    params : dict
        Key-value pairs to show.
    loc : (x, y)
        Upper-left corner in axes fraction coordinates.
    max_lines : int
        Maximum number of lines to display.
    """
    lines = []
    for k, v in params.items():
        if v is None or (isinstance(v, float) and not np.isfinite(v)):
            continue
        lines.append(f"{k}: {v}")
        if len(lines) >= max_lines:
            break

    if not lines:
        return

    text = "\n".join(lines)

    ax.text(
        loc[0], loc[1],
        text,
        transform=ax.transAxes,
        va="top",
        ha="left",
        fontsize=9,
        bbox=dict(
            boxstyle="round",
            facecolor="white",
            edgecolor="gray",
            alpha=0.85
        )
    )

# -----------------------
# Baseline (no BH) sim
# -----------------------
def build_sim_no_bh(params: dict) -> rebound.Simulation:
    epoch_iso = params["epoch"]
    kernel = params["kernel"]

    eph = load(f"{kernel}.bsp" if not str(kernel).lower().endswith(".bsp") else str(kernel))
    ts = load.timescale()
    # parse UTC ISO
    from datetime import datetime, timezone
    s = epoch_iso.strip()
    if s.endswith("Z"):
        dt = datetime.fromisoformat(s[:-1]).replace(tzinfo=timezone.utc)
    else:
        dt = datetime.fromisoformat(s)
        if dt.tzinfo is None:
            dt = dt.replace(tzinfo=timezone.utc)
        dt = dt.astimezone(timezone.utc)
    t0 = ts.from_datetime(dt)

    sun_bary = eph[SKYFIELD_KEYS["sun"]].at(t0)
    sun_r, sun_v = sun_bary.position.au, sun_bary.velocity.au_per_d

    def heliocentric_state(key):
        b = eph[key].at(t0)
        return b.position.au - sun_r, b.velocity.au_per_d - sun_v

    sim = rebound.Simulation()
    sim.G = G_AU3_Msun_day2
    sim.integrator = "whfast"
    sim.dt = float(params.get("dt_days", 0.5))

    sim.add(m=MSUN)

    order = ["mercury","venus","earth","mars","jupiter","saturn","uranus","neptune"]
    for nm in order:
        r, v = heliocentric_state(SKYFIELD_KEYS[nm])
        sim.add(m=PLANET_MASSES_MSUN[nm], x=r[0], y=r[1], z=r[2], vx=v[0], vy=v[1], vz=v[2])

    # Moon around Earth — match main script behavior (simple circular-ish)
    earth_idx = order.index("earth") + 1
    earth = sim.particles[earth_idx]
    moon_r_au = 384400.0 / AU_KM
    # NOTE: this is whatever your main script uses; if you changed it, keep in sync.
    sim.add(m=MOON_MASS_MSUN,
            x=earth.x + moon_r_au, y=earth.y, z=earth.z,
            vx=earth.vx, vy=earth.vy + 0.00257, vz=earth.vz)

    sim.move_to_com()
    return sim

def integrate_baseline_orbit_xy(sim: rebound.Simulation, body_idx: int, T_days: float, step_days: float):
    if not np.isfinite(T_days) or T_days <= 0:
        return None
    n = int(max(10, math.ceil(T_days / step_days))) + 1
    ts = np.linspace(0.0, T_days, n)
    xs = np.zeros(n); ys = np.zeros(n)
    for i, t in enumerate(ts):
        sim.integrate(t)
        p = sim.particles[body_idx]
        s = sim.particles[0]
        xs[i] = p.x - s.x
        ys[i] = p.y - s.y
    return xs, ys

def integrate_baseline_relative_orbit_xy(sim: rebound.Simulation, sat_idx: int, prim_idx: int,
                                        T_days: float, step_days: float):
    if not np.isfinite(T_days) or T_days <= 0:
        return None
    n = int(max(10, math.ceil(T_days / step_days))) + 1
    ts = np.linspace(0.0, T_days, n)
    xs = np.zeros(n); ys = np.zeros(n)
    for i, t in enumerate(ts):
        sim.integrate(t)
        sat = sim.particles[sat_idx]
        prim = sim.particles[prim_idx]
        xs[i] = sat.x - prim.x
        ys[i] = sat.y - prim.y
    return xs, ys

# -----------------------
# Plotters
# -----------------------
def plot_one_body(run_dir: Path, body: str, xb, yb, xf, yf, t_end, T_future, a_end, e_end):
    fig, ax = plt.subplots(figsize=(8,8))
    ax.set_aspect("equal")
    ax.grid(True, alpha=0.3)
    ax.set_xlabel("x [AU] (heliocentric)")
    ax.set_ylabel("y [AU] (heliocentric)")
    ax.set_title(
        f"{body}: Unperturbed @ epoch vs last orbit near t={int(round(t_end))} d\n"
        f"T_future≈{T_future:.2f} d (a_end≈{a_end:.3f} AU, e_end≈{e_end:.3f})"
    )
    c = COLOR_MAP.get(body, None)
    ax.plot(xb, yb, lw=1.6, label="Unperturbed @ epoch (1 orbit)", color=c, alpha=0.65)
    ax.plot(xf, yf, lw=1.6, label="Perturbed near end (last orbit)", color=c)
    ax.scatter([0],[0], s=70, c="gold", edgecolor="k", linewidth=0.4, label="Sun")
    ax.legend(fontsize=9)
    out = _prefix_path(run_dir, f"orbit_compare_{body}.png")
    plt.tight_layout()
    fig.savefig(out, dpi=200)
    plt.close(fig)
    print(f"Wrote {out.name}")

def plot_earth_moon_combo(run_dir: Path, earth_base, moon_base, earth_final, moon_final, t_end):
    fig, ax = plt.subplots(figsize=(9,9))
    ax.set_aspect("equal")
    ax.grid(True, alpha=0.3)
    ax.set_xlabel("x [AU] (heliocentric)")
    ax.set_ylabel("y [AU] (heliocentric)")
    ax.set_title(f"Earth + Moon: baseline @ epoch vs perturbed near t={int(round(t_end))} d")

    # baseline
    exb, eyb = earth_base
    mxb, myb = moon_base
    ax.plot(exb, eyb, lw=1.6, color=COLOR_MAP["Earth"], alpha=0.6, label="Earth baseline")
    ax.plot(mxb, myb, lw=1.2, color=COLOR_MAP["Moon"], alpha=0.6, label="Moon baseline (heliocentric w/ wiggles)")

    # final
    exf, eyf = earth_final
    mxf, myf = moon_final
    ax.plot(exf, eyf, lw=1.8, color=COLOR_MAP["Earth"], label="Earth final (last orbit)")
    ax.plot(mxf, myf, lw=1.4, color=COLOR_MAP["Moon"], label="Moon final (last orbit)")

    ax.scatter([0],[0], s=80, c="gold", edgecolor="k", linewidth=0.4, label="Sun")
    ax.legend(fontsize=9)
    out = _prefix_path(run_dir, "orbitcompare__Earth_Moon.png")
    plt.tight_layout()
    fig.savefig(out, dpi=220)
    plt.close(fig)
    print(f"Wrote {out.name}")

def plot_all_epoch(run_dir: Path, baseline_orbits: dict):
    fig, ax = plt.subplots(figsize=(9,9))
    ax.set_aspect("equal")
    ax.grid(True, alpha=0.3)
    ax.set_xlabel("x [AU] (heliocentric)")
    ax.set_ylabel("y [AU] (heliocentric)")
    ax.set_title("All bodies: baseline (unperturbed) 1 orbit @ epoch")

    for body, (x,y) in baseline_orbits.items():
        ax.plot(x, y, lw=1.1, color=COLOR_MAP.get(body, None), label=body)
    ax.scatter([0],[0], s=80, c="gold", edgecolor="k", linewidth=0.4, label="Sun")
    ax.legend(fontsize=8, ncols=2)
    out = _prefix_path(run_dir, "orbitcompare__ALL_epoch.png")
    plt.tight_layout()
    fig.savefig(out, dpi=220)
    plt.close(fig)
    print(f"Wrote {out.name}")

def plot_all_final(run_dir: Path, final_orbits: dict, t_end):
    fig, ax = plt.subplots(figsize=(9,9))
    ax.set_aspect("equal")
    ax.grid(True, alpha=0.3)
    ax.set_xlabel("x [AU] (heliocentric)")
    ax.set_ylabel("y [AU] (heliocentric)")
    ax.set_title(f"All bodies: perturbed last-orbit window near t={int(round(t_end))} d")

    for body, (x,y) in final_orbits.items():
        ax.plot(x, y, lw=1.1, color=COLOR_MAP.get(body, None), label=body)
    ax.scatter([0],[0], s=80, c="gold", edgecolor="k", linewidth=0.4, label="Sun")
    ax.legend(fontsize=8, ncols=2)
    out = _prefix_path(run_dir, "orbitcompare__ALL_final.png")
    plt.tight_layout()
    fig.savefig(out, dpi=220)
    plt.close(fig)
    print(f"Wrote {out.name}")

def plot_moon_geocentric(run_dir: Path, params: dict,
                        base_rel, final_rel,
                        t_end, T0, T_end, a0, e0, aE, eE):
    fig, ax = plt.subplots(figsize=(8, 8))
    ax.set_aspect("equal")
    ax.grid(True, alpha=0.3)
    ax.set_xlabel("Δx [AU] (Moon - Earth)")
    ax.set_ylabel("Δy [AU] (Moon - Earth)")
    ax.set_title("Moon geocentric: baseline vs perturbed")

    xb, yb = base_rel
    xf, yf = final_rel

    ax.plot(xb, yb, lw=1.5, color=COLOR_MAP["Moon"], alpha=0.65, label="Baseline (1 lunar orbit)")
    ax.plot(xf, yf, lw=1.7, color=COLOR_MAP["Moon"], label="Final (last lunar orbit)")
    ax.scatter([0], [0], s=50, c=COLOR_MAP["Earth"], edgecolor="k", linewidth=0.3, label="Earth (origin)")

    # Put parameters + derived period info in an in-plot box (instead of the title)
    info_params = dict(params or {})
    info_params["T0_days"] = float(T0)
    info_params["Tend_days"] = float(T_end)
    info_params["t_end_days"] = float(t_end)
    if np.isfinite(a0):
        info_params["a0_AU"] = float(a0)
    if np.isfinite(e0):
        info_params["e0"] = float(e0)
    if np.isfinite(aE):
        info_params["a_end_AU"] = float(aE)
    if np.isfinite(eE):
        info_params["e_end"] = float(eE)

    _add_params_box(ax, info_params, loc=(0.02, 0.98), max_lines=14)

    ax.legend(fontsize=9)
    out = _prefix_path(run_dir, "orbitcompare__Moon_geocentric.png")
    plt.tight_layout()
    fig.savefig(out, dpi=220)
    plt.close(fig)
    print(f"Wrote {out.name}")

# -----------------------
# Main
# -----------------------
def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("run_dir", help="One run subfolder (contains __orbits__*.xlsx).")
    ap.add_argument("--tend", type=float, default=36500.0,
                    help="End time (days) for the 'future orbit' window (default 36500).")
    ap.add_argument("--baseline-step", type=float, default=0.5,
                    help="Sampling step (days) for baseline orbit integration (default 0.5 day).")
    args = ap.parse_args()

    run_dir = Path(args.run_dir)
    xlsx = _find_orbits_xlsx(run_dir)

    params = {"dt_days": 0.5}
    inp = _find_prefixed_input_yaml(run_dir)
    if inp:
        with open(inp, "r", encoding="utf-8") as f:
            params.update(yaml.safe_load(f) or {})
    if "epoch" not in params or "kernel" not in params:
        raise ValueError("Could not find epoch/kernel. Ensure <prefix>__input.yaml exists in run_dir.")

    tend_req = float(args.tend)
    base_step = float(args.baseline_step)

    # Build baseline sim once; we'll rebuild per integration to avoid time carryover
    sim0 = build_sim_no_bh(params)
    mu_sun = sim0.G * MSUN
    mu_earthmoon = sim0.G * (PLANET_MASSES_MSUN["earth"] + MOON_MASS_MSUN)

    # baseline indices: 0 Sun, 1..8 planets, 9 Moon
    name_to_idx0 = {
        "Mercury": 1, "Venus": 2, "Earth": 3, "Mars": 4,
        "Jupiter": 5, "Saturn": 6, "Uranus": 7, "Neptune": 8,
        "Moon": 9,
    }
    earth_idx0 = name_to_idx0["Earth"]
    moon_idx0  = name_to_idx0["Moon"]

    baseline_orbits = {}
    final_orbits = {}

    # We'll need Earth/Moon perturbed heliocentric last-orbit arrays for the combo plot
    earth_base_xy = moon_base_xy = None
    earth_final_xy = moon_final_xy = None

    with pd.ExcelFile(xlsx) as book:
        sun_df = _read_sheet(book, "Sun")
        t_all = sun_df["time_days"].to_numpy(dtype=float)

        j_end = int(np.argmin(np.abs(t_all - tend_req)))
        t_end = float(t_all[j_end])

        sx = sun_df["x_AU"].to_numpy(dtype=float)
        sy = sun_df["y_AU"].to_numpy(dtype=float)
        sz = sun_df["z_AU"].to_numpy(dtype=float) if "z_AU" in sun_df.columns else np.zeros_like(sx)

        # Optional Sun velocities
        if {"vx_AU_per_d","vy_AU_per_d","vz_AU_per_d"}.issubset(set(sun_df.columns)):
            svx = sun_df["vx_AU_per_d"].to_numpy(dtype=float)
            svy = sun_df["vy_AU_per_d"].to_numpy(dtype=float)
            svz = sun_df["vz_AU_per_d"].to_numpy(dtype=float)
        else:
            svx = np.zeros_like(sx); svy = np.zeros_like(sx); svz = np.zeros_like(sx)

        need_vel = {"vx_AU_per_d","vy_AU_per_d","vz_AU_per_d"}

        # --- Per body orbit-compare ---
        for body in PLANET_NAMES:
            if body not in book.sheet_names:
                print(f"[skip] Missing sheet {body} in {xlsx.name}")
                continue

            df = _read_sheet(book, body)

            # positions
            bx = df["x_AU"].to_numpy(dtype=float)
            by = df["y_AU"].to_numpy(dtype=float)
            bz = df["z_AU"].to_numpy(dtype=float) if "z_AU" in df.columns else np.zeros_like(bx)

            n = min(len(t_all), len(bx))
            t = t_all[:n]
            hx = bx[:n] - sx[:n]
            hy = by[:n] - sy[:n]
            hz = bz[:n] - sz[:n]

            if not need_vel.issubset(set(df.columns)):
                raise ValueError(
                    f"Sheet '{body}' lacks velocity columns {need_vel}. "
                    "Please ensure your orbits workbook includes vx_AU_per_d, vy_AU_per_d, vz_AU_per_d."
                )

            bvx = df["vx_AU_per_d"].to_numpy(dtype=float)[:n]
            bvy = df["vy_AU_per_d"].to_numpy(dtype=float)[:n]
            bvz = df["vz_AU_per_d"].to_numpy(dtype=float)[:n]
            hvx = bvx - svx[:n]
            hvy = bvy - svy[:n]
            hvz = bvz - svz[:n]

            # end-state for period inference
            xE, yE, zE = float(hx[j_end]), float(hy[j_end]), float(hz[j_end])
            vxE, vyE, vzE = float(hvx[j_end]), float(hvy[j_end]), float(hvz[j_end])

            T_future, aE, eE = _infer_period_from_state_heliocentric(xE,yE,zE,vxE,vyE,vzE, mu_sun)
            if not np.isfinite(T_future) or T_future <= 0:
                print(f"[warn] {body}: could not compute future period; skipping.")
                continue

            t1 = t_end - T_future
            mask_future = (t >= t1) & (t <= t_end)
            if mask_future.sum() < 5:
                print(f"[warn] {body}: too few samples in last-orbit window; consider smaller output_interval_days.")
                continue

            # baseline at epoch
            sim_base = build_sim_no_bh(params)
            p0 = sim_base.particles[name_to_idx0[body]]
            s0 = sim_base.particles[0]
            x0,y0,z0 = p0.x - s0.x, p0.y - s0.y, p0.z - s0.z
            vx0,vy0,vz0 = p0.vx - s0.vx, p0.vy - s0.vy, p0.vz - s0.vz
            T0, a0, e0 = _infer_period_from_state_heliocentric(x0,y0,z0,vx0,vy0,vz0, mu_sun)
            base_xy = integrate_baseline_orbit_xy(build_sim_no_bh(params), name_to_idx0[body], T0, base_step)
            if base_xy is None:
                print(f"[warn] {body}: baseline integration failed; skipping.")
                continue

            xb, yb = base_xy
            xf, yf = hx[mask_future], hy[mask_future]

            baseline_orbits[body] = (xb, yb)
            final_orbits[body] = (xf, yf)

            if body == "Earth":
                earth_base_xy = (xb, yb)
                earth_final_xy = (xf, yf)
            if body == "Moon":
                moon_base_xy = (xb, yb)
                moon_final_xy = (xf, yf)

            plot_one_body(run_dir, body, xb, yb, xf, yf, t_end, T_future, aE, eE)

        # --- Combined plots (ALL + Earth/Moon) ---
        if baseline_orbits:
            plot_all_epoch(run_dir, baseline_orbits)
        if final_orbits:
            plot_all_final(run_dir, final_orbits, t_end)
        if earth_base_xy and moon_base_xy and earth_final_xy and moon_final_xy:
            plot_earth_moon_combo(run_dir, earth_base_xy, moon_base_xy, earth_final_xy, moon_final_xy, t_end)

        # --- Geocentric Moon plot (Moon relative to Earth) ---
        if ("Earth" in book.sheet_names) and ("Moon" in book.sheet_names):
            earth_df = _read_sheet(book, "Earth")
            moon_df  = _read_sheet(book, "Moon")

            # positions
            ex = earth_df["x_AU"].to_numpy(dtype=float)
            ey = earth_df["y_AU"].to_numpy(dtype=float)
            ez = earth_df["z_AU"].to_numpy(dtype=float) if "z_AU" in earth_df.columns else np.zeros_like(ex)

            mx = moon_df["x_AU"].to_numpy(dtype=float)
            my = moon_df["y_AU"].to_numpy(dtype=float)
            mz = moon_df["z_AU"].to_numpy(dtype=float) if "z_AU" in moon_df.columns else np.zeros_like(mx)

            n2 = min(len(t_all), len(ex), len(mx))
            t2 = t_all[:n2]

            # velocities
            if not need_vel.issubset(set(earth_df.columns)) or not need_vel.issubset(set(moon_df.columns)):
                raise ValueError("Earth/Moon sheets must include vx_AU_per_d, vy_AU_per_d, vz_AU_per_d for geocentric plot.")

            evx = earth_df["vx_AU_per_d"].to_numpy(dtype=float)[:n2]
            evy = earth_df["vy_AU_per_d"].to_numpy(dtype=float)[:n2]
            evz = earth_df["vz_AU_per_d"].to_numpy(dtype=float)[:n2]

            mvx = moon_df["vx_AU_per_d"].to_numpy(dtype=float)[:n2]
            mvy = moon_df["vy_AU_per_d"].to_numpy(dtype=float)[:n2]
            mvz = moon_df["vz_AU_per_d"].to_numpy(dtype=float)[:n2]

            # relative state (Moon - Earth), in inertial frame; Sun subtraction cancels in the difference
            dx = (mx[:n2] - ex[:n2])
            dy = (my[:n2] - ey[:n2])
            dz = (mz[:n2] - ez[:n2])

            dvx = (mvx - evx)
            dvy = (mvy - evy)
            dvz = (mvz - evz)

            # baseline geocentric period from EARLY perturbed-run relative state (Moon-Earth)
            # (This matches whatever Moon initial conditions were used in the actual run.)
            i0 = 0
            T0_geo, a0_geo, e0_geo = _infer_period_from_state_relative(
                float(dx[i0]), float(dy[i0]), float(dz[i0]),
                float(dvx[i0]), float(dvy[i0]), float(dvz[i0]),
                mu_earthmoon
            )
            if (not np.isfinite(T0_geo)) or T0_geo <= 0:
                # fallback: sidereal month
                T0_geo, a0_geo, e0_geo = 27.321661, float("nan"), float("nan")
            t0 = float(t2[i0])
            mask_base = (t2 >= t0) & (t2 <= (t0 + T0_geo))
            if mask_base.sum() < 5:
                print("[warn] Moon geocentric: too few samples in first-lunar-orbit window; skipping.")
                base_rel = None
            else:
                base_rel = (dx[mask_base], dy[mask_base])
            if base_rel is None:
                print("[warn] Moon geocentric: baseline integration failed; skipping.")
            else:
                # end-state geocentric period from final relative state at t_end index (same j_end)
                if j_end >= n2:
                    jg = n2 - 1
                else:
                    jg = j_end

                T_end_geo, aE_geo, eE_geo = _infer_period_from_state_relative(
                    float(dx[jg]), float(dy[jg]), float(dz[jg]),
                    float(dvx[jg]), float(dvy[jg]), float(dvz[jg]),
                    mu_earthmoon
                )
                if (not np.isfinite(T_end_geo)) or T_end_geo <= 0:
                    # fallback
                    T_end_geo, aE_geo, eE_geo = 27.321661, float("nan"), float("nan")

                t1g = t_end - T_end_geo
                mask_g = (t2 >= t1g) & (t2 <= t_end)
                if mask_g.sum() < 5:
                    # fallback: just take last N samples covering ~T_end_geo
                    # (works even if cadence is coarse)
                    # choose points with t within [t_end - T_end_geo, t_end]
                    mask_g = (t2 >= (t_end - T_end_geo)) & (t2 <= t_end)

                if mask_g.sum() < 5:
                    print("[warn] Moon geocentric: too few samples in last-lunar-orbit window; skipping.")
                else:
                    final_rel = (dx[mask_g], dy[mask_g])
                    plot_moon_geocentric(run_dir, params, base_rel, final_rel, t_end,
                                         T0_geo, T_end_geo, a0_geo, e0_geo, aE_geo, eE_geo)
        else:
            print("[warn] Moon geocentric: Earth or Moon sheet missing; skipping.")

if __name__ == "__main__":
    main()