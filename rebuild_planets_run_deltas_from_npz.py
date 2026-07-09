#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Rebuild version-18 style *__planets_run_deltas.csv files from snapshot_*.npz.

For each run subfolder under a parent folder (e.g. simulations/20251127_201847),
this script:
  - finds earliest and latest snapshot_*.npz
  - computes heliocentric (Sun-centered) a,e,q,Q for:
        Mercury..Neptune + Moon
  - writes:
        <subfolder_name>__planets_run_deltas.csv

This is meant to reconstruct the "full planets" CSV even if v19 only wrote inner planets.

Assumptions (true for your scripts):
  - Particle order in snapshots:
        0 Sun
        1 Mercury
        2 Venus
        3 Earth
        4 Mars
        5 Jupiter
        6 Saturn
        7 Uranus
        8 Neptune
        9 Moon
        10 BH
        11.. belt tracers...
  - Units: AU, AU/day, Msun with G = 0.0002959122082855911 AU^3/(Msun*day^2)
"""

import argparse
import math
from pathlib import Path
import numpy as np
import pandas as pd

G_AU3_Msun_day2 = 0.0002959122082855911
MSUN = 1.0

BODY_INDEX = {
    "Mercury": 1,
    "Venus": 2,
    "Earth": 3,
    "Mars": 4,
    "Jupiter": 5,
    "Saturn": 6,
    "Uranus": 7,
    "Neptune": 8,
    "Moon": 9,
}

PLANETS_FULL = ["Mercury","Venus","Earth","Mars","Jupiter","Saturn","Uranus","Neptune","Moon"]


def prefix_path(subdir: Path, name: str) -> Path:
    return subdir / f"{subdir.name}__{name}"


def find_snapshots(subdir: Path):
    snaps = sorted(subdir.glob(f"{subdir.name}__snapshot_*.npz"))
    if not snaps:
        # Some older runs might have non-prefixed snapshots:
        snaps = sorted(subdir.glob("snapshot_*.npz"))
    return snaps


def load_snapshot(path: Path):
    d = np.load(path, allow_pickle=True)
    arr = d["arr"]  # shape (N, 6) = x,y,z,vx,vy,vz
    return arr


def elements_from_state(r_vec, v_vec, mu):
    """
    Return a,e,q,Q for a heliocentric two-body osculating orbit.

    r_vec, v_vec: 3-vectors in AU and AU/day
    mu: G*(M_sun) in AU^3/day^2  (here we use Sun-only mu; consistent with your prior outputs)
    """
    r = np.linalg.norm(r_vec)
    v2 = float(np.dot(v_vec, v_vec))
    if r <= 0:
        return float("nan"), float("nan"), float("nan"), float("nan")

    h_vec = np.cross(r_vec, v_vec)
    h2 = float(np.dot(h_vec, h_vec))

    # Specific orbital energy
    eps = 0.5 * v2 - mu / r

    # a = -mu/(2 eps)
    if abs(eps) < 1e-30:
        a = float("inf")
    else:
        a = -mu / (2.0 * eps)

    # e-vector
    e_vec = np.cross(v_vec, h_vec) / mu - (r_vec / r)
    e = float(np.linalg.norm(e_vec))

    # perihelion/aphelion for bound orbits
    if math.isfinite(a) and e < 1.0:
        q = a * (1.0 - e)
        Q = a * (1.0 + e)
    else:
        # hyperbolic/unbound: Q not defined; q from h^2/(mu*(1+e))
        q = (h2 / (mu * (1.0 + e))) if (1.0 + e) > 1e-15 else float("nan")
        Q = float("nan")

    return float(a), float(e), float(q), float(Q)


def compute_table(arr, mu):
    """
    arr: snapshot array (N,6)
    Returns dict body -> (a,e,q,Q) heliocentric
    """
    sun = arr[0, :]
    sun_r = sun[0:3]
    sun_v = sun[3:6]

    out = {}
    for body in PLANETS_FULL:
        idx = BODY_INDEX[body]
        p = arr[idx, :]
        r_vec = (p[0:3] - sun_r)
        v_vec = (p[3:6] - sun_v)
        out[body] = elements_from_state(r_vec, v_vec, mu)
    return out


def process_subfolder(subdir: Path, overwrite: bool):
    snaps = find_snapshots(subdir)
    if len(snaps) < 2:
        print(f"[skip] {subdir.name}: not enough snapshots ({len(snaps)})")
        return

    before_path = snaps[0]
    after_path = snaps[-1]

    out_csv = prefix_path(subdir, "planets_run_deltas.csv")
    if out_csv.exists() and not overwrite:
        print(f"[skip] {subdir.name}: output exists")
        return

    arr_b = load_snapshot(before_path)
    arr_a = load_snapshot(after_path)

    mu = G_AU3_Msun_day2 * MSUN

    tab_b = compute_table(arr_b, mu)
    tab_a = compute_table(arr_a, mu)

    rows = []
    for body in PLANETS_FULL:
        a_b, e_b, q_b, Q_b = tab_b[body]
        a_a, e_a, q_a, Q_a = tab_a[body]
        rows.append({
            "body": body,
            "a_before_AU": a_b,
            "e_before": e_b,
            "q_before_AU": q_b,
            "Q_before_AU": Q_b,
            "a_after_AU": a_a,
            "e_after": e_a,
            "q_after_AU": q_a,
            "Q_after_AU": Q_a,
        })

    df = pd.DataFrame(rows)
    df.to_csv(out_csv, index=False)
    print(f"[ok] {subdir.name}: wrote {out_csv.name} (from {before_path.name} and {after_path.name})")


def main():
    ap = argparse.ArgumentParser(
        description="Rebuild full planets_run_deltas.csv from snapshot_*.npz in each run subfolder."
    )
    ap.add_argument(
        "parent",
        help="Simulation parent folder, e.g. simulations/20251127_201847",
    )
    ap.add_argument(
        "--overwrite",
        action="store_true",
        help="Overwrite existing __planets_run_deltas.csv files",
    )
    args = ap.parse_args()

    parent = Path(args.parent).resolve()
    if not parent.is_dir():
        raise FileNotFoundError(f"Parent folder not found: {parent}")

    subdirs = [p for p in parent.iterdir() if p.is_dir()]
    if not subdirs:
        print("[info] No subfolders found.")
        return

    for sub in sorted(subdirs, key=lambda p: p.name):
        process_subfolder(sub, overwrite=args.overwrite)


if __name__ == "__main__":
    main()
