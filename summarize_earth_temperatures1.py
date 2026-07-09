# -*- coding: utf-8 -*-
"""
summarize_earth_temperatures1.py

Post-processing helper for Solar System + BH runs.

Given a simulation parent folder such as:
    simulations/20251111_224816/

This script:
  * Reads Earth data from each *__planets_run_deltas.csv
  * Computes perihelion/aphelion before & after
  * Computes idealized temperatures (T ∝ r^-1/2)
  * Reads BH parameters from input.yaml
  * Writes a summary CSV in the parent folder
  * Prints a counter for processed subfolders
"""

import argparse
import math
from pathlib import Path

import numpy as np
import pandas as pd
import yaml


# ----------------------------------------------------------------------
# Temperature helper
# ----------------------------------------------------------------------
def effective_temp_from_distance(
    r_au: float,
    T_ref_K: float = 288.0,
    r_ref_au: float = 1.0
) -> float:
    """Idealized temperature scaling: T ∝ r^{-1/2}."""
    if r_au <= 0:
        return float("nan")
    return T_ref_K * math.sqrt(r_ref_au / r_au)


# ----------------------------------------------------------------------
# Per-parent processing
# ----------------------------------------------------------------------
def process_parent_folder(parent: Path) -> Path:
    parent = parent.resolve()
    if not parent.is_dir():
        raise FileNotFoundError(f"Parent folder not found: {parent}")

    rows = []
    run_counter = 0

    subdirs = sorted(p for p in parent.iterdir() if p.is_dir())

    for idx, sub in enumerate(subdirs, start=1):
        planets_files = list(sub.glob("*__planets_run_deltas.csv"))
        if not planets_files:
            continue

        run_counter += 1
        print(f"[{run_counter}] Processing: {sub.name}")

        planets_path = planets_files[0]

        # Load input.yaml (BH parameters)
        yaml_files = list(sub.glob("*__input.yaml"))
        params = {}
        if yaml_files:
            try:
                with open(yaml_files[0], "r", encoding="utf-8") as f:
                    params = yaml.safe_load(f) or {}
            except Exception as e:
                print(f"[warn] Could not read YAML in {sub.name}: {e}")

        bh_mass_msun = params.get("bh_mass_msun", np.nan)
        bh_rp_au     = params.get("bh_rp_au", np.nan)
        bh_vinf_kms  = params.get("bh_vinf_kms", np.nan)
        bh_inc_deg   = params.get("bh_inc_deg", np.nan)

        # Load planets deltas
        try:
            df = pd.read_csv(planets_path)
        except Exception as e:
            print(f"[warn] Could not read {planets_path}: {e}")
            continue

        if "body" not in df.columns:
            print(f"[warn] No 'body' column in {planets_path}")
            continue

        earth_rows = df[df["body"] == "Earth"]
        if earth_rows.empty:
            print(f"[warn] No Earth row in {planets_path}")
            continue

        earth = earth_rows.iloc[0]

        a_before = float(earth.get("a_before", np.nan))
        e_before = float(earth.get("e_before", np.nan))
        q_before = float(earth.get("q_before", np.nan))

        a_after = float(earth.get("a_after", np.nan))
        e_after = float(earth.get("e_after", np.nan))
        q_after = float(earth.get("q_after", np.nan))

        Q_before = a_before * (1 + e_before) if np.isfinite(a_before) and np.isfinite(e_before) else np.nan
        Q_after  = a_after  * (1 + e_after)  if np.isfinite(a_after)  and np.isfinite(e_after)  else np.nan

        T_peri_before_K = effective_temp_from_distance(q_before)
        T_aphe_before_K = effective_temp_from_distance(Q_before)
        T_peri_after_K  = effective_temp_from_distance(q_after)
        T_aphe_after_K  = effective_temp_from_distance(Q_after)

        rows.append({
            "run_folder": sub.name,

            "bh_mass_msun": bh_mass_msun,
            "bh_rp_au": bh_rp_au,
            "bh_vinf_kms": bh_vinf_kms,
            "bh_inc_deg": bh_inc_deg,

            "a_before_AU": a_before,
            "e_before": e_before,
            "q_before_AU": q_before,
            "Q_before_AU": Q_before,
            "T_peri_before_K": T_peri_before_K,
            "T_aphe_before_K": T_aphe_before_K,
            "T_peri_before_C": T_peri_before_K - 273.15,
            "T_aphe_before_C": T_aphe_before_K - 273.15,

            "a_after_AU": a_after,
            "e_after": e_after,
            "q_after_AU": q_after,
            "Q_after_AU": Q_after,
            "T_peri_after_K": T_peri_after_K,
            "T_aphe_after_K": T_aphe_after_K,
            "T_peri_after_C": T_peri_after_K - 273.15,
            "T_aphe_after_C": T_aphe_after_K - 273.15,
        })

    if not rows:
        raise RuntimeError(f"No Earth data found under {parent}")

    df_out = pd.DataFrame(rows)

    sort_cols = [c for c in ["bh_rp_au", "bh_vinf_kms", "bh_inc_deg", "run_folder"] if c in df_out.columns]
    if sort_cols:
        df_out = df_out.sort_values(sort_cols)

    out_path = parent / f"{parent.name}__earth_temp_summary.csv"
    df_out.to_csv(out_path, index=False)

    print(f"\nProcessed {run_counter} run subfolders.")
    print(f"Wrote summary to {out_path}")

    return out_path


# ----------------------------------------------------------------------
# CLI
# ----------------------------------------------------------------------
def main():
    parser = argparse.ArgumentParser(
        description="Summarize Earth perihelion/aphelion distances and temperatures "
                    "from planets_run_deltas.csv files under a simulation parent folder."
    )
    parser.add_argument(
        "parent_folder",
        help="Path to a simulation parent folder, e.g. simulations/20251111_224816"
    )
    args = parser.parse_args()

    parent = Path(args.parent_folder)
    process_parent_folder(parent)


if __name__ == "__main__":
    main()
