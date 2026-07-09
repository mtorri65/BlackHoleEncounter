#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import argparse
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from PIL import Image


def find_bh_track_from_name(filename: str, base_dir: Path) -> Path:
    timestamp = filename.split("__")[0]
    sim_root = base_dir / "simulations" / timestamp
    if not sim_root.exists():
        raise FileNotFoundError(f"Simulation folder not found: {sim_root}")

    matches = list(sim_root.rglob(filename))
    if not matches:
        raise FileNotFoundError(f"Could not find {filename} under {sim_root}")

    if len(matches) > 1:
        print("⚠️ Multiple matches found, using first:")
        for m in matches:
            print("   ", m)

    return matches[0]


def find_col(df, candidates):
    for c in candidates:
        if c in df.columns:
            return c
    raise KeyError(
        f"None of the columns {candidates} found.\nAvailable columns: {df.columns.tolist()}"
    )


def main():
    ap = argparse.ArgumentParser(description="Overlay BH sky trajectory on Vera Rubin window image")
    ap.add_argument(
        "--bh-radec",
        required=True,
        help="BH RA/Dec Excel filename (e.g. 20260121_114641__...__bh_radec__1p5__0p1.xlsx)"
    )
    ap.add_argument(
        "--image",
        default="VR-MayNov-10-mirrored-boundaries.jpg",
        help="Background sky image"
    )
    ap.add_argument(
        "--out",
        default="bh_track_overlay.png",
        help="Output image filename"
    )
    args = ap.parse_args()

    base_dir = Path.cwd()

    # Rectangle sky bounds (degrees) from your vertices
    RA_MIN  = 306.5   # 20h26m
    RA_MAX  = 311.0   # 20h44m
    DEC_MIN = -13.6   # -13d36m
    DEC_MAX = -12.1   # -12d06m

    # Locate BH track
    bh_path = find_bh_track_from_name(args.bh_radec, base_dir)
    print(f"Using BH track: {bh_path}")

    bh = pd.read_excel(bh_path)

    # Detect columns and convert RA hours->deg if needed
    ra_col  = find_col(bh, ["ra_deg", "RA_deg", "ra", "RA", "RA_hours", "ra_hours"])
    dec_col = find_col(bh, ["dec_deg", "Dec_deg", "DEC_deg", "dec", "Dec", "DEC"])

    ra_raw = bh[ra_col].astype(float).to_numpy()
    dec    = bh[dec_col].astype(float).to_numpy()

    if "hour" in ra_col.lower():
        ra = ra_raw * 15.0
        ra_units = "hours→deg"
    else:
        ra = ra_raw
        ra_units = "deg"

    print(f"Using RA column: {ra_col} ({ra_units}), Dec column: {dec_col} (deg)")

    # Load background image
    img = np.asarray(Image.open(args.image))
    ny, nx = img.shape[:2]

    # Map RA/Dec -> pixel coords
    x = (ra - RA_MIN) / (RA_MAX - RA_MIN) * nx
    y = (DEC_MAX - dec) / (DEC_MAX - DEC_MIN) * ny

    # Plot overlay
    plt.figure(figsize=(12, 6))
    plt.imshow(img)

    plt.plot(x, y, color="cyan", lw=2.8, label="BH trajectory")
    plt.scatter(x[0],  y[0],  s=100, color="lime", zorder=5, label="Start (May)")
    plt.scatter(x[-1], y[-1], s=100, color="red",  zorder=5, label="End (Nov)")

    # Direction arrows
    step = max(1, len(x) // 20)
    for i in range(0, len(x) - step, step):
        plt.arrow(
            x[i], y[i],
            x[i + 1] - x[i],
            y[i + 1] - y[i],
            color="cyan",
            alpha=0.7,
            head_width=10,
            length_includes_head=True
        )

    plt.axis("off")
    plt.legend(loc="lower left")
    plt.title("BH Sky Trajectory over Vera Rubin May–Nov Window")
    plt.tight_layout()

    # Save + show
    plt.savefig(args.out, dpi=200)
    print(f"Wrote: {args.out}")
    plt.show()


if __name__ == "__main__":
    main()
