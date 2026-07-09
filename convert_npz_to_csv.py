#!/usr/bin/env python3
import numpy as np
import pandas as pd
from pathlib import Path
import yaml

# --------- CONFIGURE THIS ---------
RUN_DIR = Path(r"simulations/20251112_230429/20251112_230429__rp10__vinf30")
OUT_DIR = RUN_DIR / "csv_frames"
# ----------------------------------

# Names/order of the main bodies in your simulation
BODY_NAMES = ["Sun","Mercury","Venus","Earth","Mars",
              "Jupiter","Saturn","Uranus","Neptune","Moon","BH"]

# Masses in solar masses (same as main script; Sun & BH handled separately)
PLANET_MASSES_MSUN = {
    "Mercury": 1.651e-7,
    "Venus":   2.447e-6,
    "Earth":   3.003e-6,
    "Mars":    3.227e-7,
    "Jupiter": 9.5479e-4,
    "Saturn":  2.8585e-4,
    "Uranus":  4.3662e-5,
    "Neptune": 5.1513e-5,
    "Moon":    3.694e-8,
}

# Arbitrary visualization sizes (roughly mass-related; unitless)
SIZES = {
    "Sun":     1.5,
    "Mercury": 0.3,
    "Venus":   0.5,
    "Earth":   0.5,
    "Mars":    0.4,
    "Jupiter": 1.0,
    "Saturn":  0.9,
    "Uranus":  0.8,
    "Neptune": 0.8,
    "Moon":    0.3,
    "BH":      2.0,
}

# RGB colours (0–255) for the 11 main objects; belt tracers will be grey
COLOURS = {
    "Sun":     (255, 220,  0),
    "Mercury": (150,150,150),
    "Venus":   (255,165,  0),
    "Earth":   ( 80,150,255),
    "Mars":    (200, 60, 40),
    "Jupiter": (200,170,120),
    "Saturn":  (230,210,120),
    "Uranus":  (140,220,220),
    "Neptune": ( 70,140,255),
    "Moon":    (220,220,220),
    "BH":      (  0,  0,  0),
}

# Belt tracer appearance
TRACER_COLOR = (180, 180, 180)
TRACER_SIZE  = 0.1   # small


def find_input_yaml(run_dir: Path) -> Path:
    """Find the copied input.yaml in the run directory."""
    candidates = list(run_dir.glob("*__input.yaml"))
    if not candidates:
        raise FileNotFoundError(f"No '*__input.yaml' found in {run_dir}")
    # If several exist, just take the first.
    return candidates[0]


def load_bh_mass(run_dir: Path) -> float:
    """Read bh_mass_msun from the run's input YAML."""
    params_path = find_input_yaml(run_dir)
    with open(params_path, "r", encoding="utf-8") as f:
        params = yaml.safe_load(f)
    return float(params["bh_mass_msun"])


def convert_one(npz_path: Path, out_path: Path, bh_mass_msun: float):
    data = np.load(npz_path, allow_pickle=True)
    arr = data["arr"]                 # shape (N, 6)
    N = arr.shape[0]

    # Base names from file, if present (Sun..BH)
    if "names" in data and len(data["names"]) > 0:
        base_names = list(data["names"])
    else:
        base_names = []

    # Build full names list of length N (main bodies + tracers)
    if len(base_names) >= N:
        names = base_names[:N]
    else:
        n_extra = N - len(base_names)
        tracer_names = [f"tracer_{i}" for i in range(n_extra)]
        names = base_names + tracer_names

    # Parse time (days) from filename: snapshot_000000.npz -> t_days = 0
    stem = npz_path.stem
    try:
        frame_str = stem.split("_")[-1]
        t_days = float(frame_str)
    except Exception:
        t_days = 0.0

    # Pre-allocate visual attributes
    mass_msun = np.zeros(N, dtype=float)
    size = np.zeros(N, dtype=float)
    color_r = np.zeros(N, dtype=int)
    color_g = np.zeros(N, dtype=int)
    color_b = np.zeros(N, dtype=int)

    for i, name in enumerate(names):
        if name in BODY_NAMES:
            # Main bodies
            if name == "Sun":
                m = 1.0
            elif name == "BH":
                m = bh_mass_msun
            else:
                m = PLANET_MASSES_MSUN.get(name, 0.0)

            mass_msun[i] = m
            size[i] = SIZES.get(name, 0.5)
            r, g, b = COLOURS.get(name, (255, 255, 255))
            color_r[i], color_g[i], color_b[i] = r, g, b
        else:
            # Belt tracers
            mass_msun[i] = 0.0
            size[i] = TRACER_SIZE
            r, g, b = TRACER_COLOR
            color_r[i], color_g[i], color_b[i] = r, g, b

    df = pd.DataFrame({
        "t_days":    np.full(N, t_days),
        "name":      names,
        "index":     np.arange(N),
        "mass_msun": mass_msun,
        "size":      size,
        "color_r":   color_r,
        "color_g":   color_g,
        "color_b":   color_b,
        "x":  arr[:, 0],
        "y":  arr[:, 1],
        "z":  arr[:, 2],
        "vx": arr[:, 3],
        "vy": arr[:, 4],
        "vz": arr[:, 5],
    })

    df.to_csv(out_path, index=False)


def main():
    OUT_DIR.mkdir(exist_ok=True)

    snapshots = sorted(RUN_DIR.glob("*snapshot_*.npz"))
    if not snapshots:
        print("No NPZ snapshots found in", RUN_DIR)
        return

    bh_mass = load_bh_mass(RUN_DIR)
    print(f"BH mass from YAML: {bh_mass} Msun")
    print(f"Found {len(snapshots)} snapshots.")

    for snap in snapshots:
        frame = snap.stem.split("_")[-1]
        out_file = OUT_DIR / f"snapshot_{frame}.csv"
        print(f"Writing {out_file}")
        convert_one(snap, out_file, bh_mass_msun=bh_mass)

    print("DONE — CSV files with size & colours are in:", OUT_DIR)


if __name__ == "__main__":
    main()
