#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Plot a 2D heliocentric snapshot of Sun + planets + BH at BH perihelion from a REBOUND SimulationArchive.

- Hash-first lookup (recommended for archives produced by solar_system_bh_rebound22.py)
- Index fallback for older archives

Run from the project root. You may pass:
- archive filename (searched under ./simulations/**/), or
- full/relative path

Examples:
  python plot_bh_perihelion.py RUNNAME__archive.bin --sun-hash Sun --bh-hash BH
  python plot_bh_perihelion.py RUNNAME__archive.bin --sun-index 0 --bh-index 9
"""

import argparse
import math
import os
from pathlib import Path
from typing import Optional, Tuple

import matplotlib.pyplot as plt
import rebound

PLANET_HASHES_DEFAULT = "Mercury,Venus,Earth,Mars,Jupiter,Saturn,Uranus,Neptune"


def load_archive(path: str):
    """Version-robust SimulationArchive loader."""
    if hasattr(rebound, "SimulationArchive"):
        return rebound.SimulationArchive(path)
    if hasattr(rebound, "Simulationarchive"):
        return rebound.Simulationarchive(path)
    raise RuntimeError("This REBOUND install does not expose SimulationArchive/Simulationarchive.")


def resolve_archive_path(user_arg: str, simulations_dir: str = "simulations") -> Path:
    """
    If `user_arg` is a file path, use it. Otherwise treat as filename and search under ./simulations/**/
    """
    p = Path(user_arg)
    if p.is_file():
        return p.resolve()

    root = Path(simulations_dir)
    if not root.is_dir():
        raise RuntimeError(f"Simulations directory not found: {root.resolve()}")

    matches = list(root.rglob(p.name))
    if not matches:
        raise RuntimeError(f"Could not find '{p.name}' under {root.resolve()}. Pass full path if needed.")
    if len(matches) > 1:
        msg = "\n".join(str(m) for m in matches[:30])
        more = "" if len(matches) <= 30 else f"\n...and {len(matches) - 30} more"
        raise RuntimeError(
            f"Multiple matches found for '{p.name}'. Please pass a full path.\nMatches:\n{msg}{more}"
        )
    return matches[0].resolve()


def infer_output_path_from_archive(archive_path: Path, suffix: str = "__perihelion_planets.png") -> Path:
    name = archive_path.name
    if name.endswith("__archive.bin"):
        base = name[:-len("__archive.bin")]
    elif name.endswith(".bin"):
        base = name[:-4]
    else:
        base = archive_path.stem
    return archive_path.parent / f"{base}{suffix}"


def get_particle(sim, hash_name: Optional[str], index: Optional[int], role: str):
    """
    Hash-first lookup; index fallback.
    """
    if hash_name:
        try:
            return sim.particles[hash_name]
        except Exception:
            if index is None:
                raise RuntimeError(f"Could not find {role} by hash '{hash_name}', and no index fallback given.")
    if index is not None:
        return sim.particles[index]
    raise RuntimeError(f"Must provide either --{role.lower()}-hash or --{role.lower()}-index.")


def heliocentric_xyz(sim, sun_p, p):
    return (p.x - sun_p.x, p.y - sun_p.y, p.z - sun_p.z)


def r_xyz(x, y, z):
    return math.sqrt(x * x + y * y + z * z)


def find_bh_perihelion_snapshot(sa, sun_hash, sun_index, bh_hash, bh_index) -> Tuple[int, float, float]:
    best_k = 0
    best_t = 0.0
    best_r = float("inf")

    for k in range(len(sa)):
        sim = sa[k]
        sun = get_particle(sim, sun_hash, sun_index, "Sun")
        bh = get_particle(sim, bh_hash, bh_index, "BH")
        x, y, z = heliocentric_xyz(sim, sun, bh)
        r = r_xyz(x, y, z)
        if r < best_r:
            best_r = r
            best_k = k
            best_t = float(sim.t)

    return best_k, best_t, float(best_r)


def parse_csv_list(s: str):
    return [x.strip() for x in s.split(",") if x.strip()]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("archive", help="Archive filename or path (e.g. RUNNAME__archive.bin)")

    # Hash-first (new world)
    ap.add_argument("--sun-hash", default=None, help="Sun particle hash (e.g. Sun)")
    ap.add_argument("--bh-hash", default=None, help="BH particle hash (e.g. BH)")
    ap.add_argument("--planet-hashes", default=PLANET_HASHES_DEFAULT,
                    help=f"Comma-separated planet hashes. Default: {PLANET_HASHES_DEFAULT}")

    # Index fallback (old world)
    ap.add_argument("--sun-index", type=int, default=None, help="Fallback Sun index (e.g. 0)")
    ap.add_argument("--bh-index", type=int, default=None, help="Fallback BH index (e.g. 9)")
    ap.add_argument("--planet-indices", default=None,
                    help="Fallback comma-separated planet indices (e.g. '1,2,3,4,5,6,7,8')")

    ap.add_argument("--lim-au", type=float, default=35.0, help="Plot half-width in AU (default 35)")
    ap.add_argument("--no-labels", action="store_true", help="Disable labels")
    args = ap.parse_args()

    archive_path = resolve_archive_path(args.archive, simulations_dir="simulations")
    if not os.path.isfile(archive_path):
        raise RuntimeError(f"File not found: {archive_path}")
    print(f"Using archive: {archive_path}")

    sa = load_archive(str(archive_path))
    if len(sa) == 0:
        raise RuntimeError("Archive is empty.")

    planet_hashes = parse_csv_list(args.planet_hashes)

    planet_indices = None
    if args.planet_indices:
        planet_indices = [int(x.strip()) for x in args.planet_indices.split(",") if x.strip()]

    best_k, best_t, best_r = find_bh_perihelion_snapshot(
        sa, args.sun_hash, args.sun_index, args.bh_hash, args.bh_index
    )
    sim = sa[best_k]
    sun = get_particle(sim, args.sun_hash, args.sun_index, "Sun")
    bh = get_particle(sim, args.bh_hash, args.bh_index, "BH")

    # Build planet particle list: try hashes; fallback to indices if needed
    planets = []
    for j, nm in enumerate(planet_hashes):
        try:
            p = get_particle(sim, nm, None, role=f"Planet({nm})")
            planets.append((nm, p))
        except Exception:
            if planet_indices is None or j >= len(planet_indices):
                raise RuntimeError(
                    f"Could not find planet '{nm}' by hash, and no usable --planet-indices fallback."
                )
            p = sim.particles[planet_indices[j]]
            planets.append((nm, p))

    # Plot
    plt.figure(figsize=(10, 10))
    plt.axhline(0, linewidth=0.8)
    plt.axvline(0, linewidth=0.8)

    # Sun at origin
    plt.scatter([0], [0], s=220, marker="*", zorder=6)
    if not args.no_labels:
        plt.text(0, 0, " Sun", va="center", ha="left")

    # Planets
    for nm, p in planets:
        x, y, z = heliocentric_xyz(sim, sun, p)
        plt.scatter([x], [y], s=90, zorder=5)
        if not args.no_labels:
            plt.text(x, y, f" {nm}", va="center", ha="left", fontsize=9)

    # BH
    x_bh, y_bh, z_bh = heliocentric_xyz(sim, sun, bh)
    plt.scatter([x_bh], [y_bh], s=240, marker="X", zorder=7)
    if not args.no_labels:
        plt.text(x_bh, y_bh, " BH", va="center", ha="left")

    plt.gca().set_aspect("equal", adjustable="box")
    plt.xlim(-args.lim_au, args.lim_au)
    plt.ylim(-args.lim_au, args.lim_au)
    plt.xlabel("x (AU, heliocentric)")
    plt.ylabel("y (AU, heliocentric)")
    plt.title(f"Planets + BH at BH perihelion | snapshot={best_k} | t={best_t:.6g} days | r_peri={best_r:.6g} AU")

    out_path = infer_output_path_from_archive(archive_path)
    plt.tight_layout()
    plt.savefig(out_path, dpi=200)

    print(f"Saved: {out_path}")
    print(f"BH perihelion snapshot: {best_k}  t={best_t}  r_peri={best_r} AU")


if __name__ == "__main__":
    main()
