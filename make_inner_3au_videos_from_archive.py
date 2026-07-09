#!/usr/bin/env python3
# -*- coding: utf-8 -*-

r"""
Generate two MP4 videos from a REBOUND SimulationArchive:
  - *__inner_xy_orbits.mp4  (2D x–y)
  - *__inner_3d_orbits.mp4  (3D)

Inner region view is limited to +/- RMAX AU around the Sun (default 3 AU),
but includes ONLY the portions of each body's trajectory that are inside that region:
- If r <= RMAX: append the real (x,y[,z]) to the trail and show the marker.
- If r >  RMAX: append NaN to break the trail and hide the marker.

Hash-first (recommended for archives from solar_system_bh_rebound24.py):
  python make_inner_3au_videos_from_archive_indicator.py RUNNAME__archive.bin --sun-hash Sun --bh-hash BH

Index fallback (older anonymous archives):
  python make_inner_3au_videos_from_archive_indicator.py RUNNAME__archive.bin --sun-index 0 --bh-index 9 --planet-start 1 --planet-count 8

Optional:
  --planet-hashes Mercury,Venus,Earth,Mars,Jupiter,Saturn,Uranus,Neptune
  --include-moon
  --rmax-au 3
  --stride 2
  --fps 30
  --trail-max 3000
  --progress-every 50
  --no-labels
"""

import argparse
import math
import time
from pathlib import Path
from typing import Dict, List, Optional, Tuple

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation, FFMpegWriter

import rebound

PLANET_NAMES = ["Mercury", "Venus", "Earth", "Mars", "Jupiter", "Saturn", "Uranus", "Neptune"]


def load_archive(path: str):
    """Version-robust SimulationArchive loader."""
    if hasattr(rebound, "SimulationArchive"):
        return rebound.SimulationArchive(path)
    if hasattr(rebound, "Simulationarchive"):
        return rebound.Simulationarchive(path)
    raise RuntimeError("This REBOUND install does not expose SimulationArchive/Simulationarchive.")


def resolve_archive_path(user_arg: str, simulations_dir: str = "simulations") -> Path:
    """
    Accept either:
      - a full/relative path to the archive, OR
      - just the archive filename (e.g. RUNNAME__archive.bin)

    If only a filename is given, search recursively under ./simulations/**/ for it.
    """
    p = Path(user_arg)

    if p.is_file():
        return p.resolve()

    root = Path(simulations_dir)
    if not root.is_dir():
        raise RuntimeError(f"Simulations directory not found: {root.resolve()}")

    matches = list(root.rglob(p.name))
    if not matches:
        raise RuntimeError(
            f"Could not find '{p.name}' under {root.resolve()}\n"
            f"Tip: verify the file name, or pass the full path."
        )
    if len(matches) > 1:
        msg = "\n".join(str(m) for m in matches[:30])
        more = "" if len(matches) <= 30 else f"\n...and {len(matches) - 30} more"
        raise RuntimeError(
            f"Multiple matches found for '{p.name}'. Please pass a full path.\n"
            f"Matches:\n{msg}{more}"
        )

    return matches[0].resolve()


def base_prefix_from_archive_name(archive_path: Path) -> str:
    """
    RUNNAME__archive.bin -> RUNNAME
    (If it doesn't match, uses stem.)
    """
    name = archive_path.name
    if name.endswith("__archive.bin"):
        return name[:-len("__archive.bin")]
    if name.endswith(".bin"):
        return name[:-4]
    return archive_path.stem


def outpaths_for_videos(archive_path: Path) -> Tuple[Path, Path]:
    """
    Output files go into the same folder as the archive:
      RUNNAME__inner_xy_orbits.mp4
      RUNNAME__inner_3d_orbits.mp4
    """
    prefix = base_prefix_from_archive_name(archive_path)
    folder = archive_path.parent
    return (
        folder / f"{prefix}__inner_xy_orbits.mp4",
        folder / f"{prefix}__inner_3d_orbits.mp4",
    )


def ensure_ffmpeg_available():
    try:
        _ = FFMpegWriter(fps=30)
    except Exception as e:
        raise RuntimeError(
            "FFmpeg is required to write MP4.\n"
            "Install ffmpeg and ensure it's on your PATH, then try again."
        ) from e


def progress_print(prefix: str, i: int, n: int, t0: float, every: int = 50):
    """Print a progress line every `every` frames (and on last frame)."""
    if (i % every) != 0 and i != (n - 1):
        return
    elapsed = time.time() - t0
    frac = (i + 1) / n
    eta = elapsed * (1 / frac - 1) if frac > 1e-12 else 0.0
    print(
        f"{prefix} [{i+1:5d}/{n}] {100*frac:6.2f}% | elapsed {elapsed:7.1f}s | ETA {eta:7.1f}s",
        flush=True,
    )


def parse_csv_list(s: str) -> List[str]:
    return [x.strip() for x in s.split(",") if x.strip()]


def get_particle(sim, hash_name: Optional[str], index: Optional[int], role: str):
    """
    Hash-first lookup; index fallback.
    Note: string hashes are case-sensitive in REBOUND.
    """
    if hash_name:
        try:
            return sim.particles[hash_name]
        except Exception:
            if index is None:
                raise RuntimeError(
                    f"Could not find {role} by hash '{hash_name}', and no index fallback provided."
                )
    if index is not None:
        return sim.particles[index]
    raise RuntimeError(f"Must provide either a hash or an index for {role}.")


def heliocentric_xyz_from_particles(sun_p, p) -> Tuple[float, float, float]:
    return (p.x - sun_p.x, p.y - sun_p.y, p.z - sun_p.z)


def r_xyz(x: float, y: float, z: float) -> float:
    return math.sqrt(x * x + y * y + z * z)


def trim_trail(arr: List[float], maxlen: int):
    if maxlen is None or maxlen <= 0:
        return
    extra = len(arr) - maxlen
    if extra > 0:
        del arr[:extra]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("archive", help="Archive filename or path (e.g. RUNNAME__archive.bin)")

    # Hash-first (new archives)
    ap.add_argument("--sun-hash", default=None, help="Hash name for Sun (e.g. Sun)")
    ap.add_argument("--bh-hash", default=None, help="Hash name for BH (e.g. BH)")
    ap.add_argument(
        "--planet-hashes",
        default=",".join(PLANET_NAMES),
        help="Comma-separated planet hash names (default Mercury,Venus,Earth,Mars,Jupiter,Saturn,Uranus,Neptune)",
    )
    ap.add_argument("--include-moon", action="store_true", help="Also render Moon (hash 'Moon' if available)")

    # Index fallback (older archives)
    ap.add_argument("--sun-index", type=int, default=0, help="Index of Sun particle (default 0)")
    ap.add_argument("--bh-index", type=int, default=9, help="Index of BH particle (default 9)")
    ap.add_argument("--planet-start", type=int, default=1, help="Starting index of planets (default 1)")
    ap.add_argument("--planet-count", type=int, default=8, help="Number of planets to include (default 8)")

    # Render options
    ap.add_argument("--rmax-au", type=float, default=3.0, help="Region radius in AU (default 3)")
    ap.add_argument("--stride", type=int, default=1, help="Use every Nth snapshot (default 1)")
    ap.add_argument("--fps", type=int, default=30, help="Output video FPS (default 30)")
    ap.add_argument(
        "--trail-max",
        type=int,
        default=3000,
        help="Max stored trail points per body (default 3000; 0 = unlimited).",
    )
    ap.add_argument(
        "--progress-every",
        type=int,
        default=50,
        help="Print console progress every N frames (default 50).",
    )
    ap.add_argument("--no-labels", action="store_true", help="Disable labels")
    args = ap.parse_args()

    archive_path = resolve_archive_path(args.archive, simulations_dir="simulations")
    print(f"Using archive: {archive_path}")

    ensure_ffmpeg_available()

    sa = load_archive(str(archive_path))
    if len(sa) == 0:
        raise RuntimeError("Archive is empty.")

    # Frames
    stride = max(1, int(args.stride))
    frame_indices = list(range(0, len(sa), stride))
    total_frames = len(frame_indices)
    print(f"Snapshots in archive: {len(sa)} | Using {total_frames} frames (stride={stride})")

    out_xy, out_3d = outpaths_for_videos(archive_path)
    print(f"Will write:\n  {out_xy}\n  {out_3d}")

    # Decide which mode to use:
    # If user provided --sun-hash/--bh-hash, we assume hash mode for bodies.
    use_hash_mode = bool(args.sun_hash or args.bh_hash)

    if use_hash_mode:
        planet_hashes = parse_csv_list(args.planet_hashes)
        bodies_hash: List[Tuple[str, str]] = [(nm, nm) for nm in planet_hashes]
        if args.include_moon:
            bodies_hash.append(("Moon", "Moon"))
        # BH included last
        bh_hash = args.bh_hash or "BH"
        bodies_hash.append(("BH", bh_hash))
        bodies_index: List[Tuple[str, int]] = []
    else:
        # Index mode (legacy)
        sim0 = sa[0]
        needed_max = max(args.sun_index, args.bh_index, args.planet_start + args.planet_count - 1)
        if sim0.N <= needed_max:
            raise RuntimeError(
                f"Archive snapshots have N={sim0.N} particles, but you requested index {needed_max}. "
                "Adjust --planet-start/--planet-count/--bh-index accordingly."
            )
        planet_indices = list(range(args.planet_start, args.planet_start + args.planet_count))
        planet_labels = [
            PLANET_NAMES[j] if j < len(PLANET_NAMES) else f"Planet{j+1}"
            for j in range(args.planet_count)
        ]
        bodies_index = list(zip(planet_labels, planet_indices)) + [("BH", args.bh_index)]
        bodies_hash = []

    # Trail storage (NaNs break line segments)
    max_trail = None if args.trail_max == 0 else max(10, int(args.trail_max))
    trails_xy: Dict[str, Dict[str, List[float]]] = {}
    trails_3d: Dict[str, Dict[str, List[float]]] = {}
    body_names = [nm for (nm, _) in (bodies_hash if use_hash_mode else bodies_index)]
    for name in body_names:
        trails_xy[name] = {"x": [], "y": []}
        trails_3d[name] = {"x": [], "y": [], "z": []}

    writer = FFMpegWriter(fps=args.fps, codec="libx264", bitrate=2500)

    # -----------------
    # 2D XY animation
    # -----------------
    print("\n=== Rendering XY video ===")
    t0_xy = time.time()

    fig2, ax2 = plt.subplots(figsize=(8, 8))
    ax2.set_aspect("equal", adjustable="box")
    ax2.set_xlim(-args.rmax_au, args.rmax_au)
    ax2.set_ylim(-args.rmax_au, args.rmax_au)
    ax2.set_xlabel("x (AU, heliocentric)")
    ax2.set_ylabel("y (AU, heliocentric)")

    # Sun marker (origin)
    ax2.scatter([0], [0], s=180, marker="*", zorder=5)

    line_artists_2d = {}
    point_artists_2d = {}
    text_artists_2d = {}

    for name in body_names:
        (ln,) = ax2.plot([], [], linewidth=1.5)
        pt = ax2.scatter([], [], s=50, marker="o", zorder=6)
        line_artists_2d[name] = ln
        point_artists_2d[name] = pt
        if not args.no_labels:
            text_artists_2d[name] = ax2.text(0, 0, "", fontsize=9)

    def init_xy():
        for name in body_names:
            line_artists_2d[name].set_data([], [])
            point_artists_2d[name].set_offsets(np.empty((0, 2)))
            if not args.no_labels:
                text_artists_2d[name].set_text("")
        return list(line_artists_2d.values()) + list(point_artists_2d.values()) + list(text_artists_2d.values())

    def update_xy(frame_idx_pos: int):
        k = frame_indices[frame_idx_pos]
        sim = sa[k]

        # sun particle
        sun = get_particle(sim, args.sun_hash or "Sun", args.sun_index, "Sun") if use_hash_mode else sim.particles[args.sun_index]

        if frame_idx_pos == 0:
            # quick visibility debug (first 20 particles by distance)
            print("\n--- first-frame heliocentric distances (first 20 indices) ---")
            for idx in range(min(20, sim.N)):
                p = sim.particles[idx]
                dx, dy, dz = p.x - sun.x, p.y - sun.y, p.z - sun.z
                rr = (dx*dx + dy*dy + dz*dz) ** 0.5
                print(f"idx {idx:2d}  r(AU)={rr:.6f}")
            print("----------------------------------------------------------\n", flush=True)

        if use_hash_mode:
            iterator = bodies_hash
            for name, h in iterator:
                p = get_particle(sim, h, None, role=name)
                x, y, z = heliocentric_xyz_from_particles(sun, p)
                rr = r_xyz(x, y, z)
                if rr <= args.rmax_au:
                    trails_xy[name]["x"].append(x)
                    trails_xy[name]["y"].append(y)
                else:
                    trails_xy[name]["x"].append(np.nan)
                    trails_xy[name]["y"].append(np.nan)

                if max_trail is not None:
                    trim_trail(trails_xy[name]["x"], max_trail)
                    trim_trail(trails_xy[name]["y"], max_trail)

                line_artists_2d[name].set_data(trails_xy[name]["x"], trails_xy[name]["y"])

                if rr <= args.rmax_au:
                    point_artists_2d[name].set_offsets(np.array([[x, y]]))
                    if not args.no_labels:
                        text_artists_2d[name].set_position((x, y))
                        text_artists_2d[name].set_text(f" {name}")
                else:
                    point_artists_2d[name].set_offsets(np.empty((0, 2)))
                    if not args.no_labels:
                        text_artists_2d[name].set_text("")
        else:
            for name, i in bodies_index:
                p = sim.particles[i]
                x, y, z = heliocentric_xyz_from_particles(sun, p)
                rr = r_xyz(x, y, z)

                if rr <= args.rmax_au:
                    trails_xy[name]["x"].append(x)
                    trails_xy[name]["y"].append(y)
                else:
                    trails_xy[name]["x"].append(np.nan)
                    trails_xy[name]["y"].append(np.nan)

                if max_trail is not None:
                    trim_trail(trails_xy[name]["x"], max_trail)
                    trim_trail(trails_xy[name]["y"], max_trail)

                line_artists_2d[name].set_data(trails_xy[name]["x"], trails_xy[name]["y"])

                if rr <= args.rmax_au:
                    point_artists_2d[name].set_offsets(np.array([[x, y]]))
                    if not args.no_labels:
                        text_artists_2d[name].set_position((x, y))
                        text_artists_2d[name].set_text(f" {name}")
                else:
                    point_artists_2d[name].set_offsets(np.empty((0, 2)))
                    if not args.no_labels:
                        text_artists_2d[name].set_text("")

        ax2.set_title(
            f"Inner <= {args.rmax_au:g} AU | frame {frame_idx_pos+1}/{total_frames} | "
            f"snapshot {k} | t={sim.t:.6g} days"
        )
        progress_print("XY ", frame_idx_pos, total_frames, t0_xy, every=max(1, args.progress_every))
        return list(line_artists_2d.values()) + list(point_artists_2d.values()) + list(text_artists_2d.values())

    anim2 = FuncAnimation(
        fig2,
        update_xy,
        frames=total_frames,
        init_func=init_xy,
        blit=False,
        interval=1000 / max(1, args.fps),
    )

    anim2.save(str(out_xy), writer=writer)
    plt.close(fig2)
    print(f"Wrote: {out_xy}")
    print(f"XY done in {time.time() - t0_xy:.1f}s")

    # -----------------
    # 3D animation
    # -----------------
    print("\n=== Rendering 3D video ===")
    t0_3d = time.time()
    from mpl_toolkits.mplot3d import Axes3D  # noqa: F401

    fig3 = plt.figure(figsize=(8, 8))
    ax3 = fig3.add_subplot(111, projection="3d")
    ax3.set_xlim(-args.rmax_au, args.rmax_au)
    ax3.set_ylim(-args.rmax_au, args.rmax_au)
    ax3.set_zlim(-args.rmax_au, args.rmax_au)
    ax3.set_xlabel("x (AU)")
    ax3.set_ylabel("y (AU)")
    ax3.set_zlabel("z (AU)")

    ax3.scatter([0], [0], [0], s=180, marker="*")  # Sun

    line_artists_3d = {}
    point_artists_3d = {}
    text_artists_3d = {}

    for name in body_names:
        ln, = ax3.plot([], [], [], linewidth=1.2)
        pt = ax3.scatter([], [], [], s=45, marker="o")
        line_artists_3d[name] = ln
        point_artists_3d[name] = pt
        if not args.no_labels:
            text_artists_3d[name] = ax3.text(0, 0, 0, "", fontsize=9)

    def init_3d():
        for name in body_names:
            line_artists_3d[name].set_data([], [])
            line_artists_3d[name].set_3d_properties([])
            point_artists_3d[name]._offsets3d = ([], [], [])
            if not args.no_labels:
                text_artists_3d[name].set_text("")
        return list(line_artists_3d.values()) + list(point_artists_3d.values()) + list(text_artists_3d.values())

    def update_3d(frame_idx_pos: int):
        k = frame_indices[frame_idx_pos]
        sim = sa[k]
        sun = get_particle(sim, args.sun_hash or "Sun", args.sun_index, "Sun") if use_hash_mode else sim.particles[args.sun_index]

        if use_hash_mode:
            for name, h in bodies_hash:
                p = get_particle(sim, h, None, role=name)
                x, y, z = heliocentric_xyz_from_particles(sun, p)
                rr = r_xyz(x, y, z)

                if rr <= args.rmax_au:
                    trails_3d[name]["x"].append(x)
                    trails_3d[name]["y"].append(y)
                    trails_3d[name]["z"].append(z)
                else:
                    trails_3d[name]["x"].append(np.nan)
                    trails_3d[name]["y"].append(np.nan)
                    trails_3d[name]["z"].append(np.nan)

                if max_trail is not None:
                    trim_trail(trails_3d[name]["x"], max_trail)
                    trim_trail(trails_3d[name]["y"], max_trail)
                    trim_trail(trails_3d[name]["z"], max_trail)

                xs, ys, zs = trails_3d[name]["x"], trails_3d[name]["y"], trails_3d[name]["z"]
                line_artists_3d[name].set_data(xs, ys)
                line_artists_3d[name].set_3d_properties(zs)

                if rr <= args.rmax_au:
                    point_artists_3d[name]._offsets3d = ([x], [y], [z])
                    if not args.no_labels:
                        text_artists_3d[name].set_position((x, y))
                        text_artists_3d[name].set_3d_properties(z, zdir="z")
                        text_artists_3d[name].set_text(f" {name}")
                else:
                    point_artists_3d[name]._offsets3d = ([], [], [])
                    if not args.no_labels:
                        text_artists_3d[name].set_text("")
        else:
            for name, i in bodies_index:
                p = sim.particles[i]
                x, y, z = heliocentric_xyz_from_particles(sun, p)
                rr = r_xyz(x, y, z)

                if rr <= args.rmax_au:
                    trails_3d[name]["x"].append(x)
                    trails_3d[name]["y"].append(y)
                    trails_3d[name]["z"].append(z)
                else:
                    trails_3d[name]["x"].append(np.nan)
                    trails_3d[name]["y"].append(np.nan)
                    trails_3d[name]["z"].append(np.nan)

                if max_trail is not None:
                    trim_trail(trails_3d[name]["x"], max_trail)
                    trim_trail(trails_3d[name]["y"], max_trail)
                    trim_trail(trails_3d[name]["z"], max_trail)

                xs, ys, zs = trails_3d[name]["x"], trails_3d[name]["y"], trails_3d[name]["z"]
                line_artists_3d[name].set_data(xs, ys)
                line_artists_3d[name].set_3d_properties(zs)

                if rr <= args.rmax_au:
                    point_artists_3d[name]._offsets3d = ([x], [y], [z])
                    if not args.no_labels:
                        text_artists_3d[name].set_position((x, y))
                        text_artists_3d[name].set_3d_properties(z, zdir="z")
                        text_artists_3d[name].set_text(f" {name}")
                else:
                    point_artists_3d[name]._offsets3d = ([], [], [])
                    if not args.no_labels:
                        text_artists_3d[name].set_text("")

        ax3.set_title(
            f"Inner <= {args.rmax_au:g} AU | frame {frame_idx_pos+1}/{total_frames} | "
            f"snapshot {k} | t={sim.t:.6g} days"
        )
        progress_print("3D ", frame_idx_pos, total_frames, t0_3d, every=max(1, args.progress_every))
        return list(line_artists_3d.values()) + list(point_artists_3d.values()) + list(text_artists_3d.values())

    anim3 = FuncAnimation(
        fig3,
        update_3d,
        frames=total_frames,
        init_func=init_3d,
        blit=False,
        interval=1000 / max(1, args.fps),
    )

    anim3.save(str(out_3d), writer=writer)
    plt.close(fig3)
    print(f"Wrote: {out_3d}")
    print(f"3D done in {time.time() - t0_3d:.1f}s")


if __name__ == "__main__":
    main()
