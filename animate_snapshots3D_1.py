#!/usr/bin/env python3
"""
3D animation of Solar System + BH + belt from snapshot_*.npz.

- Assumes NPZ format from solar_system_bh_reboundXX.py:
    arr: (N, 6) -> x,y,z,vx,vy,vz in AU and AU/day
    names: array of strings for the first 11 bodies (Sun..BH), or empty.

- First 11 particles are:
    0: Sun
    1: Mercury
    2: Venus
    3: Earth
    4: Mars
    5: Jupiter
    6: Saturn
    7: Uranus
    8: Neptune
    9: Moon
    10: BH
  Remaining particles are belt tracers.

Usage:
    python animate_snapshots3D.py path/to/run_dir --mp4 orbit3D.mp4
    python animate_snapshots3D.py path/to/run_dir --gif orbit3D.gif
"""

import argparse
from pathlib import Path
import re

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation, FFMpegWriter, PillowWriter
from mpl_toolkits.mplot3d import Axes3D  # noqa: F401 (needed for 3D)


# ------------------------------------------------------------
# Configuration: body names, colors, sizes
# ------------------------------------------------------------
BODY_NAMES = [
    "Sun", "Mercury", "Venus", "Earth", "Mars",
    "Jupiter", "Saturn", "Uranus", "Neptune", "Moon", "BH"
]

# Distinct, readable colors for the 11 main bodies
BODY_COLORS = {
    "Sun":     "gold",
    "Mercury": "red",
    "Venus":   "limegreen",
    "Earth":   "dodgerblue",
    "Mars":    "orangered",
    "Jupiter": "sandybrown",
    "Saturn":  "khaki",
    "Uranus":  "turquoise",
    "Neptune": "mediumblue",
    "Moon":    "lightgray",
    "BH":      "black",
}

# Visually nice sizes for plotting (these are marker areas, not radii)
BODY_SIZES = {
    "Sun":     80.0,
    "Mercury": 20.0,
    "Venus":   26.0,
    "Earth":   28.0,
    "Mars":    24.0,
    "Jupiter": 40.0,
    "Saturn":  36.0,
    "Uranus":  32.0,
    "Neptune": 32.0,
    "Moon":    18.0,
    "BH":      50.0,
}
BELT_SIZE = 4.0
BELT_COLOR = "darkgray"
BELT_ALPHA = 0.4

# ------------------------------------------------------------
# Helpers
# ------------------------------------------------------------

def infer_time_days(path: Path) -> float:
    """Extract integer time in days from 'snapshot_XXXXXX' in filename."""
    m = re.search(r"snapshot_(\d+)", path.name)
    if m:
        return float(m.group(1))
    return 0.0


def load_frames(run_dir: Path, stride: int = 1):
    """
    Load all snapshot_*.npz in run_dir with given stride.

    Returns:
        times: list[float]
        frames: list[np.ndarray], each shape (N, 3), heliocentric (Sun at ~0)
        n_main: 11 (main bodies)
    """
    snaps = sorted(run_dir.glob("*snapshot_*.npz"))
    if not snaps:
        raise FileNotFoundError(f"No snapshot_*.npz found in {run_dir}")

    times = []
    frames = []

    for snap in snaps[::stride]:
        data = np.load(snap, allow_pickle=True)
        arr = data["arr"]              # (N, 6)
        pos = arr[:, :3]               # (N, 3) -> x,y,z in AU

        # Use particle 0 (Sun) as origin → heliocentric coords
        sun_pos = pos[0].copy()
        pos_rel = pos - sun_pos        # (N,3)

        t = infer_time_days(snap)
        times.append(t)
        frames.append(pos_rel)

    n_main = min(11, frames[0].shape[0])
    return times, frames, n_main


# ------------------------------------------------------------
# Animation setup
# ------------------------------------------------------------

def make_animation(run_dir: Path,
                   output: Path | None = None,
                   make_gif: bool = False,
                   stride: int = 1,
                   fps: int = 30):
    times, frames, n_main = load_frames(run_dir, stride=stride)
    n_frames = len(frames)
    print(f"Loaded {n_frames} frames from {run_dir} (stride={stride}).")
    print(f"Each frame has {frames[0].shape[0]} particles (11 main + belt).")

    # --- create figure ---
    fig = plt.figure(figsize=(8, 8))
    ax = fig.add_subplot(111, projection="3d")

    # Limits: try to choose something decent from all frames
    all_pos = np.vstack(frames)  # (n_frames*N, 3)
    max_range = np.max(np.linalg.norm(all_pos[:, :2], axis=1))
    if max_range == 0:
        max_range = 1.0
    lim = max_range * 1.1

    ax.set_xlim(-lim, lim)
    ax.set_ylim(-lim, lim)
    ax.set_zlim(-lim, lim)

    ax.set_xlabel("x [AU]")
    ax.set_ylabel("y [AU]")
    ax.set_zlabel("z [AU]")
    ax.set_title("Solar System + BH flyby (heliocentric 3D)")

    # We’ll create one scatter per main body, plus one for belt
    main_scatters = []
    for i in range(n_main):
        name = BODY_NAMES[i] if i < len(BODY_NAMES) else f"Body{i}"
        color = BODY_COLORS.get(name, "white")
        size = BODY_SIZES.get(name, 30.0)
        scat = ax.scatter([], [], [], s=size, color=color, label=name)
        main_scatters.append(scat)

    # Belt tracers
    belt_scatter = ax.scatter([], [], [], s=BELT_SIZE,
                              color=BELT_COLOR, alpha=BELT_ALPHA,
                              linewidths=0)

    # Legend for main bodies
    ax.legend(loc="upper right", fontsize=8)

    # Text overlay for time
    time_text = ax.text2D(0.02, 0.95, "", transform=ax.transAxes)

    # --- init and update functions ---

    def init():
        for scat in main_scatters:
            scat._offsets3d = ([], [], [])
        belt_scatter._offsets3d = ([], [], [])
        time_text.set_text("")
        return main_scatters + [belt_scatter, time_text]

    def update(frame_idx):
        pos = frames[frame_idx]
        N = pos.shape[0]

        # main bodies
        for i, scat in enumerate(main_scatters):
            if i < N:
                x = pos[i, 0]
                y = pos[i, 1]
                z = pos[i, 2]
                scat._offsets3d = ([x], [y], [z])
            else:
                scat._offsets3d = ([], [], [])

        # belt
        if N > n_main:
            belt = pos[n_main:, :]
            belt_scatter._offsets3d = (belt[:, 0], belt[:, 1], belt[:, 2])
        else:
            belt_scatter._offsets3d = ([], [], [])

        time_text.set_text(f"t = {times[frame_idx]:.1f} days")
        return main_scatters + [belt_scatter, time_text]

    anim = FuncAnimation(
        fig, update, init_func=init,
        frames=n_frames, interval=1000.0 / fps, blit=False
    )

    if output is None:
        # Just show interactively
        plt.show()
    else:
        output = Path(output)
        if make_gif:
            print(f"Saving GIF to {output} ...")
            writer = PillowWriter(fps=fps)
            anim.save(str(output), writer=writer)
        else:
            print(f"Saving MP4 to {output} ...")
            writer = FFMpegWriter(fps=fps, bitrate=2000)
            anim.save(str(output), writer=writer)
        print("Done.")
    plt.close(fig)


# ------------------------------------------------------------
# CLI
# ------------------------------------------------------------

def main():
    ap = argparse.ArgumentParser(
        description="Animate 3D Solar System + BH snapshots (NPZ)."
    )
    ap.add_argument(
        "run_dir",
        help="Directory containing snapshot_*.npz (one simulation run)."
    )
    ap.add_argument(
        "--mp4",
        metavar="FILE",
        help="Write MP4 via ffmpeg to FILE (e.g. orbit3D.mp4)."
    )
    ap.add_argument(
        "--gif",
        metavar="FILE",
        help="Write GIF via Pillow to FILE (e.g. orbit3D.gif)."
    )
    ap.add_argument(
        "--stride",
        type=int,
        default=1,
        help="Use every Nth snapshot (default=1)."
    )
    ap.add_argument(
        "--fps",
        type=int,
        default=30,
        help="Frames per second for output (default=30)."
    )
    args = ap.parse_args()

    run_dir = Path(args.run_dir)
    if not run_dir.is_dir():
        raise SystemExit(f"{run_dir} is not a directory.")

    if args.gif and args.mp4:
        raise SystemExit("Choose either --gif OR --mp4, not both.")

    if args.gif:
        make_animation(run_dir, output=Path(args.gif), make_gif=True,
                       stride=args.stride, fps=args.fps)
    elif args.mp4:
        make_animation(run_dir, output=Path(args.mp4), make_gif=False,
                       stride=args.stride, fps=args.fps)
    else:
        # No output -> just show interactively
        make_animation(run_dir, output=None,
                       stride=args.stride, fps=args.fps)


if __name__ == "__main__":
    main()
