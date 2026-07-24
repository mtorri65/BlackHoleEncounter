#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
reconstruct_trajectories_from_archive.py

Reconstruct heliocentric body trajectories from REBOUND SimulationArchive
files and browse the resulting images interactively.

Given a folder, the script decides between two modes:

  - Single-run mode: an archive (*__archive.bin / *.bin) sits directly
    inside the given folder -> exactly one image is generated.
  - Batch mode: no archive sits directly inside the given folder, so each
    immediate subfolder that contains its own archive is treated as one
    run -> one image per run subfolder.

For each run, all snapshots in its SimulationArchive are walked once to
build heliocentric (Sun-relative) x/y trajectories for the requested
bodies, which are then rendered to a PNG at the requested AU field of
view. A parameter box on each image reports the six BH encounter
parameters for that run (bh_rp_au, bh_vinf_kms, bh_inc_deg,
bh_tperi_offset_days, bh_Omega_deg, bh_omega_deg), read from the run's
subfolder-name tokens (the same "__rp10__vinf30__inc20__toff2920__Om90__om180"
convention solar_system_bh_rebound26.py encodes them with), falling
back to the run's copied input.yaml if the folder name can't be parsed.

Once generated, the images are displayed one at a time in a matplotlib
window with keyboard navigation.

Usage
-----
    python reconstruct_trajectories_from_archive.py simulations/20260314_234712 --fov-x 40 --fov-y 40
    python reconstruct_trajectories_from_archive.py simulations/20260314_234712/20260314_234712__rp10__vinf30__inc20__toff2920__Om90__om180 --fov 10
    python reconstruct_trajectories_from_archive.py simulations/20260314_234712 --bodies Mercury,Venus,Earth,Mars,BH

Controls (image viewer)
------------------------
Left / Right : previous / next image
R            : prompt for new fov_x / fov_y and regenerate all images
H            : toggle help banner
Q or Esc     : quit
"""

import argparse
import re
import sys
from pathlib import Path
from typing import Optional

import matplotlib
# Qt's screen-enumeration can fail with "Unable to open monitor interface"
# on some Windows VDI/virtual-desktop setups even though a GUI is otherwise
# usable. Tk is more tolerant of that, and ships with a standard Python
# install, so prefer it before pyplot picks a default backend.
try:
    matplotlib.use("TkAgg")
except Exception:
    pass

import matplotlib.pyplot as plt
import numpy as np
import yaml
import rebound

# ---------------------------------------------------------------------------
# Constants
# ---------------------------------------------------------------------------

BODY_NAMES = [
    "Sun", "Mercury", "Venus", "Earth", "Mars",
    "Jupiter", "Saturn", "Uranus", "Neptune", "Moon", "BH",
]
BODY_INDEX_FALLBACK = {name: i for i, name in enumerate(BODY_NAMES)}

DEFAULT_BODIES = [
    "Mercury", "Venus", "Earth", "Mars",
    "Jupiter", "Saturn", "Uranus", "Neptune", "Moon", "BH",
]

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
    "BH": "black",
}

BH_PARAM_KEYS = [
    ("rp", "bh_rp_au"),
    ("vinf", "bh_vinf_kms"),
    ("inc", "bh_inc_deg"),
    ("toff", "bh_tperi_offset_days"),
    ("Om", "bh_Omega_deg"),
    ("om", "bh_omega_deg"),
]

FOLDER_TOKEN_RE = re.compile(
    r"__rp(?P<rp>[^_]+)"
    r"__vinf(?P<vinf>[^_]+)"
    r"__inc(?P<inc>[^_]+)"
    r"__toff(?P<toff>[^_]+)"
    r"__Om(?P<Om>[^_]+)"
    r"__om(?P<om>[^_]+)"
)

HELP_TEXT = "Left/Right: prev/next   R: regenerate at new FOV   H: help   Q/Esc: quit"

# ---------------------------------------------------------------------------
# Archive helpers (hash-first, index-fallback -- same pattern used elsewhere
# in this project, e.g. plot_bh_perihelion.py / make_inner_3au_videos_from_archive.py)
# ---------------------------------------------------------------------------

def load_archive(path: str):
    """Version-robust SimulationArchive loader."""
    if hasattr(rebound, "SimulationArchive"):
        return rebound.SimulationArchive(path)
    if hasattr(rebound, "Simulationarchive"):
        return rebound.Simulationarchive(path)
    raise RuntimeError("This REBOUND install does not expose SimulationArchive/Simulationarchive.")


def archive_uses_hashes(sim) -> bool:
    """
    Decide once per archive whether particle hashes are usable at all, by
    checking if "Sun" resolves. Mixing hash lookup and index-fallback on a
    per-body basis is unsafe: if hashes work for some bodies but not others,
    a missing hash almost certainly means that body isn't in this particular
    archive (e.g. a stripped-down/custom archive) rather than that it needs
    the canonical index -- falling back to the canonical index in that case
    would silently plot the wrong particle under the wrong name.
    """
    try:
        _ = sim.particles["Sun"]
        return True
    except Exception:
        return False


def get_particle(sim, name: str, use_hash: bool):
    """
    use_hash=True:  hash-only lookup (raises if this specific body isn't
                    present in the archive -- it is NOT guessed by index).
    use_hash=False: index-only lookup via the canonical BODY_NAMES order,
                    for genuinely anonymous/old-style archives with no
                    hashes at all.
    """
    if use_hash:
        return sim.particles[name]
    idx = BODY_INDEX_FALLBACK.get(name)
    if idx is None or idx >= sim.N:
        raise KeyError(name)
    return sim.particles[idx]


def find_archive_in_dir(d: Path) -> Optional[Path]:
    for pattern in ("*__archive.bin", "*archive*.bin", "*.bin"):
        matches = sorted(d.glob(pattern))
        if matches:
            if len(matches) > 1:
                print(f"[warn] Multiple archive candidates in {d}, using: {matches[0].name}")
            return matches[0]
    return None


def discover_runs(root: Path):
    """
    Returns a list of (run_dir, archive_path) pairs.

    If an archive sits directly in `root`, root is treated as a single run.
    Otherwise, every immediate subfolder of `root` that contains its own
    archive is treated as one run.
    """
    direct = find_archive_in_dir(root)
    if direct is not None:
        return [(root, direct)]

    runs = []
    for sub in sorted(p for p in root.iterdir() if p.is_dir()):
        arch = find_archive_in_dir(sub)
        if arch is None:
            continue
        runs.append((sub, arch))

    if not runs:
        raise FileNotFoundError(
            f"No SimulationArchive (*.bin) found directly in {root} "
            f"or in any of its immediate subfolders."
        )
    return runs


# ---------------------------------------------------------------------------
# BH parameter extraction
# ---------------------------------------------------------------------------

def _decode_token(tok: str) -> float:
    return float(tok.replace("m", "-").replace("p", "."))


def extract_bh_params(run_dir: Path) -> dict:
    """
    Resolve the six BH encounter parameters for one run.

    Primary source: the run subfolder name tokens (rp/vinf/inc/toff/Om/om),
    which the core engine writes with the *resolved* per-run value. The
    copied "__input.yaml" in a run folder is a verbatim copy of the master
    sweep config and therefore may still contain range strings -- it is
    only used as a best-effort fallback when the folder name can't be
    parsed (e.g. a single-value, non-swept run).
    """
    m = FOLDER_TOKEN_RE.search(run_dir.name)
    if m:
        return {label: _decode_token(m.group(tok)) for tok, label in BH_PARAM_KEYS}

    params = {label: None for _, label in BH_PARAM_KEYS}
    yaml_candidates = list(run_dir.glob("*__input.yaml")) + list(run_dir.glob("input.yaml"))
    if yaml_candidates:
        try:
            with open(yaml_candidates[0], "r") as f:
                cfg = yaml.safe_load(f) or {}
            for _, label in BH_PARAM_KEYS:
                v = cfg.get(label)
                if isinstance(v, (int, float)):
                    params[label] = float(v)
                elif isinstance(v, str) and "," not in v:
                    try:
                        params[label] = float(v)
                    except ValueError:
                        pass
        except Exception as e:
            print(f"[warn] Could not read {yaml_candidates[0]}: {e}")
    return params


def format_bh_params(params: dict) -> str:
    lines = []
    for _, label in BH_PARAM_KEYS:
        v = params.get(label)
        lines.append(f"{label}: {v:g}" if v is not None else f"{label}: N/A")
    return "\n".join(lines)


# ---------------------------------------------------------------------------
# Trajectory reconstruction (single pass over the archive)
# ---------------------------------------------------------------------------

def reconstruct_tracks(archive_path: Path, bodies, show_progress=True):
    """
    Walk every snapshot in the archive once, returning a dict:
        body_name -> (x_rel_AU ndarray, y_rel_AU ndarray)
    positions are heliocentric (Sun subtracted each snapshot).
    """
    sa = load_archive(str(archive_path))
    n = len(sa)
    if n == 0:
        raise ValueError(f"Archive is empty: {archive_path}")

    wanted = [b for b in bodies if b != "Sun"]
    raw = {b: {"x": [], "y": []} for b in wanted}
    missing_warned = set()

    use_hash = archive_uses_hashes(sa[0])
    print(f"  [{archive_path.name}] particle lookup mode: {'hash' if use_hash else 'index (anonymous archive)'}")

    for k in range(n):
        if show_progress and (k % max(1, n // 20) == 0 or k == n - 1):
            print(f"\r  [{archive_path.name}] snapshot {k + 1}/{n}", end="", flush=True)
        sim = sa[k]
        try:
            sun = get_particle(sim, "Sun", use_hash)
        except Exception:
            continue
        for b in wanted:
            try:
                p = get_particle(sim, b, use_hash)
            except Exception:
                if b not in missing_warned:
                    print(f"\n[warn] Body '{b}' not found in {archive_path.name}; skipping it for this run.")
                    missing_warned.add(b)
                continue
            raw[b]["x"].append(p.x - sun.x)
            raw[b]["y"].append(p.y - sun.y)

    if show_progress:
        print()

    tracks = {}
    for b, d in raw.items():
        if b in missing_warned:
            continue
        tracks[b] = (np.array(d["x"]), np.array(d["y"]))
    return tracks


# ---------------------------------------------------------------------------
# Rendering
# ---------------------------------------------------------------------------

def infer_image_path(archive_path: Path) -> Path:
    name = archive_path.name
    if name.endswith("__archive.bin"):
        base = name[: -len("__archive.bin")]
    elif name.endswith(".bin"):
        base = name[:-4]
    else:
        base = archive_path.stem
    return archive_path.parent / f"{base}__reconstructed_trajectories.png"


def render_run_image(tracks: dict, bh_params: dict, bodies, fov_x: float, fov_y: float, out_path: Path):
    fig, ax = plt.subplots(figsize=(8, 8))
    ax.set_aspect("equal", adjustable="box")
    ax.grid(True, alpha=0.3)
    ax.set_xlabel("x [AU] heliocentric")
    ax.set_ylabel("y [AU] heliocentric")
    ax.set_title(f"FOV = {fov_x:g} AU x {fov_y:g} AU", fontsize=10)

    ax.scatter([0], [0], s=70, c="gold", edgecolor="k", linewidth=0.5, label="Sun", zorder=5)

    for b in bodies:
        if b == "Sun" or b not in tracks:
            continue
        x, y = tracks[b]
        if x.size == 0:
            continue
        color = COLOR_MAP.get(b)
        lw = 1.7 if b == "BH" else 1.0
        marker = "x" if b == "BH" else "o"
        ax.plot(x, y, lw=lw, color=color, label=b, alpha=0.95)
        ax.scatter([x[-1]], [y[-1]], s=30 if b == "BH" else 12, color=color, marker=marker, zorder=6)

    ax.set_xlim(-fov_x / 2.0, fov_x / 2.0)
    ax.set_ylim(-fov_y / 2.0, fov_y / 2.0)
    ax.legend(fontsize=8, ncols=2, loc="upper right")

    ax.text(
        0.02, 0.98,
        format_bh_params(bh_params),
        transform=ax.transAxes,
        va="top", ha="left", fontsize=9, family="monospace",
        bbox=dict(boxstyle="round", facecolor="white", edgecolor="gray", alpha=0.85),
    )

    fig.tight_layout()
    fig.savefig(out_path, dpi=180)
    plt.close(fig)


# ---------------------------------------------------------------------------
# Interactive viewer
# ---------------------------------------------------------------------------

class TrajectoryViewer:
    def __init__(self, runs, bodies, fov_x, fov_y, show_progress=True):
        """
        runs: list of dicts with keys: run_dir, archive_path, tracks, bh_params, image_path
        """
        self.runs = runs
        self.bodies = bodies
        self.fov_x = float(fov_x)
        self.fov_y = float(fov_y)
        self.i = 0
        self.show_help = True
        self.show_progress = show_progress

        self.fig, self.ax = plt.subplots(figsize=(9, 9))
        self.cid = self.fig.canvas.mpl_connect("key_press_event", self._on_key)
        try:
            self.fig.canvas.manager.set_window_title("Reconstructed Trajectories")
        except Exception:
            pass

        self.generate_all(self.fov_x, self.fov_y)
        self._draw()

    def generate_all(self, fov_x, fov_y):
        total = len(self.runs)
        print(f"Generating {total} image(s) at FOV {fov_x:g} x {fov_y:g} AU...")
        for idx, run in enumerate(self.runs, start=1):
            print(f"  [{idx}/{total}] Rendering folder: {run['run_dir'].name}")
            render_run_image(
                run["tracks"], run["bh_params"],
                self.bodies, fov_x, fov_y, run["image_path"],
            )
            print(f"    Saved: {run['image_path']}")
        self.fov_x, self.fov_y = fov_x, fov_y

    def _draw(self):
        self.ax.clear()
        self.ax.axis("off")
        run = self.runs[self.i]
        img = plt.imread(run["image_path"])
        self.ax.imshow(img)
        self.ax.set_title(f"[{self.i + 1}/{len(self.runs)}] {run['run_dir'].name}", fontsize=10)
        if self.show_help:
            self.ax.text(
                0.5, -0.02, HELP_TEXT, transform=self.ax.transAxes,
                ha="center", va="top", fontsize=9,
            )
        self.fig.tight_layout()
        self.fig.canvas.draw_idle()

    def _prompt_regenerate(self):
        print(f"\nCurrent FOV: {self.fov_x:g} x {self.fov_y:g} AU")
        try:
            raw_x = input(f"New fov_x (AU) [blank = keep {self.fov_x:g}]: ").strip()
            raw_y = input(f"New fov_y (AU) [blank = keep {self.fov_y:g}]: ").strip()
        except EOFError:
            print("No console input available; regeneration cancelled.")
            return
        try:
            new_x = float(raw_x) if raw_x else self.fov_x
            new_y = float(raw_y) if raw_y else self.fov_y
        except ValueError:
            print("Invalid number entered; regeneration cancelled.")
            return
        if new_x <= 0 or new_y <= 0:
            print("FOV values must be positive; regeneration cancelled.")
            return
        self.generate_all(new_x, new_y)
        self._draw()

    def _on_key(self, event):
        key = (event.key or "").lower()
        if key == "right":
            self.i = (self.i + 1) % len(self.runs)
            self._draw()
        elif key == "left":
            self.i = (self.i - 1) % len(self.runs)
            self._draw()
        elif key == "r":
            self._prompt_regenerate()
        elif key == "h":
            self.show_help = not self.show_help
            self._draw()
        elif key in ("escape", "q"):
            plt.close(self.fig)


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

def parse_bodies(s: str):
    if not s:
        return list(DEFAULT_BODIES)
    return [b.strip() for b in s.split(",") if b.strip()]


def main():
    ap = argparse.ArgumentParser(
        description="Reconstruct body trajectories from SimulationArchive file(s) and browse the resulting images.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    ap.add_argument("folder", help="Run folder (holds an archive directly) or parent folder (holds run subfolders).")
    ap.add_argument("--fov", type=float, default=None, help="Square field of view in AU. Sets both --fov-x and --fov-y.")
    ap.add_argument("--fov-x", type=float, default=80.0, help="Horizontal field of view in AU. Default: 80.")
    ap.add_argument("--fov-y", type=float, default=80.0, help="Vertical field of view in AU. Default: 80.")
    ap.add_argument(
        "--bodies", type=str, default=",".join(DEFAULT_BODIES),
        help="Comma-separated body list to plot. Default: all planets, Moon, BH.",
    )
    ap.add_argument("--no-progress", dest="progress", action="store_false", default=True,
                     help="Disable console progress messages while reconstructing archives.")
    args = ap.parse_args()

    fov_x = args.fov if args.fov is not None else args.fov_x
    fov_y = args.fov if args.fov is not None else args.fov_y
    if fov_x <= 0 or fov_y <= 0:
        raise ValueError("FOV values must be positive.")

    root = Path(args.folder)
    if not root.is_dir():
        raise FileNotFoundError(f"Folder not found: {root}")

    bodies = parse_bodies(args.bodies)

    run_pairs = discover_runs(root)
    print(f"Found {len(run_pairs)} run(s) under {root}")

    runs = []
    total = len(run_pairs)
    for idx, (run_dir, archive_path) in enumerate(run_pairs, start=1):
        print(f"[{idx}/{total}] Processing folder: {run_dir.name}")
        tracks = reconstruct_tracks(archive_path, bodies, show_progress=args.progress)
        bh_params = extract_bh_params(run_dir)
        runs.append({
            "run_dir": run_dir,
            "archive_path": archive_path,
            "tracks": tracks,
            "bh_params": bh_params,
            "image_path": infer_image_path(archive_path),
        })

    try:
        TrajectoryViewer(runs, bodies, fov_x, fov_y, show_progress=args.progress)
        plt.show()
    except Exception as e:
        print(f"\n[warn] Could not open the interactive viewer ({e}).")
        print("All images were still generated and saved next to each run's archive:")
        for run in runs:
            print(f"  {run['image_path']}")


if __name__ == "__main__":
    main()
