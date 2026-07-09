#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
gallery_heliocentric_replot_viewer.py

Interactive gallery for Solar System + BH simulation runs.

Unlike gallery_heliocentric_viewer.py, this script does NOT zoom into an
existing PNG. Instead, it reads each run's __orbits__*.xlsx workbook and
redraws the heliocentric plot live with a user-specified AU field of view.

This avoids pixelated zooming when the original PNG was autoscaled to include
a very distant black hole.

Usage examples
--------------
Show all runs under a sweep root, redrawing each with an 80 AU x 80 AU view:

    python gallery_heliocentric_replot_viewer.py simulations/20260314_234712 --fov-x 80 --fov-y 80

Use a 120 AU square view:

    python gallery_heliocentric_replot_viewer.py simulations/20260314_234712 --fov 120

Show only selected bodies:

    python gallery_heliocentric_replot_viewer.py simulations/20260314_234712 --fov 80 --bodies Mercury,Venus,Earth,Mars,Jupiter,Saturn,Uranus,Neptune

Controls
--------
Left / Right arrows : previous / next run
A/D                 : pan left/right by pan fraction
W/S                 : pan up/down by pan fraction
+ / -               : shrink/enlarge field of view by zoom-step
0                   : reset pan and field of view to initial values
H                   : toggle help
Q or Esc            : quit

Notes
-----
- The original PNGs and Excel files are not modified.
- The plot is Sun-centered by default.
- Field of view is in AU.
- The script searches recursively for *__orbits__*.xlsx files under the root.
"""

import argparse
import os
import glob
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt


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

HELP_TEXT = (
    "←/→: prev/next   W/A/S/D: pan   +/-: FOV in/out   "
    "0: reset view   H: help   Esc/Q: quit"
)


def _progress_line(enabled: bool, workbook: str, label: str, body: str = "", time_days=None):
    """Print a same-line progress message in the console."""
    if not enabled:
        return
    wb = str(workbook)
    if len(wb) > 48:
        wb = "..." + wb[-45:]
    if time_days is None:
        msg = f"{wb:<48} | {label:<10} | {body:<10}"
    else:
        try:
            msg = f"{wb:<48} | {label:<10} | {body:<10} | time_days = {float(time_days):12.3f}"
        except Exception:
            msg = f"{wb:<48} | {label:<10} | {body:<10} | time_days = {time_days}"
    print("\r" + msg.ljust(120), end="", flush=True)


def _progress_done(enabled: bool):
    """Finish the current same-line progress message."""
    if enabled:
        print("", flush=True)


def find_orbit_workbooks(root: str):
    patterns = [
        os.path.join(root, "**", "*__orbits__*.xlsx"),
        os.path.join(root, "**", "orbits__*.xlsx"),
    ]
    files = []
    for pattern in patterns:
        files.extend(glob.glob(pattern, recursive=True))
    # De-duplicate while preserving sorted order.
    files = sorted(set(files))
    return files


def _read_sheet(book: pd.ExcelFile, sheet_name: str):
    df = pd.read_excel(book, sheet_name=sheet_name)
    df.columns = [str(c).strip() for c in df.columns]
    return df


def _safe_array(df: pd.DataFrame, col: str):
    if col not in df.columns:
        return np.array([], dtype=float)
    return df[col].to_numpy(dtype=float)


def _load_run_data(xlsx_path: str, bodies, show_progress=False):
    """
    Load heliocentric x/y tracks from one run workbook.

    Returns
    -------
    dict
        body -> (x_rel_AU, y_rel_AU)
    """
    tracks = {}
    workbook_label = f"{Path(xlsx_path).parent.name} / {Path(xlsx_path).name}"
    _progress_line(show_progress, workbook_label, "Opening", "", "")
    with pd.ExcelFile(xlsx_path) as book:
        if "Sun" not in book.sheet_names:
            raise ValueError(f"Workbook lacks Sun sheet: {xlsx_path}")

        sun_df = _read_sheet(book, "Sun")
        sx = _safe_array(sun_df, "x_AU")
        sy = _safe_array(sun_df, "y_AU")
        if sx.size == 0 or sy.size == 0:
            raise ValueError(f"Sun sheet lacks x_AU/y_AU columns: {xlsx_path}")

        for body in bodies:
            if body not in book.sheet_names:
                continue
            _progress_line(show_progress, workbook_label, "Loading", body, "")
            df = _read_sheet(book, body)
            x = _safe_array(df, "x_AU")
            y = _safe_array(df, "y_AU")
            t = _safe_array(df, "time_days")
            n = min(len(x), len(y), len(sx), len(sy))
            if n == 0:
                continue
            if len(t) >= n:
                _progress_line(show_progress, workbook_label, "Loaded", body, t[n - 1])
            tracks[body] = (x[:n] - sx[:n], y[:n] - sy[:n])

    return tracks


def _short_title(path: str):
    p = Path(path)
    return f"{p.parent.name} / {p.name}"


class ReplotViewer:
    def __init__(
        self,
        files,
        bodies,
        fov_x,
        fov_y,
        pan_fraction=0.15,
        zoom_step=1.25,
        show_legend=True,
        show_bh=True,
        linewidth=1.0,
        show_progress=True,
    ):
        if not files:
            raise FileNotFoundError("No *__orbits__*.xlsx files found.")

        self.files = files
        self.n = len(files)
        self.i = 0
        self.bodies = list(bodies)

        if not show_bh and "BH" in self.bodies:
            self.bodies.remove("BH")

        self.initial_fov_x = float(fov_x)
        self.initial_fov_y = float(fov_y)
        self.fov_x = float(fov_x)
        self.fov_y = float(fov_y)
        self.cx = 0.0
        self.cy = 0.0

        self.pan_fraction = float(pan_fraction)
        self.zoom_step = float(zoom_step)
        self.show_help = True
        self.show_legend = bool(show_legend)
        self.linewidth = float(linewidth)
        self.show_progress = bool(show_progress)

        self.fig, self.ax = plt.subplots(figsize=(9, 9))
        self.cid = self.fig.canvas.mpl_connect("key_press_event", self._on_key)

        try:
            self.fig.canvas.manager.set_window_title("Heliocentric Replot Gallery")
        except Exception:
            pass

        self._draw()

    def _reset_view(self):
        self.fov_x = self.initial_fov_x
        self.fov_y = self.initial_fov_y
        self.cx = 0.0
        self.cy = 0.0

    def _draw(self):
        self.ax.clear()
        self.ax.set_aspect("equal", adjustable="box")
        self.ax.grid(True, alpha=0.3)
        self.ax.set_xlabel("x [AU] heliocentric")
        self.ax.set_ylabel("y [AU] heliocentric")

        path = self.files[self.i]
        title = f"[{self.i + 1}/{self.n}] {_short_title(path)}\nFOV = {self.fov_x:g} AU × {self.fov_y:g} AU, center=({self.cx:g}, {self.cy:g}) AU"
        self.ax.set_title(title, fontsize=10)

        try:
            tracks = _load_run_data(path, self.bodies, show_progress=self.show_progress)
        except Exception as e:
            self.ax.text(
                0.5, 0.5,
                f"Could not load workbook:\n{path}\n\n{e}",
                transform=self.ax.transAxes,
                ha="center",
                va="center",
                fontsize=10,
                bbox=dict(facecolor="white", alpha=0.9, edgecolor="gray"),
            )
            self.fig.canvas.draw_idle()
            return

        # Sun at origin in heliocentric coordinates.
        self.ax.scatter([0], [0], s=70, c="gold", edgecolor="k", linewidth=0.5, label="Sun", zorder=5)

        workbook_label = f"{Path(path).parent.name} / {Path(path).name}"
        for body in self.bodies:
            if body not in tracks:
                continue
            x, y = tracks[body]
            _progress_line(self.show_progress, workbook_label, "Plotting", body, len(x))
            c = COLOR_MAP.get(body, None)
            lw = self.linewidth * (1.7 if body == "BH" else 1.0)
            marker = "x" if body == "BH" else "o"
            msize = 35 if body == "BH" else 12
            self.ax.plot(x, y, lw=lw, color=c, label=body, alpha=0.95)
            if len(x) > 0:
                self.ax.scatter([x[-1]], [y[-1]], s=msize, color=c, marker=marker, zorder=6)

        xmin = self.cx - self.fov_x / 2.0
        xmax = self.cx + self.fov_x / 2.0
        ymin = self.cy - self.fov_y / 2.0
        ymax = self.cy + self.fov_y / 2.0
        self.ax.set_xlim(xmin, xmax)
        self.ax.set_ylim(ymin, ymax)

        if self.show_legend:
            self.ax.legend(fontsize=8, ncols=2, loc="upper right")

        if self.show_help:
            self.ax.text(
                0.01, 0.01,
                HELP_TEXT,
                transform=self.ax.transAxes,
                fontsize=9,
                color="w",
                bbox=dict(facecolor="black", alpha=0.45, pad=4),
            )

        self.fig.tight_layout()
        _progress_done(self.show_progress)
        self.fig.canvas.draw_idle()

    def _on_key(self, event):
        key = (event.key or "").lower()

        if key == "right":
            self.i = (self.i + 1) % self.n
            self._reset_view()
            self._draw()
            return
        if key == "left":
            self.i = (self.i - 1) % self.n
            self._reset_view()
            self._draw()
            return

        # Pan with WASD. Leave arrow keys for prev/next to keep original gallery behavior.
        dx = self.pan_fraction * self.fov_x
        dy = self.pan_fraction * self.fov_y
        if key == "a":
            self.cx -= dx
        elif key == "d":
            self.cx += dx
        elif key == "w":
            self.cy += dy
        elif key == "s":
            self.cy -= dy
        elif key in ("+", "="):
            # Zoom in by shrinking the AU field of view.
            self.fov_x /= self.zoom_step
            self.fov_y /= self.zoom_step
        elif key in ("-", "_"):
            # Zoom out by enlarging the AU field of view.
            self.fov_x *= self.zoom_step
            self.fov_y *= self.zoom_step
        elif key == "0":
            self._reset_view()
        elif key == "h":
            self.show_help = not self.show_help
        elif key == "l":
            self.show_legend = not self.show_legend
        elif key in ("escape", "q"):
            plt.close(self.fig)
            return
        else:
            return

        self._draw()


def parse_bodies(s: str):
    if not s:
        return DEFAULT_BODIES
    return [b.strip() for b in s.split(",") if b.strip()]


def main():
    ap = argparse.ArgumentParser(
        description="Interactively replot heliocentric simulation workbooks with AU field of view."
    )
    ap.add_argument(
        "root",
        help="Root folder containing run subfolders, or one specific run folder."
    )
    ap.add_argument(
        "--fov",
        type=float,
        default=None,
        help="Square field of view in AU. Sets both --fov-x and --fov-y."
    )
    ap.add_argument(
        "--fov-x",
        type=float,
        default=80.0,
        help="Horizontal field of view in AU. Default: 80."
    )
    ap.add_argument(
        "--fov-y",
        type=float,
        default=80.0,
        help="Vertical field of view in AU. Default: 80."
    )
    ap.add_argument(
        "--bodies",
        type=str,
        default=",".join(DEFAULT_BODIES),
        help="Comma-separated body list to plot. Default: all planets, Moon, BH."
    )
    ap.add_argument(
        "--no-bh",
        action="store_true",
        help="Do not plot the BH track."
    )
    ap.add_argument(
        "--pan-fraction",
        type=float,
        default=0.15,
        help="Fraction of current FOV moved per W/A/S/D key press. Default: 0.15."
    )
    ap.add_argument(
        "--zoom-step",
        type=float,
        default=1.25,
        help="FOV multiplier/divider for +/- keys. Default: 1.25."
    )
    ap.add_argument(
        "--linewidth",
        type=float,
        default=1.0,
        help="Base line width. Default: 1.0."
    )
    ap.add_argument(
        "--no-legend",
        action="store_true",
        help="Hide legend by default. Press L in the viewer to toggle."
    )
    ap.add_argument(
        "--progress",
        action="store_true",
        default=True,
        help="Show same-line console progress with workbook name while loading/plotting each body. Default: on."
    )
    ap.add_argument(
        "--no-progress",
        dest="progress",
        action="store_false",
        help="Disable console progress messages."
    )

    args = ap.parse_args()

    fov_x = args.fov if args.fov is not None else args.fov_x
    fov_y = args.fov if args.fov is not None else args.fov_y

    if fov_x <= 0 or fov_y <= 0:
        raise ValueError("FOV values must be positive.")
    if args.zoom_step <= 1.0:
        raise ValueError("--zoom-step must be > 1.")

    files = find_orbit_workbooks(args.root)
    if not files:
        print(f"No orbit workbooks found under: {args.root}")
        print("Expected files matching *__orbits__*.xlsx")
        sys.exit(1)

    viewer = ReplotViewer(
        files=files,
        bodies=parse_bodies(args.bodies),
        fov_x=fov_x,
        fov_y=fov_y,
        pan_fraction=args.pan_fraction,
        zoom_step=args.zoom_step,
        show_legend=not args.no_legend,
        show_bh=not args.no_bh,
        linewidth=args.linewidth,
        show_progress=args.progress,
    )
    plt.show()


if __name__ == "__main__":
    main()
