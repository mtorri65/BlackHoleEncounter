#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Keyboard gallery for heliocentric plots.

- Shows one "*__heliocentric_plot.png" at a time from all subfolders of a root run folder
- Left / Right arrows to move backward / forward
- Esc or Q to quit
- H to toggle help banner

Usage (example):
    python gallery_heliocentric_viewer.py "simulations/20251107_223214"

Notes:
- Works recursively: finds all matches under the given root.
- Sorts images by folder name then filename.
"""

import argparse
import os
import glob
import sys
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.image import imread

HELP_TEXT = "←/→: prev/next   H: help   Esc/Q: quit"

def find_images(root: str):
    pattern = os.path.join(root, "**", "*__heliocentric_plot.png")
    files = sorted(glob.glob(pattern, recursive=True))
    return files

class ImageViewer:
    def __init__(self, files):
        if not files:
            raise FileNotFoundError("No *__heliocentric_plot.png images found.")
        self.files = files
        self.n = len(files)
        self.i = 0
        self.show_help = True

        self.fig, self.ax = plt.subplots(figsize=(9, 9))
        self.ax.axis("off")
        self.img_artist = None
        self.text_artist = None
        # Draw first image
        self._draw()

        # Connect events
        self.cid = self.fig.canvas.mpl_connect("key_press_event", self._on_key)

        # Window title (best-effort; some backends ignore)
        try:
            self.fig.canvas.manager.set_window_title("Heliocentric Plot Gallery")
        except Exception:
            pass

    def _read(self, path):
        # imread supports PNG natively
        img = imread(path)
        # Ensure range 0..1 for floats
        if img.dtype.kind == 'f' and (img.max() > 1.0 or img.min() < 0.0):
            img = np.clip(img, 0.0, 255.0) / 255.0
        return img

    def _title(self):
        fn = os.path.basename(self.files[self.i])
        parent = os.path.basename(os.path.dirname(self.files[self.i]))
        return f"[{self.i+1}/{self.n}] {parent} / {fn}"

    def _draw(self):
        self.ax.clear()
        self.ax.axis("off")
        path = self.files[self.i]
        img = self._read(path)
        self.img_artist = self.ax.imshow(img)
        # Title on top
        self.ax.set_title(self._title(), fontsize=11)
        # Help banner (toggle with H)
        if self.show_help:
            self.text_artist = self.ax.text(
                0.01, 0.01, HELP_TEXT, transform=self.ax.transAxes,
                fontsize=10, color="w", bbox=dict(facecolor="black", alpha=0.4, pad=4)
            )
        self.fig.tight_layout()
        self.fig.canvas.draw_idle()

    def _on_key(self, event):
        key = (event.key or "").lower()
        if key in ("right",):
            self.i = (self.i + 1) % self.n
            self._draw()
        elif key in ("left",):
            self.i = (self.i - 1) % self.n
            self._draw()
        elif key in ("h",):
            self.show_help = not self.show_help
            self._draw()
        elif key in ("escape", "q"):
            plt.close(self.fig)

def main():
    ap = argparse.ArgumentParser(description="Flip through *__heliocentric_plot.png images with ←/→.")
    ap.add_argument("root", help="Root folder containing run subfolders (e.g. simulations/20251107_223214)")
    args = ap.parse_args()

    files = find_images(args.root)
    if not files:
        print(f"No images found under: {args.root}")
        sys.exit(1)

    viewer = ImageViewer(files)
    plt.show()

if __name__ == "__main__":
    main()
