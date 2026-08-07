"""Sky tracks (RA/Dec) of the black hole and planets during a flyby.

Plots the apparent path of every major body across the sky as seen from two
observers -- by default Earth and Mars -- over a chosen window.

Why the shipped ``*_radec*.xlsx`` file is not used
--------------------------------------------------
That file contains the black hole only. It is not needed: the REBOUND runs are
integrated in the **equatorial J2000 frame**, so right ascension and declination
for any target from any observer come straight out of the state vectors with no
frame rotation at all:

    v   = r_target - r_observer
    RA  = atan2(v_y, v_x)
    Dec = asin(v_z / |v|)

Verified against the shipped BH file: geocentric RA/Dec agree to six decimal
places at t = 0, t = 100 and at periapsis.

Positions are **geometric**, not apparent: no light-time correction, no
aberration, no refraction. Light time is minutes at these ranges and aberration
is ~20 arcsec, both far below anything visible at this plot scale, but they mean
these are not ephemeris-grade coordinates.

Usage
-----
    python plot_sky_tracks.py <parquet_dir>/<run>/orbits.parquet \
        --start 2027-07-01 --end 2027-10-31 --out sky.png
"""

from __future__ import annotations

import argparse
import datetime as dt
from pathlib import Path

import numpy as np
import pandas as pd

EPOCH = dt.datetime(1873, 9, 1)

INK, INK2, MUTED = "#0b0b0b", "#52514e", "#898781"
GRID, AXIS, SURFACE = "#e1e0d9", "#c3c2b7", "#fcfcfb"

# Draw order. The BH leads because it is the subject of the plot.
BODIES = ["BH", "Sun", "Mercury", "Venus", "Earth", "Mars", "Jupiter", "Saturn"]

# Colour is bound to the *body*, not to its position in the drawing order.
# Each panel omits its own observer, so an index-based assignment would shift
# every colour after that body and make the same hue mean Mars in one panel and
# Earth in the other. The palette's categorical slots are taken in their
# validated order, which is what the adjacent-pair guarantee covers -- and line
# tracks are exactly the adjacent-pair case.
BODY_COLOUR = {
    "BH":      "#2a78d6",   # slot 1  blue
    "Sun":     "#eb6834",   # slot 2  orange
    "Mercury": "#1baf7a",   # slot 3  aqua
    "Venus":   "#eda100",   # slot 4  yellow
    "Earth":   "#e87ba4",   # slot 5  magenta
    "Mars":    "#008300",   # slot 6  green
    "Jupiter": "#4a3aa7",   # slot 7  violet
    "Saturn":  "#e34948",   # slot 8  red
}


def to_t(datestr: str) -> int:
    return (dt.datetime.strptime(datestr, "%Y-%m-%d") - EPOCH).days


def to_date(t: float) -> dt.datetime:
    return EPOCH + dt.timedelta(days=float(t))


def sky_track(frames: dict, observer: str, target: str):
    """Geometric RA [deg], Dec [deg] and range [AU] of ``target`` from ``observer``."""
    o, s = frames[observer], frames[target]
    v = np.c_[s.x_au.values - o.x_au.values,
              s.y_au.values - o.y_au.values,
              s.z_au.values - o.z_au.values]
    r = np.linalg.norm(v, axis=1)
    ra = np.degrees(np.arctan2(v[:, 1], v[:, 0])) % 360.0
    dec = np.degrees(np.arcsin(v[:, 2] / r))
    return ra, dec, r


def break_wraps(ra: np.ndarray, *others):
    """Insert NaN where RA crosses 0/360 so the line breaks instead of sweeping.

    Without this, a body near RA = 0 draws a spurious horizontal streak across
    the whole plot every time it wraps.
    """
    ra = ra.copy()
    jump = np.where(np.abs(np.diff(ra)) > 180.0)[0]
    out = [ra] + [np.asarray(o, dtype=float).copy() for o in others]
    for arr in out:
        arr[jump] = np.nan
    return out


def _style(ax):
    ax.set_facecolor(SURFACE)
    for name, sp in ax.spines.items():
        sp.set_visible(name not in ("top", "right"))
        sp.set_color(AXIS)
    ax.tick_params(colors=MUTED, labelsize=9)
    ax.xaxis.label.set_color(INK2)
    ax.yaxis.label.set_color(INK2)
    ax.grid(True, color=GRID, linewidth=0.8)
    ax.set_axisbelow(True)


def main() -> int:
    p = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    p.add_argument("parquet", type=Path, help="A run's orbits.parquet")
    p.add_argument("--start", default="2027-07-01")
    p.add_argument("--end", default="2027-10-31")
    p.add_argument("--observers", nargs=2, default=["Earth", "Mars"])
    p.add_argument("--out", type=Path, default=Path("sky_tracks.png"))
    args = p.parse_args()

    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    t0, t1 = to_t(args.start), to_t(args.end)
    d = pd.read_parquet(args.parquet)
    d = d[(d.t_days >= t0) & (d.t_days <= t1)]
    frames = {b: g.sort_values("t_days") for b, g in d.groupby("body", observed=True)}
    times = frames[args.observers[0]].t_days.values

    # Month boundaries, for reading timing off the tracks.
    months = [i for i in range(1, len(times))
              if to_date(times[i]).month != to_date(times[i - 1]).month]

    fig, axes = plt.subplots(1, 2, figsize=(15, 6.4))
    fig.patch.set_facecolor(SURFACE)

    for ax, obs in zip(axes, args.observers):
        _style(ax)
        for body in [b for b in BODIES if b != obs and b in frames]:
            ra, dec, rng = sky_track(frames, obs, body)
            ra, dec = break_wraps(ra, dec)
            c = BODY_COLOUR[body]
            is_bh = body == "BH"
            ax.plot(ra, dec, color=c, lw=2.6 if is_bh else 1.6,
                    alpha=1.0 if is_bh else 0.85, zorder=5 if is_bh else 3,
                    label=f"{body}" if not is_bh else "Black hole",
                    solid_capstyle="round")
            if is_bh:
                # Month ticks and endpoints, so the sweep can be read in time.
                ax.scatter(ra[months], dec[months], s=26, color=c,
                           edgecolors=SURFACE, linewidths=0.8, zorder=6)
                for i, lab in ((0, args.start[:7]), (len(ra) - 1, args.end[:7])):
                    if np.isfinite(ra[i]):
                        ax.annotate(lab, (ra[i], dec[i]), textcoords="offset points",
                                    xytext=(6, 6), fontsize=8, color=c, weight="bold")

        ax.invert_xaxis()          # RA increases eastward, i.e. to the left
        ax.set_xlabel("Right ascension [deg]")
        ax.set_ylabel("Declination [deg]")
        ax.set_title(f"Sky as seen from {obs}", color=INK, fontsize=12, loc="left")
        # Each panel gets its own legend: the body sets differ by one (each
        # observer cannot see itself), so a single shared legend would be wrong.
        leg = ax.legend(loc="best", frameon=False, fontsize=8.5, ncol=2)
        for t in leg.get_texts():
            t.set_color(INK2)

    fig.suptitle(
        f"Apparent tracks, {args.start} to {args.end} — black-hole periapsis July 2027",
        color=INK, fontsize=13, x=0.01, ha="left")
    fig.tight_layout(rect=(0, 0, 1, 0.95))
    fig.savefig(args.out, dpi=150, facecolor=SURFACE)
    plt.close(fig)
    print(f"Wrote {args.out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
