"""Animate the black hole's path across the sky, as MP4 or GIF.

Pacing is the design problem this script exists to solve. Real time is useless
here: the object is effectively motionless for 160 years and then crosses the
sky in two, so a clock-paced animation is almost all stillness followed by one
blurred frame.

So by default the animation advances at constant *apparent speed* -- equal arc
along the drawn track per frame -- and lets the calendar run at whatever rate
that implies. The readout sprints through years per second early on and crawls
to days per second at perihelion, which is the story told as a rate rather than
as geometry. ``--pace time`` gives the honest-but-dull linear-in-time version
for comparison.

ffmpeg comes from the ``imageio-ffmpeg`` package if it is installed, so no
system-wide install is needed; a system ffmpeg on PATH is used otherwise. With
neither, pass ``--gif``.

Usage
-----
    python animate_sky_track.py <run_dir>
    python animate_sky_track.py <run_dir> --seconds 45 --width 1920
    python animate_sky_track.py <run_dir> --gif --width 900
"""

from __future__ import annotations

import argparse
import datetime as dt
import glob
from pathlib import Path

import numpy as np
import pandas as pd
import yaml
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.animation as manim
from matplotlib.colors import LinearSegmentedColormap, Normalize
from matplotlib.collections import LineCollection

SURFACE = "#050a1e"
INK, MUTED, FAINT = "#e8e8ea", "#9a9aa2", "#232a48"
RAMP = LinearSegmentedColormap.from_list("amber", ["#6b2f16", "#d95926", "#ffc79a"])
LABEL, HEAD = "#c9b96a", "#fff0dd"
MON = ['Jan','Feb','Mar','Apr','May','Jun','Jul','Aug','Sep','Oct','Nov','Dec']

# Chart orientation. RA_LEFT_EDGE is the RA hour at the left edge, which is also
# where the track is allowed to break; RA_RIGHT says whether RA increases
# rightward.
#
# The default is the star-atlas convention: RA increases to the *left*, so the
# axis reads 24h at the left edge down to 0h at the right. That is how the sky
# looks from inside the celestial sphere, i.e. looking up, and it means the
# chart can be compared directly against any printed atlas. Setting
# RA_RIGHT = True mirrors it into a left-to-right graph, which reads more
# naturally as a plot but is backwards as a sky map.
#
# Seaming at 0h/24h matters independently of orientation: this run spans RA
# 5.05h-23.17h and never crosses 0h, so the track draws as one unbroken curve.
# The original 8h seam cut it into two disconnected pieces.
RA_LEFT_EDGE = 24.0
RA_RIGHT = False


def ra_to_x(ra_hours):
    """Map right ascension in hours to the chart's x coordinate."""
    if RA_RIGHT:
        return (np.asarray(ra_hours) - RA_LEFT_EDGE) % 24.0
    return (RA_LEFT_EDGE - np.asarray(ra_hours)) % 24.0


def tick_label(h: int) -> str:
    """RA label for the tick at chart x = h."""
    raw = RA_LEFT_EDGE + h if RA_RIGHT else RA_LEFT_EDGE - h
    ra = raw % 24.0
    # Wrap into [0, 24] rather than [0, 24): at the far edge the modulo lands
    # back on the seam hour, and showing 0h at both ends makes the axis look
    # like it restarts instead of spanning the sky once.
    if ra == 0.0 and raw != 0.0:
        ra = 24.0
    return f"{int(round(ra))}ʰ"


CONST = [(5.6,40,"Auriga"),(4.5,17,"Taurus"),(5.5,0,"Orion"),(7.0,24,"Gemini"),
    (8.5,20,"Cancer"),(10.5,15,"Leo"),(13.3,-2,"Virgo"),(15.2,-16,"Libra"),
    (16.8,-30,"Scorpius"),(19.0,-28,"Sagittarius"),(21.0,-19,"Capricornus"),
    (22.5,-10,"Aquarius"),(0.8,12,"Pisces"),(2.6,21,"Aries"),(1.7,-10,"Cetus"),
    (22.5,20,"Pegasus"),(20.5,42,"Cygnus"),(18.8,37,"Lyra"),(19.7,5,"Aquila"),
    (17.2,-5,"Ophiuchus"),(17.2,30,"Hercules"),(14.7,30,"Bootes"),
    (11.5,55,"Ursa Major"),(3.4,45,"Perseus"),(0.8,38,"Andromeda"),
    (1.0,60,"Cassiopeia"),(6.8,-24,"Canis Major"),(3.5,-20,"Eridanus"),
    (10.0,-20,"Hydra"),(12.4,-18,"Corvus"),(12.0,40,"Canes Ven.")]


def _fallback_labels(ax) -> None:
    """Approximate constellation names, used when no star catalogue is present."""
    for rah, dcd, nm in CONST:
        ax.text(float(ra_to_x(rah)), dcd, nm, color=LABEL, fontsize=9,
                ha="center", va="center", zorder=3, alpha=0.8)


def find_ffmpeg() -> str | None:
    try:
        import imageio_ffmpeg
        return imageio_ffmpeg.get_ffmpeg_exe()
    except Exception:
        from shutil import which
        return which("ffmpeg")


def main() -> int:
    p = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    p.add_argument("run_dir", type=Path)
    p.add_argument("--seconds", type=float, default=30.0, help="Playback length.")
    p.add_argument("--fps", type=int, default=30)
    p.add_argument("--width", type=int, default=1600, help="Pixels; height follows.")
    p.add_argument("--pace", choices=["arc", "time"], default="arc",
                   help="'arc' = constant apparent speed; 'time' = linear in time.")
    p.add_argument("--from-year", type=int, default=None,
                   help="Restrict the animation to this year onward.")
    p.add_argument("--to-year", type=int, default=None,
                   help="Restrict the animation to this year and earlier.")
    p.add_argument("--stars", type=float, default=5.0, metavar="MAG",
                   help="Faintest star drawn in the backdrop; 0 disables it "
                        "and falls back to plain name labels.")
    p.add_argument("--gif", action="store_true")
    p.add_argument("--out", type=Path, default=None)
    args = p.parse_args()

    cfg = yaml.safe_load(open(glob.glob(str(args.run_dir / "*input.yaml"))[0],
                              encoding="utf-8"))
    epoch = dt.datetime.strptime(str(cfg["epoch"])[:19], "%Y-%m-%dT%H:%M:%S")
    toff = float(cfg["bh_tperi_offset_days"])
    win = float(cfg.get("output_dense_window_days", 0) or 0)

    W = args.width - (args.width % 2)              # H.264 needs even dimensions
    FIG_W, DPI = 16.0, W / 16.0
    FIG_H = 8.5
    H = int(round(FIG_H * DPI)); H -= H % 2
    FIG_H = H / DPI
    X_SPAN, Y_SPAN = 24.0, 150.0

    d = pd.read_parquet(glob.glob(str(args.run_dir / "*orbits*.parquet"))[0])
    g = {b: v.sort_values("t_days") for b, v in d.groupby("body", observed=True)}
    tn = g["BH"].t_days.values
    lo, hi = toff - win, toff + win
    t = (np.unique(np.concatenate([tn[(tn < lo) | (tn > hi)],
                                   np.arange(lo, hi + .05, .05)])) if win else tn)
    P = {b: np.array([np.interp(t, g[b].t_days, g[b][c].values)
                      for c in ("x_au", "y_au", "z_au")]) for b in ("BH", "Earth", "Sun")}
    v = P["BH"] - P["Earth"]
    r = np.linalg.norm(v, axis=0)
    ra = np.degrees(np.arctan2(v[1], v[0])) % 360 / 15.0
    dec = np.degrees(np.arcsin(v[2] / r))
    rs = np.linalg.norm(P["BH"] - P["Sun"], axis=0)
    yrs = np.array([epoch.year + (x + epoch.timetuple().tm_yday) / 365.25 for x in t])
    X, Y = ra_to_x(ra), dec

    dX, dY = np.diff(X), np.diff(Y)
    seam = np.abs(dX) > 12
    seg = np.hypot(dX * FIG_W / X_SPAN, dY * FIG_H / Y_SPAN); seg[seam] = 0.0
    arc = np.concatenate([[0.0], np.cumsum(seg)])

    # Optional date window. Arc is re-accumulated inside it so the pacing is
    # even across the window rather than across the whole track.
    yr_of = np.array([(epoch + dt.timedelta(days=float(x))).year for x in t])
    lo_i = int(np.argmax(yr_of >= args.from_year)) if args.from_year else 0
    hi_i = (len(t) - 1 - int(np.argmax(yr_of[::-1] <= args.to_year))
            if args.to_year else len(t) - 1)
    if hi_i <= lo_i:
        raise SystemExit("Empty date window.")

    n_frames = int(args.seconds * args.fps)
    if args.pace == "arc":
        a0, a1 = arc[lo_i], arc[hi_i]
        frames = np.searchsorted(arc, np.linspace(a0, a1, n_frames))
    else:
        frames = np.searchsorted(t, np.linspace(t[lo_i], t[hi_i], n_frames))
    frames = frames.clip(lo_i, hi_i)

    fig, ax = plt.subplots(figsize=(FIG_W, FIG_H), dpi=DPI)
    fig.patch.set_facecolor(SURFACE); ax.set_facecolor(SURFACE)
    # Fix the limits up front: the backdrop needs them to decide which
    # constellation names fall on the chart, and adding its collections would
    # otherwise drag the autoscale around.
    ax.set_xlim(0, 24); ax.set_ylim(-70, 80)
    for h in range(0, 25, 2): ax.axvline(h, color=FAINT, lw=0.8, zorder=1)
    for dd in range(-60, 81, 20): ax.axhline(dd, color=FAINT, lw=0.8, zorder=1)
    ax.axhline(0, color="#3a4468", lw=1.3, zorder=1)
    if args.stars > 0:
        import sky_backdrop
        if not sky_backdrop.draw_backdrop(ax, ra_to_x, max_mag=args.stars,
                                          label_size=8.5):
            print("  no star data; run fetch_constellation_data.py once")
            _fallback_labels(ax)
    else:
        _fallback_labels(ax)
    bounds = [0] + list(np.where(seam)[0] + 1) + [len(X)]
    for a, b in zip(bounds[:-1], bounds[1:]):
        if b - a > 1:
            ax.plot(X[a:b], Y[a:b], color="#3a2a2a", lw=1.4, zorder=4)

    if lo_i > 0:      # what the object already did, so the window has context
        for a, b in zip(bounds[:-1], bounds[1:]):
            e = min(b, lo_i + 1)
            if e - a > 1:
                ax.plot(X[a:e], Y[a:e], color="#6b4a3a", lw=2.0, zorder=5)

    norm = Normalize(yrs[lo_i:hi_i + 1].min(), yrs[lo_i:hi_i + 1].max())
    trail = LineCollection([], cmap=RAMP, norm=norm, linewidth=3.0, zorder=6)
    ax.add_collection(trail)
    glow, = ax.plot([], [], "o", ms=28, color="#ffb27a", alpha=0.22, zorder=8)
    head, = ax.plot([], [], "o", ms=12, color=HEAD, mec="#d95926", mew=2.2, zorder=9)

    ax.set_xlim(0, 24); ax.set_ylim(-70, 80)   # re-assert: add_collection autoscales
    ax.set_xticks(range(0, 25, 2))
    ax.set_xticklabels([tick_label(h) for h in range(0, 25, 2)])
    ax.set_yticks(range(-60, 81, 20))
    ax.set_yticklabels([f"{q:+d}°" for q in range(-60, 81, 20)])
    ax.tick_params(colors=MUTED, labelsize=10)
    for sp in ax.spines.values(): sp.set_color(MUTED)
    ax.set_xlabel("right ascension", color=MUTED, fontsize=11)
    ax.set_ylabel("declination", color=MUTED, fontsize=11)
    pace_txt = ("at constant apparent speed" if args.pace == "arc" else "in real time")
    y_lo = (epoch + dt.timedelta(days=float(t[lo_i]))).year
    y_hi = (epoch + dt.timedelta(days=float(t[hi_i]))).year
    ax.set_title(f"A black hole crosses the sky — {y_lo} to {y_hi}, {pace_txt}",
                 color=INK, fontsize=15, loc="left", pad=12)
    readout = ax.text(0.985, 0.055, "", transform=ax.transAxes, ha="right",
                      va="bottom", color=INK, fontsize=16, family="monospace",
                      zorder=12, linespacing=1.6,
                      bbox=dict(boxstyle="round,pad=0.5", fc="#0d1430",
                                ec="#2c3660", alpha=0.94))
    ax.text(0.012, 0.03,
            ("Equal distance along the track per frame — the calendar, not the motion, is what accelerates."
             if args.pace == "arc" else
             "Linear in time: 160 years of stillness, then the whole sky in two."),
            transform=ax.transAxes, color=MUTED, fontsize=9.5, va="bottom", zorder=12)
    fig.tight_layout()

    def draw(n):
        i = frames[n]
        j = max(lo_i + 1, i)
        ok = ~seam[lo_i:j]
        pts = np.array([X[lo_i:j + 1], Y[lo_i:j + 1]]).T.reshape(-1, 1, 2)
        trail.set_segments(np.concatenate([pts[:-1], pts[1:]], axis=1)[ok])
        trail.set_array(yrs[lo_i:j][ok])
        head.set_data([X[i]], [Y[i]]); glow.set_data([X[i]], [Y[i]])
        w = epoch + dt.timedelta(days=float(t[i]))
        nxt = frames[min(n + 1, len(frames) - 1)]
        rate = abs(t[nxt] - t[i]) * args.fps / 365.25
        rs_txt = f"{rate:,.1f} yr/s" if rate >= 1 else f"{rate * 365.25:.1f} d/s"
        readout.set_text(f"{w.day:2d} {MON[w.month - 1]} {w.year}\n"
                         f"{rs[i]:8.2f} AU from Sun\n{rs_txt:>15}")
        return trail, head, glow, readout

    ext = "gif" if args.gif else "mp4"
    dest = args.out or Path(f"bh_sky_track.{ext}")
    ani = manim.FuncAnimation(fig, draw, frames=len(frames), interval=1000 / args.fps,
                              blit=False)
    if args.gif:
        ani.save(dest, writer=manim.PillowWriter(fps=args.fps), dpi=DPI)
    else:
        exe = find_ffmpeg()
        if not exe:
            raise SystemExit("No ffmpeg found. `pip install imageio-ffmpeg`, or use --gif.")
        matplotlib.rcParams["animation.ffmpeg_path"] = exe
        writer = manim.FFMpegWriter(
            fps=args.fps, codec="libx264",
            extra_args=["-pix_fmt", "yuv420p", "-crf", "18", "-preset", "slow"])
        ani.save(dest, writer=writer, dpi=DPI)
    import os
    print(f"wrote {dest}  {os.path.getsize(dest) / 1e6:.1f} MB  "
          f"{W}x{H}  {len(frames)} frames  {len(frames) / args.fps:.1f} s")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
