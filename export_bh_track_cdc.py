"""Write the black hole's sky track as a Cartes du Ciel user-defined-object file.

Format, reverse-engineered from CdC's own shipped example row:

    NGP 12h51m26s +27d07m42s 0 14 0.00 30.00 0 Example of user defined object: ...
    |   |         |          | |  |    |     | |
    |   RA        Dec        ? ?  magn size' ? free-text description
    name

Only five of the nine fields could be matched against the dialog columns (name,
RA, Dec, magnitude, size in arcminutes). The two unlabelled integers after Dec
and the one before the description are almost certainly type and colour codes,
but rather than guess they are **copied verbatim from the example**, so every
generated point renders exactly like the row CdC ships and is known to draw.

Coordinates are J2000. That is not an assumption: CdC's own example is the North
Galactic Pole at 12h51m26s +27d07m42s, which is its J2000 position (the standard
value is 12h51m26.28s, +27d07'41.7"), so the file is read as J2000 and precessed
for display. Using equinox-of-date values here would misplace the track by 1.6
degrees at the 1885 epoch.

Sampling is at constant angular spacing along the heliocentric path -- see
``export_bh_track.py``, which this shares its logic with.

Usage
-----
    python export_bh_track_cdc.py <run_dir> --step-deg 1.5
    python export_bh_track_cdc.py <run_dir> --step-deg 0.2 --from 2046 --to 2049
"""

from __future__ import annotations

import argparse
import datetime as dt
import glob
from pathlib import Path

import numpy as np
import pandas as pd
import yaml

# Copied verbatim from CdC's shipped example so the markers are known to render.
F4, F5, F8 = "0", "14", "0"


def fmt_ra(deg: float) -> str:
    h = deg / 15.0
    hh = int(h); m = (h - hh) * 60
    mm = int(m); ss = int(round((m - mm) * 60))
    if ss == 60: ss, mm = 0, mm + 1
    if mm == 60: mm, hh = 0, hh + 1
    return f"{hh % 24:02d}h{mm:02d}m{ss:02d}s"


def fmt_dec(deg: float) -> str:
    sign = "+" if deg >= 0 else "-"
    a = abs(deg)
    dd = int(a); m = (a - dd) * 60
    mm = int(m); ss = int(round((m - mm) * 60))
    if ss == 60: ss, mm = 0, mm + 1
    if mm == 60: mm, dd = 0, dd + 1
    return f"{sign}{dd:02d}d{mm:02d}m{ss:02d}s"


def main() -> int:
    p = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    p.add_argument("run_dir", type=Path)
    p.add_argument("--step-deg", type=float, default=1.5,
                   help="Angular spacing along the track [deg].")
    p.add_argument("--size-arcmin", type=float, default=5.0,
                   help="Marker size; CdC's example uses 30.")
    p.add_argument("--label-every", type=int, default=10,
                   help="Give every Nth point a year label; others get a dot.")
    p.add_argument("--from", dest="y0", type=int, default=None)
    p.add_argument("--to", dest="y1", type=int, default=None)
    p.add_argument("--out", type=Path, default=Path("bh_track_cdc.txt"))
    p.add_argument("--fine-days", type=float, default=0.05)
    args = p.parse_args()

    cfgf = glob.glob(str(args.run_dir / "*input.yaml"))[0]
    cfg = yaml.safe_load(open(cfgf, encoding="utf-8"))
    epoch = dt.datetime.strptime(str(cfg["epoch"])[:19], "%Y-%m-%dT%H:%M:%S")
    toff = float(cfg["bh_tperi_offset_days"])
    win = float(cfg.get("output_dense_window_days", 0) or 0)

    d = pd.read_parquet(glob.glob(str(args.run_dir / "*orbits*.parquet"))[0])
    g = {b: v.sort_values("t_days") for b, v in d.groupby("body", observed=True)}
    tn = g["BH"].t_days.values
    lo, hi = toff - win, toff + win
    fine = np.arange(lo, hi + args.fine_days, args.fine_days)
    times = np.unique(np.concatenate([tn[(tn < lo) | (tn > hi)], fine])) if win else tn

    def pos(b):
        return np.array([np.interp(times, g[b].t_days, g[b][c].values)
                         for c in ("x_au", "y_au", "z_au")])

    B, S, E = pos("BH"), pos("Sun"), pos("Earth")
    vh, vg = B - S, B - E
    rh = np.linalg.norm(vh, axis=0)
    uh = vh / rh
    step = np.degrees(np.arccos(np.clip((uh[:, :-1] * uh[:, 1:]).sum(axis=0), -1, 1)))
    arc = np.concatenate([[0.0], np.cumsum(step)])
    want = np.arange(0.0, arc[-1] + args.step_deg * 0.5, args.step_deg)
    idx = np.unique(np.concatenate([[0], np.searchsorted(arc, want).clip(0, len(times) - 1),
                                    [len(times) - 1]]))

    dates = np.array([epoch + dt.timedelta(days=float(times[i])) for i in idx])
    if args.y0 or args.y1:
        keep = np.array([(args.y0 or 0) <= x.year <= (args.y1 or 9999) for x in dates])
        idx, dates = idx[keep], dates[keep]

    r = np.linalg.norm(vg[:, idx], axis=0)
    ra = np.degrees(np.arctan2(vg[1, idx], vg[0, idx])) % 360
    dec = np.degrees(np.arcsin(vg[2, idx] / r))
    rs = rh[idx]

    lines = []
    for k, (when, a, dd, dist) in enumerate(zip(dates, ra, dec, rs)):
        label = f"BH{when.year}" if k % args.label_every == 0 else "."
        desc = (f"black hole {when:%Y-%m-%d}  {dist:.2f} AU from the Sun")
        lines.append(f"{label} {fmt_ra(a)} {fmt_dec(dd)} {F4} {F5} 0.00 "
                     f"{args.size_arcmin:.2f} {F8} {desc}")
    args.out.write_text("\n".join(lines) + "\n", encoding="utf-8")

    print(f"{len(lines)} points -> {args.out}")
    print(f"  spacing {args.step_deg} deg, marker {args.size_arcmin}'")
    print(f"  span {dates[0]:%Y-%m-%d} to {dates[-1]:%Y-%m-%d}")
    print(f"  labelled points: {sum(1 for l in lines if not l.startswith('. '))}")
    print("\nfirst two lines:")
    for l in lines[:2]:
        print("   " + l)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
