"""Export the black hole's sky track as a CSV, evenly spaced along the sky.

Written for plotting the track in external planetarium software (SkyChart /
Cartes du Ciel and similar), which knows its own projection and so registers the
points exactly -- something an image overlay cannot do reliably.

Sampling
--------
Points are spaced at a roughly constant *angular* step along the track rather
than at constant time, because the object's apparent speed varies by orders of
magnitude: it drifts a few arcseconds a year at 860 AU and sweeps degrees a day
around perihelion.

The step is measured along the **heliocentric** direction -- the true path --
not the geocentric one. Geocentrically the track also carries an annual parallax
loop, which is a few arcminutes when the hole is far (invisible at chart scale)
but several degrees once it is inside ~40 AU. Spacing by geocentric arc would
spend most of the points chasing those loops instead of following the path.

For the default configuration this means 98% of the samples fall within about
five years of perihelion, because that is where 98% of the motion happens. The
long quiet centuries genuinely need only a handful of points.

Interpolation
-------------
Only inside the engine's dense-output window, where the trajectory is logged
daily and sub-day interpolation of Earth's position is harmless. Outside it the
log is 30-day and Earth moves ~30 deg between samples, so interpolating there
would produce badly wrong *geocentric* positions; native sample times are used
instead. This costs nothing, since the heliocentric motion out there is under
0.4 deg per 30-day step.

Coordinates
-----------
The simulation runs in the equatorial J2000 frame, so J2000 columns are exact.
Equinox-of-date columns are also written, precessed per point, for charts drawn
in "apparent" coordinates. **Set your chart to whichever you use** -- at these
epochs the two differ by up to ~1.6 deg, which is far larger than any plotting
tolerance.

Usage
-----
    python export_bh_track.py simulations/<STAMP>/<RUN> --step-deg 0.5
    python export_bh_track.py simulations/<STAMP>/<RUN> --step-deg 0.25 --out track.csv
"""

from __future__ import annotations

import argparse
import datetime as dt
import glob
from pathlib import Path

import numpy as np
import pandas as pd
import yaml


def read_epoch(run_dir: Path) -> dt.datetime:
    """Simulation epoch from the run's own copied config."""
    for cand in glob.glob(str(run_dir / "*input.yaml")):
        raw = yaml.safe_load(open(cand, encoding="utf-8")).get("epoch")
        if raw:
            return dt.datetime.strptime(str(raw)[:19], "%Y-%m-%dT%H:%M:%S")
    raise SystemExit(f"No epoch found in {run_dir}")


def hms(ra_deg: float) -> str:
    h = ra_deg / 15.0
    m = (h % 1) * 60
    return f"{int(h):02d} {int(m):02d} {(m % 1) * 60:05.2f}"


def dms(dec_deg: float) -> str:
    s = "+" if dec_deg >= 0 else "-"
    a = abs(dec_deg)
    m = (a % 1) * 60
    return f"{s}{int(a):02d} {int(m):02d} {(m % 1) * 60:04.1f}"


def main() -> int:
    p = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    p.add_argument("run_dir", type=Path, help="A single run folder.")
    p.add_argument("--step-deg", type=float, default=0.5,
                   help="Target angular spacing along the track [deg].")
    p.add_argument("--out", type=Path, default=None)
    p.add_argument("--fine-days", type=float, default=0.05,
                   help="Interpolation step inside the dense window [days].")
    args = p.parse_args()

    epoch = read_epoch(args.run_dir)
    cfg = yaml.safe_load(open(glob.glob(str(args.run_dir / "*input.yaml"))[0],
                              encoding="utf-8"))
    toff = float(cfg["bh_tperi_offset_days"])
    win = float(cfg.get("output_dense_window_days", 0) or 0)

    d = pd.read_parquet(glob.glob(str(args.run_dir / "*orbits*.parquet"))[0])
    g = {b: v.sort_values("t_days") for b, v in d.groupby("body", observed=True)}
    t_native = g["BH"].t_days.values

    # Candidate times: native everywhere, refined inside the dense window only.
    if win > 0:
        lo, hi = toff - win, toff + win
        fine = np.arange(lo, hi + args.fine_days, args.fine_days)
        times = np.unique(np.concatenate([t_native[(t_native < lo) | (t_native > hi)],
                                          fine]))
    else:
        times = t_native

    def pos(body):
        return np.array([np.interp(times, g[body].t_days, g[body][c].values)
                         for c in ("x_au", "y_au", "z_au")])

    B, S, E = pos("BH"), pos("Sun"), pos("Earth")
    vh, vg = B - S, B - E
    rh, rg = np.linalg.norm(vh, axis=0), np.linalg.norm(vg, axis=0)
    uh = vh / rh

    # Cumulative arc along the heliocentric direction, then resample uniformly.
    step = np.degrees(np.arccos(np.clip((uh[:, :-1] * uh[:, 1:]).sum(axis=0), -1, 1)))
    arc = np.concatenate([[0.0], np.cumsum(step)])
    want = np.arange(0.0, arc[-1] + args.step_deg * 0.5, args.step_deg)
    idx = np.searchsorted(arc, want).clip(0, len(times) - 1)
    # Always carry the true endpoints: the object stops moving long before the
    # integration ends, so the last uniform-arc step lands well short of it and
    # the final stretch would otherwise be missing entirely.
    idx = np.unique(np.concatenate([[0], idx, [len(times) - 1]]))

    def radec(v):
        r = np.linalg.norm(v, axis=0)
        return np.degrees(np.arctan2(v[1], v[0])) % 360, np.degrees(np.arcsin(v[2] / r))

    ra_g, dec_g = radec(vg[:, idx])
    ra_h, dec_h = radec(vh[:, idx])
    ts = times[idx]
    dates = [epoch + dt.timedelta(days=float(x)) for x in ts]

    # Equinox of date, per point, via skyfield's precession.
    from skyfield.api import load
    from skyfield.positionlib import ICRF
    sky_ts = load.timescale()
    ra_d, dec_d = np.zeros(len(idx)), np.zeros(len(idx))
    for i, (when, vv) in enumerate(zip(dates, vg[:, idx].T)):
        tt = sky_ts.utc(when.year, when.month, when.day,
                        when.hour, when.minute, when.second)
        r, dd, _ = ICRF(vv, t=tt).radec(epoch=tt)
        ra_d[i], dec_d[i] = r.hours * 15.0, dd.degrees

    astep = np.concatenate([[0.0], np.degrees(np.arccos(np.clip(
        (uh[:, idx][:, :-1] * uh[:, idx][:, 1:]).sum(axis=0), -1, 1)))])

    out = pd.DataFrame({
        "date_utc": [x.strftime("%Y-%m-%d %H:%M") for x in dates],
        "t_days": np.round(ts, 3),
        "ra_deg_j2000": np.round(ra_g, 5),
        "dec_deg_j2000": np.round(dec_g, 5),
        "ra_hms_j2000": [hms(v) for v in ra_g],
        "dec_dms_j2000": [dms(v) for v in dec_g],
        "ra_deg_of_date": np.round(ra_d, 5),
        "dec_deg_of_date": np.round(dec_d, 5),
        "ra_hms_of_date": [hms(v) for v in ra_d],
        "dec_dms_of_date": [dms(v) for v in dec_d],
        "ra_deg_j2000_helio": np.round(ra_h, 5),
        "dec_deg_j2000_helio": np.round(dec_h, 5),
        "dist_earth_au": np.round(rg[idx], 4),
        "dist_sun_au": np.round(rh[idx], 4),
        # annual parallax semi-amplitude = 1 AU / distance
        "parallax_arcmin": np.round(3437.75 / rh[idx], 2),
        "step_from_prev_deg": np.round(astep, 4),
    })
    dest = args.out or Path(f"bh_track_{args.run_dir.name[-40:]}.csv")
    out.to_csv(dest, index=False)

    print(f"{len(out)} points, target spacing {args.step_deg}°")
    print(f"  achieved step: median {astep[1:].mean():.3f}°, "
          f"max {astep[1:].max():.3f}°, min {astep[1:].min():.3f}°")
    print(f"  span {dates[0]:%Y-%m-%d} to {dates[-1]:%Y-%m-%d}")
    q = np.array([ (x - epoch).days for x in dates ])
    print(f"  {(np.abs(q - toff) <= win).sum()} of {len(out)} points "
          f"({100*(np.abs(q-toff)<=win).mean():.0f}%) lie within the dense window")
    print(f"\nWrote {dest}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
