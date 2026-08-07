"""Topocentric sky view of the flyby from a fixed site on Earth.

Converts the simulated geocentric positions into **altitude and azimuth** for an
observer at a given latitude and longitude, answering "where in my sky is it, and
when can I see it" rather than the geocentric RA/Dec of ``plot_sky_tracks.py``.

Where Earth's rotation comes from
---------------------------------
It is *not* in the simulation. REBOUND integrates Earth as a point mass, so the
runs contain no spin phase at all. Earth's rotation is, however, completely
decoupled from the flyby dynamics -- a passing mass changes the orbit, not the
day -- so it is supplied analytically from Greenwich Mean Sidereal Time at the
UTC of each sample, and combined with the simulated positions.

This is a real assumption, not a formality: it means the *date* of an event is
simulated but the *time of night* is reconstructed. It would break for a body
close enough to raise meaningful tides, which is not the case here.

Calendar drift: check which sweep your run came from
----------------------------------------------------
The date a plot is labelled with and the season it actually shows can disagree,
because the simulated year is not exactly the real one. How badly depends on when
the run was produced.

**Runs from before the Moon fix (engine commit be8cba8) -- includes the whole
``simulations/20260724_230314`` sweep.** The Moon was added at a fixed offset from
Earth with an orbital velocity, but Earth's velocity was never adjusted to
compensate, so the Earth-Moon pair started with ~12.4 m/s of spurious momentum.
That put Earth's semi-major axis at 1.00059 AU instead of 1.000018 and its year at
**365.570 days**. Over the 154 years from the 1873 epoch the simulated northern
solstice walks from 20 June in 1874 (correct) to 29 June by 1900, 31 July by 2000
and **9 August by 2026** -- about **+0.33 days/year, reaching ~49 days**. A plot
labelled "2027-09-15" from such a run shows the solar geometry of roughly 28 July.

**Runs from after the fix**, which reads the real Moon from the ephemeris, have a
correct orbit. The measured tropical year is then 365.2563 d against the true
365.2422: a residual of **+0.0141 days/year, about +4.4 days over the full 308-year
span** -- roughly a fifth of a solar-disc diameter's worth of seasonal phase, and
below the resolution of anything plotted here.

That residual is not an error left over from the fix. 0.0141 d/yr is 20.4
minutes/yr, which is exactly the 50.3 arcsec/yr precession of the equinoxes. A
point-mass integration applies no torque to Earth's spin axis, so the axis never
precesses and "seasons" repeat on the *sidereal* year rather than the tropical
one. It is irreducible in this model, and small enough not to matter.

Either way the plot is internally self-consistent -- the Sun position used for the
day/night shading is the simulated one. The caveat is only about mapping it onto a
real-world calendar date.

None of this affects ``plot_sky_tracks.py``: RA/Dec are inertial and carry no
calendar dependence. It bites only once a site on a rotating Earth is involved,
because that ties the inertial sky to a date.

What is and is not included
---------------------------
* **Topocentric parallax** is applied -- the observer sits on Earth's surface,
  not at its centre. At 1.1 AU this is only ~8 arcsec, but it costs nothing.
* **Light-time, aberration, refraction, nutation and polar motion** are not.
  Refraction alone reaches ~35 arcmin at the horizon, so altitudes near 0 deg
  are approximate.
* Positions between the simulation's samples are linearly interpolated. Fine at
  the 1-day cadence, except for a body moving fast at closest approach where the
  cadence is itself the limit. **Check the cadence of your run first**: sweeps
  using ``output_dense_window_days`` log every 30 days away from the encounter,
  and linear interpolation across 30 days of orbital motion is meaningless. This
  script is only valid inside the dense window.

Usage
-----
    python plot_local_sky.py <run>/orbits.parquet --lat 45.5 --lon 9.2 \
        --date 2027-09-15 --days 5 --out local_sky.png
"""

from __future__ import annotations

import argparse
import datetime as dt
from pathlib import Path

import numpy as np
import pandas as pd

EPOCH = dt.datetime(1873, 9, 1)
R_EARTH_AU = 6378.137 / 1.495978707e8      # equatorial radius in AU

INK, INK2, MUTED = "#0b0b0b", "#52514e", "#898781"
GRID, AXIS, SURFACE = "#e1e0d9", "#c3c2b7", "#fcfcfb"
NIGHT = "#dbe6f5"

BODIES = ["BH", "Sun", "Mercury", "Venus", "Mars", "Jupiter", "Saturn", "Moon"]
BODY_COLOUR = {
    "BH": "#2a78d6", "Sun": "#eb6834", "Mercury": "#1baf7a", "Venus": "#eda100",
    "Earth": "#e87ba4", "Mars": "#008300", "Jupiter": "#4a3aa7", "Saturn": "#e34948",
    "Moon": "#898781",
}


def julian_date(when: dt.datetime) -> float:
    """Julian Date from a UTC datetime (Gregorian, valid for our range)."""
    y, m = when.year, when.month
    if m <= 2:
        y, m = y - 1, m + 12
    a = y // 100
    b = 2 - a + a // 4
    day = (when.day + when.hour / 24.0 + when.minute / 1440.0
           + when.second / 86400.0)
    return (int(365.25 * (y + 4716)) + int(30.6001 * (m + 1))
            + day + b - 1524.5)


def gmst_deg(when: dt.datetime) -> float:
    """Greenwich Mean Sidereal Time [deg]. Standard polynomial in JD."""
    jd = julian_date(when)
    d = jd - 2451545.0
    t = d / 36525.0
    g = (280.46061837 + 360.98564736629 * d
         + 0.000387933 * t * t - t ** 3 / 38710000.0)
    return g % 360.0


def altaz(vec_eq: np.ndarray, lat_deg: float, lst_deg: float):
    """Altitude and azimuth [deg] of an equatorial vector for a site.

    ``vec_eq`` is the topocentric vector in the equatorial frame; azimuth is
    measured from north through east.
    """
    r = np.linalg.norm(vec_eq, axis=-1)
    ra = np.degrees(np.arctan2(vec_eq[..., 1], vec_eq[..., 0])) % 360.0
    dec = np.degrees(np.arcsin(vec_eq[..., 2] / r))
    H = np.radians((lst_deg - ra) % 360.0)          # hour angle
    phi, d = np.radians(lat_deg), np.radians(dec)
    alt = np.arcsin(np.sin(phi) * np.sin(d) + np.cos(phi) * np.cos(d) * np.cos(H))
    az = np.arctan2(-np.sin(H) * np.cos(d),
                    np.cos(phi) * np.sin(d) - np.sin(phi) * np.cos(d) * np.cos(H))
    return np.degrees(alt), np.degrees(az) % 360.0


def observer_offset(lat_deg: float, lst_deg: np.ndarray) -> np.ndarray:
    """Geocentric position of the observing site [AU], equatorial frame."""
    phi = np.radians(lat_deg)
    th = np.radians(lst_deg)
    return R_EARTH_AU * np.c_[np.cos(phi) * np.cos(th),
                              np.cos(phi) * np.sin(th),
                              np.full_like(th, np.sin(phi))]


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
    p.add_argument("parquet", type=Path)
    p.add_argument("--lat", type=float, required=True, help="Latitude [deg, +N]")
    p.add_argument("--lon", type=float, required=True, help="Longitude [deg, +E]")
    p.add_argument("--date", default="2027-09-15", help="Centre date (UTC)")
    p.add_argument("--days", type=float, default=5.0, help="Window length [days]")
    p.add_argument("--site", default=None, help="Label for the site")
    p.add_argument("--out", type=Path, default=Path("local_sky.png"))
    args = p.parse_args()

    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    centre = dt.datetime.strptime(args.date, "%Y-%m-%d")
    t_mid = (centre - EPOCH).days
    t0, t1 = t_mid - args.days / 2, t_mid + args.days / 2

    d = pd.read_parquet(args.parquet)
    d = d[(d.t_days >= t0 - 2) & (d.t_days <= t1 + 2)]
    frames = {b: g.sort_values("t_days") for b, g in d.groupby("body", observed=True)}

    # Fine grid: rotation dominates, so sample every 5 minutes.
    tt = np.arange(t0, t1, 5.0 / 1440.0)
    times = [EPOCH + dt.timedelta(days=float(t)) for t in tt]
    lst = np.array([gmst_deg(w) for w in times]) + args.lon
    site_vec = observer_offset(args.lat, lst)

    def topo(body):
        e, s = frames["Earth"], frames[body]
        geo = np.c_[np.interp(tt, e.t_days, s.x_au.values - e.x_au.values),
                    np.interp(tt, e.t_days, s.y_au.values - e.y_au.values),
                    np.interp(tt, e.t_days, s.z_au.values - e.z_au.values)]
        return geo - site_vec          # centre -> surface

    alt = {b: altaz(topo(b), args.lat, lst)[0]
           for b in BODIES if b in frames and b != "Earth"}
    az = {b: altaz(topo(b), args.lat, lst)[1]
          for b in BODIES if b in frames and b != "Earth"}

    site = args.site or f"{args.lat:+.1f}, {args.lon:+.1f}"
    fig, (ax1, ax2) = plt.subplots(
        1, 2, figsize=(15, 6.2),
        gridspec_kw={"width_ratios": [1.35, 1]})
    fig.patch.set_facecolor(SURFACE)

    # --- Panel A: altitude vs time, with night shaded --------------------
    _style(ax1)
    sun = alt["Sun"]
    night = sun < -6.0                      # civil twilight
    ax1.fill_between(tt - t_mid, -90, 90, where=night, color=NIGHT,
                     zorder=0, linewidth=0, label="night (Sun < −6°)")
    for b, a in alt.items():
        ax1.plot(tt - t_mid, a, color=BODY_COLOUR[b],
                 lw=2.4 if b == "BH" else 1.4,
                 alpha=1.0 if b == "BH" else 0.85,
                 zorder=5 if b == "BH" else 3,
                 label="Black hole" if b == "BH" else b)
    ax1.axhline(0, color=INK2, lw=1.2, zorder=4)
    ax1.set_ylim(-90, 90)
    ax1.set_xlabel(f"Days from {args.date}")
    ax1.set_ylabel("Altitude [deg]   (0 = horizon)")
    ax1.set_title(f"Visibility from {site}", color=INK, fontsize=12, loc="left")
    leg = ax1.legend(loc="lower left", frameon=False, fontsize=8, ncol=3)
    for t in leg.get_texts():
        t.set_color(INK2)

    # --- Panel B: polar alt/az for the darkest hours of the centre night --
    ax2.remove()
    ax2 = fig.add_subplot(1, 2, 2, projection="polar")
    ax2.set_facecolor(SURFACE)
    mask = night & (np.abs(tt - t_mid) < 0.5)
    for b in alt:
        a, z = alt[b][mask], az[b][mask]
        vis = a > 0
        if not vis.any():
            continue
        ax2.plot(np.radians(z[vis]), 90 - a[vis], color=BODY_COLOUR[b],
                 lw=2.4 if b == "BH" else 1.4, zorder=5 if b == "BH" else 3)
        i = int(np.argmax(a))
        if a[i] > 0:
            ax2.annotate(b, (np.radians(z[i]), 90 - a[i]),
                         fontsize=8, color=BODY_COLOUR[b], weight="bold")
    ax2.set_theta_zero_location("N")
    ax2.set_theta_direction(-1)
    ax2.set_rlim(0, 90)
    ax2.set_rgrids([30, 60, 90], labels=["60°", "30°", "0°"], color=MUTED, fontsize=8)
    ax2.set_xticks(np.radians([0, 90, 180, 270]))
    ax2.set_xticklabels(["N", "E", "S", "W"], color=INK2, fontsize=9)
    ax2.grid(color=GRID, linewidth=0.8)
    ax2.set_title(f"Night sky, {args.date}\n(horizon at the rim)",
                  color=INK, fontsize=12, loc="left")

    fig.suptitle(
        f"Local sky during the flyby — {site} — black-hole periapsis July 2027",
        color=INK, fontsize=13, x=0.01, ha="left")
    fig.tight_layout(rect=(0, 0, 1, 0.94))
    fig.savefig(args.out, dpi=150, facecolor=SURFACE)
    plt.close(fig)
    print(f"Wrote {args.out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
