"""Recover Mars's post-flyby orbital elements from a BH sweep.

Companion to ``extract_earth_elements.py``. The method is identical --- the
simulation runs in the equatorial J2000 frame, and REBOUND applies no torque to
a point mass, so each planet's spin axis is fixed in inertial space while the
black hole tilts its orbital plane. Given that axis:

    obliquity   eps   = angle(spin axis, orbit normal)
    equinox dir e_eq  = h_hat x s_perp
    lambda_p          = angle from e_eq to the eccentricity vector, about h

The one difference is the axis itself. For Earth the frame's z-axis *is* the
spin axis, which is what made the original recovery so clean. Mars's pole points
elsewhere, so it is constructed from the IAU right ascension and declination.

Validation: applied at t = 0 this returns a = 1.5237 AU, e = 0.0934,
obliquity = 25.18 deg against Mars's true 1.524 / 0.0934 / 25.19.

Reads the Parquet tree written by ``convert_orbits_to_parquet.py`` rather than
the original workbooks, which makes the whole sweep a few minutes.

Usage
-----
    python extract_mars_elements.py simulations/<STAMP>_parquet --workers 5
"""

from __future__ import annotations

import argparse
import glob
import os
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path

import numpy as np
import pandas as pd

GM_SUN = 2.959122082855911e-4      # AU^3 / day^2
DAYS_PER_YEAR_EARTH = 365.256363

# IAU 2015 north-pole orientations, equatorial J2000 (right ascension, declination).
IAU_POLES = {
    "Mars": (317.68143, 52.88650),
    "Mercury": (281.0103, 61.4155),
    "Venus": (272.76, 67.16),
    "Jupiter": (268.056595, 64.495303),
    "Saturn": (40.589, 83.537),
    "Uranus": (257.311, -15.175),
    "Neptune": (299.36, 43.46),
}
# Sidereal year today [days], used to scale the post-flyby year via Kepler III.
YEAR_TODAY = {"Mars": 686.98, "Mercury": 87.97, "Venus": 224.70,
              "Jupiter": 4332.6, "Saturn": 10759.2,
              "Uranus": 30685.4, "Neptune": 60189.0}
A_TODAY = {"Mars": 1.5237, "Mercury": 0.3871, "Venus": 0.7233,
           "Jupiter": 5.2044, "Saturn": 9.5826,
           "Uranus": 19.2184, "Neptune": 30.1104}


def spin_axis(body: str) -> np.ndarray:
    """Unit spin-axis vector in the equatorial J2000 frame, from the IAU pole."""
    ra, dec = np.radians(IAU_POLES[body])
    return np.array([np.cos(dec) * np.cos(ra), np.cos(dec) * np.sin(ra), np.sin(dec)])


def elements_from_state(r: np.ndarray, v: np.ndarray, spin: np.ndarray) -> dict:
    """Heliocentric a, e, obliquity and longitude of perihelion."""
    h = np.cross(r, v)
    h_norm, r_norm = np.linalg.norm(h), np.linalg.norm(r)
    if h_norm == 0 or r_norm == 0:
        return dict(a_au=np.nan, ecc=np.nan, obliquity_deg=np.nan,
                    lon_perihelion_deg=np.nan)
    h_hat = h / h_norm
    e_vec = np.cross(v, h) / GM_SUN - r / r_norm
    ecc = float(np.linalg.norm(e_vec))
    inv_a = 2.0 / r_norm - v.dot(v) / GM_SUN
    a = float(1.0 / inv_a) if abs(inv_a) > 1e-30 else np.inf

    cos_eps = float(np.clip(h_hat.dot(spin), -1.0, 1.0))
    eps = float(np.degrees(np.arccos(cos_eps)))

    s_perp = spin - cos_eps * h_hat
    s_norm = np.linalg.norm(s_perp)
    if s_norm < 1e-12 or ecc < 1e-12:
        lam_p = np.nan                     # degenerate: zero obliquity or circular
    else:
        e_eq = np.cross(h_hat, s_perp / s_norm)
        e_hat = e_vec / ecc
        lam_p = float(np.degrees(np.arctan2(
            np.cross(e_eq, e_hat).dot(h_hat), e_eq.dot(e_hat))) % 360.0)
    return dict(a_au=a, ecc=ecc, obliquity_deg=eps, lon_perihelion_deg=lam_p)


def process_run(args) -> dict | None:
    path, body = args
    spin = spin_axis(body)
    try:
        d = pd.read_parquet(path)
    except Exception as exc:                       # noqa: BLE001
        return {"run": os.path.basename(os.path.dirname(path)),
                "error": f"{type(exc).__name__}: {exc}"}
    row = {"run": os.path.basename(os.path.dirname(path))}
    for when, t in (("before", d.t_days.min()), ("after", d.t_days.max())):
        e = d[d.t_days == t].set_index("body")
        if body not in e.index or "Sun" not in e.index:
            return {**row, "error": f"missing {body} or Sun"}
        P = lambda b: e.loc[b, ["x_au", "y_au", "z_au"]].values.astype(float)
        V = lambda b: e.loc[b, ["vx", "vy", "vz"]].values.astype(float)
        el = elements_from_state(P(body) - P("Sun"), V(body) - V("Sun"), spin)
        row.update({f"{k}_{when}": v for k, v in el.items()})

    a = row["a_au_after"]
    bound = bool(np.isfinite(a) and a > 0 and row["ecc_after"] < 1.0)
    row["bound_after"] = bound
    # Kepler III, scaled from the body's present-day year.
    row["days_per_year_after"] = (
        YEAR_TODAY[body] * (a / A_TODAY[body]) ** 1.5 if bound else np.nan)
    return row


def main() -> int:
    p = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    p.add_argument("parquet_dir", type=Path,
                   help="Parquet tree, e.g. simulations/<STAMP>_parquet")
    p.add_argument("--body", default="Mars", choices=sorted(IAU_POLES))
    p.add_argument("--out", type=Path, default=None)
    p.add_argument("--workers", type=int, default=1)
    args = p.parse_args()

    files = sorted(glob.glob(str(args.parquet_dir / "*" / "orbits.parquet")))
    if not files:
        raise SystemExit(f"No orbits.parquet found under {args.parquet_dir}")
    jobs = [(f, args.body) for f in files]

    rows = []
    if args.workers <= 1:
        for i, j in enumerate(jobs, 1):
            r = process_run(j)
            if r:
                rows.append(r)
            if i % 50 == 0 or i == len(jobs):
                print(f"  {i}/{len(jobs)}", flush=True)
    else:
        with ProcessPoolExecutor(max_workers=args.workers) as ex:
            futs = [ex.submit(process_run, j) for j in jobs]
            for i, f in enumerate(as_completed(futs), 1):
                r = f.result()
                if r:
                    rows.append(r)
                if i % 50 == 0 or i == len(jobs):
                    print(f"  {i}/{len(jobs)}", flush=True)

    df = pd.DataFrame(rows).sort_values("run").reset_index(drop=True)
    out = args.out or args.parquet_dir.parent / (
        args.parquet_dir.name.replace("_parquet", "") + f"_{args.body.lower()}_elements.csv")
    df.to_csv(out, index=False)
    n_err = int(df["error"].notna().sum()) if "error" in df.columns else 0
    print(f"\n{args.body}: {len(df)} runs ({n_err} errors); "
          f"bound after flyby: {int(df['bound_after'].sum())}")
    print(f"Wrote {out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
