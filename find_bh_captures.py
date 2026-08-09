"""Which planets leave with the black hole, rather than merely leaving?

An ejection and a capture look identical in the per-run deltas: both show a
planet no longer bound to the Sun. They are very different outcomes. A capture
means the planet is on a closed orbit around the black hole and departs the
system as its satellite; a free ejection means it drifts off alone.

Telling them apart needs the final state vectors, which the deltas CSV does not
carry, so this reads the last logged timestep of each run's trajectory Parquet
and evaluates the two-body energy of every body against both attractors:

    eps = v^2/2 - mu/r        bound when negative

Classification, per body:

    solar      bound to the Sun (the normal case)
    CAPTURED   bound to the black hole and not to the Sun
    free       bound to neither -- ejected alone into interstellar space
    ambiguous  bound to both by this test; the two-body approximation has
               broken down and the run needs looking at individually

A caveat worth keeping in mind: "bound" here is a two-body energy at one
instant, 154 years after the encounter in the default configuration. A wide,
weakly bound capture may not survive the next passing star, and this test cannot
see that. It is a statement about the end of the integration, not about
permanence.

Usage
-----
    python find_bh_captures.py simulations/<STAMP> --workers 5
    python find_bh_captures.py simulations/<STAMP> --only-captures
"""

from __future__ import annotations

import argparse
import glob
import os
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path

import numpy as np
import pandas as pd

G = 2.959122082855911e-4          # AU^3 / (Msun day^2)
M_SUN = 1.0

# Masses [Msun]; negligible against Sun and BH but free to include.
MASSES = {
    "Mercury": 1.651e-7, "Venus": 2.447e-6, "Earth": 3.003e-6, "Mars": 3.227e-7,
    "Jupiter": 9.5479e-4, "Saturn": 2.8585e-4, "Uranus": 4.3662e-5,
    "Neptune": 5.1513e-5, "Moon": 3.6923e-8,
}
COLS = ["body", "t_days", "x_au", "y_au", "z_au", "vx", "vy", "vz"]


def _elements(dr, dv, mu):
    """(bound, a, e) of a relative state under standard gravitational parameter mu."""
    r = float(np.linalg.norm(dr))
    v2 = float(dv.dot(dv))
    if r == 0:
        return False, np.nan, np.nan
    eps = 0.5 * v2 - mu / r
    if eps >= 0:
        return False, np.nan, np.nan
    a = -mu / (2.0 * eps)
    h = np.cross(dr, dv)
    e = float(np.sqrt(max(0.0, 1.0 - h.dot(h) / (mu * a))))
    return True, a, e


def process_run(args) -> list[dict]:
    path, bh_mass = args
    run = os.path.basename(os.path.dirname(path))
    try:
        d = pd.read_parquet(path, columns=COLS)
    except Exception as exc:                       # noqa: BLE001 - record, don't abort
        return [{"run": run, "body": None, "error": f"{type(exc).__name__}: {exc}"}]

    last = d[d.t_days == d.t_days.max()].set_index("body")
    if "Sun" not in last.index or "BH" not in last.index:
        return [{"run": run, "body": None, "error": "missing Sun or BH"}]

    def state(b):
        row = last.loc[b]
        return (np.array([row.x_au, row.y_au, row.z_au], dtype=float),
                np.array([row.vx, row.vy, row.vz], dtype=float))

    r_sun, v_sun = state("Sun")
    r_bh, v_bh = state("BH")
    out = []
    for body in MASSES:
        if body not in last.index:
            continue
        r, v = state(body)
        m = MASSES[body]
        b_sun, a_sun, e_sun = _elements(r - r_sun, v - v_sun, G * (M_SUN + m))
        b_bh, a_bh, e_bh = _elements(r - r_bh, v - v_bh, G * (bh_mass + m))
        if b_sun and not b_bh:
            state_lbl = "solar"
        elif b_bh and not b_sun:
            state_lbl = "CAPTURED"
        elif not b_sun and not b_bh:
            state_lbl = "free"
        else:
            state_lbl = "ambiguous"
        rec = {"run": run, "body": body, "state": state_lbl,
               "a_sun_au": a_sun, "e_sun": e_sun,
               "a_bh_au": a_bh, "e_bh": e_bh,
               "r_from_bh_au": float(np.linalg.norm(r - r_bh))}
        if b_bh:
            # Kepler III about the black hole.
            rec["period_bh_days"] = 2 * np.pi * np.sqrt(
                a_bh ** 3 / (G * (bh_mass + m)))
        out.append(rec)
    return out


def main() -> int:
    p = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    p.add_argument("parent_dir", type=Path, help="Sweep folder, e.g. simulations/<STAMP>")
    p.add_argument("--bh-mass", type=float, default=0.1, help="BH mass [Msun].")
    p.add_argument("--out", type=Path, default=None)
    p.add_argument("--workers", type=int, default=1)
    p.add_argument("--only-captures", action="store_true",
                   help="Print only the runs with a capture.")
    args = p.parse_args()

    files = sorted(set(glob.glob(str(args.parent_dir / "*" / "orbits.parquet")))
                   | set(glob.glob(str(args.parent_dir / "*" / "*orbits*.parquet"))))
    if not files:
        raise SystemExit(f"No orbits Parquet found under {args.parent_dir}")
    jobs = [(f, args.bh_mass) for f in files]

    rows = []
    if args.workers <= 1:
        for i, j in enumerate(jobs, 1):
            rows += process_run(j)
            if i % 50 == 0 or i == len(jobs):
                print(f"  {i}/{len(jobs)}", flush=True)
    else:
        with ProcessPoolExecutor(max_workers=args.workers) as ex:
            futs = [ex.submit(process_run, j) for j in jobs]
            for i, f in enumerate(as_completed(futs), 1):
                rows += f.result()
                if i % 50 == 0 or i == len(jobs):
                    print(f"  {i}/{len(jobs)}", flush=True)

    df = pd.DataFrame(rows).sort_values(["run", "body"]).reset_index(drop=True)
    out = args.out or args.parent_dir.parent / f"{args.parent_dir.name}_bh_captures.csv"
    df.to_csv(out, index=False)

    ok = df[df.get("state").notna()] if "state" in df.columns else df.iloc[:0]
    print(f"\n{len(files)} runs, {len(ok)} body-outcomes\n")
    for lbl in ("solar", "free", "CAPTURED", "ambiguous"):
        n = int((ok.state == lbl).sum())
        print(f"  {lbl:10s} {n:>6}")
    cap = ok[ok.state == "CAPTURED"]
    if len(cap):
        print(f"\nCaptures: {len(cap)} across "
              f"{cap.run.nunique()} runs ({100 * cap.run.nunique() / len(files):.1f}%)")
        print(f"\n{'run':<58}{'body':>9}{'a_bh':>9}{'e_bh':>7}{'period':>11}")
        for _, r in cap.sort_values("a_bh_au").iterrows():
            per = r.get("period_bh_days", np.nan)
            per_s = f"{per/365.25:.1f} yr" if np.isfinite(per) else "-"
            print(f"{r.run[-58:]:<58}{r.body:>9}{r.a_bh_au:>9.3f}{r.e_bh:>7.3f}{per_s:>11}")
    print(f"\nWrote {out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
