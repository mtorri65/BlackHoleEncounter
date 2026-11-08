"""Earth -> Mars launch windows from one run's orbit log (a porkchop search).

For every departure day and every time of flight, solves Lambert's problem
(Sun-centred, single revolution, prograde with respect to Earth's orbit normal)
between Earth on the departure day and Mars on the arrival day, using the
positions and velocities the core engine logged. A transfer is costed as

    burn from a 300 km circular Earth orbit to the departure v_inf
  + burn from the arrival v_inf into a 300 km circular Mars orbit

which is the usual all-propulsive figure; aerobraking at Mars would cut the
second term. For each departure day the cheapest time of flight is kept, and a
"window" is a local minimum of that curve, reported with the span of dates
within --margin km/s of it.

The transfer arc feels only the Sun. That is exact enough once the hole is far
away (it is 14 AU out by late 2049 in the adopted run), but a departure in
2045-2047 would fly through the encounter, and its cost here is only indicative.

Needs a run logged at a uniform 1-day cadence: set ``output_interval_days: 1``
and ``output_dense_window_days: 0`` in the YAML. The engine's default two-rate
log is too coarse away from periapsis.

Usage:
    python launch_windows.py <run folder | orbits .parquet> --from 2048-01-01 --to 2200-01-01
    python launch_windows.py <run> --from 2053-07-01 --to 2053-08-01 --dep-step 1 --tof-step 1
    python launch_windows.py <run> --from 2025-06-01 --to 2047-01-01 --tof-max 600 --max-dv 7.5

Writes <out>_daily.csv (best transfer for each departure day) and
<out>_windows.csv, and prints the window table.
"""
import argparse
from datetime import datetime, timedelta, timezone
from pathlib import Path

import numpy as np
import pandas as pd
import yaml

MU = 0.0002959122082855911            # G*Msun [AU^3/day^2], as in the engine
AUD = 149597870.7 / 86400.0           # AU/day -> km/s
MU_E, R_E = 398600.4418, 6378.137 + 300.0   # km^3/s^2, km (300 km parking orbit)
MU_M, R_M = 42828.37, 3389.5 + 300.0


def load_run(arg):
    """Return (parquet path, epoch) from a run folder or a parquet file in one."""
    p = Path(arg)
    pq = p if p.suffix == ".parquet" else next(p.glob("*__orbits__*.parquet"))
    cfg = next(pq.parent.glob("*__input.yaml"))
    epoch_s = str(yaml.safe_load(open(cfg, encoding="utf-8"))["epoch"]).strip()
    epoch = datetime.fromisoformat(epoch_s.replace("Z", "+00:00"))
    if epoch.tzinfo is None:
        epoch = epoch.replace(tzinfo=timezone.utc)
    return pq, epoch


def stumpff(z):
    C = np.empty_like(z); S = np.empty_like(z)
    p, n, s = z > 1e-8, z < -1e-8, np.abs(z) <= 1e-8
    sp, sn = np.sqrt(z[p]), np.sqrt(-z[n])
    C[p] = (1 - np.cos(sp)) / z[p];   S[p] = (sp - np.sin(sp)) / sp**3
    C[n] = (np.cosh(sn) - 1) / -z[n]; S[n] = (np.sinh(sn) - sn) / sn**3
    C[s] = 0.5 - z[s] / 24;           S[s] = 1 / 6 - z[s] / 120
    return C, S


def lambert(r1, r2, dt, hdir):
    """Vectorised universal-variable Lambert solver (Curtis, Alg. 5.2).

    Solved by bisection on z rather than Newton: slower, but it cannot diverge,
    which matters when millions of cases are solved blind. Returns v1, v2 and a
    mask of cases that did not converge.
    """
    n1, n2 = np.linalg.norm(r1, axis=1), np.linalg.norm(r2, axis=1)
    dth = np.arccos(np.clip(np.einsum("ij,ij->i", r1, r2) / (n1 * n2), -1, 1))
    retro = np.einsum("ij,ij->i", np.cross(r1, r2), hdir) < 0
    dth = np.where(retro, 2 * np.pi - dth, dth)
    A = np.sin(dth) * np.sqrt(n1 * n2 / (1 - np.cos(dth)))
    lo = np.full_like(dt, -400.0)
    hi = np.full_like(dt, 4 * np.pi**2 - 1e-6)
    for _ in range(70):
        z = 0.5 * (lo + hi)
        C, S = stumpff(z)
        y = n1 + n2 + A * (z * S - 1) / np.sqrt(C)
        yp = np.maximum(y, 1e-300)
        F = (yp / C) ** 1.5 * S + A * np.sqrt(yp) - np.sqrt(MU) * dt
        F = np.where(y < 0, -1.0, F)            # y < 0: z too small
        lo = np.where(F < 0, z, lo)
        hi = np.where(F >= 0, z, hi)
    z = 0.5 * (lo + hi)
    C, S = stumpff(z)
    y = n1 + n2 + A * (z * S - 1) / np.sqrt(C)
    f = 1 - y / n1; g = A * np.sqrt(y / MU); gd = 1 - y / n2
    v1 = (r2 - f[:, None] * r1) / g[:, None]
    v2 = (gd[:, None] * r2 - r1) / g[:, None]
    bad = (y <= 0) | (hi - lo > 1e-6) | ~np.isfinite(g)
    return v1, v2, bad


def burn(vinf, mu, r):
    """Burn [km/s] between a circular orbit of radius r and a hyperbola of excess vinf."""
    return np.sqrt(vinf**2 + 2 * mu / r) - np.sqrt(mu / r)


def search(states, epoch, d_from, d_to, dep_step, tof):
    t, rE, vE, rM, vM, rB = states
    t0 = int(t[0])
    deps = np.arange(d_from, d_to, dep_step)
    deps = deps[deps + tof.max() <= t[-1]]
    D, T = np.meshgrid(deps, tof, indexing="ij")
    D, T = D.ravel(), T.ravel()
    i1, i2 = (D - t0).astype(int), (D + T - t0).astype(int)
    h = np.cross(rE[i1], vE[i1]); h /= np.linalg.norm(h, axis=1)[:, None]
    res = np.empty((len(D), 3))
    for s in range(0, len(D), 400_000):
        sl = slice(s, s + 400_000)
        v1, v2, bad = lambert(rE[i1[sl]], rM[i2[sl]], T[sl], h[sl])
        vd = np.linalg.norm(v1 - vE[i1[sl]], axis=1) * AUD
        va = np.linalg.norm(v2 - vM[i2[sl]], axis=1) * AUD
        tot = burn(vd, MU_E, R_E) + burn(va, MU_M, R_M)
        tot[bad] = np.inf
        res[sl] = np.c_[vd, va, tot]
    shape = (len(deps), len(tof))
    tot = res[:, 2].reshape(shape)
    k = np.argmin(tot, axis=1); j = np.arange(len(deps))
    best = pd.DataFrame({
        "t_dep": deps, "tof_days": tof[k], "dv_total_kms": tot[j, k],
        "vinf_dep_kms": res[:, 0].reshape(shape)[j, k],
        "vinf_arr_kms": res[:, 1].reshape(shape)[j, k],
    })
    i = (best.t_dep - t0).astype(int).to_numpy()
    lonE = np.degrees(np.arctan2(rE[i, 1], rE[i, 0]))
    lonM = np.degrees(np.arctan2(rM[i, 1], rM[i, 0]))
    # Mars's heliocentric longitude minus Earth's: > 0 means Mars leads Earth.
    best["mars_lead_deg"] = (lonM - lonE + 180) % 360 - 180
    best["r_earth_au"] = np.linalg.norm(rE[i], axis=1)
    best["r_mars_au"] = np.linalg.norm(rM[i], axis=1)
    best["r_bh_au"] = np.linalg.norm(rB[i], axis=1)
    day = lambda d: (epoch + timedelta(days=float(d))).date().isoformat()
    best["depart"] = [day(d) for d in best.t_dep]
    best["arrive"] = [day(d + T) for d, T in zip(best.t_dep, best.tof_days)]
    return best


def windows(best, margin, dep_step, half_width_days=120):
    """Local minima of the best-cost curve, each with the span within `margin`."""
    v = best.dv_total_kms.to_numpy()
    w = max(1, half_width_days // dep_step)
    rows = []
    for i in range(1, len(v) - 1):
        lo, hi = max(0, i - w), min(len(v), i + w + 1)
        if np.isfinite(v[i]) and v[i] == v[lo:hi].min():
            a = i
            while a > 0 and v[a - 1] <= v[i] + margin:
                a -= 1
            b = i
            while b < len(v) - 1 and v[b + 1] <= v[i] + margin:
                b += 1
            r = best.iloc[i].to_dict()
            r["opens"], r["closes"] = best.depart.iloc[a], best.depart.iloc[b]
            r["width_days"] = (b - a) * dep_step
            rows.append(r)
    return pd.DataFrame(rows)


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("run", help="run folder, or its __orbits__*.parquet")
    ap.add_argument("--from", dest="d_from", required=True, help="first departure date, YYYY-MM-DD")
    ap.add_argument("--to", dest="d_to", required=True, help="last departure date, YYYY-MM-DD")
    ap.add_argument("--dep-step", type=int, default=2, help="departure grid step [days]")
    ap.add_argument("--tof-min", type=int, default=30)
    ap.add_argument("--tof-max", type=int, default=900, help="longest flight considered [days]")
    ap.add_argument("--tof-step", type=int, default=4)
    ap.add_argument("--margin", type=float, default=0.5,
                    help="window span = dates within this many km/s of the minimum")
    ap.add_argument("--max-dv", type=float, default=np.inf,
                    help="drop minima costlier than this (the shallow dips between windows)")
    ap.add_argument("--out", default=None, help="output prefix (default: next to the parquet)")
    a = ap.parse_args()

    pq, epoch = load_run(a.run)
    df = pd.read_parquet(pq)

    def body(nm):
        b = df[df.body == nm].sort_values("t_days")
        return (b.t_days.to_numpy(), b[["x_au", "y_au", "z_au"]].to_numpy(float),
                b[["vx", "vy", "vz"]].to_numpy(float))

    t, rS, vS = body("Sun")
    _, rE, vE = body("Earth")
    _, rM, vM = body("Mars")
    _, rB, _ = body("BH")
    if not np.allclose(np.diff(t), 1.0):
        raise SystemExit("This run is not logged at a uniform 1-day cadence; see the docstring.")
    # The engine logs in its own frame; the Sun drifts in it after the encounter.
    states = (t, rE - rS, vE - vS, rM - rS, vM - vS, rB - rS)

    day = lambda s: (datetime.fromisoformat(s).replace(tzinfo=timezone.utc) - epoch).days
    tof = np.arange(a.tof_min, a.tof_max + 1, a.tof_step, dtype=float)
    best = search(states, epoch, day(a.d_from), day(a.d_to), a.dep_step, tof)
    W = windows(best, a.margin, a.dep_step)
    if len(W):
        W = W[W.dv_total_kms <= a.max_dv]

    out = a.out or str(pq.parent / f"launch_windows_{a.d_from}_{a.d_to}")
    best.drop(columns="t_dep").to_csv(out + "_daily.csv", index=False)
    W.drop(columns="t_dep", errors="ignore").to_csv(out + "_windows.csv", index=False)

    cols = ["depart", "opens", "closes", "width_days", "tof_days", "arrive", "dv_total_kms",
            "vinf_dep_kms", "vinf_arr_kms", "mars_lead_deg", "r_bh_au"]
    with pd.option_context("display.width", 250, "display.max_rows", 1000):
        print(W[cols].round(3).to_string(index=False) if len(W) else "(no windows)")
    if len(W) > 1:
        print("spacing [yr]:", np.round(np.diff(W.t_dep.to_numpy()) / 365.25, 2))
    print(f"\nwrote {out}_daily.csv and {out}_windows.csv")


if __name__ == "__main__":
    main()
