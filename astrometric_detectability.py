"""Could 19th-century astrometry have seen the black hole coming?

Answers two linked questions for the flyby scenario:

  horizon  How large an angular residual does the BH induce in each planet's
           position, as a function of how far away it is? Sets the distance at
           which period astrometry would have noticed it.

  vinf     Given a detection threshold on a given date, and a fixed periapsis
           date, what v_infinity does the BH need? The observability constraint
           and the arrival date together pin the approach speed.

  rates    Is the required speed astrophysically plausible? Compares dynamical
           friction and gas-drag timescales against the alternative explanation
           (selection from the low tail of the halo velocity distribution).

Method
------
Two N-body integrations, identical but for the BH, differenced. The observable
is the **angular** residual in a planet's heliocentric direction: 19th-century
meridian work measures directions and nothing else, so the radial part of the
perturbation is invisible and is excluded.

The residual is then detrended -- best-fit constant and linear terms removed.
An astronomer facing unexplained residuals refits the orbit, which absorbs a
constant into the mean longitude at epoch and a linear trend into the mean
motion; neither is evidence of a perturber. Only what survives is. This is a
conservative proxy: a full six-element refit would absorb more still, so the
detrended figures are an upper bound on what a real analysis would have kept.

Both detrended and raw values are reported.

Consistency with the engine
---------------------------
The BH's hyperbolic orbit is constructed exactly as in
``solar_system_bh_rebound26.py``, including its sign convention: that engine
writes ``r = a(e cosh F - 1)`` with ``a < 0``, so r is negative and the position
comes out antipodal, then negates the velocity ("inbound branch") to make the
state consistent. The net effect is the whole orbit rotated 180 degrees in its
plane -- the true argument of periapsis is ``bh_omega_deg + 180``. Reproduced
here so these numbers apply to the actual sweeps. Verified: toff=56210,
vinf=25 gives r_BH-Sun = 820.473 AU at epoch, matching the engine's printout.

Usage
-----
    python astrometric_detectability.py horizon
    python astrometric_detectability.py vinf --peri 2047-07-26 --when 1885
    python astrometric_detectability.py rates
"""

from __future__ import annotations

import argparse
import datetime as dt
import math

import numpy as np
import rebound
from skyfield.api import load
from skyfield.constants import AU_KM

G = 0.0002959122082855911          # AU^3 / (Msun day^2)
EPOCH = dt.datetime(1873, 9, 1)
RAD2AS = 206264.806
AU_PC = 4.84814e-6
G_PC = 4.30091e-3                  # pc/Msun (km/s)^2
KMS_PER_PCMYR = 1.02271

PLANET_MASSES = {
    "mercury": 1.651e-7, "venus": 2.447e-6, "earth": 3.003e-6, "mars": 3.227e-7,
    "jupiter": 9.5479e-4, "saturn": 2.8585e-4, "uranus": 4.3662e-5,
    "neptune": 5.1513e-5,
}
KEYS = {"mercury": "mercury", "venus": "venus", "earth": "earth",
        "mars": "mars barycenter", "jupiter": "jupiter barycenter",
        "saturn": "saturn barycenter", "uranus": "uranus barycenter",
        "neptune": "neptune barycenter"}
NAMES = list(PLANET_MASSES)

_STATE: dict | None = None


def states():
    """Heliocentric states of the planets at the 1873 epoch (cached)."""
    global _STATE
    if _STATE is None:
        eph = load("de440s.bsp")
        ts = load.timescale()
        t0 = ts.from_datetime(EPOCH.replace(tzinfo=dt.timezone.utc))
        sr = eph["sun"].at(t0).position.au
        sv = eph["sun"].at(t0).velocity.au_per_d
        _STATE = {n: (eph[KEYS[n]].at(t0).position.au - sr,
                      eph[KEYS[n]].at(t0).velocity.au_per_d - sv) for n in NAMES}
    return _STATE


def bh_state(toff_days, vinf_kms, rp_au, bh_mass=0.1):
    """BH heliocentric state at epoch for periapsis toff_days later.

    Omega = inc = omega = 0, so the perifocal frame is the reference frame.
    """
    mu = G * (1.0 + bh_mass)
    vinf = vinf_kms / (AU_KM / 86400.0)
    a = -mu / vinf ** 2
    e = 1.0 + rp_au / abs(a)
    n = math.sqrt(mu / abs(a) ** 3)
    M = -n * toff_days
    F = math.asinh(M / e)
    for _ in range(300):
        f = e * math.sinh(F) - F - M
        dF = -f / (e * math.cosh(F) - 1.0)
        F += math.copysign(min(abs(dF), 1.0), dF)
        if abs(dF) < 1e-14:
            break
    cF, sF = math.cosh(F), math.sinh(F)
    nu = math.atan2(math.sqrt(e * e - 1.0) * sF / (e * cF - 1.0),
                    (e - cF) / (e * cF - 1.0))
    r = a * (e * cF - 1.0)                       # negative, as in the engine
    h = math.sqrt(mu * abs(a) * (e * e - 1.0))
    x = r * np.array([math.cos(nu), math.sin(nu), 0.0])
    v = (mu / h) * np.array([-math.sin(nu), e + math.cos(nu), 0.0])
    return x, -v                                 # engine's "inbound branch"


def r_at(tau_days, vinf_kms, rp_au, bh_mass=0.1):
    """Heliocentric distance [AU] tau_days BEFORE periapsis."""
    return float(np.linalg.norm(bh_state(tau_days, vinf_kms, rp_au, bh_mass)[0]))


def _sim(toff, vinf, rp, with_bh):
    sim = rebound.Simulation()
    sim.G = G
    sim.integrator = "ias15"
    sim.add(m=1.0, x=0, y=0, z=0, vx=0, vy=0, vz=0)
    for n in NAMES:
        r, v = states()[n]
        sim.add(m=PLANET_MASSES[n], x=r[0], y=r[1], z=r[2],
                vx=v[0], vy=v[1], vz=v[2])
    if with_bh:
        x, v = bh_state(toff, vinf, rp)
        sim.add(m=0.1, x=x[0], y=x[1], z=x[2], vx=v[0], vy=v[1], vz=v[2])
    return sim


def _track(toff, vinf, rp, with_bh, times):
    sim = _sim(toff, vinf, rp, with_bh)
    pos = np.zeros((len(times), len(NAMES), 3))
    dbh = np.zeros(len(times))
    for i, t in enumerate(times):
        sim.integrate(float(t))
        s = sim.particles[0]
        for j in range(len(NAMES)):
            p = sim.particles[j + 1]
            pos[i, j] = (p.x - s.x, p.y - s.y, p.z - s.z)
        if with_bh:
            b = sim.particles[-1]
            dbh[i] = math.dist((b.x, b.y, b.z), (s.x, s.y, s.z))
    return pos, dbh


def _detrend(y, t):
    """Remove best-fit constant + linear -- what an orbit refit absorbs."""
    A = np.c_[np.ones_like(t), t]
    return y - A @ np.linalg.lstsq(A, y, rcond=None)[0]


def residuals(toff, vinf, rp, times, base=None):
    """Peak angular residual [arcsec] per planet: (detrended, raw)."""
    if base is None:
        base, _ = _track(0, vinf, rp, False, times)
    pert, dbh = _track(toff, vinf, rp, True, times)
    det, raw = {}, {}
    for j, n in enumerate(NAMES):
        b = base[:, j]
        u0 = b / np.linalg.norm(b, axis=1)[:, None]
        u1 = pert[:, j] / np.linalg.norm(pert[:, j], axis=1)[:, None]
        d = u1 - u0
        # Resolve the offset onto two sky-plane axes so each is detrended alone.
        e1 = np.cross(u0, np.cross(b, np.gradient(b, axis=0)))
        e1 /= np.linalg.norm(e1, axis=1)[:, None]
        e2 = np.cross(u0, e1)
        c1, c2 = (d * e1).sum(1), (d * e2).sum(1)
        raw[n] = float(np.hypot(c1, c2).max() * RAD2AS)
        det[n] = float(np.hypot(_detrend(c1, times),
                                _detrend(c2, times)).max() * RAD2AS)
    return det, raw, dbh


def cmd_horizon(args):
    times = np.linspace(0.0, args.window * 365.25, 361)
    base, _ = _track(0, args.vinf, args.rp, False, times)
    print(f"Peak DETRENDED angular residual over a {args.window:.0f}-year window "
          f"[arcsec], v_inf={args.vinf} km/s\n")
    hdr = f"{'D_BH mid':>9}  " + "  ".join(f"{n[:7]:>8}" for n in NAMES)
    print(hdr + "\n" + f"{'[AU]':>9}  " + "  ".join(f"{'[as]':>8}" for _ in NAMES))
    rows = []
    for toff in args.toff:
        det, raw, dbh = residuals(toff, args.vinf, args.rp, times, base)
        dmid = dbh[len(dbh) // 2]
        rows.append((dmid, det, raw))
        print(f"{dmid:>9.1f}  " + "  ".join(f"{det[n]:>8.2e}" for n in NAMES))
    print("\nSame, RAW (no refit absorbed; upper bound):")
    print(hdr)
    for dmid, det, raw in rows:
        print(f"{dmid:>9.1f}  " + "  ".join(f"{raw[n]:>8.2e}" for n in NAMES))


def cmd_vinf(args):
    peri = dt.datetime.strptime(args.peri, "%Y-%m-%d")
    end = dt.datetime(args.when, 1, 1)
    t0 = (dt.datetime(args.when - int(args.window), 1, 1) - EPOCH).days
    t1 = (end - EPOCH).days
    times = np.linspace(t0, t1, 361)
    toff = (peri - EPOCH).days
    print(f"periapsis {peri:%Y-%m-%d}  ->  bh_tperi_offset_days = {toff}")
    print(f"{args.window:.0f}-year window ending {args.when}; "
          f"body = {args.body}; rp = {args.rp} AU\n")
    base, _ = _track(0, 25.0, args.rp, False, times)
    print(f"{'v_inf':>7} {'r@start':>9} {'r@end':>9} {'detrended':>11} {'raw':>9}")
    for v in args.vinf_grid:
        det, raw, dbh = residuals(toff, v, args.rp, times, base)
        print(f"{v:>7.1f} {dbh[0]:>9.1f} {dbh[-1]:>9.1f} "
              f"{det[args.body]:>10.2f}\" {raw[args.body]:>8.2f}\"")
    print(f"\nr_BH at the {EPOCH:%Y-%m-%d} epoch:")
    for v in args.vinf_grid:
        print(f"  v_inf={v:5.1f} km/s -> {r_at(toff, v, args.rp):7.1f} AU")
    print("\nSensitivity of the constraint to periapsis distance "
          "(rp is negligible against ~400 AU):")
    tau = (peri - end).days
    for rp in (0.25, 0.5, 1.0, 1.5):
        print(f"  rp={rp:4.2f} AU -> r({args.when}) = "
              f"{r_at(tau, 10.0, rp):.1f} AU at v_inf=10")


def cmd_rates(args):
    M = args.mass

    def t_df(v, rho, lnL=10.0):
        a = 4 * np.pi * G_PC ** 2 * M * rho * lnL / v ** 2
        return (v / a) / KMS_PER_PCMYR                     # Myr

    hub = 13.8e3
    print(f"Timescale to change v by ~{args.v:.0f} km/s for a {M} Msun object:\n")
    for lab, rho, lnL in (("dynamical friction, halo (0.01 Msun/pc3)", 0.01, 10.0),
                          ("dynamical friction, disk midplane (0.1)", 0.1, 10.0),
                          ("dynamical friction, dense GMC (100)", 100.0, 10.0),
                          ("Bondi-Hoyle gas drag, ISM (0.03)", 0.03, 1.0)):
        t = t_df(args.v, rho, lnL)
        print(f"  {lab:42s} {t:.2e} Myr = {t / hub:.1e} Hubble times")
    print("  (friction scales as 1/M: ~1e6 Msun in the disk for a Hubble-time brake)")

    sig, vsun = args.sigma, args.vsun
    v = np.linspace(0, args.vcut, 4000)
    f = (v / (np.sqrt(2 * np.pi) * sig * vsun)
         * (np.exp(-(v - vsun) ** 2 / (2 * sig ** 2))
            - np.exp(-(v + vsun) ** 2 / (2 * sig ** 2))))
    P = float(np.trapezoid(f, v))
    print(f"\nSelection instead of deceleration (sigma={sig}, v_sun={vsun} km/s):")
    print(f"  P(|v_rel| < {args.vcut:.0f} km/s) = {P:.2e}  (1 in {1 / P:,.0f})")

    def vsigma(vv, rp_au):
        vesc = np.sqrt(2 * G_PC / (rp_au * AU_PC))
        return vv * np.pi * (rp_au * AU_PC) ** 2 * (1.0 + vesc ** 2 / vv ** 2)

    print("\n  Slow objects are rare but focused, so the RATE is barely suppressed:")
    for rp in (0.5, 1.5):
        vesc = np.sqrt(2 * G_PC / (rp * AU_PC))
        print(f"    rp={rp} AU: v_esc={vesc:.1f} km/s, focusing x"
              f"{1 + vesc ** 2 / args.vcut ** 2:.0f}, "
              f"rate({args.vcut:.0f})/rate(270) = "
              f"{vsigma(args.vcut, rp) / vsigma(270.0, rp):.2f}")

    n = args.rho_dm / M
    rp = args.rp
    slow = n * P * vsigma(args.vcut, rp) * KMS_PER_PCMYR
    allv = n * vsigma(270.0, rp) * KMS_PER_PCMYR
    print(f"\n  Encounters within {rp} AU over 10 Gyr, if these are 100% of dark matter:")
    print(f"    any speed          : {allv * 10e3:.1e}")
    print(f"    and v_inf < {args.vcut:.0f} km/s : {slow * 10e3:.1e}")


def main() -> int:
    p = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    sub = p.add_subparsers(dest="cmd", required=True)

    h = sub.add_parser("horizon", help="Residual vs BH distance, all planets.")
    h.add_argument("--window", type=float, default=30.0, help="Baseline [yr].")
    h.add_argument("--vinf", type=float, default=25.0)
    h.add_argument("--rp", type=float, default=0.5)
    h.add_argument("--toff", type=float, nargs="+",
                   default=[110000, 90000, 70000, 56210, 45000, 35000, 28000,
                            22000, 17000])
    h.set_defaults(func=cmd_horizon)

    v = sub.add_parser("vinf", help="Solve v_inf from an observability constraint.")
    v.add_argument("--peri", default="2047-07-26", help="Periapsis date.")
    v.add_argument("--when", type=int, default=1885, help="Year the limit applies.")
    v.add_argument("--window", type=float, default=30.0, help="Baseline [yr].")
    v.add_argument("--body", default="uranus", choices=NAMES)
    v.add_argument("--rp", type=float, default=0.5)
    v.add_argument("--vinf-grid", type=float, nargs="+",
                   default=[8, 10, 12, 14, 16, 19, 22, 25])
    v.set_defaults(func=cmd_vinf)

    r = sub.add_parser("rates", help="Astrophysical plausibility of a slow BH.")
    r.add_argument("--mass", type=float, default=0.1, help="BH mass [Msun].")
    r.add_argument("--v", type=float, default=200.0, help="Halo speed [km/s].")
    r.add_argument("--vcut", type=float, default=10.0, help="Required v_inf [km/s].")
    r.add_argument("--sigma", type=float, default=156.0, help="Halo dispersion.")
    r.add_argument("--vsun", type=float, default=232.0, help="Sun through halo.")
    r.add_argument("--rho-dm", type=float, default=0.01, help="Local DM [Msun/pc3].")
    r.add_argument("--rp", type=float, default=1.5, help="Encounter radius [AU].")
    r.set_defaults(func=cmd_rates)

    args = p.parse_args()
    args.func(args)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
