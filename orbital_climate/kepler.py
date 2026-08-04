"""Kepler's equation and orbital geometry.

Given the mean anomaly ``M`` and eccentricity ``e`` (elliptical, 0 <= e < 1)
this module solves for the eccentric anomaly ``E`` via Newton iteration on

    E - e * sin(E) = M

and derives the true anomaly ``nu`` and heliocentric radius ``r``:

    r = a (1 - e cos E)
    tan(nu/2) = sqrt((1+e)/(1-e)) tan(E/2)

All angles are in radians. Functions are vectorised over ``M`` via numpy.
"""

from __future__ import annotations

import numpy as np

TWO_PI = 2.0 * np.pi


def solve_eccentric_anomaly(M, ecc: float, tol: float = 1e-13, max_iter: int = 100):
    """Solve ``E - e sin E = M`` for the eccentric anomaly ``E`` [rad].

    Parameters
    ----------
    M : float or array_like
        Mean anomaly [rad]. May be any real value; it is wrapped to [-pi, pi]
        for a well-conditioned Newton start, and the wrapping is undone so the
        returned ``E`` stays on the same branch as ``M``.
    ecc : float
        Eccentricity, 0 <= ecc < 1.
    tol : float
        Convergence tolerance on the update ``|dE|`` [rad].
    max_iter : int
        Maximum Newton iterations.

    Returns
    -------
    ndarray
        Eccentric anomaly ``E`` [rad], same shape as ``M``.
    """
    if not (0.0 <= ecc < 1.0):
        raise ValueError(f"Eccentricity must satisfy 0 <= e < 1; got {ecc}.")

    M = np.asarray(M, dtype=float)

    # Wrap to [-pi, pi] for a stable start, remembering the removed full turns.
    turns = np.round(M / TWO_PI)
    M_wrapped = M - turns * TWO_PI

    # Standard robust initial guess.
    E = M_wrapped + ecc * np.sin(M_wrapped)

    for _ in range(max_iter):
        f = E - ecc * np.sin(E) - M_wrapped
        fp = 1.0 - ecc * np.cos(E)
        dE = -f / fp
        E = E + dE
        if np.all(np.abs(dE) < tol):
            break

    return E + turns * TWO_PI


def true_anomaly(E, ecc: float):
    """True anomaly ``nu`` [rad] from eccentric anomaly ``E`` [rad].

    Uses the numerically stable half-angle ``atan2`` form so the result lands
    in the correct quadrant and stays continuous across perihelion.
    """
    E = np.asarray(E, dtype=float)
    beta = np.sqrt((1.0 + ecc) / (1.0 - ecc))
    return 2.0 * np.arctan2(beta * np.sin(E / 2.0), np.cos(E / 2.0))


def radius(E, a: float, ecc: float):
    """Heliocentric radius ``r = a (1 - e cos E)`` in the units of ``a``."""
    E = np.asarray(E, dtype=float)
    return a * (1.0 - ecc * np.cos(E))


def mean_anomaly(nu, ecc: float):
    """Inverse map: mean anomaly ``M`` [rad] from true anomaly ``nu`` [rad].

    Convenience for tests (round-trip M -> E -> nu -> M).
    """
    nu = np.asarray(nu, dtype=float)
    # Eccentric anomaly from true anomaly.
    E = 2.0 * np.arctan2(
        np.sqrt(1.0 - ecc) * np.sin(nu / 2.0),
        np.sqrt(1.0 + ecc) * np.cos(nu / 2.0),
    )
    return E - ecc * np.sin(E)


def kepler_state(M, a: float, ecc: float):
    """Convenience bundle: return ``(E, nu, r)`` for mean anomaly ``M``.

    ``r`` is in the units of ``a``.
    """
    E = solve_eccentric_anomaly(M, ecc)
    nu = true_anomaly(E, ecc)
    r = radius(E, a, ecc)
    return E, nu, r
