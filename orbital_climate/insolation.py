"""Top-of-atmosphere insolation.

Given the orbital state (from :mod:`orbital_climate.kepler`) this module
computes the **daily-mean** insolation on a horizontal surface as a function of
latitude and orbital position, following the standard astronomical formulation.

Pipeline (all angles in radians internally):

    solar longitude   lambda = nu + lambda_p
    declination       delta  = asin(sin(eps) * sin(lambda))
    sunset hour angle H0     = acos(-tan(phi) * tan(delta))     [clamped]
    daily-mean flux   Q      = (S0 / pi) / r^2
                               * (H0 sin(phi) sin(delta) + cos(phi) cos(delta) sin(H0))

Here ``r`` is the heliocentric distance in AU and ``S0`` is the solar constant
at 1 AU, so ``S0 / r^2`` is the instantaneous flux at distance ``r``.

Time-averaging over the orbit is done by sampling the **mean anomaly** uniformly
(``M`` advances linearly in time), which is what makes the annual-mean identity

    <S> = S0 / (4 a^2 sqrt(1 - e^2))

exact and usable as an end-to-end unit test.
"""

from __future__ import annotations

import numpy as np

from .config import Config
from .kepler import kepler_state

TWO_PI = 2.0 * np.pi


def solar_longitude(nu, config: Config):
    """Solar longitude ``lambda = nu + lambda_p`` [rad]."""
    return np.asarray(nu, dtype=float) + config.lon_perihelion_rad


def declination(lam, config: Config):
    """Solar declination ``delta = asin(sin(eps) sin(lambda))`` [rad]."""
    lam = np.asarray(lam, dtype=float)
    return np.arcsin(np.sin(config.obliquity_rad) * np.sin(lam))


def sunset_hour_angle(phi, delta):
    """Sunset hour angle ``H0`` [rad], with polar day/night clamping.

    ``cos(H0) = -tan(phi) tan(delta)``. When the argument leaves [-1, 1] the
    Sun either never sets (polar day, ``H0 = pi``) or never rises
    (polar night, ``H0 = 0``); the clamp handles both branches.

    Broadcasts over ``phi`` and ``delta``.
    """
    phi = np.asarray(phi, dtype=float)
    delta = np.asarray(delta, dtype=float)
    cos_H = -np.tan(phi) * np.tan(delta)
    cos_H = np.clip(cos_H, -1.0, 1.0)
    return np.arccos(cos_H)


def daily_mean_insolation(phi, M, config: Config):
    """Daily-mean insolation [W m^-2] at latitude(s) ``phi`` and mean anomaly ``M``.

    Parameters
    ----------
    phi : float or array_like
        Latitude(s) [rad], in [-pi/2, pi/2].
    M : float or array_like
        Mean anomaly [rad]; ``M`` advances linearly in time (``M = 2*pi*t/P``).
    config : Config
        Orbital / stellar parameters.

    Returns
    -------
    ndarray
        Daily-mean insolation, broadcast over ``phi`` x ``M``. If both inputs
        are 1-D the result is 2-D with shape ``(len(M), len(phi))``.
    """
    phi = np.asarray(phi, dtype=float)
    M = np.asarray(M, dtype=float)

    # Orbital position: radius in AU and declination for each M.
    _, nu, r = kepler_state(M, config.a_au, config.ecc)
    lam = solar_longitude(nu, config)
    delta = declination(lam, config)

    # Arrange so that M varies along axis 0 and phi along axis 1 when both are
    # vectors; scalars broadcast naturally.
    if phi.ndim and M.ndim:
        phi_b = phi[np.newaxis, :]
        delta_b = delta[:, np.newaxis]
        r_b = r[:, np.newaxis]
    else:
        phi_b, delta_b, r_b = phi, delta, r

    H0 = sunset_hour_angle(phi_b, delta_b)
    flux = config.S0 / (r_b * r_b)  # instantaneous flux at distance r [W m^-2]

    return (flux / np.pi) * (
        H0 * np.sin(phi_b) * np.sin(delta_b)
        + np.cos(phi_b) * np.cos(delta_b) * np.sin(H0)
    )


def annual_mean_insolation(phi, config: Config, n_time: int = 2000):
    """Time-mean daily insolation [W m^-2] at latitude(s) ``phi``.

    The average is taken over one orbit by sampling the mean anomaly uniformly
    (equivalent to a uniform-in-time average, since ``M`` is linear in time).
    """
    M = np.linspace(0.0, TWO_PI, n_time, endpoint=False)
    Q = daily_mean_insolation(phi, M, config)
    # Average over the time axis (axis 0 when phi is a vector; else the whole array).
    axis = 0 if np.ndim(phi) else None
    return np.mean(Q, axis=axis)


def global_annual_mean_insolation(config: Config, n_lat: int = 400, n_time: int = 2000):
    """Global-and-annual mean insolation [W m^-2].

    Area-weighted (``cos phi``) latitude average of the annual-mean insolation.
    Should equal the analytic result ``S0 / (4 a^2 sqrt(1 - e^2))``; this is the
    primary end-to-end validation of the Kepler + insolation pipeline.
    """
    # Gauss-like uniform grid in sin(phi) gives equal-area weighting for free:
    # d(area) ∝ cos(phi) d(phi) = d(sin phi).
    sin_phi = np.linspace(-1.0, 1.0, n_lat)
    phi = np.arcsin(np.clip(sin_phi, -1.0, 1.0))
    annual = annual_mean_insolation(phi, config, n_time=n_time)
    # Equal-area mean over latitude == simple mean in sin(phi) via trapezoid.
    return np.trapezoid(annual, sin_phi) / (sin_phi[-1] - sin_phi[0])


def analytic_global_annual_mean(config: Config) -> float:
    """Closed-form global annual mean ``S0 / (4 a^2 sqrt(1 - e^2))`` [W m^-2]."""
    return config.S0 / (4.0 * config.a_au ** 2 * np.sqrt(1.0 - config.ecc ** 2))
