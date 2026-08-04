"""Tests for the insolation module.

The two headline tests are the physics validation targets from the context doc:
  * the analytic global-annual-mean identity  <S> = S0 / (4 a^2 sqrt(1-e^2))
  * the 65 deg N June-peak insolation dropping ~480 -> ~400 W/m^2 for
    e = 0.117, lambda_p = 283 deg.
"""

import numpy as np
import pytest

from orbital_climate.config import Config
from orbital_climate.insolation import (
    declination,
    sunset_hour_angle,
    daily_mean_insolation,
    global_annual_mean_insolation,
    analytic_global_annual_mean,
)

TWO_PI = 2.0 * np.pi


def _peak_at_latitude(lat_deg, config, n_time=4000):
    phi = np.radians(lat_deg)
    M = np.linspace(0.0, TWO_PI, n_time, endpoint=False)
    Q = daily_mean_insolation(phi, M, config)
    return float(np.max(Q))


# ---------------------------------------------------------------------------
# Analytic annual-mean identity  (end-to-end pipeline check)
# ---------------------------------------------------------------------------
@pytest.mark.parametrize("ecc", [0.0, 0.0167, 0.117, 0.3])
def test_global_annual_mean_matches_analytic(ecc):
    cfg = Config(ecc=ecc)
    numerical = global_annual_mean_insolation(cfg, n_lat=800, n_time=4000)
    analytic = analytic_global_annual_mean(cfg)
    assert numerical == pytest.approx(analytic, rel=2e-4)


def test_annual_mean_at_e0_is_S0_over_4():
    cfg = Config(ecc=0.0)
    assert analytic_global_annual_mean(cfg) == pytest.approx(cfg.S0 / 4.0, rel=1e-12)


def test_annual_mean_increases_with_eccentricity():
    """Jensen's inequality: <1/r^2> grows with e at fixed a."""
    low = analytic_global_annual_mean(Config(ecc=0.0167))
    high = analytic_global_annual_mean(Config(ecc=0.117))
    assert high > low


def test_semimajor_axis_scaling():
    """Global annual mean scales as 1/a^2."""
    base = global_annual_mean_insolation(Config(a_au=1.0, ecc=0.05), n_lat=400, n_time=2000)
    far = global_annual_mean_insolation(Config(a_au=2.0, ecc=0.05), n_lat=400, n_time=2000)
    assert far == pytest.approx(base / 4.0, rel=1e-3)


# ---------------------------------------------------------------------------
# 65 deg N June-peak Milankovitch diagnostic
# ---------------------------------------------------------------------------
def test_65N_peak_present_day():
    """Present-day 65 N summer peak ~480 W/m^2."""
    cfg = Config(ecc=0.0167, lon_perihelion_deg=283.0)
    peak = _peak_at_latitude(65.0, cfg)
    assert 465.0 < peak < 490.0


def test_65N_peak_perturbed():
    """Perturbed (e=0.117) 65 N summer peak ~400 W/m^2."""
    cfg = Config(ecc=0.117, lon_perihelion_deg=283.0)
    peak = _peak_at_latitude(65.0, cfg)
    assert 385.0 < peak < 410.0


def test_65N_peak_drop_fraction():
    """The perturbation cuts 65 N summer peak insolation by ~17%."""
    present = _peak_at_latitude(65.0, Config(ecc=0.0167, lon_perihelion_deg=283.0))
    perturbed = _peak_at_latitude(65.0, Config(ecc=0.117, lon_perihelion_deg=283.0))
    drop = (present - perturbed) / present
    assert 0.14 < drop < 0.20


# ---------------------------------------------------------------------------
# Geometry sanity checks
# ---------------------------------------------------------------------------
def test_declination_bounded_by_obliquity():
    cfg = Config(obliquity_deg=23.44)
    lam = np.linspace(0.0, TWO_PI, 200)
    delta = declination(lam, cfg)
    assert np.max(np.abs(delta)) <= cfg.obliquity_rad + 1e-12


def test_polar_night_zero_insolation():
    """Winter pole in polar night receives zero daily-mean insolation."""
    cfg = Config(ecc=0.0, lon_perihelion_deg=0.0)
    # At M=0 (perihelion, lambda=0) declination is 0; push to northern winter:
    # choose an orbital position with strongly negative declination, i.e.
    # lambda ~ 270 deg (southern summer / northern winter).
    # Find M giving lambda ~ 270 by scanning.
    M = np.linspace(0.0, TWO_PI, 2000, endpoint=False)
    phi_north_pole = np.radians(85.0)
    Q = daily_mean_insolation(phi_north_pole, M, cfg)
    # There must be a stretch of the year with exactly zero insolation (polar night).
    assert np.min(Q) == pytest.approx(0.0, abs=1e-9)


def test_polar_day_positive_insolation():
    """Summer pole in polar day receives positive insolation."""
    cfg = Config(ecc=0.0, lon_perihelion_deg=0.0)
    M = np.linspace(0.0, TWO_PI, 2000, endpoint=False)
    Q = daily_mean_insolation(np.radians(85.0), M, cfg)
    assert np.max(Q) > 0.0


def test_hour_angle_clamping():
    """Hour angle saturates at 0 (polar night) and pi (polar day)."""
    # High latitude, high declination -> polar day -> H0 = pi.
    H_day = sunset_hour_angle(np.radians(80.0), np.radians(23.0))
    assert float(H_day) == pytest.approx(np.pi, abs=1e-12)
    # High latitude, opposite-sign declination -> polar night -> H0 = 0.
    H_night = sunset_hour_angle(np.radians(80.0), np.radians(-23.0))
    assert float(H_night) == pytest.approx(0.0, abs=1e-12)


def test_equator_equinox_insolation():
    """At the equator with delta=0, daily-mean insolation is S0/pi (times 1/r^2)."""
    cfg = Config(ecc=0.0, lon_perihelion_deg=0.0)
    # M=0 -> nu=0 -> lambda=0 -> delta=0, r=a=1.
    Q = daily_mean_insolation(0.0, 0.0, cfg)
    assert float(Q) == pytest.approx(cfg.S0 / np.pi, rel=1e-9)
