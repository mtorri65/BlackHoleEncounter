"""Tests for the Kepler solver and orbital geometry."""

import numpy as np
import pytest

from orbital_climate.kepler import (
    solve_eccentric_anomaly,
    true_anomaly,
    radius,
    mean_anomaly,
    kepler_state,
)

TWO_PI = 2.0 * np.pi


@pytest.mark.parametrize("ecc", [0.0, 0.0167, 0.117, 0.5, 0.9])
def test_keplers_equation_residual(ecc):
    """Returned E must satisfy E - e sin E = M to tight tolerance."""
    M = np.linspace(-2 * TWO_PI, 2 * TWO_PI, 501)
    E = solve_eccentric_anomaly(M, ecc)
    residual = E - ecc * np.sin(E) - M
    assert np.max(np.abs(residual)) < 1e-11


def test_circular_limit():
    """For e = 0: E = M and true anomaly nu = M (mod 2pi)."""
    M = np.linspace(0.0, TWO_PI, 50, endpoint=False)
    E = solve_eccentric_anomaly(M, 0.0)
    np.testing.assert_allclose(E, M, atol=1e-12)
    nu = true_anomaly(E, 0.0)
    # nu and M agree modulo 2pi.
    diff = (nu - M + np.pi) % TWO_PI - np.pi
    np.testing.assert_allclose(diff, 0.0, atol=1e-12)


@pytest.mark.parametrize("ecc", [0.0167, 0.117, 0.5])
def test_roundtrip_M_to_nu_to_M(ecc):
    """M -> E -> nu -> M must be the identity (modulo 2pi)."""
    M = np.linspace(0.01, TWO_PI - 0.01, 200)
    E, nu, _ = kepler_state(M, a=1.0, ecc=ecc)
    M_back = mean_anomaly(nu, ecc)
    diff = (M_back - M + np.pi) % TWO_PI - np.pi
    np.testing.assert_allclose(diff, 0.0, atol=1e-10)


@pytest.mark.parametrize("ecc", [0.0167, 0.117, 0.5])
def test_perihelion_aphelion_radii(ecc):
    """r = a(1-e) at perihelion (M=0), a(1+e) at aphelion (M=pi)."""
    a = 1.0
    _, _, r_peri = kepler_state(0.0, a, ecc)
    _, _, r_apo = kepler_state(np.pi, a, ecc)
    assert r_peri == pytest.approx(a * (1.0 - ecc), rel=1e-12)
    assert r_apo == pytest.approx(a * (1.0 + ecc), rel=1e-12)


def test_true_anomaly_quadrants():
    """nu is 0 at perihelion, pi at aphelion, and monotonic in between."""
    ecc = 0.3
    M = np.linspace(0.0, TWO_PI, 361)
    _, nu, _ = kepler_state(M, 1.0, ecc)
    nu_unwrapped = np.unwrap(nu)
    assert nu_unwrapped[0] == pytest.approx(0.0, abs=1e-9)
    # Monotonically increasing over one orbit.
    assert np.all(np.diff(nu_unwrapped) > 0)


def test_known_solution():
    """Spot-check against an independently computed value.

    For e = 0.2, M = 1.0 rad, the eccentric anomaly solving E - 0.2 sin E = 1
    is approximately 1.1853242 rad.
    """
    E = solve_eccentric_anomaly(1.0, 0.2)
    assert float(E) == pytest.approx(1.1853242, abs=1e-6)


def test_vectorization_matches_scalar():
    """Vectorised call equals element-wise scalar calls."""
    ecc = 0.117
    M = np.array([0.3, 1.1, 2.7, 4.9])
    E_vec = solve_eccentric_anomaly(M, ecc)
    E_scalar = np.array([float(solve_eccentric_anomaly(m, ecc)) for m in M])
    np.testing.assert_allclose(E_vec, E_scalar, atol=1e-13)


def test_invalid_eccentricity_raises():
    with pytest.raises(ValueError):
        solve_eccentric_anomaly(1.0, 1.0)
    with pytest.raises(ValueError):
        solve_eccentric_anomaly(1.0, -0.1)
