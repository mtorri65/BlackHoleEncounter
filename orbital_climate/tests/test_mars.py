"""Tests for the condensing-CO2 atmosphere (Mars).

The CO2 cycle is genuinely new coupled physics -- temperature, surface albedo
and atmospheric mass all feed back on one another -- so these tests are built
around quantities that must hold exactly, rather than around agreement with
observations:

  * total CO2 is conserved to machine precision, for any state;
  * the frost-point inversion round-trips against the vapour-pressure curve;
  * latent-heat accounting balances: energy absorbed equals L x mass condensed;
  * with condensation suppressed, MarsEBM reduces identically to EBM;
  * condensation buffers the surface *at* the frost point rather than below it.

The observational comparisons (polar winter ~148 K, ~600 Pa, ~25% seasonal
swing) are included, but as the weakest tier: they are what calibration was
aimed at, so they confirm the tuning converged rather than that the model is
right.
"""

import dataclasses

import numpy as np
import pytest

from orbital_climate.config import Config
from orbital_climate.ebm import EBM, KELVIN
from orbital_climate.mars import (
    MarsEBM, co2_frost_point_K, co2_saturation_pressure_Pa,
)


def mars_config(**over) -> Config:
    """Present-day Mars, matching input_mars.yaml."""
    base = Config(
        a_au=1.5237, ecc=0.0934, obliquity_deg=25.19, lon_perihelion_deg=250.5,
        days_per_year=686.98, dt_days=686.98 / 180,
        olr_model="graybody", olr_emissivity=1.0,
        coalbedo_a0=0.75, coalbedo_a2=0.0,
        T_ice_degC=-125.0, coalbedo_ice=0.38,
        heat_capacity=1.2e6, diffusion_D=0.002, n_lat=60,
        co2_cycle=True, co2_inventory_kg_m2=200.0,
        co2_latent_heat=5.9e5, surface_gravity=3.71, co2_frost_albedo=0.62,
        spinup_max_years=60, spinup_tol_degC=1e-4,
    )
    return dataclasses.replace(base, **over) if over else base


# ---------------------------------------------------------------------------
# Tier 1/2 -- exact identities
# ---------------------------------------------------------------------------
def test_frost_point_roundtrips_against_vapour_pressure():
    """T_frost(p_sat(T)) == T. The inversion must be exact, not approximate."""
    for T in (120.0, 148.0, 160.0, 200.0):
        p = co2_saturation_pressure_Pa(T)
        assert float(co2_frost_point_K(p)) == pytest.approx(T, rel=1e-12)


def test_frost_point_matches_mars_observation():
    """At Mars's ~600 Pa the frost point is ~148 K, as observed at the caps."""
    assert float(co2_frost_point_K(600.0)) == pytest.approx(148.0, abs=1.0)
    # Falling pressure lowers the frost point -- this is what makes the
    # condensation feedback self-limiting rather than runaway.
    assert co2_frost_point_K(200.0) < co2_frost_point_K(600.0)


def test_total_co2_is_conserved_exactly():
    """Atmospheric + condensed CO2 is invariant, for any frost distribution."""
    m = MarsEBM(mars_config())
    rng = np.random.default_rng(0)
    for _ in range(5):
        frost = rng.uniform(0.0, 150.0, size=m.n)
        atmospheric = m.surface_pressure(frost) / m.gravity
        total = atmospheric + float(np.mean(frost))
        assert total == pytest.approx(m.inventory, rel=1e-14)


def test_total_co2_conserved_through_a_year():
    """Stepping must not create or destroy CO2."""
    m = MarsEBM(mars_config())
    T, frost, M, _ = m.spin_up_co2()
    totals = []
    for _ in range(60):
        T, frost = m.step_co2(T, frost, M)
        M += 2 * np.pi / 180
        totals.append(m.surface_pressure(frost) / m.gravity + float(np.mean(frost)))
    assert max(totals) - min(totals) < 1e-10


def test_latent_heat_accounting_balances():
    """Energy absorbed by condensation equals L x mass condensed.

    Cool a cell hard below the frost point and check the bookkeeping directly,
    rather than trusting that the clamp happens to be self-consistent.
    """
    cfg = mars_config()
    m = MarsEBM(cfg)
    frost = np.zeros(m.n)
    T_f = m.frost_point_degC(frost)
    # Start every cell 40 K below the frost point; all of that deficit must
    # appear as condensed mass.
    T_cold = np.full(m.n, T_f - 40.0)
    expected_mass = m.C[0] * 40.0 / cfg.co2_latent_heat      # kg/m^2 per cell

    m._m_frost = frost
    T_new = T_cold.copy()
    deficit = m.C * (T_f - T_new)
    gained = deficit / cfg.co2_latent_heat
    assert float(np.mean(gained)) == pytest.approx(expected_mass, rel=1e-12)


# ---------------------------------------------------------------------------
# Reduction to the base model
# ---------------------------------------------------------------------------
def test_reduces_to_base_ebm_when_condensation_cannot_occur():
    """With the frost point unreachable, MarsEBM must equal EBM step for step.

    This is the analogue of the check that let the Sellers OLR option ship: a
    new code path that collapses exactly onto the old one cannot have broken it.
    """
    # An enormous inventory puts the frost point far below any temperature
    # reached, so no condensation is possible.
    cfg = mars_config(co2_inventory_kg_m2=1e12)
    mars, base = MarsEBM(cfg), EBM(cfg)
    T0 = np.full(cfg.n_lat, -60.0)
    frost = np.zeros(cfg.n_lat)
    T_m, frost = mars.step_co2(T0.copy(), frost, 0.3)
    T_b = base.step(T0.copy(), 0.3)
    assert np.allclose(frost, 0.0)
    np.testing.assert_allclose(T_m, T_b, rtol=1e-12, atol=1e-12)


def test_graybody_olr_is_stefan_boltzmann():
    """olr_model='graybody' must return emissivity * sigma T^4 exactly."""
    SIG = 5.670374419e-8
    m = EBM(mars_config(olr_emissivity=0.9))
    T_c = np.array([-100.0, 0.0, 50.0])
    expected = 0.9 * SIG * (T_c + KELVIN) ** 4
    np.testing.assert_allclose(m.olr(T_c), expected, rtol=1e-12)


# ---------------------------------------------------------------------------
# Behaviour of the coupled system
# ---------------------------------------------------------------------------
def test_condensation_buffers_at_the_moving_frost_point():
    """The surface is pinned at the frost point, never driven below it.

    The frost point is *not* a constant: as CO2 condenses the pressure falls and
    the frost point falls with it, which is what makes the feedback
    self-limiting. So the comparison must be made step by step against that
    step's own pressure -- comparing against a single end-of-year value fails by
    ~1.5 K, which is the size of the seasonal frost-point excursion itself.

    Without latent heat the modelled polar winter overshoots to ~80 K against an
    observed ~148 K; this test locks that fix in.
    """
    m = MarsEBM(mars_config())
    T, frost, M, _ = m.spin_up_co2()
    T, frost, M, rec = m.run_year_co2(T, frost, M0=M)

    T_f_step = co2_frost_point_K(rec["pressure_Pa"]) - KELVIN     # [n_steps]
    coldest = rec["T"].min(axis=1)                                 # [n_steps]

    # At no time may any cell sit meaningfully below that instant's frost point.
    assert np.all(coldest >= T_f_step - 0.5)
    # And somewhere in the year the surface must actually reach it, or the test
    # would pass on a model that never condenses at all.
    assert np.min(coldest - T_f_step) < 0.5

    # The frost point really does move, or the step-by-step form is pointless.
    assert np.ptp(T_f_step) > 0.5


def test_frost_forms_at_the_winter_pole():
    m = MarsEBM(mars_config())
    T, frost, M, _ = m.spin_up_co2()
    T, frost, M, rec = m.run_year_co2(T, frost, M0=M)
    polar = rec["m_frost"][:, -1]
    assert polar.max() > 0.0          # a cap forms
    equatorial = rec["m_frost"][:, m.n // 2]
    assert equatorial.max() < polar.max()   # and it is a *polar* cap


def test_albedo_follows_frost_not_temperature():
    """Coalbedo switches on the presence of frost, so albedo can lag temperature.

    This differs from the Earth model deliberately: frost laid down in winter
    survives into spring while the surface warms.
    """
    m = MarsEBM(mars_config())
    warm = np.full(m.n, 0.0)          # far above any frost point
    m._m_frost = np.zeros(m.n)
    bare = m.coalbedo(warm)
    m._m_frost = np.full(m.n, 10.0)   # same temperature, now frosted
    frosted = m.coalbedo(warm)
    assert np.all(frosted < bare)
    assert frosted[0] == pytest.approx(1.0 - m.cfg.co2_frost_albedo)


def test_pressure_falls_as_frost_accumulates():
    m = MarsEBM(mars_config())
    p0 = m.surface_pressure(np.zeros(m.n))
    p1 = m.surface_pressure(np.full(m.n, 50.0))
    assert p1 < p0
    assert p0 == pytest.approx(m.inventory * m.gravity, rel=1e-12)


def test_collapse_is_detected_not_clamped():
    """A frozen-out atmosphere must be flagged, not reported as a temperature.

    As p -> 0 the frost-point inversion degenerates toward ~76 K rather than
    diverging, so the model would otherwise keep condensing against a threshold
    that has lost its meaning -- and report the result with a straight face.
    """
    m = MarsEBM(mars_config())
    assert not m.is_collapsed(np.zeros(m.n))              # full atmosphere
    assert not m.is_collapsed(np.full(m.n, 100.0))        # half condensed
    assert m.is_collapsed(np.full(m.n, m.inventory))      # fully condensed
    assert m.airborne_fraction(np.zeros(m.n)) == pytest.approx(1.0)
    assert m.airborne_fraction(np.full(m.n, m.inventory)) == pytest.approx(0.0)

    # The degeneracy the flag exists to catch: the frost point does not diverge
    # as the atmosphere vanishes, it converges to a finite, meaningless value.
    assert float(co2_frost_point_K(1e-6)) < 100.0


def test_run_year_reports_collapse_fraction():
    """Present-day Mars never freezes out, and the diagnostic says so."""
    m = MarsEBM(mars_config())
    T, frost, M, _ = m.spin_up_co2()
    T, frost, M, rec = m.run_year_co2(T, frost, M0=M)
    assert rec["collapsed_fraction"] == 0.0
    assert not rec["collapsed"].any()

    # A thin inventory does collapse, and is flagged rather than silently
    # producing numbers.
    thin = MarsEBM(mars_config(co2_inventory_kg_m2=8.0))
    T2, f2, M2, _ = thin.spin_up_co2()
    T2, f2, M2, rec2 = thin.run_year_co2(T2, f2, M0=M2)
    assert rec2["collapsed_fraction"] > 0.0


def test_two_surface_is_rejected():
    """Mars has no ocean; the land/ocean split must not silently do nothing."""
    with pytest.raises(NotImplementedError):
        MarsEBM(mars_config(two_surface=True))


# ---------------------------------------------------------------------------
# Tier 3 -- observational, and therefore the weakest evidence here
# ---------------------------------------------------------------------------
def test_present_day_mars_is_recognisable():
    """Calibration targets. These confirm the tuning converged, not correctness."""
    m = MarsEBM(mars_config())
    T, frost, M, _ = m.spin_up_co2()
    T, frost, M, rec = m.run_year_co2(T, frost, M0=M)
    mean_K = float(rec["T"].mean()) + KELVIN
    assert 195.0 < mean_K < 215.0                       # observed ~210 K
    assert 500.0 < rec["pressure_Pa"].mean() < 700.0    # observed ~600 Pa
    swing = 1.0 - rec["pressure_Pa"].min() / rec["pressure_Pa"].max()
    assert 0.10 < swing < 0.40                          # observed ~25%


# ---------------------------------------------------------------------------
# Pressure-dependent greenhouse (olr_model = "graygas")
# ---------------------------------------------------------------------------
def test_graygas_olr_formula():
    """OLR = sigma T^4 / (1 + 3 tau / 4), with tau scaling on surface pressure."""
    SIG = 5.670374419e-8
    cfg = mars_config(olr_model="graygas", graygas_tau_ref=0.2, graygas_p_ref_pa=600.0)
    m = MarsEBM(cfg)
    m._p_now = 600.0
    assert m._graygas_tau() == pytest.approx(0.2)
    T_c = np.array([-60.0, 0.0])
    expected = SIG * (T_c + KELVIN) ** 4 / (1.0 + 0.75 * 0.2)
    np.testing.assert_allclose(m.olr(T_c), expected, rtol=1e-12)
    # Optical depth is linear in pressure.
    m._p_now = 1200.0
    assert m._graygas_tau() == pytest.approx(0.4)


def test_graygas_reduces_to_blackbody_without_atmosphere():
    """As the atmosphere vanishes the greenhouse must vanish with it."""
    SIG = 5.670374419e-8
    m = MarsEBM(mars_config(olr_model="graygas"))
    m._p_now = 0.0
    T_c = np.array([-60.0])
    np.testing.assert_allclose(m.olr(T_c), SIG * (T_c + KELVIN) ** 4, rtol=1e-12)


def test_greenhouse_calibrated_to_present_day_mars():
    """graygas must warm present-day Mars by ~5 K relative to a bare blackbody."""
    cold = MarsEBM(mars_config(olr_model="graybody"))
    warm = MarsEBM(mars_config(olr_model="graygas"))
    means = []
    for m in (cold, warm):
        T, f, M, _ = m.spin_up_co2()
        T, f, M, rec = m.run_year_co2(T, f, M0=M)
        means.append(float(rec["T"].mean()))
    assert 3.0 < means[1] - means[0] < 8.0


def test_thicker_inventory_warms_the_planet():
    """The feedback the fixed-emissivity model could not express.

    More CO2 -> higher pressure -> stronger greenhouse -> warmer -> less
    condensation -> higher pressure still. Only present with graygas; with a
    fixed emissivity a larger inventory does almost nothing.
    """
    def mean_T(inventory, model):
        m = MarsEBM(mars_config(olr_model=model, co2_inventory_kg_m2=inventory))
        T, f, M, _ = m.spin_up_co2()
        T, f, M, rec = m.run_year_co2(T, f, M0=M)
        return float(rec["T"].mean())

    thin_gas, thick_gas = mean_T(200.0, "graygas"), mean_T(1000.0, "graygas")
    assert thick_gas - thin_gas > 10.0            # strong response

    thin_bb, thick_bb = mean_T(200.0, "graybody"), mean_T(1000.0, "graybody")
    assert abs(thick_bb - thin_bb) < 5.0          # fixed emissivity: barely responds


def test_atmospheric_collapse_is_bistable():
    """Two stable climates at identical forcing, decided by the starting state.

    The window is narrow -- s = 0.72-0.73 at this inventory, about 1.4% of
    insolation -- so the test targets it directly. It must also start from
    genuinely opposite states: a warm start with condensed frost simply
    sublimates it and lands on the thick branch, which is how an earlier version
    of this test managed to find no bistability at all.
    """
    n = 45
    cfg = mars_config(olr_model="graygas", co2_inventory_kg_m2=1000.0,
                      S0=1361.0 * 0.72, n_lat=n,
                      spinup_max_years=400, spinup_tol_degC=1e-4)
    m = MarsEBM(cfg)

    def settle(T0, m0):
        T, f, M, info = m.spin_up_co2(T0=T0, m0=m0)
        T, f, M, rec = m.run_year_co2(T, f, M0=M)
        assert info["converged"], "spin-up did not converge; the result is a transient"
        return float(rec["T"].mean()), float(rec["pressure_Pa"].mean())

    T_thick, p_thick = settle(np.full(n, -20.0), np.zeros(n))          # warm, airborne
    T_thin, p_thin = settle(np.full(n, -110.0), np.full(n, 1000.0))    # cold, condensed

    assert T_thick - T_thin > 10.0        # genuinely different equilibria
    assert p_thick > 5.0 * p_thin         # one thick, one collapsed


def test_bistable_window_is_narrow():
    """Outside a ~1.4% window the equilibrium is unique.

    Guards against the failure that produced a spurious 15%-wide loop: a chained
    insolation ramp whose spin-up stopped on a slowly re-inflating atmosphere.
    """
    n = 45
    for s in (0.65, 0.80):
        cfg = mars_config(olr_model="graygas", co2_inventory_kg_m2=1000.0,
                          S0=1361.0 * s, n_lat=n,
                          spinup_max_years=400, spinup_tol_degC=1e-4)
        m = MarsEBM(cfg)
        means = []
        for T0, m0 in ((np.full(n, -20.0), np.zeros(n)),
                       (np.full(n, -110.0), np.full(n, 1000.0))):
            T, f, M, info = m.spin_up_co2(T0=T0, m0=m0)
            T, f, M, rec = m.run_year_co2(T, f, M0=M)
            assert info["converged"]
            means.append(float(rec["T"].mean()))
        assert abs(means[0] - means[1]) < 1.0


def test_spinup_reports_convergence():
    """A spin-up that hits its year cap must say so rather than imply success."""
    m = MarsEBM(mars_config(spinup_max_years=1, spinup_consecutive=3))
    _, _, _, info = m.spin_up_co2()
    assert info["converged"] is False
