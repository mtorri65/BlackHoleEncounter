"""Tests for the energy-balance model core, experiment, and sweep harness.

Physics-validation tests are built to be independent of seasonal detail where
possible, so they pin down the discretisation rather than a tuned number:

  * the diffusion operator conserves energy and annihilates a constant field;
  * with D=0 and no ice each latitude reaches the analytic local radiative
    equilibrium  <T> = (<Q> a - A) / B  (exercises the whole stepper);
  * at equilibrium the global-and-annual energy budget closes;
  * the tuned default reproduces the ~288 K global mean;
  * the semi-implicit scheme is stable at large timesteps.
"""

import numpy as np
import pytest

from orbital_climate.config import Config
from orbital_climate.ebm import EBM
from orbital_climate.experiment import run_equilibrium, run_perturbation
from orbital_climate.insolation import annual_mean_insolation


# ---------------------------------------------------------------------------
# Diffusion operator
# ---------------------------------------------------------------------------
def test_diffusion_annihilates_constant():
    """A spatially constant temperature has zero diffusive tendency."""
    m = EBM(Config(n_lat=120))
    ones = np.ones(m.n)
    np.testing.assert_allclose(m.L @ ones, 0.0, atol=1e-10)


def test_diffusion_conserves_energy():
    """Equal-area sum of the diffusive tendency is zero for any field."""
    m = EBM(Config(n_lat=120))
    rng = np.random.default_rng(0)
    for _ in range(5):
        T = rng.normal(size=m.n) * 20.0
        # Equal-area grid => area-weighted sum is the plain sum.
        assert abs(np.sum(m.L @ T)) < 1e-8


def test_diffusion_smooths_variance():
    """Pure diffusion (no forcing) reduces the spatial variance each step."""
    m = EBM(Config(n_lat=120, diffusion_D=0.6))
    T = np.cos(3 * m.phi) * 10.0
    v0 = np.var(T)
    # One explicit diffusion sub-step (small) must not increase variance.
    T1 = T + 1e-3 * (m.L @ T)
    assert np.var(T1) <= v0 + 1e-12


# ---------------------------------------------------------------------------
# Analytic local radiative equilibrium (D = 0, no ice)
# ---------------------------------------------------------------------------
def test_local_radiative_equilibrium_no_diffusion():
    """With D=0 and no ice, <T>_i = (<Q>_i a0 - A) / B at every latitude."""
    cfg = Config(
        n_lat=120,
        diffusion_D=0.0,
        coalbedo_a2=0.0,          # uniform coalbedo a0
        T_ice_degC=-1.0e9,        # never freeze -> no ice feedback
        dt_days=1.0,
        spinup_max_years=120,
        spinup_tol_degC=1e-6,
    )
    m = EBM(cfg)
    T, M, info = m.spin_up()
    _, _, rec = m.run_year(T, M0=M)
    T_time_mean = rec["T"].mean(axis=0)

    Q_annual = annual_mean_insolation(m.phi, cfg)
    T_analytic = (Q_annual * cfg.coalbedo_a0 - cfg.olr_A) / cfg.olr_B
    np.testing.assert_allclose(T_time_mean, T_analytic, atol=0.05)


# ---------------------------------------------------------------------------
# Global energy budget closure
# ---------------------------------------------------------------------------
def test_global_energy_budget_closes():
    """At equilibrium, <absorbed> == <OLR> in the global-and-annual mean.

    Diffusion conserves energy globally, so the only balance is radiative:
    mean(Q * coalbedo) == A + B * mean(T).
    """
    cfg = Config(n_lat=120)
    m = EBM(cfg)
    T, M, info = m.spin_up()
    _, _, rec = m.run_year(T, M0=M)

    Tt = rec["T"]                                  # [n_steps, n_lat]
    absorbed = 0.0
    for k in range(Tt.shape[0]):
        Q = m.insolation(rec["M"][k])
        absorbed += np.mean(Q * m.coalbedo(Tt[k]))
    absorbed /= Tt.shape[0]

    olr = cfg.olr_A + cfg.olr_B * np.mean(Tt)
    assert absorbed == pytest.approx(olr, abs=0.05)


# ---------------------------------------------------------------------------
# Benchmark: present-day global mean ~288 K
# ---------------------------------------------------------------------------
def test_present_day_global_mean_288K():
    eq = run_equilibrium(Config())
    T_K = eq.global_mean + 273.15
    assert 287.0 < T_K < 289.0


def test_equilibrium_profile_is_warm_at_equator_cold_at_poles():
    eq = run_equilibrium(Config())
    i_eq = int(np.argmin(np.abs(eq.lat_deg - 0.0)))
    assert eq.T[i_eq] > eq.T[0] and eq.T[i_eq] > eq.T[-1]


# ---------------------------------------------------------------------------
# Eccentricity injection warms the annual mean (Jensen's inequality via EBM)
# ---------------------------------------------------------------------------
def test_higher_eccentricity_warms_equilibrium():
    low = run_equilibrium(Config(ecc=0.0167)).global_mean
    high = run_equilibrium(Config(ecc=0.117)).global_mean
    assert high > low


# ---------------------------------------------------------------------------
# Ice-albedo feedback: dimmer sun pushes the ice edge equatorward
# ---------------------------------------------------------------------------
def test_ice_edge_moves_equatorward_when_cooled():
    warm = run_equilibrium(Config(S0=1361.0))
    cold = run_equilibrium(Config(S0=1290.0))
    assert cold.iceline_lat_nh is not None and warm.iceline_lat_nh is not None
    assert cold.iceline_lat_nh < warm.iceline_lat_nh


# ---------------------------------------------------------------------------
# Semi-implicit stability at large timestep
# ---------------------------------------------------------------------------
def test_semi_implicit_stable_large_dt():
    small = run_equilibrium(Config(dt_days=1.0)).global_mean
    big = run_equilibrium(Config(dt_days=30.0)).global_mean
    assert np.isfinite(big)
    assert abs(big - small) < 1.0   # same equilibrium within discretisation error


# ---------------------------------------------------------------------------
# Perturbation experiment plumbing
# ---------------------------------------------------------------------------
def test_perturbation_runs_and_records():
    base = Config()
    pert = Config(ecc=0.117)
    res = run_perturbation(base, pert, n_years=5)
    assert res.years.tolist() == [1, 2, 3, 4, 5]
    for key in ("global_mean", "nh_mean", "sh_mean", "Tdiag_summer_max"):
        assert res.diagnostics[key].shape == (5,)
        assert np.all(np.isfinite(res.diagnostics[key]))


# ---------------------------------------------------------------------------
# Sweep harness
# ---------------------------------------------------------------------------
def test_sweep_writes_table(tmp_path):
    from orbital_climate.sweep import run_sweep

    run_dir = run_sweep(
        Config(n_lat=90, dt_days=5.0),
        {"ecc": [0.0167, 0.117]},
        out_dir=tmp_path,
        workers=1,
        stamp="testrun",
    )
    files = list(run_dir.glob("sweep_results.*"))
    assert len(files) == 1

    import pandas as pd
    path = files[0]
    df = pd.read_parquet(path) if path.suffix == ".parquet" else pd.read_csv(path)
    assert len(df) == 2
    assert "ecc" in df.columns and "global_mean_K" in df.columns
    # Higher eccentricity -> warmer global mean.
    df = df.sort_values("ecc").reset_index(drop=True)
    assert df.loc[1, "global_mean_K"] > df.loc[0, "global_mean_K"]


def test_diag_lat_is_configurable():
    """The tracked diagnostic latitude follows Config.diag_lat_deg."""
    eq_n = run_equilibrium(Config(diag_lat_deg=65.0))
    eq_s = run_equilibrium(Config(diag_lat_deg=-65.0))
    assert eq_n.diag_lat_deg == pytest.approx(65.0, abs=1.0)
    assert eq_s.diag_lat_deg == pytest.approx(-65.0, abs=1.0)

    # A perturbation run tracks the same configured latitude.
    res = run_perturbation(Config(), Config(ecc=0.117, diag_lat_deg=45.0), n_years=3)
    assert res.diag_lat_deg == pytest.approx(45.0, abs=1.0)
    assert "Tdiag_summer_max" in res.diagnostics


def test_equilibrium_reports_seasonal_extremes():
    """Seasonal extremes must bracket the annual mean and be self-consistent."""
    eq = run_equilibrium(Config())
    # Summer peak above winter trough, and a strictly positive seasonal range.
    assert eq.diag_summer_max > eq.diag_winter_min
    assert eq.diag_seasonal_range == pytest.approx(
        eq.diag_summer_max - eq.diag_winter_min, abs=1e-12)
    # Hemispheric means bracket (or equal) the global mean.
    assert min(eq.nh_mean, eq.sh_mean) <= eq.global_mean <= max(eq.nh_mean, eq.sh_mean)


def test_obliquity_drives_seasonal_range():
    """A larger axial tilt produces a larger seasonal swing at mid-latitudes."""
    small = run_equilibrium(Config(obliquity_deg=10.0)).diag_seasonal_range
    large = run_equilibrium(Config(obliquity_deg=35.0)).diag_seasonal_range
    assert large > small


# ---------------------------------------------------------------------------
# Nonlinear OLR (Sellers 1969) and the runaway-greenhouse diagnostic
# ---------------------------------------------------------------------------
SIGMA_SB = 5.670374419e-8


def test_olr_forms_agree_at_present_day():
    """Sellers and the linear fit are both calibrated near 288 K."""
    m = EBM(Config(olr_model="sellers"))
    lin = EBM(Config(olr_model="linear"))
    T = np.array([15.0])                      # 288.15 K
    assert m.olr(T).item() == pytest.approx(lin.olr(T).item(), rel=0.02)


def test_sellers_approaches_blackbody_when_frozen():
    """The key fix: a cold, dry atmosphere must emit close to sigma T^4.

    The linear form instead collapses toward zero (it reaches OLR = 0 at
    175.9 K), which is qualitatively wrong.
    """
    m = EBM(Config(olr_model="sellers"))
    lin = EBM(Config(olr_model="linear"))
    for T_K in (180.0, 200.0, 220.0):
        T_c = np.array([T_K - 273.15])
        bb = SIGMA_SB * T_K ** 4
        sellers = m.olr(T_c).item()
        linear = lin.olr(T_c).item()
        assert sellers > linear                    # Sellers radiates more when cold
        assert 0.6 * bb < sellers <= bb            # and stays a sane fraction of blackbody
    # At the coldest point the linear form is nearly dead while Sellers is not.
    assert lin.olr(np.array([180.0 - 273.15])).item() < 20.0
    assert m.olr(np.array([180.0 - 273.15])).item() > 40.0


def test_linear_olr_never_exceeds_its_zero_point():
    """Documents the linear form's hard floor: OLR = 0 at 175.9 K."""
    lin = EBM(Config(olr_model="linear"))
    T_zero = -Config().olr_A / Config().olr_B          # degC
    assert lin.olr(np.array([T_zero])).item() == pytest.approx(0.0, abs=1e-9)
    assert lin.olr(np.array([T_zero - 10.0])).item() < 0.0     # unphysical below


def test_sellers_reproduces_288K():
    eq = run_equilibrium(Config(olr_model="sellers", n_lat=90))
    assert 287.0 < eq.global_mean + 273.15 < 289.5


def test_unknown_olr_model_raises():
    m = EBM(Config(olr_model="not_a_model"))
    with pytest.raises(ValueError):
        m.olr(np.array([0.0]))


def test_runaway_flag_fires_only_when_too_hot():
    """Earth is stable; moving it well inside the runaway limit is not."""
    earth = run_equilibrium(Config(olr_model="sellers", n_lat=90))
    assert not earth.runaway
    assert earth.absorbed_mean < Config().olr_runaway_limit

    hot = run_equilibrium(Config(olr_model="sellers", a_au=0.80, n_lat=90))
    assert hot.runaway
    assert hot.absorbed_mean > Config().olr_runaway_limit


def test_absorbed_matches_olr_at_equilibrium():
    """Energy balance: global-mean absorbed flux equals global-mean OLR."""
    cfg = Config(olr_model="sellers", n_lat=90)
    eq = run_equilibrium(cfg)
    m = EBM(cfg)
    olr_mean = float(np.mean(m.olr(eq.T)))
    assert eq.absorbed_mean == pytest.approx(olr_mean, abs=2.0)


# ---------------------------------------------------------------------------
# Two-surface land/ocean mode
# ---------------------------------------------------------------------------
def test_land_fraction_matches_earth():
    """The zonal land-fraction profile integrates to Earth's ~29% land cover."""
    from orbital_climate.ebm import earth_land_fraction
    m = EBM(Config(two_surface=True, n_lat=180))
    # Equal-area grid, so the plain cell mean is the area-weighted global mean.
    assert m.land_frac.mean() == pytest.approx(0.29, abs=0.02)
    # Antarctica is land; the southern mid-latitudes are nearly all ocean.
    assert earth_land_fraction(-85.0) > 0.9
    assert earth_land_fraction(-55.0) < 0.1


def test_two_surface_state_layout():
    """Two-surface state is [land, ocean]; blending returns one value per latitude."""
    m = EBM(Config(two_surface=True, n_lat=90))
    assert m.n_surf == 2
    T = np.concatenate([np.full(90, 10.0), np.full(90, 20.0)])
    T_l, T_o = m.split(T)
    assert T_l.shape == (90,) and T_o.shape == (90,)
    blended = m.blend(T)
    assert blended.shape == (90,)
    # Blend is the land-fraction-weighted average of the two surfaces.
    np.testing.assert_allclose(blended, m.land_frac * 10.0 + (1 - m.land_frac) * 20.0)


def test_blend_is_idempotent():
    """Blending an already-blended state leaves it unchanged."""
    m = EBM(Config(two_surface=True, n_lat=90))
    T = np.concatenate([np.full(90, 10.0), np.full(90, 20.0)])
    once = m.blend(T)
    np.testing.assert_allclose(m.blend(once), once)


def test_land_swings_more_than_ocean():
    """The whole point: land's small heat capacity gives a far larger seasonal swing."""
    m = EBM(Config(two_surface=True, n_lat=90))
    T, M, _ = m.spin_up()
    _, _, rec = m.run_year(T, M0=M)
    T_land, T_ocean = m.split(rec["T"])
    i65 = int(np.argmin(np.abs(m.lat_deg - 65.0)))
    land_range = np.ptp(T_land[:, i65])
    ocean_range = np.ptp(T_ocean[:, i65])
    assert land_range > 3.0 * ocean_range
    # Calibrated against observed 65 N seasonality (~40 K land, ~9 K ocean).
    assert 25.0 < land_range < 55.0
    assert 4.0 < ocean_range < 15.0


def test_two_surface_preserves_global_mean():
    """Splitting the surface changes seasonality, not the radiative equilibrium."""
    single = run_equilibrium(Config(two_surface=False, n_lat=90)).global_mean
    double = run_equilibrium(Config(two_surface=True, n_lat=90)).global_mean
    assert abs(double - single) < 1.5


def test_two_surface_reports_land_ocean_extremes():
    eq = run_equilibrium(Config(two_surface=True, n_lat=90))
    assert eq.diag_land_seasonal_range > eq.diag_ocean_seasonal_range
    # The blended range must fall between the two surfaces' ranges.
    assert (eq.diag_ocean_seasonal_range <= eq.diag_seasonal_range
            <= eq.diag_land_seasonal_range)
    # Single-surface mode leaves these unset.
    assert run_equilibrium(Config(n_lat=90)).diag_land_seasonal_range is None


def test_stronger_coupling_damps_land_swing():
    """Raising the land<->ocean exchange pulls land toward the ocean's mild cycle."""
    weak = run_equilibrium(Config(two_surface=True, n_lat=90,
                                  land_ocean_coupling=1.0)).diag_land_seasonal_range
    strong = run_equilibrium(Config(two_surface=True, n_lat=90,
                                    land_ocean_coupling=8.0)).diag_land_seasonal_range
    assert strong < weak


def test_sweep_rejects_unknown_field(tmp_path):
    from orbital_climate.sweep import run_sweep
    with pytest.raises(ValueError):
        run_sweep(Config(), {"not_a_field": [1, 2]}, out_dir=tmp_path, workers=1)
