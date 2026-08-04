"""Orbital Perturbation Climate Model.

A 1-D (latitude-resolved) seasonal energy-balance model with orbital forcing,
in the North / Budyko / Sellers lineage.

This first version delivers the physics foundation only:
    * kepler     -- solve Kepler's equation; true anomaly and heliocentric radius
    * insolation -- daily-mean top-of-atmosphere insolation vs. latitude and season
    * config     -- YAML-backed configuration with sensible defaults

The energy-balance core and sweep harness build on top of these later.
"""

from .config import Config, load_config
from .kepler import solve_eccentric_anomaly, true_anomaly, radius, kepler_state
from .insolation import (
    solar_longitude,
    declination,
    sunset_hour_angle,
    daily_mean_insolation,
    annual_mean_insolation,
)
from .ebm import EBM
from .experiment import (
    run_equilibrium,
    run_perturbation,
    perturbed_config_from_scenario,
    save_perturbation_outputs,
)
from .sweep import run_sweep

__all__ = [
    "Config",
    "load_config",
    "solve_eccentric_anomaly",
    "true_anomaly",
    "radius",
    "kepler_state",
    "solar_longitude",
    "declination",
    "sunset_hour_angle",
    "daily_mean_insolation",
    "annual_mean_insolation",
    "EBM",
    "run_equilibrium",
    "run_perturbation",
    "perturbed_config_from_scenario",
    "save_perturbation_outputs",
    "run_sweep",
]
