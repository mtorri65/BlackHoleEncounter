"""Configuration for the orbital climate model.

Physical constants and orbital/model parameters live here. Values may be
overridden via a YAML file (see ``input_climate.yaml``), mirroring the
``input.yaml`` convention used elsewhere in this repository.

Angles in the public config are expressed in **degrees** (as a human writes
them in YAML); the physics modules convert to radians internally.
"""

from __future__ import annotations

import math
from dataclasses import dataclass, asdict, fields
from pathlib import Path

import yaml


# ---------------------------------------------------------------------------
# Physical constants
# ---------------------------------------------------------------------------
S0_DEFAULT = 1361.0             # Solar constant at 1 AU [W m^-2]
AU_M = 1.495978707e11           # 1 astronomical unit [m]
DAYS_PER_YEAR = 365.256363      # Sidereal year [days]


@dataclass
class Config:
    """Orbital + model configuration.

    Orbital elements describe Earth's heliocentric orbit. The default values
    reproduce present-day Earth; the perturbed scenario from the context doc
    corresponds to ``e = 0.117`` with ``a`` essentially unchanged.
    """

    # --- Stellar / orbital forcing ---
    S0: float = S0_DEFAULT          # solar constant at 1 AU [W m^-2]
    a_au: float = 1.0               # semi-major axis [AU]
    ecc: float = 0.0167             # orbital eccentricity [-]
    obliquity_deg: float = 23.44    # axial tilt epsilon [deg]
    lon_perihelion_deg: float = 283.0  # longitude of perihelion lambda_p [deg]
                                    #   ~283 deg places perihelion near Jan 3

    # --- Numerics for the insolation daily-mean sampling ---
    days_per_year: float = DAYS_PER_YEAR

    # --- EBM core (North / Budyko / Sellers) ---
    # Outgoing longwave radiation linearisation OLR = A + B*T, with T in degC.
    olr_A: float = 203.3            # OLR intercept [W m^-2]
    olr_B: float = 2.09             # OLR slope [W m^-2 degC^-1]

    # --- OLR parameterisation ---
    # "linear"  : Budyko (1969) form OLR = A + B*T. Calibrated near 288 K; it
    #             reaches OLR = 0 at 175.9 K and behaves qualitatively wrongly
    #             below ~230 K, where a drying atmosphere should push OLR *up*
    #             toward blackbody, not down toward zero.
    # "sellers" : Sellers (1969) OLR = sigma T^4 [1 - m tanh(19 T^6 1e-16)],
    #             T in Kelvin. Agrees with the linear form at 288 K to 0.3% but
    #             correctly approaches blackbody emission as the planet freezes.
    olr_model: str = "linear"
    # Atmospheric attenuation parameter. 0.51 is calibrated so the Sellers form
    # reproduces the same 288 K present-day equilibrium as the tuned linear one
    # (the raw Sellers value 0.5 gives 285.7 K once the ice-albedo feedback and
    # the full latitude profile are accounted for).
    sellers_m: float = 0.51

    # Emissivity for olr_model = "graybody" (OLR = emissivity * sigma T^4).
    # 1.0 is correct for a body with negligible greenhouse effect; Mars needs
    # 1.00, Earth would need ~0.60.
    olr_emissivity: float = 1.0

    # --- olr_model = "graygas": pressure-dependent greenhouse ---
    # Grey-gas atmosphere in radiative equilibrium:
    #     OLR = sigma T^4 / (1 + 3 tau / 4),    tau = tau_ref * (p / p_ref)
    # Unlike "graybody" this responds to atmospheric mass, so a planet whose
    # caps sublimate develops a real greenhouse. tau_ref = 0.1316 reproduces
    # Mars's observed ~5 K greenhouse at 600 Pa.
    #
    # Only meaningful together with a CO2 inventory large enough for the
    # pressure to move: at the default 200 kg/m^2 the ceiling is 742 Pa (1.24x
    # today) and the greenhouse varies by under 1 K, so the two settings are
    # one physical assumption rather than two.
    graygas_tau_ref: float = 0.1316
    graygas_p_ref_pa: float = 600.0

    # Simpson-Nakajima ceiling: the maximum OLR a moist atmosphere can sustain
    # (Nakajima, Hayashi & Abe 1992). If the absorbed flux exceeds it there is
    # *no* equilibrium -- the planet enters a runaway greenhouse. This is used as
    # a diagnostic flag, not as a cap, because capping would invent a stable
    # state that physically does not exist.
    olr_runaway_limit: float = 300.0    # [W m^-2]
    # Meridional heat transport coefficient in D * d/dx[(1-x^2) dT/dx].
    diffusion_D: float = 0.58       # [W m^-2 degC^-1]
    # Coalbedo (absorbed fraction) = a0 + a2*P2(x) where ice-free; a_ice where
    # T < T_ice. Tuned so the annual-mean global temperature is ~15 degC / 288 K.
    coalbedo_a0: float = 0.676      # mean ice-free coalbedo [-] (tuned to 288 K)
    coalbedo_a2: float = -0.200     # P2 latitude structure of coalbedo [-]
    coalbedo_ice: float = 0.38      # coalbedo over ice (albedo 0.62) [-]
    T_ice_degC: float = -10.0       # ice-formation temperature threshold [degC]
    # Mixed-layer heat capacity. tau = C/B sets the thermal inertia; the doc's
    # ~1.6 yr ocean time constant corresponds to C ~ 1.05e8 J m^-2 degC^-1.
    heat_capacity: float = 1.05e8   # C [J m^-2 degC^-1] (single-surface mode)

    # --- Two-surface land/ocean mode (North & Coakley 1979) ---
    # When enabled, each latitude carries *separate* land and ocean temperatures
    # coupled by a zonal exchange term. An area-weighted blend of heat capacities
    # would be dominated by the ocean value and barely change the seasonal cycle;
    # separate temperatures preserve the ~175x contrast in thermal inertia, which
    # is what produces continental-scale seasonal extremes.
    two_surface: bool = False       # opt-in; False reproduces the single-surface model
    heat_capacity_land: float = 1.2e6    # C_land  [J m^-2 degC^-1] (~1 m soil, tau ~ 7 d)
    heat_capacity_ocean: float = 2.1e8   # C_ocean [J m^-2 degC^-1] (50 m mixed layer)
    # nu [W m^-2 degC^-1]: zonal land<->ocean exchange. Calibrated so the 65 N
    # seasonal range is ~40 K over land and ~9 K over ocean, matching observed
    # continental vs. maritime seasonality.
    land_ocean_coupling: float = 3.5
    land_fraction_override: float | None = None   # force a uniform land fraction if set

    # --- Condensable-atmosphere (CO2) cycle: see orbital_climate/mars.py ---
    # On Mars the atmosphere itself condenses onto the winter pole, removing
    # mass from the atmosphere and buffering the surface at the frost point via
    # latent heat. Without it, modelled polar winters overshoot to ~80 K against
    # an observed ~148 K.
    co2_cycle: bool = False
    co2_inventory_kg_m2: float = 200.0    # total CO2 (atmosphere + caps) [kg/m^2]
    co2_latent_heat: float = 5.9e5        # sublimation enthalpy [J/kg]
    surface_gravity: float = 3.71         # [m/s^2]; Mars
    co2_frost_albedo: float = 0.62        # fresh CO2 frost; coalbedo = 1 - this
    # Fraction of the inventory that must remain airborne for the model to be
    # meaningful. Below it the atmosphere has essentially all condensed: the
    # frost point degenerates (it tends to ~76 K as p -> 0), so temperatures and
    # pressures reported from that state describe nothing. Detected and flagged
    # rather than clamped, for the same reason as the Simpson-Nakajima runaway.
    co2_collapse_threshold: float = 0.01
    # Spin-up convergence for the frost field, as a *fraction of the inventory*
    # per year. An absolute threshold is meaningless when the inventory ranges
    # over orders of magnitude: 1e-3 kg/m^2/yr is tight against 200 kg/m^2 and
    # useless against 1000, where it will declare convergence on an atmosphere
    # that is still visibly re-inflating.
    # 1e-5 of the inventory per year. Tighter than this sits below the
    # year-to-year numerical noise on the cap mass (~1e-3 kg/m^2 at a 1000
    # kg/m^2 inventory), making convergence unreachable on the collapsed branch
    # -- which reads as "did not converge" when the state is in fact stable to
    # five significant figures over 2000 years.
    co2_spinup_tol_frac: float = 1e-5
    # Consecutive years that must satisfy both criteria. A slowly creeping state
    # can meet a per-year tolerance indefinitely; requiring several in a row
    # catches that.
    spinup_consecutive: int = 3
    n_lat: int = 180                # number of latitude cells (cell-centred in x=sin(phi))

    # Latitude at which the seasonal peak/trough temperature and peak insolation
    # are tracked (the Milankovitch diagnostic; 65 N by convention). May be
    # negative for a Southern-Hemisphere diagnostic.
    diag_lat_deg: float = 65.0

    # --- Time stepping ---
    dt_days: float = 2.0            # timestep [days] (semi-implicit -> stable at large dt)
    spinup_max_years: int = 200     # cap on spin-up length [years]
    spinup_tol_degC: float = 1e-4   # year-over-year global-mean change to call equilibrium

    # ------------------------------------------------------------------
    # Derived / convenience properties
    # ------------------------------------------------------------------
    @property
    def obliquity_rad(self) -> float:
        return math.radians(self.obliquity_deg)

    @property
    def lon_perihelion_rad(self) -> float:
        return math.radians(self.lon_perihelion_deg)

    @property
    def a_m(self) -> float:
        """Semi-major axis in metres."""
        return self.a_au * AU_M

    def to_dict(self) -> dict:
        return asdict(self)


def load_config(path: str | Path | None = None, **overrides) -> Config:
    """Load a :class:`Config`, optionally from a YAML file, with keyword overrides.

    Unknown keys in the YAML file are ignored with no error so that a shared
    config file may carry parameters for model stages not yet implemented
    (EBM core, sweep harness).
    """
    data: dict = {}
    if path is not None:
        with open(path, "r", encoding="utf-8") as fh:
            loaded = yaml.safe_load(fh) or {}
        if not isinstance(loaded, dict):
            raise ValueError(f"Config file {path!r} must contain a mapping at top level.")
        data.update(loaded)
    data.update(overrides)

    valid = {f.name for f in fields(Config)}
    filtered = {k: v for k, v in data.items() if k in valid}
    return Config(**filtered)
