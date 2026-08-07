"""Condensable-atmosphere (CO2) cycle for a Mars-like planet.

Earth's atmosphere does not condense. Mars's does: roughly a quarter of it
freezes onto the winter pole each year and sublimates back in spring. That
single fact couples three things an Earth EBM keeps separate --- surface
temperature, surface albedo, and atmospheric mass --- and it is the reason a
recalibrated-but-otherwise-unchanged Earth model gets Martian winters wrong by
about 70 K.

The mechanism
-------------
Where the surface would cool below the CO2 frost point, CO2 condenses instead.
The latent heat released holds the surface *at* the frost point no matter how
much more energy it radiates away. So the polar winter temperature is not set
by the energy balance at all --- it is pinned by thermodynamics.

Three couplings follow, and all three are represented here:

* **Latent buffering.** Condensation absorbs the energy deficit, clamping T to
  the frost point rather than letting it run away.
* **Mass exchange.** Condensed CO2 leaves the atmosphere, lowering the surface
  pressure. The frost point falls with it (Clausius-Clapeyron), so condensation
  is self-limiting --- a *negative* feedback, unlike ice-albedo.
* **Albedo.** Fresh CO2 frost is bright (albedo ~0.62). Unlike the Earth model,
  albedo here follows the *presence of frost*, not a temperature threshold ---
  frost laid down in winter survives into spring while the surface is warming.

Frost point
-----------
CO2 saturation vapour pressure, p(T) = 1.382e12 exp(-3182.48 / T) Pa, inverted:

    T_frost(p) = 3182.48 / (ln(1.382e12) - ln p)

which gives 148 K at Mars's ~600 Pa, matching the observed winter cap
temperature, and falls to 141 K at 200 Pa.

State
-----
``MarsEBM`` carries a second state variable alongside temperature: the frost
mass per unit area at each latitude. Because it changes the shape of the state,
it is a subclass rather than a flag on ``EBM`` --- the Earth model is untouched.
"""

from __future__ import annotations

import numpy as np

from .config import Config
from .ebm import EBM, KELVIN

# CO2 saturation vapour pressure constants: p(Pa) = A_SVP * exp(-B_SVP / T_K)
_A_SVP = 1.382e12
_B_SVP = 3182.48
_LN_A = float(np.log(_A_SVP))


def co2_frost_point_K(p_pa) -> np.ndarray | float:
    """CO2 frost-point temperature [K] at pressure ``p_pa`` [Pa]."""
    p = np.maximum(np.asarray(p_pa, dtype=float), 1e-6)
    return _B_SVP / (_LN_A - np.log(p))


def co2_saturation_pressure_Pa(T_K) -> np.ndarray | float:
    """CO2 saturation vapour pressure [Pa] at temperature ``T_K`` [K]."""
    T = np.maximum(np.asarray(T_K, dtype=float), 1.0)
    return _A_SVP * np.exp(-_B_SVP / T)


class MarsEBM(EBM):
    """EBM with a condensing CO2 atmosphere.

    The state is ``(T, m_frost)``: temperature [degC] and frost mass per unit
    area [kg/m^2] at each latitude. Total CO2 --- atmospheric plus condensed ---
    is conserved exactly by construction.
    """

    def __init__(self, config: Config):
        super().__init__(config)
        if self.two_surface:
            raise NotImplementedError(
                "MarsEBM does not support the two-surface land/ocean split; "
                "Mars has no ocean, so a single surface is the right model.")
        self.inventory = float(config.co2_inventory_kg_m2)
        self.latent = float(config.co2_latent_heat)
        self.gravity = float(config.surface_gravity)
        self.frost_coalbedo = 1.0 - float(config.co2_frost_albedo)
        # Frost present at the current step; consulted by coalbedo().
        self._m_frost = np.zeros(self.n)
        # Surface pressure at the current step; consulted by _graygas_tau().
        self._p_now = self.surface_pressure(self._m_frost)

    # ------------------------------------------------------------------
    # Atmospheric state
    # ------------------------------------------------------------------
    def surface_pressure(self, m_frost: np.ndarray) -> float:
        """Global surface pressure [Pa] given the condensed mass distribution.

        Equal-area grid, so the plain cell mean is the global mean column of
        condensed CO2. What remains in the atmosphere is the inventory minus
        that, and pressure is its weight.
        """
        condensed = float(np.mean(m_frost))
        column = max(self.inventory - condensed, 0.0)
        return column * self.gravity

    def frost_point_degC(self, m_frost: np.ndarray) -> float:
        """Frost-point temperature [degC] at the current surface pressure."""
        return float(co2_frost_point_K(self.surface_pressure(m_frost))) - KELVIN

    def total_co2(self, m_frost: np.ndarray) -> float:
        """Total CO2 inventory [kg/m^2]; must be conserved. Diagnostic."""
        return float(np.mean(m_frost)) + max(
            self.inventory - float(np.mean(m_frost)), 0.0)

    def airborne_fraction(self, m_frost: np.ndarray) -> float:
        """Fraction of the CO2 inventory still in the atmosphere."""
        if self.inventory <= 0:
            return 0.0
        return max(self.inventory - float(np.mean(m_frost)), 0.0) / self.inventory

    def is_collapsed(self, m_frost: np.ndarray) -> bool:
        """Has the atmosphere essentially all condensed onto the surface?

        Beyond this point the model reports numbers it is not entitled to. The
        frost point is obtained by inverting the vapour-pressure curve, and as
        ``p -> 0`` that inversion degenerates -- it tends toward ~76 K rather
        than diverging, so the model will happily keep condensing against a
        frost point that has lost its meaning.

        This is flagged rather than clamped, on the same reasoning as the
        Simpson-Nakajima runaway in the Earth model: a cap would manufacture a
        state that does not physically exist, whereas a flag says the model has
        left its domain and the honest output is "the atmosphere freezes out"
        rather than a temperature.
        """
        return self.airborne_fraction(m_frost) < self.cfg.co2_collapse_threshold

    # ------------------------------------------------------------------
    # Physics overrides
    # ------------------------------------------------------------------
    def coalbedo(self, T: np.ndarray) -> np.ndarray:
        """Absorbed fraction, bright wherever CO2 frost is lying.

        Earth's version switches on temperature. Here it switches on the
        *presence of frost*, which is not the same thing: frost deposited in
        winter persists into spring while the surface warms, so albedo lags
        temperature. That lag is physical and affects the cap retreat.
        """
        P2 = 0.5 * (3.0 * self._x_eff ** 2 - 1.0)
        bare = self.cfg.coalbedo_a0 + self.cfg.coalbedo_a2 * P2
        return np.where(self._m_frost > 0.0, self.frost_coalbedo, bare)

    def _graygas_tau(self) -> float:
        """Grey-gas optical depth at the current surface pressure.

        tau = tau_ref * (p / p_ref), so a planet whose caps sublimate develops a
        real greenhouse. This closes a feedback the fixed-emissivity model could
        not express:

            warmer -> caps sublimate -> higher pressure -> stronger greenhouse
                   -> warmer

        which is positive, and therefore capable of producing a tipping point in
        the same way the ice-albedo feedback does. Only meaningful when the CO2
        inventory is large enough for pressure to move appreciably.
        """
        p_ref = max(float(self.cfg.graygas_p_ref_pa), 1e-9)
        return float(self.cfg.graygas_tau_ref) * self._p_now / p_ref

    # ------------------------------------------------------------------
    # Time integration
    # ------------------------------------------------------------------
    def step_co2(self, T: np.ndarray, m_frost: np.ndarray, M: float):
        """Advance ``(T, m_frost)`` one timestep.

        The radiative/dynamical step runs first, producing the temperature the
        surface *would* reach with no phase change. Condensation or sublimation
        then reconciles that against the frost point, converting the surplus or
        deficit into latent heat and mass.
        """
        self._m_frost = m_frost                    # seen by coalbedo()
        self._p_now = self.surface_pressure(m_frost)   # seen by _graygas_tau()
        T_pred = super().step(T, M)                # radiation + transport + diffusion
        m = m_frost.copy()

        T_f = self.frost_point_degC(m)
        C = self.C                                  # J/m^2/degC

        # --- Condensation: surface wants to fall below the frost point ---
        cond = T_pred < T_f
        if np.any(cond):
            deficit = C[cond] * (T_f - T_pred[cond])       # J/m^2 to be absorbed
            m[cond] += deficit / self.latent
            T_pred[cond] = T_f

        # --- Sublimation: frost present and the surface wants to warm ---
        subl = (~cond) & (m > 0.0) & (T_pred > T_f)
        if np.any(subl):
            surplus = C[subl] * (T_pred[subl] - T_f)       # J/m^2 available
            need = surplus / self.latent                    # kg/m^2 to remove all of it
            have = m[subl]
            removed = np.minimum(have, need)
            m[subl] = have - removed
            # Fully buffered where frost remains; leftover energy warms the rest.
            leftover = (need - removed) * self.latent
            T_pred[subl] = T_f + leftover / C[subl]

        self._m_frost = m
        return T_pred, m

    def spin_up_co2(self, T0=None, m0=None, M0: float = 0.0):
        """Integrate whole years until the seasonal cycle repeats.

        Convergence is judged on the global mean temperature, as in the base
        class, but the frost field must also settle --- a model can reach a
        steady mean while its caps are still growing year on year, so the
        cap mass is checked too.
        """
        T = np.full(self.n, 15.0) if T0 is None else np.array(T0, dtype=float)
        m = np.zeros(self.n) if m0 is None else np.array(m0, dtype=float)

        n_steps, dM = self._year_steps()
        M = float(M0)
        prev_T = self.global_mean(T)
        prev_m = float(np.mean(m))
        drift_T = drift_m = np.inf
        # Frost tolerance scales with the inventory; an absolute one would be
        # tight for a thin atmosphere and meaningless for a thick one.
        m_tol = max(self.cfg.co2_spinup_tol_frac * self.inventory, 1e-12)
        streak = 0
        years = 0
        for years in range(1, self.cfg.spinup_max_years + 1):
            for _ in range(n_steps):
                T, m = self.step_co2(T, m, M)
                M += dM
            cur_T, cur_m = self.global_mean(T), float(np.mean(m))
            drift_T, drift_m = abs(cur_T - prev_T), abs(cur_m - prev_m)
            prev_T, prev_m = cur_T, cur_m
            if drift_T < self.cfg.spinup_tol_degC and drift_m < m_tol:
                streak += 1
                if streak >= self.cfg.spinup_consecutive:
                    break
            else:
                streak = 0          # a creeping state never accumulates a streak
        converged = streak >= self.cfg.spinup_consecutive
        return T, m, M, {"years": years, "drift_T": drift_T, "drift_m": drift_m,
                         "converged": converged}

    def run_year_co2(self, T: np.ndarray, m: np.ndarray, M0: float = 0.0):
        """Integrate one orbital year, recording temperature, frost and pressure."""
        n_steps, dM = self._year_steps()
        M = float(M0)
        rec_T = np.empty((n_steps, self.n))
        rec_m = np.empty((n_steps, self.n))
        rec_p = np.empty(n_steps)
        rec_day = np.empty(n_steps)
        rec_collapsed = np.zeros(n_steps, dtype=bool)
        for k in range(n_steps):
            T, m = self.step_co2(T, m, M)
            M += dM
            rec_T[k], rec_m[k] = T, m
            rec_p[k] = self.surface_pressure(m)
            rec_collapsed[k] = self.is_collapsed(m)
            rec_day[k] = (M - M0) / (2.0 * np.pi) * self.cfg.days_per_year
        return T, m, M, {"T": rec_T, "m_frost": rec_m,
                         "pressure_Pa": rec_p, "day": rec_day,
                         "collapsed": rec_collapsed,
                         # Fraction of the year spent with the atmosphere frozen
                         # out. Any value above zero means the run's reported
                         # temperatures and pressures are outside the model's
                         # domain for at least part of the year.
                         "collapsed_fraction": float(np.mean(rec_collapsed))}
