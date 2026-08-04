"""One-dimensional (latitude-resolved) seasonal energy-balance model.

North / Budyko / Sellers formulation:

    C dT/dt = Q(x,t) * a(x,T) - (A + B*T) + D d/dx[(1 - x^2) dT/dx]

with

    x   = sin(latitude)        independent variable (equal-area in x)
    T   = surface temperature  [degC]
    Q   = daily-mean insolation from :mod:`orbital_climate.insolation`
    a   = coalbedo (absorbed fraction), with an ice-albedo step
    A,B = linear OLR coefficients
    D   = meridional diffusion coefficient
    C   = mixed-layer heat capacity

Discretisation
--------------
The domain x in [-1, 1] is split into ``n`` equal cells (cell-centred). The
diffusion term is a flux divergence with weight ``w = 1 - x^2`` evaluated at
cell interfaces; because ``w = 0`` at the poles (x = +/-1), the no-flux
boundary condition is satisfied automatically and the operator conserves
energy exactly (its area-weighted sum is zero for any field).

Time stepping is **semi-implicit (IMEX)**: the linear diffusion and OLR terms
are treated implicitly (a constant tridiagonal solve, factorised once), while
the nonlinear insolation*coalbedo source is evaluated explicitly at the old
temperature. This is unconditionally stable for the stiff diffusion/relaxation
terms, so the timestep is set by accuracy, not stability.
"""

from __future__ import annotations

import numpy as np
from scipy.linalg import lu_factor, lu_solve

from .config import Config
from .insolation import daily_mean_insolation, annual_mean_insolation

TWO_PI = 2.0 * np.pi
SECONDS_PER_DAY = 86400.0
SIGMA_SB = 5.670374419e-8      # Stefan-Boltzmann constant [W m^-2 K^-4]
KELVIN = 273.15

# Earth's zonal land fraction, by 10-degree band centre (south -> north).
# Antarctica dominates the far south; the northern mid-latitudes are the most
# continental; the southern mid-latitudes are almost pure ocean.
_LAND_FRAC_LAT = np.array([-85., -75., -65., -55., -45., -35., -25., -15., -5.,
                           5., 15., 25., 35., 45., 55., 65., 75., 85.])
_LAND_FRAC_VAL = np.array([1.00, 0.80, 0.10, 0.01, 0.03, 0.11, 0.23, 0.23, 0.24,
                           0.23, 0.28, 0.38, 0.45, 0.52, 0.57, 0.53, 0.25, 0.05])


def earth_land_fraction(lat_deg: np.ndarray) -> np.ndarray:
    """Land fraction at the given latitudes, interpolated from Earth's zonal profile."""
    return np.interp(np.asarray(lat_deg, dtype=float), _LAND_FRAC_LAT, _LAND_FRAC_VAL)


class EBM:
    """Seasonal 1-D energy-balance model on an equal-area latitude grid."""

    def __init__(self, config: Config):
        self.cfg = config
        n = int(config.n_lat)
        self.n = n

        # Cell-centred grid in x = sin(latitude); equal cell width => equal area.
        self.dx = 2.0 / n
        self.x = -1.0 + (np.arange(n) + 0.5) * self.dx          # cell centres
        self.phi = np.arcsin(np.clip(self.x, -1.0, 1.0))        # latitude [rad]
        self.lat_deg = np.degrees(self.phi)

        # Interface positions and diffusion weights w = 1 - x^2 at interfaces.
        x_iface = -1.0 + np.arange(n + 1) * self.dx             # n+1 interfaces
        w_iface = 1.0 - x_iface ** 2                            # zero at both poles

        # Tridiagonal diffusion operator L: (L T)_i = [w_{i+1/2}(T_{i+1}-T_i)
        #   - w_{i-1/2}(T_i - T_{i-1})] / dx^2. Boundary fluxes vanish (w=0).
        self.L = self._build_diffusion_matrix(w_iface)

        self.two_surface = bool(config.two_surface)
        self.n_surf = 2 if self.two_surface else 1
        self._dt_s = config.dt_days * SECONDS_PER_DAY

        if self.two_surface:
            # Land fraction per cell, clamped away from 0/1 so the energy-conserving
            # exchange term (which scales by f/(1-f)) stays finite at the poles.
            if config.land_fraction_override is not None:
                f = np.full(n, float(config.land_fraction_override))
            else:
                f = earth_land_fraction(self.lat_deg)
            self.land_frac = np.clip(f, 1e-3, 1.0 - 1e-3)
            # State vector is [T_land (n), T_ocean (n)].
            self.C = np.concatenate([
                np.full(n, float(config.heat_capacity_land)),
                np.full(n, float(config.heat_capacity_ocean)),
            ])
            self._x_eff = np.tile(self.x, 2)
        else:
            self.land_frac = None
            self.C = np.full(n, float(config.heat_capacity))
            self._x_eff = self.x

        # Precompute the constant implicit operator for the semi-implicit step
        # and factorise it once; it is reused every step.
        self._lu = lu_factor(self._build_implicit_operator())

    # ------------------------------------------------------------------
    # Operator construction
    # ------------------------------------------------------------------
    def _build_implicit_operator(self) -> np.ndarray:
        """Assemble the implicit matrix solved each timestep.

        Single surface:  M = diag(C/dt + B) - D L

        Two surfaces, with land fraction f and exchange coefficient nu:

            C_l dT_l/dt = ... + nu (T_o - T_l)
            C_o dT_o/dt = ... - nu (T_o - T_l) f/(1-f)

        The f/(1-f) factor makes the exchange conserve energy: the flux per unit
        *land* area entering the land column equals the flux per unit *ocean*
        area leaving the ocean column once both are weighted by their areas.
        """
        n, cfg = self.n, self.cfg
        diff = cfg.diffusion_D * self.L
        if not self.two_surface:
            return np.diag(self.C / self._dt_s + cfg.olr_B) - diff

        nu = float(cfg.land_ocean_coupling)
        g = nu * self.land_frac / (1.0 - self.land_frac)   # ocean-side scaling
        C_l, C_o = self.C[:n], self.C[n:]

        M = np.zeros((2 * n, 2 * n))
        M[:n, :n] = np.diag(C_l / self._dt_s + cfg.olr_B + nu) - diff
        M[n:, n:] = np.diag(C_o / self._dt_s + cfg.olr_B + g) - diff
        M[:n, n:] = -nu * np.eye(n)
        M[n:, :n] = -np.diag(g)
        return M

    def _build_diffusion_matrix(self, w_iface: np.ndarray) -> np.ndarray:
        n = self.n
        dx2 = self.dx ** 2
        L = np.zeros((n, n))
        for i in range(n):
            wl = w_iface[i]        # left interface  (i-1/2)
            wr = w_iface[i + 1]    # right interface (i+1/2)
            if i > 0:
                L[i, i - 1] = wl / dx2
            if i < n - 1:
                L[i, i + 1] = wr / dx2
            L[i, i] = -(wl + wr) / dx2
        return L

    # ------------------------------------------------------------------
    # Physics helpers
    # ------------------------------------------------------------------
    def coalbedo(self, T: np.ndarray) -> np.ndarray:
        """Absorbed fraction a(x,T) with an ice-albedo step at ``T_ice``.

        In two-surface mode this is evaluated independently on each surface, so
        land can be frozen while the ocean at the same latitude is not.
        """
        P2 = 0.5 * (3.0 * self._x_eff ** 2 - 1.0)
        a_noice = self.cfg.coalbedo_a0 + self.cfg.coalbedo_a2 * P2
        return np.where(T > self.cfg.T_ice_degC, a_noice, self.cfg.coalbedo_ice)

    # ------------------------------------------------------------------
    # State-vector helpers (single vs two surface)
    # ------------------------------------------------------------------
    def split(self, T: np.ndarray):
        """Return ``(T_land, T_ocean)``; both are ``T`` itself in single-surface mode."""
        T = np.asarray(T, dtype=float)
        if not self.two_surface:
            return T, T
        return T[..., :self.n], T[..., self.n:]

    def blend(self, T: np.ndarray) -> np.ndarray:
        """Area-weighted land/ocean blend, giving one temperature per latitude.

        Idempotent: a state that is already blended (last axis of length ``n``)
        is returned unchanged, so callers may blend defensively without having to
        track whether an upstream step already did.
        """
        T = np.asarray(T, dtype=float)
        if not self.two_surface or T.shape[-1] == self.n:
            return T
        T_l, T_o = self.split(T)
        return self.land_frac * T_l + (1.0 - self.land_frac) * T_o

    def olr(self, T: np.ndarray) -> np.ndarray:
        """Outgoing longwave radiation [W m^-2] for temperature ``T`` [degC].

        ``linear``  -- Budyko (1969):  A + B*T
        ``sellers`` -- Sellers (1969): sigma T_K^4 [1 - m tanh(19 T_K^6 1e-16)]

        The Sellers form matches the linear one at present-day Earth (235.1 vs
        234.3 W/m^2 at 288 K) but, unlike it, tends to blackbody emission as the
        planet freezes and its atmosphere dries -- the physically correct limit.
        """
        model = str(self.cfg.olr_model).lower()
        if model == "linear":
            return self.cfg.olr_A + self.cfg.olr_B * np.asarray(T, dtype=float)
        if model == "sellers":
            # Clamp to a small positive absolute temperature: the explicit
            # source can transiently probe unphysical values during spin-up.
            T_K = np.maximum(np.asarray(T, dtype=float) + KELVIN, 1.0)
            return SIGMA_SB * T_K ** 4 * (
                1.0 - self.cfg.sellers_m * np.tanh(19.0 * T_K ** 6 * 1e-16))
        raise ValueError(f"Unknown olr_model {self.cfg.olr_model!r}; "
                         "expected 'linear' or 'sellers'.")

    def insolation(self, M: float) -> np.ndarray:
        """Daily-mean insolation Q(x) at orbital mean anomaly ``M`` [rad]."""
        return daily_mean_insolation(self.phi, float(M), self.cfg)

    def annual_mean_insolation(self) -> np.ndarray:
        """Annual-mean insolation profile Q(x) (for annual-mean equilibria)."""
        return annual_mean_insolation(self.phi, self.cfg)

    # ------------------------------------------------------------------
    # Diagnostics (equal-area => simple means over cells)
    # ------------------------------------------------------------------
    def global_mean(self, T: np.ndarray) -> float:
        """Area-weighted global mean. Equal-area grid => plain cell mean.

        In two-surface mode the land/ocean blend is taken first, so each surface
        contributes in proportion to the area it actually covers.
        """
        return float(np.mean(self.blend(T)))

    def ice_line_lat(self, T: np.ndarray):
        """Latitude [deg] of the northern-hemisphere ice edge, or None if none.

        Ice edge = poleward-most latitude in the NH where T crosses ``T_ice``.
        Returns +90 if the NH is ice-free to the pole, or None if fully glaciated.
        Two-surface states are blended first, matching :meth:`global_mean`.
        """
        T = self.blend(T)
        north = self.lat_deg >= 0
        Tn = T[north]
        latn = self.lat_deg[north]
        icy = Tn <= self.cfg.T_ice_degC
        if not icy.any():
            return 90.0
        if icy.all():
            return None
        # First icy cell moving poleward.
        idx = np.argmax(icy)  # lowest-latitude icy cell in NH
        return float(latn[idx])

    # ------------------------------------------------------------------
    # Time integration
    # ------------------------------------------------------------------
    def step(self, T: np.ndarray, M: float) -> np.ndarray:
        """Advance temperature one timestep with the semi-implicit scheme.

        The implicit operator carries a linear ``B*T`` relaxation term so it can
        stay constant and be factorised once. Any nonlinear OLR is handled by
        adding ``B*T - OLR(T)`` to the explicit source: the ``B*T`` contributions
        cancel exactly at convergence, so ``olr_B`` acts purely as a numerical
        preconditioner and the converged state satisfies the true balance

            D d/dx[(1-x^2) dT/dx] + Q a(T) - OLR(T) = 0

        With ``olr_model = "linear"`` this reduces identically to the original
        scheme (the source collapses to ``Q a(T) - A``).
        """
        Q = self.insolation(M)
        if self.two_surface:
            Q = np.tile(Q, 2)          # both surfaces see the same insolation
        source = Q * self.coalbedo(T) - self.olr(T) + self.cfg.olr_B * T
        rhs = (self.C / self._dt_s) * T + source
        return lu_solve(self._lu, rhs)

    def _year_steps(self):
        """Number of timesteps per orbital year and the per-step dM."""
        n_steps = int(round(self.cfg.days_per_year / self.cfg.dt_days))
        dM = TWO_PI / n_steps
        return n_steps, dM

    def spin_up(self, T0: np.ndarray | None = None, M0: float = 0.0):
        """Integrate whole years until the seasonal cycle repeats.

        Returns ``(T, M, info)`` where ``T`` is the state at the end of the
        last completed year, ``M`` the corresponding mean anomaly, and ``info``
        a dict with the number of years run and the final year-over-year drift.
        """
        if T0 is None:
            # Start from the annual-mean radiative-diffusive guess.
            T = np.full(self.n_surf * self.n, 15.0)
        else:
            T = np.array(T0, dtype=float)
            if T.size == self.n and self.two_surface:
                T = np.tile(T, 2)      # promote a single-surface profile

        n_steps, dM = self._year_steps()
        M = float(M0)
        prev_mean = self.global_mean(T)
        drift = np.inf
        years = 0
        for years in range(1, self.cfg.spinup_max_years + 1):
            for _ in range(n_steps):
                T = self.step(T, M)
                M += dM
            cur_mean = self.global_mean(T)
            drift = abs(cur_mean - prev_mean)
            prev_mean = cur_mean
            if drift < self.cfg.spinup_tol_degC:
                break
        return T, M, {"years": years, "final_drift_degC": drift}

    def run_year(self, T: np.ndarray, M0: float = 0.0, record: bool = True):
        """Integrate one orbital year, optionally recording per-step diagnostics.

        Returns ``(T_end, M_end, record_dict)``. When ``record`` is True the
        record contains arrays over the year: mean anomaly ``M``, day-of-year,
        global-mean T, and full temperature field ``T`` (shape [n_steps, n_lat]).
        """
        n_steps, dM = self._year_steps()
        M = float(M0)
        rec_M = np.empty(n_steps)
        rec_gm = np.empty(n_steps)
        rec_T = np.empty((n_steps, self.n_surf * self.n)) if record else None
        for k in range(n_steps):
            T = self.step(T, M)
            M += dM
            if record:
                rec_M[k] = M
                rec_gm[k] = self.global_mean(T)
                rec_T[k] = T
        out = None
        if record:
            day = (rec_M - M0) / TWO_PI * self.cfg.days_per_year
            out = {"M": rec_M, "day": day, "global_mean": rec_gm, "T": rec_T}
        return T, M, out
