"""Experiments on the EBM: equilibrium and sudden-orbital-change transients.

Two entry points:

* :func:`run_equilibrium` -- spin a configuration up to its repeating seasonal
  cycle and report equilibrium diagnostics.
* :func:`run_perturbation` -- spin up a *baseline* orbit, then switch the
  orbital elements at t=0 and integrate the *perturbed* orbit forward,
  recording year-by-year transient diagnostics (global/hemispheric mean
  temperature, ice-edge latitude, and the 65 N summer peak temperature and
  insolation that drive Milankovitch glacial inception).

Nothing here changes the physics in :mod:`orbital_climate.ebm`; the EBM's heat
capacity provides the thermal inertia that lags and damps the transient.
"""

from __future__ import annotations

import dataclasses
from dataclasses import dataclass, field
from pathlib import Path

import numpy as np

from .config import Config
from .ebm import EBM
from .insolation import daily_mean_insolation

TWO_PI = 2.0 * np.pi


def _nearest_cell(ebm: EBM, lat_deg: float) -> int:
    return int(np.argmin(np.abs(ebm.lat_deg - lat_deg)))


def _annual_diags(ebm: EBM, rec: dict, i_diag: int) -> dict:
    """Reduce a recorded year to scalar diagnostics.

    ``i_diag`` is the grid cell of the diagnostic latitude (``Config.diag_lat_deg``);
    its seasonal max/min are reported as ``Tdiag_summer_max`` / ``Tdiag_winter_min``.
    """
    T_raw = rec["T"]                               # [n_steps, n_surf*n_lat]
    T = ebm.blend(T_raw)                            # [n_steps, n_lat] area-blended
    T_annual = T.mean(axis=0)                       # annual-mean profile
    north = ebm.lat_deg >= 0
    south = ebm.lat_deg < 0
    out = {
        "global_mean": float(T_annual.mean()),
        "nh_mean": float(T_annual[north].mean()),
        "sh_mean": float(T_annual[south].mean()),
        "iceline_lat_nh": ebm.ice_line_lat(T_annual),
        "Tdiag_summer_max": float(T[:, i_diag].max()),
        "Tdiag_winter_min": float(T[:, i_diag].min()),
    }
    if ebm.two_surface:
        # Land and ocean seasonal extremes separately -- the blend above damps
        # the continental swing that motivates the two-surface model.
        T_land, T_ocean = ebm.split(T_raw)
        out.update({
            "Tdiag_land_summer_max": float(T_land[:, i_diag].max()),
            "Tdiag_land_winter_min": float(T_land[:, i_diag].min()),
            "Tdiag_ocean_summer_max": float(T_ocean[:, i_diag].max()),
            "Tdiag_ocean_winter_min": float(T_ocean[:, i_diag].min()),
        })
    return out


@dataclass
class EquilibriumResult:
    config: Config
    T: np.ndarray                      # equilibrium temperature profile [degC]
    lat_deg: np.ndarray
    spinup_years: int
    global_mean: float
    iceline_lat_nh: float | None
    diag_lat_deg: float                # diagnostic latitude actually used [deg]
    peak_diag_insol: float             # peak daily-mean insolation at diag_lat [W/m^2]
    # Seasonal structure at equilibrium. The annual mean above averages these
    # away, so a world with a mild mean can still have lethal seasonal extremes.
    nh_mean: float                     # annual-mean T, northern hemisphere [degC]
    sh_mean: float                     # annual-mean T, southern hemisphere [degC]
    diag_summer_max: float             # warmest daily-mean T at diag_lat [degC]
    diag_winter_min: float             # coldest daily-mean T at diag_lat [degC]
    # Two-surface mode only (None otherwise): the land/ocean split that the
    # area-weighted blend above averages away.
    diag_land_summer_max: float | None = None
    diag_land_winter_min: float | None = None
    diag_ocean_summer_max: float | None = None
    diag_ocean_winter_min: float | None = None
    # Runaway-greenhouse diagnostic (Simpson-Nakajima). When the absorbed flux
    # exceeds what a moist atmosphere can radiate there is no equilibrium at
    # all, so the reported temperature is meaningless rather than merely
    # imprecise.
    absorbed_mean: float = float("nan")   # global annual-mean absorbed flux [W/m^2]
    runaway: bool = False

    @property
    def diag_seasonal_range(self) -> float:
        """Peak-to-trough seasonal temperature range at the diagnostic latitude [K]."""
        return self.diag_summer_max - self.diag_winter_min

    @property
    def diag_land_seasonal_range(self) -> float | None:
        """Seasonal range over land at the diagnostic latitude [K], if available."""
        if self.diag_land_summer_max is None:
            return None
        return self.diag_land_summer_max - self.diag_land_winter_min

    @property
    def diag_ocean_seasonal_range(self) -> float | None:
        """Seasonal range over ocean at the diagnostic latitude [K], if available."""
        if self.diag_ocean_summer_max is None:
            return None
        return self.diag_ocean_summer_max - self.diag_ocean_winter_min


def run_equilibrium(config: Config) -> EquilibriumResult:
    """Spin ``config`` up to equilibrium and return diagnostics.

    The recorded final year yields the seasonal extremes for free -- the same
    reduction :func:`_annual_diags` applies to each year of a transient, so
    equilibrium and transient diagnostics are defined identically.
    """
    ebm = EBM(config)
    T, M, info = ebm.spin_up()
    # One more recorded year for annual-mean profile + diagnostic peak insolation.
    _, _, rec = ebm.run_year(T, M0=M)
    T_annual = ebm.blend(rec["T"]).mean(axis=0)
    i_diag = _nearest_cell(ebm, config.diag_lat_deg)
    phi_diag = ebm.phi[i_diag]
    M_year = np.linspace(0.0, TWO_PI, 2000, endpoint=False)
    peak_insol = float(daily_mean_insolation(phi_diag, M_year, config).max())

    # Runaway check: at equilibrium the global-mean OLR equals the global-mean
    # absorbed flux, so comparing the latter to the Simpson-Nakajima ceiling
    # tells us whether a moist atmosphere could have balanced at all.
    absorbed = 0.0
    T_steps = rec["T"]
    for k in range(T_steps.shape[0]):
        Q = ebm.insolation(rec["M"][k])
        if ebm.two_surface:
            Q = np.tile(Q, 2)
        absorbed += float(np.mean(ebm.blend(Q * ebm.coalbedo(T_steps[k]))))
    absorbed /= T_steps.shape[0]

    diags = _annual_diags(ebm, rec, i_diag)
    return EquilibriumResult(
        config=config,
        T=T_annual,
        lat_deg=ebm.lat_deg,
        spinup_years=info["years"],
        global_mean=diags["global_mean"],
        iceline_lat_nh=diags["iceline_lat_nh"],
        diag_lat_deg=float(ebm.lat_deg[i_diag]),
        peak_diag_insol=peak_insol,
        nh_mean=diags["nh_mean"],
        sh_mean=diags["sh_mean"],
        diag_summer_max=diags["Tdiag_summer_max"],
        diag_winter_min=diags["Tdiag_winter_min"],
        diag_land_summer_max=diags.get("Tdiag_land_summer_max"),
        diag_land_winter_min=diags.get("Tdiag_land_winter_min"),
        diag_ocean_summer_max=diags.get("Tdiag_ocean_summer_max"),
        diag_ocean_winter_min=diags.get("Tdiag_ocean_winter_min"),
        absorbed_mean=absorbed,
        runaway=bool(absorbed > config.olr_runaway_limit),
    )


@dataclass
class PerturbationResult:
    baseline: Config
    perturbed: Config
    lat_deg: np.ndarray
    years: np.ndarray                  # transient year index (1..n_years)
    diagnostics: dict                  # name -> array over years
    T_baseline_eq: np.ndarray          # baseline equilibrium annual profile
    T_final: np.ndarray                # perturbed profile after n_years
    baseline_spinup_years: int
    diag_lat_deg: float                # diagnostic latitude actually used [deg]


def run_perturbation(
    baseline: Config,
    perturbed: Config,
    n_years: int = 60,
) -> PerturbationResult:
    """Spin up ``baseline``, switch to ``perturbed`` at t=0, track the transient.

    The perturbed run reuses the baseline equilibrium temperature field as its
    initial condition and the same grid/heat-capacity, so only the orbital
    forcing (and any other changed parameters) differ.
    """
    base_ebm = EBM(baseline)
    T0, M0, base_info = base_ebm.spin_up()
    _, _, base_rec = base_ebm.run_year(T0, M0=M0)
    T_base_annual = base_ebm.blend(base_rec["T"]).mean(axis=0)

    pert_ebm = EBM(perturbed)
    i_diag = _nearest_cell(pert_ebm, perturbed.diag_lat_deg)

    diag_names = [
        "global_mean", "nh_mean", "sh_mean",
        "iceline_lat_nh", "Tdiag_summer_max", "Tdiag_winter_min",
    ]
    if pert_ebm.two_surface:
        diag_names += [
            "Tdiag_land_summer_max", "Tdiag_land_winter_min",
            "Tdiag_ocean_summer_max", "Tdiag_ocean_winter_min",
        ]
    diags = {k: [] for k in diag_names}

    T = np.array(T0, dtype=float)
    M = 0.0
    for _ in range(n_years):
        T, M, rec = pert_ebm.run_year(T, M0=M)
        yd = _annual_diags(pert_ebm, rec, i_diag)
        for k in diag_names:
            diags[k].append(yd[k])

    diags = {k: np.array(v, dtype=float) for k, v in diags.items()}
    return PerturbationResult(
        baseline=baseline,
        perturbed=perturbed,
        lat_deg=pert_ebm.lat_deg,
        years=np.arange(1, n_years + 1),
        diagnostics=diags,
        T_baseline_eq=T_base_annual,
        T_final=pert_ebm.blend(T),
        baseline_spinup_years=base_info["years"],
        diag_lat_deg=float(pert_ebm.lat_deg[i_diag]),
    )


# ---------------------------------------------------------------------------
# Output helpers
# ---------------------------------------------------------------------------
def save_perturbation_outputs(result: PerturbationResult, out_dir: str | Path) -> Path:
    """Write transient CSV + summary + plots to ``out_dir`` (created if needed)."""
    out = Path(out_dir)
    out.mkdir(parents=True, exist_ok=True)

    # Transient CSV.
    import csv
    csv_path = out / "transient.csv"
    cols = ["year"] + list(result.diagnostics.keys())
    with open(csv_path, "w", newline="", encoding="utf-8") as fh:
        writer = csv.writer(fh)
        writer.writerow(cols)
        for j, yr in enumerate(result.years):
            row = [int(yr)] + [
                ("" if result.diagnostics[k][j] is None or np.isnan(result.diagnostics[k][j])
                 else f"{result.diagnostics[k][j]:.5f}")
                for k in result.diagnostics
            ]
            writer.writerow(row)

    _plot_perturbation(result, out)
    _write_summary(result, out)
    return out


def _plot_perturbation(result: PerturbationResult, out: Path) -> None:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    d = result.diagnostics
    yrs = result.years

    fig, axs = plt.subplots(2, 2, figsize=(11, 8))
    axs[0, 0].plot(yrs, d["global_mean"], label="global")
    axs[0, 0].plot(yrs, d["nh_mean"], label="N. hemisphere")
    axs[0, 0].plot(yrs, d["sh_mean"], label="S. hemisphere")
    axs[0, 0].set_xlabel("year after perturbation")
    axs[0, 0].set_ylabel("mean T [degC]")
    axs[0, 0].set_title("Mean temperature transient")
    axs[0, 0].legend()

    diag_label = f"{result.diag_lat_deg:.0f} deg"
    axs[0, 1].plot(yrs, d["Tdiag_summer_max"], color="tab:red")
    axs[0, 1].set_xlabel("year after perturbation")
    axs[0, 1].set_ylabel(f"{diag_label} summer peak T [degC]")
    axs[0, 1].set_title(f"{diag_label} summer peak temperature")

    axs[1, 0].plot(yrs, d["iceline_lat_nh"], color="tab:blue")
    axs[1, 0].set_xlabel("year after perturbation")
    axs[1, 0].set_ylabel("NH ice-edge latitude [deg]")
    axs[1, 0].set_title("Northern ice-edge migration")
    axs[1, 0].invert_yaxis()

    axs[1, 1].plot(result.lat_deg, result.T_baseline_eq, label="baseline eq.")
    axs[1, 1].plot(result.lat_deg, result.T_final, label=f"year {int(yrs[-1])}")
    axs[1, 1].axhline(result.perturbed.T_ice_degC, ls="--", color="0.6", lw=1)
    axs[1, 1].set_xlabel("latitude [deg]")
    axs[1, 1].set_ylabel("annual-mean T [degC]")
    axs[1, 1].set_title("Annual-mean temperature profile")
    axs[1, 1].legend()

    fig.tight_layout()
    fig.savefig(out / "transient.png", dpi=130)
    plt.close(fig)


def _write_summary(result: PerturbationResult, out: Path) -> None:
    d = result.diagnostics
    b, p = result.baseline, result.perturbed
    lines = [
        "Sudden orbital-change experiment",
        "================================",
        f"baseline:  a={b.a_au:.4f} AU, e={b.ecc:.4f}, lon_peri={b.lon_perihelion_deg:.1f} deg",
        f"perturbed: a={p.a_au:.4f} AU, e={p.ecc:.4f}, lon_peri={p.lon_perihelion_deg:.1f} deg",
        f"baseline spin-up: {result.baseline_spinup_years} years",
        f"transient length: {int(result.years[-1])} years",
        "",
        f"global mean T:  {d['global_mean'][0]:.3f} -> {d['global_mean'][-1]:.3f} degC "
        f"(delta {d['global_mean'][-1] - d['global_mean'][0]:+.3f})",
        f"NH mean T:      {d['nh_mean'][0]:.3f} -> {d['nh_mean'][-1]:.3f} degC "
        f"(delta {d['nh_mean'][-1] - d['nh_mean'][0]:+.3f})",
        f"{result.diag_lat_deg:.0f}deg summer max: "
        f"{d['Tdiag_summer_max'][0]:.3f} -> {d['Tdiag_summer_max'][-1]:.3f} degC "
        f"(delta {d['Tdiag_summer_max'][-1] - d['Tdiag_summer_max'][0]:+.3f})",
        f"NH ice edge:    {d['iceline_lat_nh'][0]:.2f} -> {d['iceline_lat_nh'][-1]:.2f} deg",
    ]
    (out / "summary.txt").write_text("\n".join(lines) + "\n", encoding="utf-8")


def perturbed_config_from_scenario(baseline: Config, ecc: float,
                                   lon_perihelion_deg: float | None = None) -> Config:
    """Return a copy of ``baseline`` with eccentricity (and optionally lon_peri) changed.

    The context-doc scenario ``r_p -> 0.885, r_a -> 1.118`` maps to
    ``e = 0.117`` with ``a`` essentially unchanged, so the default perturbation
    is a pure eccentricity injection.
    """
    changes = {"ecc": ecc}
    if lon_perihelion_deg is not None:
        changes["lon_perihelion_deg"] = lon_perihelion_deg
    return dataclasses.replace(baseline, **changes)
