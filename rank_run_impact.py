"""Rank BH-flyby sweep runs by how much they disturbed the Solar System.

Walks a parent sweep folder (e.g. simulations/<STAMP>/), reads each run's
*__planets_run_deltas.csv (a/e/q before vs. after, written by
solar_system_bh_rebound26.py), and computes a single interestingness score per
run so the most disruptive parameter combinations can be found without
eyeballing 672 CSVs by hand.

Per-body raw impact combines three dimensionless, scale-free terms:

    E_i = |1 - a_before / a_after|      relative orbital-energy change
                                         (1/a is proportional to -energy; this
                                         stays continuous through the
                                         bound -> hyperbolic transition, so an
                                         ejected planet naturally scores >= 1
                                         with no special-casing)
    X_i = |de| / (1 - e_before)         eccentricity change as a fraction of
                                         the "gap to escape" (e=1)
    Q_i = |dq| / q_before                relative perihelion change

    r_i = E_i + X_i + Q_i

Raw scores span many orders of magnitude (a barely-perturbed Neptune vs. an
ejected Mercury), so each body's contribution is log-compressed against a
data-driven noise floor r0 (the median r_i across the whole sweep) before
being combined into the run's total score:

    s_i   = log10(1 + r_i / r0)
    Score = sum_i  w_i * s_i

Earth is up-weighted (default 5x) per the user's stated interest; the Moon is
excluded by default since its *heliocentric* elements swing with its monthly
orbit regardless of any BH, so before/after snapshots mostly capture lunar
phase rather than perturbation.

Usage
-----
    python rank_run_impact.py simulations/20260724_230314
    python rank_run_impact.py simulations/20260724_230314 --top 30 --out ranked.csv
"""

from __future__ import annotations

import argparse
import math
import re
from pathlib import Path

import numpy as np
import pandas as pd

DEFAULT_WEIGHTS = {
    "Mercury": 1.0, "Venus": 1.0, "Earth": 5.0, "Mars": 1.0,
    "Jupiter": 1.0, "Saturn": 1.0, "Uranus": 1.0, "Neptune": 1.0,
    "Moon": 0.0,
}

_RUN_NAME_RE = re.compile(
    r"__rp(?P<rp>[^_]+)__vinf(?P<vinf>[^_]+)__inc(?P<inc>[^_]+)"
    r"__toff(?P<toff>[^_]+)__Om(?P<Om>[^_]+)__om(?P<om>[^_]+)$"
)


def _token_to_float(token: str) -> float:
    """Invert solar_system_bh_rebound26.safe_num_for_filename: 'm' -> '-', 'p' -> '.'."""
    return float(token.replace("p", ".").replace("m", "-"))


def _parse_run_params(run_name: str) -> dict:
    """Recover the swept BH parameters from a run's subfolder name, or {} if it
    doesn't match the expected '..._rp<>__vinf<>__inc<>__toff<>__Om<>__om<>' pattern."""
    m = _RUN_NAME_RE.search(run_name)
    if not m:
        return {}
    return {
        "bh_rp_au": _token_to_float(m["rp"]),
        "bh_vinf_kms": _token_to_float(m["vinf"]),
        "bh_inc_deg": _token_to_float(m["inc"]),
        "bh_tperi_offset_days": _token_to_float(m["toff"]),
        "bh_Omega_deg": _token_to_float(m["Om"]),
        "bh_omega_deg": _token_to_float(m["om"]),
    }


def _per_body_raw_score(row: pd.Series) -> float:
    a_before, e_before, q_before = row["a_before"], row["e_before"], row["q_before"]
    a_after, e_after, q_after = row["a_after"], row["e_after"], row["q_after"]

    if a_after == 0.0:
        E = float("inf")
    else:
        E = abs(1.0 - a_before / a_after)

    gap = max(1.0 - e_before, 1e-6)   # avoid div-by-zero for near-parabolic starts
    X = abs(e_after - e_before) / gap

    Q = abs(q_after - q_before) / max(q_before, 1e-12)

    return E + X + Q


def _collect_run_files(parent_dir: Path) -> list[Path]:
    return sorted(parent_dir.glob("*/*__planets_run_deltas.csv"))


def rank_runs(
    parent_dir: Path,
    weights: dict[str, float] = DEFAULT_WEIGHTS,
    r0: float | None = None,
) -> pd.DataFrame:
    files = _collect_run_files(parent_dir)
    if not files:
        raise FileNotFoundError(
            f"No *__planets_run_deltas.csv files found under {parent_dir}/*/"
        )

    per_run_tables = []
    for f in files:
        df = pd.read_csv(f)
        df["run"] = f.parent.name
        df["raw_r"] = df.apply(_per_body_raw_score, axis=1)
        per_run_tables.append(df)
    all_rows = pd.concat(per_run_tables, ignore_index=True)

    # Data-driven noise floor: median raw impact across every (run, body) pair
    # that has any weight assigned (i.e. excluding the Moon by default).
    weighted_bodies = {b for b, w in weights.items() if w > 0}
    finite = all_rows[
        all_rows["body"].isin(weighted_bodies) & np.isfinite(all_rows["raw_r"])
    ]["raw_r"]
    if r0 is None:
        r0 = float(finite.median()) if len(finite) else 1e-4
        if r0 <= 0:
            r0 = 1e-4

    all_rows["s"] = all_rows["raw_r"].apply(
        lambda r: math.log10(1.0 + r / r0) if math.isfinite(r) else math.log10(1.0 + 1e12)
    )
    all_rows["weight"] = all_rows["body"].map(weights).fillna(0.0)
    all_rows["ejected"] = all_rows["e_after"] >= 1.0

    rows = []
    for run_name, g in all_rows.groupby("run"):
        score = float((g["s"] * g["weight"]).sum())
        earth_row = g[g["body"] == "Earth"]
        earth_s = float(earth_row["s"].iloc[0]) if len(earth_row) else float("nan")
        earth_da = float(earth_row["da"].iloc[0]) if len(earth_row) else float("nan")
        earth_de = float(earth_row["de"].iloc[0]) if len(earth_row) else float("nan")
        ejected = g.loc[g["ejected"], "body"].tolist()
        top_body_idx = g["raw_r"].idxmax()
        row = {
            "run": run_name,
            "Score": score,
            "Earth_s": earth_s,
            "Earth_da_AU": earth_da,
            "Earth_de": earth_de,
            "n_ejected": len(ejected),
            "ejected_bodies": ",".join(ejected) if ejected else "",
            "top_impact_body": g.loc[top_body_idx, "body"],
            "top_impact_raw_r": float(g.loc[top_body_idx, "raw_r"]),
        }
        row.update(_parse_run_params(run_name))
        rows.append(row)

    out = pd.DataFrame(rows).sort_values("Score", ascending=False).reset_index(drop=True)
    out.attrs["r0"] = r0
    return out


# ---------------------------------------------------------------------------
# Plot
# ---------------------------------------------------------------------------
# Sequential single-hue blue ramp (see the dataviz skill's reference palette).
# Ordinal steps (n_ejected = 0..4): lightest clears the 2:1 on-surface floor.
_ORDINAL_BLUE = ["#86b6ef", "#5598e7", "#2a78d6", "#1c5cab", "#104281"]
# Full continuous ramp (light -> dark) for the heatmap's magnitude encoding.
_SEQUENTIAL_BLUE_STEPS = [
    "#cde2fb", "#b7d3f6", "#9ec5f4", "#86b6ef", "#6da7ec", "#5598e7",
    "#3987e5", "#2a78d6", "#256abf", "#1c5cab", "#184f95", "#104281", "#0d366b",
]

_INK_PRIMARY = "#0b0b0b"
_INK_SECONDARY = "#52514e"
_INK_MUTED = "#898781"
_GRIDLINE = "#e1e0d9"
_AXIS = "#c3c2b7"
_SURFACE = "#fcfcfb"


def _style_axes(ax) -> None:
    ax.set_facecolor(_SURFACE)
    for spine_name, spine in ax.spines.items():
        if spine_name in ("top", "right"):
            spine.set_visible(False)
        else:
            spine.set_color(_AXIS)
    ax.tick_params(colors=_INK_MUTED, labelsize=9)
    ax.xaxis.label.set_color(_INK_SECONDARY)
    ax.yaxis.label.set_color(_INK_SECONDARY)
    ax.grid(True, color=_GRIDLINE, linewidth=0.8, zorder=0)
    ax.set_axisbelow(True)


def make_impact_plot(ranked: pd.DataFrame, out_path: Path) -> None:
    """Write a 2-panel figure: Score vs. rp (colored by ejection count) and a
    max-Score heatmap over the (rp, inclination) plane."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.colors import LinearSegmentedColormap, BoundaryNorm
    from matplotlib.patches import Patch

    if "bh_rp_au" not in ranked.columns or ranked["bh_rp_au"].isna().all():
        raise ValueError(
            "No run names matched the expected '..._rp<>__vinf<>__inc<>__toff<>"
            "__Om<>__om<>' pattern -- cannot plot vs. swept parameters."
        )

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(13, 5.5))
    fig.patch.set_facecolor(_SURFACE)

    # --- Panel A: Score vs periapsis, colored by ejection count (ordinal) ---
    _style_axes(ax1)
    n_ej = ranked["n_ejected"].clip(upper=4).astype(int)
    cmap_ord = LinearSegmentedColormap.from_list("ordinal_blue", _ORDINAL_BLUE, N=5)
    # rp only takes a handful of discrete swept values, so points at each one
    # stack in an exact vertical line; a small deterministic jitter separates
    # them without disturbing the x-axis reading.
    rng = np.random.default_rng(0)
    rp_span = ranked["bh_rp_au"].max() - ranked["bh_rp_au"].min()
    jitter = rng.uniform(-0.012, 0.012, size=len(ranked)) * max(rp_span, 1.0)
    ax1.scatter(
        ranked["bh_rp_au"] + jitter, ranked["Score"],
        c=n_ej, cmap=cmap_ord, vmin=-0.5, vmax=4.5,
        s=42, alpha=0.85, linewidths=0.5, edgecolors=_SURFACE, zorder=3,
    )
    ax1.set_xlabel("BH periapsis distance $r_p$ [AU]")
    ax1.set_ylabel("Impact score")
    ax1.set_title("Run impact vs. periapsis distance", color=_INK_PRIMARY,
                  fontsize=11, loc="left")
    legend_handles = [
        Patch(facecolor=_ORDINAL_BLUE[k], edgecolor="none",
              label=("4+" if k == 4 else str(k)))
        for k in range(5)
    ]
    leg = ax1.legend(
        handles=legend_handles, title="bodies ejected", loc="upper right",
        frameon=False, fontsize=8, title_fontsize=8.5,
    )
    leg.get_title().set_color(_INK_SECONDARY)
    for text in leg.get_texts():
        text.set_color(_INK_SECONDARY)

    # --- Panel B: max Score heatmap over (rp, inclination) ---
    _style_axes(ax2)
    pivot = ranked.pivot_table(
        index="bh_inc_deg", columns="bh_rp_au", values="Score", aggfunc="max"
    ).sort_index(ascending=True)
    cmap_seq = LinearSegmentedColormap.from_list(
        "sequential_blue", _SEQUENTIAL_BLUE_STEPS, N=256
    )
    im = ax2.imshow(
        pivot.values, aspect="auto", origin="lower", cmap=cmap_seq,
        extent=[
            pivot.columns.min() - 0.1, pivot.columns.max() + 0.1,
            pivot.index.min() - 5, pivot.index.max() + 5,
        ],
    )
    ax2.set_xticks(sorted(pivot.columns))
    ax2.set_yticks(sorted(pivot.index))
    ax2.set_xlabel("BH periapsis distance $r_p$ [AU]")
    ax2.set_ylabel("Inclination $i$ [deg]")
    ax2.set_title(
        "Worst-case score by ($r_p$, $i$)\n(max over $\\Omega,\\omega$)",
        color=_INK_PRIMARY, fontsize=11, loc="left",
    )
    ax2.grid(False)
    cbar = fig.colorbar(im, ax=ax2, fraction=0.046, pad=0.04)
    cbar.set_label("Max impact score", color=_INK_SECONDARY, fontsize=9)
    cbar.ax.tick_params(colors=_INK_MUTED, labelsize=8)
    cbar.outline.set_edgecolor(_AXIS)

    fig.tight_layout()
    fig.savefig(out_path, dpi=150, facecolor=_SURFACE)
    plt.close(fig)


def main() -> int:
    p = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    p.add_argument("parent_dir", type=Path,
                   help="Sweep parent folder, e.g. simulations/20260724_230314")
    p.add_argument("--out", type=Path, default=None,
                   help="Write the full ranked table here (CSV). "
                        "Default: <parent_dir>_impact_ranking.csv next to parent_dir.")
    p.add_argument("--top", type=int, default=20, help="Rows to print to console.")
    p.add_argument("--earth-weight", type=float, default=5.0)
    p.add_argument("--moon-weight", type=float, default=0.0)
    p.add_argument("--r0", type=float, default=None,
                   help="Override the noise-floor r0 instead of using the sweep median.")
    p.add_argument("--plot", type=Path, default=None,
                   help="Write a 2-panel impact plot (Score vs rp; rp x inc heatmap) here.")
    args = p.parse_args()

    weights = dict(DEFAULT_WEIGHTS)
    weights["Earth"] = args.earth_weight
    weights["Moon"] = args.moon_weight

    ranked = rank_runs(args.parent_dir, weights=weights, r0=args.r0)

    out_path = args.out or args.parent_dir.parent / f"{args.parent_dir.name}_impact_ranking.csv"
    ranked.to_csv(out_path, index=False)

    print(f"Ranked {len(ranked)} runs from {args.parent_dir}")
    print(f"Noise floor r0 = {ranked.attrs['r0']:.4g}")
    print(f"Wrote full ranking to {out_path}\n")

    with pd.option_context("display.width", 160, "display.max_colwidth", 40):
        cols = ["run", "Score", "Earth_s", "Earth_de", "n_ejected", "ejected_bodies",
                "top_impact_body"]
        print(ranked[cols].head(args.top).to_string(index=False))

    if args.plot is not None:
        make_impact_plot(ranked, args.plot)
        print(f"\nWrote impact plot to {args.plot}")

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
