"""Recover Earth's full post-flyby orbital elements from a BH sweep.

The per-run ``*__planets_run_deltas.csv`` files carry only a, e and q -- the
orbit's *size and shape*. Driving the seasonal climate model also needs its
*orientation*: the longitude of perihelion (lambda_p) and the obliquity. Both
are recoverable from the state vectors logged in ``*__orbits__*.xlsx``.

Why this works
--------------
The REBOUND simulation runs in the **equatorial J2000 frame** (verified: Earth's
orbit is inclined 23.46 deg to the z-axis at t=0, i.e. exactly the obliquity), so
the frame's z-axis *is* Earth's spin axis. A point-mass integration applies no
torque to that axis, so it stays fixed in inertial space while the BH tilts the
*orbital plane*. Hence:

    obliquity   eps    = angle between the spin axis (z) and the orbit normal h
    equinox dir  e_eq  = h_hat x s_perp        (s_perp = spin axis projected
                                                into the orbital plane)
    lambda_p           = angle from e_eq to the eccentricity vector, about h

The equinox convention follows the insolation module: declination obeys
sin(delta) = sin(eps) sin(lambda) with lambda = 0 at the northern vernal
equinox, which places lambda = 0 at 90 deg from s_perp along the direction of
orbital motion.

Validation: applied to t=0 this recovers lambda_p = 282.3 deg at the 1873 epoch
(present-day value is ~283 deg) and reproduces the CSV's a_before / e_before to
five decimals.

Input formats
-------------
Reads a run's Parquet trajectory log if it has one -- native output since engine
v26, or generated from an older sweep by ``convert_orbits_to_parquet.py`` -- and
falls back to the legacy ``*__orbits__*.xlsx`` workbook otherwise.

The workbook path is slow by nature: ~220 MB each, of which only the first and
last rows are needed, so the sheet XML is stream-decompressed straight out of the
xlsx zip and only the head/tail bytes are kept -- ~30x faster than a full
openpyxl scan, and still far slower than reading two columns of Parquet.

Usage
-----
    python extract_earth_elements.py simulations/20260724_230314
    python extract_earth_elements.py simulations/20260724_230314 --workers 5
"""

from __future__ import annotations

import argparse
import re
import zipfile
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path

import numpy as np
import pandas as pd

GM_SUN = 2.959122082855911e-4      # AU^3 / day^2
SPIN_AXIS = np.array([0.0, 0.0, 1.0])   # equatorial J2000 -> z is Earth's spin axis
DAYS_PER_YEAR = 365.256363

_ROW_RE = re.compile(rb"<row[^>]*>.*?</row>", re.S)
_CELL_RE = re.compile(rb'<c r="([A-Z]+)\d+"[^>]*>\s*<v>([^<]*)</v>', re.S)

_NS_REL = "{http://schemas.openxmlformats.org/officeDocument/2006/relationships}"

_COLS = ["A", "B", "C", "D", "E", "F", "G"]   # time, x, y, z, vx, vy, vz


def _sheet_map(z: zipfile.ZipFile) -> dict[str, str]:
    """Map sheet display name -> zip member path, via the workbook rels.

    Parsed as XML rather than by regex: attribute order is not guaranteed by the
    spec (this writer emits Target before Id), so positional patterns break.
    """
    import xml.etree.ElementTree as ET

    rels = ET.fromstring(z.read("xl/_rels/workbook.xml.rels"))
    rid_to_target = {
        rel.get("Id"): rel.get("Target")
        for rel in rels
        if rel.get("Id") and rel.get("Target")
    }

    wb = ET.fromstring(z.read("xl/workbook.xml"))
    out = {}
    for sheet in wb.iter():
        if not sheet.tag.endswith("}sheet") and sheet.tag != "sheet":
            continue
        name = sheet.get("name")
        rid = sheet.get(f"{_NS_REL}id") or sheet.get("id")
        target = rid_to_target.get(rid)
        if name and target:
            out[name] = "xl/" + target.lstrip("/").removeprefix("xl/")
    return out


def _head_tail_rows(z: zipfile.ZipFile, member: str, tail_bytes: int = 65536):
    """Stream a sheet, returning (first_row_values, last_row_values) as float lists."""
    head = b""
    tail = b""
    with z.open(member) as fh:
        while True:
            chunk = fh.read(1 << 20)
            if not chunk:
                break
            if len(head) < tail_bytes:
                head += chunk
            tail = (tail + chunk)[-tail_bytes:]

    def parse(buf: bytes, which: str):
        rows = _ROW_RE.findall(buf)
        if not rows:
            return None
        # Skip the header row (row 1) when taking the first data row.
        for raw in (rows if which == "first" else reversed(rows)):
            cells = dict(_CELL_RE.findall(raw))
            try:
                vals = [float(cells[c.encode()]) for c in _COLS]
            except (KeyError, ValueError):
                continue
            return vals
        return None

    return parse(head, "first"), parse(tail, "last")


def elements_from_state(r: np.ndarray, v: np.ndarray) -> dict:
    """Heliocentric a, e, obliquity and lambda_p from a state vector."""
    h = np.cross(r, v)
    h_norm = np.linalg.norm(h)
    r_norm = np.linalg.norm(r)
    if h_norm == 0 or r_norm == 0:
        return {"a_au": np.nan, "ecc": np.nan,
                "obliquity_deg": np.nan, "lon_perihelion_deg": np.nan}

    h_hat = h / h_norm
    e_vec = np.cross(v, h) / GM_SUN - r / r_norm
    ecc = float(np.linalg.norm(e_vec))

    inv_a = 2.0 / r_norm - v.dot(v) / GM_SUN
    a = float(1.0 / inv_a) if abs(inv_a) > 1e-30 else np.inf

    # Obliquity: angle between the (fixed) spin axis and the orbit normal.
    cos_eps = float(np.clip(h_hat.dot(SPIN_AXIS), -1.0, 1.0))
    eps = float(np.degrees(np.arccos(cos_eps)))

    # Vernal-equinox direction in the orbital plane, then lambda_p to perihelion.
    s_perp = SPIN_AXIS - cos_eps * h_hat
    s_norm = np.linalg.norm(s_perp)
    if s_norm < 1e-12 or ecc < 1e-12:
        lam_p = np.nan          # degenerate: zero obliquity or circular orbit
    else:
        e_eq = np.cross(h_hat, s_perp / s_norm)
        e_peri = e_vec / ecc
        lam_p = float(np.degrees(np.arctan2(
            np.cross(e_eq, e_peri).dot(h_hat), e_eq.dot(e_peri))) % 360.0)

    return {"a_au": a, "ecc": ecc, "obliquity_deg": eps, "lon_perihelion_deg": lam_p}


def _head_tail_parquet(path: Path):
    """(first, last) [t,x,y,z,vx,vy,vz] rows for Earth and Sun from a Parquet log.

    Since engine v26 the trajectory log is written as long-format Parquet
    (``body, t_days, x_au, ...``) rather than a workbook; ``convert_orbits_to_-
    parquet.py`` produces the identical schema from older sweeps. Only the
    endpoints are needed, so the frame is filtered to the two bodies first.
    """
    d = pd.read_parquet(
        path, columns=["body", "t_days", "x_au", "y_au", "z_au", "vx", "vy", "vz"])
    out = {}
    for name in ("Earth", "Sun"):
        g = d[d["body"] == name].sort_values("t_days")
        if g.empty:
            return None
        v = g[["t_days", "x_au", "y_au", "z_au", "vx", "vy", "vz"]].to_numpy(dtype=float)
        out[name] = (v[0].tolist(), v[-1].tolist())
    return out


def process_run(run_dir: Path) -> dict | None:
    """Extract Earth's before/after elements for one run folder.

    Reads whichever trajectory log the run has: Parquet if present (native since
    engine v26, or produced by the converter), otherwise the legacy workbook.
    """
    pq = list(run_dir.glob("*orbits*.parquet")) or list(run_dir.glob("orbits.parquet"))
    xl = list(run_dir.glob("*__orbits__*.xlsx"))
    if not pq and not xl:
        return None
    try:
        if pq:
            got = _head_tail_parquet(pq[0])
            if got is None:
                return None
            (e_first, e_last), (s_first, s_last) = got["Earth"], got["Sun"]
        else:
            with zipfile.ZipFile(xl[0]) as z:
                smap = _sheet_map(z)
                if "Earth" not in smap or "Sun" not in smap:
                    return None
                e_first, e_last = _head_tail_rows(z, smap["Earth"])
                s_first, s_last = _head_tail_rows(z, smap["Sun"])
    except Exception as exc:                      # noqa: BLE001 - report, don't abort sweep
        return {"run": run_dir.name, "error": f"{type(exc).__name__}: {exc}"}

    if not all((e_first, e_last, s_first, s_last)):
        return {"run": run_dir.name, "error": "could not parse head/tail rows"}

    ef, el = np.array(e_first), np.array(e_last)
    sf, sl = np.array(s_first), np.array(s_last)

    before = elements_from_state(ef[1:4] - sf[1:4], ef[4:7] - sf[4:7])
    after = elements_from_state(el[1:4] - sl[1:4], el[4:7] - sl[4:7])

    row = {"run": run_dir.name, "t_end_days": float(el[0])}
    row.update({f"{k}_before": v for k, v in before.items()})
    row.update({f"{k}_after": v for k, v in after.items()})
    # Kepler's third law: the year length changes with the new semi-major axis.
    a_after = after["a_au"]
    row["days_per_year_after"] = (
        DAYS_PER_YEAR * a_after ** 1.5 if np.isfinite(a_after) and a_after > 0 else np.nan
    )
    row["bound_after"] = bool(np.isfinite(a_after) and a_after > 0 and after["ecc"] < 1.0)
    return row


def main() -> int:
    p = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    p.add_argument("parent_dir", type=Path, help="Sweep folder, e.g. simulations/20260724_230314")
    p.add_argument("--out", type=Path, default=None, help="Output CSV path.")
    p.add_argument("--workers", type=int, default=1, help="Worker processes.")
    p.add_argument("--limit", type=int, default=None, help="Process only the first N runs.")
    args = p.parse_args()

    run_dirs = sorted(d for d in args.parent_dir.iterdir() if d.is_dir())
    if args.limit:
        run_dirs = run_dirs[: args.limit]
    if not run_dirs:
        raise SystemExit(f"No run subfolders under {args.parent_dir}")

    rows = []
    if args.workers <= 1:
        for i, d in enumerate(run_dirs, 1):
            r = process_run(d)
            if r:
                rows.append(r)
            if i % 25 == 0 or i == len(run_dirs):
                print(f"  {i}/{len(run_dirs)} runs processed", flush=True)
    else:
        with ProcessPoolExecutor(max_workers=args.workers) as ex:
            futs = {ex.submit(process_run, d): d for d in run_dirs}
            for i, fut in enumerate(as_completed(futs), 1):
                r = fut.result()
                if r:
                    rows.append(r)
                if i % 25 == 0 or i == len(run_dirs):
                    print(f"  {i}/{len(run_dirs)} runs processed", flush=True)

    df = pd.DataFrame(rows).sort_values("run").reset_index(drop=True)
    out = args.out or args.parent_dir.parent / f"{args.parent_dir.name}_earth_elements.csv"
    df.to_csv(out, index=False)

    n_err = df["error"].notna().sum() if "error" in df.columns else 0
    print(f"\nProcessed {len(df)} runs ({n_err} errors). Wrote {out}")
    if "bound_after" in df.columns:
        print(f"Earth bound after flyby: {int(df['bound_after'].sum())} / {len(df)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
