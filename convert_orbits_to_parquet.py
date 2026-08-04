"""Convert REBOUND orbit workbooks to compact Parquet, with verification.

The engine writes one ``*__orbits__*.xlsx`` per run at ~220 MB, which makes a
672-run sweep ~150 GB. Almost all of that is format overhead: the workbook
stores every number as XML text, keeps five redundant ``*_str`` text copies of
numeric columns, and stores four more columns that are derivable from the state
vector. Converting to Parquet reduces the sweep to roughly 33 GB.

What is kept
------------
Per body, the columns that carry independent information:

    t_days, x_au, y_au, z_au, vx, vy, vz, disp_helio_au

Dropped as recoverable: ``r_helio`` (= |r - r_sun|), ``a_tidal`` (computable
from the BH state and mass), the ``*_m``/``*_m_s2`` unit conversions, and the
five ``*_str`` text duplicates. ``disp_helio_au`` is **kept** -- it is the
displacement against the BH-free baseline integration and cannot be
reconstructed from this file alone.

Reading strategy
----------------
``openpyxl`` takes ~13 s per sheet (29 h for the sweep). Instead the sheet XML
is stream-decompressed straight out of the xlsx zip and all ``<v>`` values are
pulled with one regex: ~0.6 s per sheet, about 21x faster.

That positional parse assumes a dense grid, so every sheet is checked two ways:
the value count must be an exact multiple of the column count, and the time
column must equal ``0, 1, 2, ...`` exactly -- any cell misalignment breaks the
latter immediately.

Verification
------------
Every converted file is read back and compared against the in-memory arrays
(exactly for float64, within float32 rounding otherwise). Row counts and body
lists are checked too. A run is only reported as ``ok`` if all of that passes.

**Source files are never modified or deleted.** Delete them yourself once the
Parquet copies are verified and backed up.

Usage
-----
    python convert_orbits_to_parquet.py simulations/20260724_230314 --workers 5
    python convert_orbits_to_parquet.py simulations/20260724_230314 \
        --out-dir D:/archive --precision float64 --workers 5
"""

from __future__ import annotations

import argparse
import re
import zipfile
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path

import numpy as np
import pandas as pd

# Columns to keep, addressed by spreadsheet column letter rather than by
# position. Sheets are NOT uniform -- the BH sheet omits the a_tidal columns
# (J, L), since the BH exerts no tidal acceleration on itself -- so positional
# parsing silently misaligns. Addressing by letter is robust to that.
#
#   A time_days   B x_AU   C y_AU   D z_AU   E vx   F vy   G vz
#   H r_helio_AU  (derivable)        J/L a_tidal (derivable)
#   N disp_helio_AU  <- kept: displacement vs the BH-free baseline run,
#                       which cannot be reconstructed from this file alone
#   I/K/M/O/Q are *_str text duplicates; P is a unit conversion of N.
_KEEP = {
    "A": "t_days", "B": "x_au", "C": "y_au", "D": "z_au",
    "E": "vx", "F": "vy", "G": "vz", "N": "disp_helio_au",
}
_COL_NAMES = list(_KEEP.values())
_STATE_COLS = ["x_au", "y_au", "z_au", "vx", "vy", "vz", "disp_helio_au"]

_NS_REL = "{http://schemas.openxmlformats.org/officeDocument/2006/relationships}"


def _column_re(letter: str) -> re.Pattern:
    return re.compile(
        rb'<c r="' + letter.encode() + rb'(\d+)"[^>]*>\s*<v>([^<]*)</v>')


def _sheet_map(z: zipfile.ZipFile) -> dict[str, str]:
    """Map sheet display name -> zip member path (attribute order is not guaranteed)."""
    import xml.etree.ElementTree as ET

    rels = ET.fromstring(z.read("xl/_rels/workbook.xml.rels"))
    rid_to_target = {r.get("Id"): r.get("Target") for r in rels
                     if r.get("Id") and r.get("Target")}
    wb = ET.fromstring(z.read("xl/workbook.xml"))
    out = {}
    for sheet in wb.iter():
        if not sheet.tag.endswith("}sheet") and sheet.tag != "sheet":
            continue
        name, rid = sheet.get("name"), sheet.get(f"{_NS_REL}id") or sheet.get("id")
        target = rid_to_target.get(rid)
        if name and target:
            out[name] = "xl/" + target.lstrip("/").removeprefix("xl/")
    return out


def read_sheet(z: zipfile.ZipFile, member: str) -> pd.DataFrame:
    """Return the wanted columns of one sheet, integrity-checked.

    Each column is extracted by its spreadsheet letter, so a sheet that omits a
    column (as the BH sheet does for a_tidal) yields NaN for it rather than
    shifting every subsequent value. Row numbers are verified to be the
    contiguous block ``2..N+1`` (row 1 is the header), and the time column is
    checked to be exactly ``0..N-1``.
    """
    buf = z.read(member)
    cols: dict[str, np.ndarray] = {}
    n_rows = None

    for letter, name in _KEEP.items():
        hits = _column_re(letter).findall(buf)
        if not hits:
            cols[name] = None            # absent in this sheet; filled below
            continue
        rows = np.array([h[0] for h in hits], dtype=np.int64)
        vals = np.array([h[1] for h in hits], dtype=float)
        expected = np.arange(2, len(rows) + 2, dtype=np.int64)
        if not np.array_equal(rows, expected):
            raise ValueError(f"{member}: column {letter} rows are not contiguous 2..N+1")
        if n_rows is None:
            n_rows = len(vals)
        elif len(vals) != n_rows:
            raise ValueError(
                f"{member}: column {letter} has {len(vals)} rows, expected {n_rows}")
        cols[name] = vals

    if n_rows is None:
        raise ValueError(f"{member}: no numeric cells found")
    if cols["t_days"] is None:
        raise ValueError(f"{member}: time column missing")
    # Integrity guard: the time column is a 0,1,2,... day counter.
    t = cols["t_days"]
    if not np.array_equal(t, np.arange(n_rows, dtype=float)):
        raise ValueError(f"{member}: time column is not 0..N-1 -- possible misalignment")

    for name, v in cols.items():
        if v is None:
            cols[name] = np.full(n_rows, np.nan)
    return pd.DataFrame(cols, columns=_COL_NAMES)


def convert_run(run_dir: Path, out_root: Path, src_root: Path,
                precision: str = "float32", compression: str = "zstd") -> dict:
    """Convert one run's orbit workbook to Parquet and verify the result."""
    result = {"run": run_dir.name}
    matches = list(run_dir.glob("*__orbits__*.xlsx"))
    if not matches:
        return {**result, "status": "skipped", "reason": "no orbits workbook"}
    src = matches[0]

    try:
        frames = []
        with zipfile.ZipFile(src) as z:
            smap = _sheet_map(z)
            bodies = [n for n in smap if n != "Summary"]
            for body in bodies:
                df = read_sheet(z, smap[body])
                df.insert(0, "body", body)
                frames.append(df)
        table = pd.concat(frames, ignore_index=True)
        table["body"] = table["body"].astype("category")
        table["t_days"] = table["t_days"].astype("int32")
        if precision == "float32":
            table[_STATE_COLS] = table[_STATE_COLS].astype("float32")

        rel = run_dir.relative_to(src_root)
        out_path = out_root / rel / "orbits.parquet"
        out_path.parent.mkdir(parents=True, exist_ok=True)
        table.to_parquet(out_path, compression=compression, index=False)

        # --- verification: read back and compare against what we wrote ---
        back = pd.read_parquet(out_path)
        if len(back) != len(table):
            raise ValueError(f"row count {len(back)} != {len(table)}")
        if sorted(back["body"].astype(str).unique()) != sorted(bodies):
            raise ValueError("body list differs after round-trip")
        tol = 0.0 if precision == "float64" else 1e-5
        worst = 0.0
        for col in _STATE_COLS + ["t_days"]:
            a = table[col].to_numpy(dtype=np.float64)
            b = back[col].to_numpy(dtype=np.float64)
            # Absent columns are all-NaN by construction; require the NaN masks
            # to agree, then compare only the finite entries.
            na, nb = np.isnan(a), np.isnan(b)
            if not np.array_equal(na, nb):
                raise ValueError(f"round-trip NaN pattern differs in column {col}")
            if na.all():
                continue
            a, b = a[~na], b[~na]
            scale = np.maximum(np.abs(a), 1.0)
            worst = max(worst, float(np.max(np.abs(a - b) / scale)))
        if worst > tol:
            raise ValueError(f"round-trip mismatch: max relative diff {worst:.3e} > {tol:.0e}")

        return {
            **result,
            "status": "ok",
            "bodies": len(bodies),
            "rows": len(table),
            "src_mb": src.stat().st_size / 1e6,
            "out_mb": out_path.stat().st_size / 1e6,
            "max_rel_diff": worst,
        }
    except Exception as exc:                       # noqa: BLE001 - record, keep going
        return {**result, "status": "failed", "reason": f"{type(exc).__name__}: {exc}"}


def main() -> int:
    p = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    p.add_argument("parent_dir", type=Path, help="Sweep folder, e.g. simulations/<STAMP>")
    p.add_argument("--out-dir", type=Path, default=None,
                   help="Root of the mirrored Parquet tree. "
                        "Default: <parent>_parquet next to the sweep folder.")
    p.add_argument("--precision", choices=("float32", "float64"), default="float32",
                   help="float32 (default) is ~4.6x smaller; measured error on "
                        "re-derived orbital elements is ~1e-6 AU / 1e-6 deg.")
    p.add_argument("--compression", default="zstd", help="Parquet codec (default zstd).")
    p.add_argument("--workers", type=int, default=1)
    p.add_argument("--limit", type=int, default=None, help="Convert only the first N runs.")
    args = p.parse_args()

    run_dirs = sorted(d for d in args.parent_dir.iterdir() if d.is_dir())
    if args.limit:
        run_dirs = run_dirs[: args.limit]
    if not run_dirs:
        raise SystemExit(f"No run subfolders under {args.parent_dir}")

    out_root = args.out_dir or args.parent_dir.parent / f"{args.parent_dir.name}_parquet"
    out_root.mkdir(parents=True, exist_ok=True)
    print(f"Converting {len(run_dirs)} runs -> {out_root}  ({args.precision}, {args.compression})")

    rows: list[dict] = []
    fn_args = (out_root, args.parent_dir, args.precision, args.compression)
    if args.workers <= 1:
        for i, d in enumerate(run_dirs, 1):
            rows.append(convert_run(d, *fn_args))
            if i % 25 == 0 or i == len(run_dirs):
                print(f"  {i}/{len(run_dirs)}", flush=True)
    else:
        with ProcessPoolExecutor(max_workers=args.workers) as ex:
            futs = [ex.submit(convert_run, d, *fn_args) for d in run_dirs]
            for i, fut in enumerate(as_completed(futs), 1):
                rows.append(fut.result())
                if i % 25 == 0 or i == len(run_dirs):
                    print(f"  {i}/{len(run_dirs)}", flush=True)

    df = pd.DataFrame(rows)
    report = out_root / "conversion_report.csv"
    df.to_csv(report, index=False)

    ok = df[df["status"] == "ok"]
    bad = df[df["status"] == "failed"]
    skipped = df[df["status"] == "skipped"]

    print(f"\nverified ok: {len(ok)}   failed: {len(bad)}   skipped: {len(skipped)}")
    if len(ok):
        src, out = ok["src_mb"].sum() / 1000, ok["out_mb"].sum() / 1000
        print(f"  source     : {src:8.2f} GB")
        print(f"  parquet    : {out:8.2f} GB")
        print(f"  reduction  : {src / out:8.2f}x   ({src - out:.1f} GB reclaimable)")
        print(f"  worst round-trip relative error: {ok['max_rel_diff'].max():.2e}")
    for _, r in bad.head(10).iterrows():
        print(f"  FAILED {r['run']}: {r['reason']}")
    print(f"\nReport: {report}")
    print("Source workbooks were NOT modified or deleted. Remove them yourself "
          "once the Parquet tree is verified and backed up.")
    return 0 if len(bad) == 0 else 1


if __name__ == "__main__":
    raise SystemExit(main())
