# -*- coding: utf-8 -*-
"""
Count rows in *__planets_run_deltas.csv files under a simulation parent folder.

Subfolders are processed in *chronological creation order* (not alphabetical).
Chronological order is taken from the subfolder's filesystem creation time,
with a fallback to modification time if creation time is unavailable.

Writes:
  <PARENT_NAME>__planets_run_row_counts_chrono.csv
"""

import argparse
from pathlib import Path
from datetime import datetime
import pandas as pd


def get_created_ts(p: Path) -> float:
    """
    Return a best-effort 'creation time' timestamp for sorting.

    On Windows, st_ctime is creation time.
    On some Unix filesystems, st_ctime is metadata-change time, not true creation.
    For robustness, we use:
      - st_ctime (Windows creation time)
      - fallback to st_mtime if needed (modification time)
    """
    st = p.stat()
    # On Windows, this is creation time. Else it's still a usable ordering signal.
    t = getattr(st, "st_ctime", None)
    if t is None:
        t = st.st_mtime
    return float(t)


def fmt_ts(ts: float) -> str:
    return datetime.fromtimestamp(ts).strftime("%Y-%m-%d %H:%M:%S")


def process_parent_folder(parent: Path) -> Path:
    parent = parent.resolve()
    if not parent.is_dir():
        raise FileNotFoundError(f"Parent folder not found: {parent}")

    subfolders = [p for p in parent.iterdir() if p.is_dir()]
    # Sort by creation time (chronological)
    subfolders.sort(key=get_created_ts)

    rows = []

    for sub in subfolders:
        created = get_created_ts(sub)
        modified = float(sub.stat().st_mtime)

        matches = list(sub.glob("*__planets_run_deltas.csv"))
        if not matches:
            # still record that it had no file (optional); comment out if you prefer skipping
            rows.append({
                "subfolder": sub.name,
                "created_time": fmt_ts(created),
                "modified_time": fmt_ts(modified),
                "csv_file": "",
                "row_count": None,
                "status": "no_csv_found",
            })
            continue

        csv_path = matches[0]

        try:
            df = pd.read_csv(csv_path)
            n_rows = len(df)
            status = "ok"
        except Exception as e:
            print(f"[WARN] Could not read {csv_path}: {e}")
            n_rows = None
            status = "read_error"

        rows.append({
            "subfolder": sub.name,
            "created_time": fmt_ts(created),
            "modified_time": fmt_ts(modified),
            "csv_file": csv_path.name,
            "row_count": n_rows,
            "status": status,
        })

    if not rows:
        raise RuntimeError(f"No subfolders found under {parent}")

    out_df = pd.DataFrame(rows)
    out_path = parent / f"{parent.name}__planets_run_row_counts_chrono.csv"
    out_df.to_csv(out_path, index=False)

    print(f"[OK] Wrote {out_path}")
    return out_path


def main():
    parser = argparse.ArgumentParser(
        description="Count rows in *__planets_run_deltas.csv files in chronological subfolder creation order."
    )
    parser.add_argument(
        "parent_folder",
        help="Simulation parent folder (e.g. simulations/20251127_201847)"
    )
    args = parser.parse_args()

    process_parent_folder(Path(args.parent_folder))


if __name__ == "__main__":
    main()
