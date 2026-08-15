"""Build the compact star/constellation files the sky-track charts draw as backdrop.

Run once. It downloads the Hipparcos catalogue (~53 MB) and a Stellarium
constellation-figure file, then writes two small CSVs next to this script:

    sky_stars.csv             bright stars:  ra_deg, dec_deg, vmag
    sky_constellation_lines.csv   figure segments: const, ra1,dec1, ra2,dec2

The big catalogue is only ever a build input -- it goes to a temp directory and
is not kept, because the distilled files are ~150 KB together and that is all
the charts need. Re-run with --keep-raw to hold on to it.

Coordinates
-----------
Hipparcos positions are ICRS, which for plotting purposes is J2000. The
simulation runs in equatorial J2000 too, so stars and black-hole track share a
frame and no precession is applied anywhere. At chart scale the catalogue's
J1991.25 proper-motion epoch is irrelevant -- a fast star moves a few
arcseconds, far under a pixel.

Constellation figures
---------------------
From Stellarium's ``western`` sky culture, pinned to tag v0.21.3. The figures
are not standardised -- the IAU defines constellation *boundaries*, not the
stick figures -- so this is one conventional choice among several. The file is
pinned rather than tracked to master because Stellarium reorganised its
skyculture layout afterwards and the current path 404s.

Usage
-----
    python fetch_constellation_data.py
    python fetch_constellation_data.py --max-mag 6.0
"""

from __future__ import annotations

import argparse
import ssl
import tempfile
import urllib.request
from pathlib import Path

import numpy as np
import pandas as pd

HIP_URL = "https://cdsarc.cds.unistra.fr/ftp/cats/I/239/hip_main.dat"
FAB_URL = ("https://raw.githubusercontent.com/Stellarium/stellarium/"
           "v0.21.3/skycultures/western/constellationship.fab")

HERE = Path(__file__).resolve().parent
STARS_CSV = HERE / "sky_stars.csv"
LINES_CSV = HERE / "sky_constellation_lines.csv"


def _ctx() -> ssl.SSLContext:
    """CA bundle from certifi; the system store lacks the CDS issuer here."""
    try:
        import certifi
        return ssl.create_default_context(cafile=certifi.where())
    except ImportError:
        return ssl.create_default_context()


def download(url: str, dest: Path) -> Path:
    print(f"  fetching {url}")
    with urllib.request.urlopen(url, timeout=120, context=_ctx()) as r:
        dest.write_bytes(r.read())
    print(f"    -> {dest.name}  {dest.stat().st_size/1e6:.1f} MB")
    return dest


def main() -> int:
    p = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    p.add_argument("--max-mag", type=float, default=6.5,
                   help="Faintest star to keep [V mag]. 6.5 is naked-eye limit.")
    p.add_argument("--keep-raw", action="store_true",
                   help="Keep the downloaded Hipparcos catalogue.")
    args = p.parse_args()

    tmp = Path(args.keep_raw and HERE or tempfile.mkdtemp(prefix="skydata_"))
    hip_file = tmp / "hip_main.dat"
    fab_file = tmp / "constellationship.fab"

    print("Downloading source catalogues:")
    if not hip_file.exists():
        download(HIP_URL, hip_file)
    download(FAB_URL, fab_file)

    from skyfield.data import hipparcos, stellarium

    print("\nParsing Hipparcos...")
    with hip_file.open("rb") as f:
        df = hipparcos.load_dataframe(f)
    # load_dataframe already drops rows without a parallax/position solution.
    df = df[np.isfinite(df.ra_degrees) & np.isfinite(df.dec_degrees)]
    print(f"  {len(df):,} stars with positions")

    with fab_file.open("rb") as f:
        cons = stellarium.parse_constellations(f)
    print(f"  {len(cons)} constellation figures")

    # --- figure segments -----------------------------------------------------
    # Look each HIP number up in the catalogue. A few figure stars are missing
    # from the position-solved subset; those segments are dropped rather than
    # drawn to a wrong place.
    ra, dec = df.ra_degrees, df.dec_degrees
    rows, dropped = [], 0
    for abbrev, segments in cons:
        for h1, h2 in segments:
            if h1 not in ra.index or h2 not in ra.index:
                dropped += 1
                continue
            rows.append((abbrev, ra[h1], dec[h1], ra[h2], dec[h2]))
    lines = pd.DataFrame(rows, columns=["const", "ra1", "dec1", "ra2", "dec2"])
    lines.round(5).to_csv(LINES_CSV, index=False)
    print(f"\n  {len(lines)} segments -> {LINES_CSV.name}"
          + (f"  ({dropped} dropped, star not in catalogue)" if dropped else ""))

    # --- bright stars --------------------------------------------------------
    bright = df[df.magnitude <= args.max_mag]
    out = pd.DataFrame({
        "ra_deg": bright.ra_degrees.round(5),
        "dec_deg": bright.dec_degrees.round(5),
        "vmag": bright.magnitude.round(2),
    }).sort_values("vmag")
    out.to_csv(STARS_CSV, index=False)
    print(f"  {len(out)} stars <= mag {args.max_mag} -> {STARS_CSV.name}"
          f"  ({STARS_CSV.stat().st_size/1e3:.0f} KB)")

    if not args.keep_raw:
        hip_file.unlink(missing_ok=True)
        fab_file.unlink(missing_ok=True)
        print("\n  removed the raw catalogue (pass --keep-raw to retain it)")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
