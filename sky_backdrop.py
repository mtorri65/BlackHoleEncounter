"""Draw a star field and constellation figures behind a sky chart.

Reads the compact CSVs built by ``fetch_constellation_data.py``. Run that once
first; if the files are missing, ``draw_backdrop`` degrades to name labels only
rather than failing, so the charts still render on a machine that has never
fetched them.

Everything here is deliberately recessive. The backdrop exists so the viewer can
place the black hole against something familiar -- it is context, not data, and
must never compete with the track drawn over it. Hence dim blue-greys, thin
lines, and a magnitude cut well brighter than the catalogue limit: at mag 6.5
the 8,870 stars turn the panel into noise, while ~1,600 at mag 5.0 read as
constellations at a glance.

Both the star and figure coordinates are J2000, matching the simulation frame,
so no precession is applied. See ``fetch_constellation_data.py``.
"""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent
STARS_CSV = HERE / "sky_stars.csv"
LINES_CSV = HERE / "sky_constellation_lines.csv"

STAR_COLOR = "#b9c9e8"
FIG_COLOR = "#37527f"
NAME_COLOR = "#c9b96a"

# IAU abbreviation -> name, for the labels. Genitive forms are not used.
NAMES = {
    "And": "Andromeda", "Ant": "Antlia", "Aps": "Apus", "Aql": "Aquila",
    "Aqr": "Aquarius", "Ara": "Ara", "Ari": "Aries", "Aur": "Auriga",
    "Boo": "Bootes", "CMa": "Canis Major", "CMi": "Canis Minor",
    "CVn": "Canes Venatici", "Cae": "Caelum", "Cam": "Camelopardalis",
    "Cap": "Capricornus", "Car": "Carina", "Cas": "Cassiopeia",
    "Cen": "Centaurus", "Cep": "Cepheus", "Cet": "Cetus", "Cha": "Chamaeleon",
    "Cir": "Circinus", "Cnc": "Cancer", "Col": "Columba", "Com": "Coma Ber.",
    "CrA": "Corona Aus.", "CrB": "Corona Bor.", "Crt": "Crater",
    "Cru": "Crux", "Crv": "Corvus", "Cyg": "Cygnus", "Del": "Delphinus",
    "Dor": "Dorado", "Dra": "Draco", "Equ": "Equuleus", "Eri": "Eridanus",
    "For": "Fornax", "Gem": "Gemini", "Gru": "Grus", "Her": "Hercules",
    "Hor": "Horologium", "Hya": "Hydra", "Hyi": "Hydrus", "Ind": "Indus",
    "LMi": "Leo Minor", "Lac": "Lacerta", "Leo": "Leo", "Lep": "Lepus",
    "Lib": "Libra", "Lup": "Lupus", "Lyn": "Lynx", "Lyr": "Lyra",
    "Men": "Mensa", "Mic": "Microscopium", "Mon": "Monoceros",
    "Mus": "Musca", "Nor": "Norma", "Oct": "Octans", "Oph": "Ophiuchus",
    "Ori": "Orion", "Pav": "Pavo", "Peg": "Pegasus", "Per": "Perseus",
    "Phe": "Phoenix", "Pic": "Pictor", "PsA": "Piscis Aus.",
    "Psc": "Pisces", "Pup": "Puppis", "Pyx": "Pyxis", "Ret": "Reticulum",
    "Scl": "Sculptor", "Sco": "Scorpius", "Sct": "Scutum",
    "Ser": "Serpens", "Sex": "Sextans", "Sge": "Sagitta", "Sgr": "Sagittarius",
    "Tau": "Taurus", "Tel": "Telescopium", "TrA": "Triangulum Aus.",
    "Tri": "Triangulum", "Tuc": "Tucana", "UMa": "Ursa Major",
    "UMi": "Ursa Minor", "Vel": "Vela", "Vir": "Virgo", "Vol": "Volans",
    "Vul": "Vulpecula",
}


def load() -> tuple[pd.DataFrame | None, pd.DataFrame | None]:
    stars = pd.read_csv(STARS_CSV) if STARS_CSV.exists() else None
    lines = pd.read_csv(LINES_CSV) if LINES_CSV.exists() else None
    return stars, lines


def _unwrap(x1: np.ndarray, x2: np.ndarray) -> np.ndarray:
    """Shift x1 so each segment takes the short way round the 24h wrap."""
    d = x2 - x1
    return x1 + 24.0 * (d > 12.0) - 24.0 * (d < -12.0)


def star_sizes(vmag: np.ndarray, faintest: float) -> np.ndarray:
    """Marker areas that fall off steeply enough to show the bright pattern.

    Linear-in-magnitude sizing makes every star look alike; true flux scaling
    (2.512**-mag) makes Sirius a blob and everything else invisible. This sits
    between the two.
    """
    return 0.75 * np.clip(faintest - vmag + 0.45, 0.05, None) ** 1.9


def draw_backdrop(ax, ra_to_x, *, max_mag: float = 5.0, labels: bool = True,
                  label_size: float = 8.5, alpha: float = 1.0,
                  x_period: float = 24.0) -> bool:
    """Draw stars, constellation figures and names onto ``ax``.

    ``ra_to_x`` maps right ascension in hours to the chart's x coordinate, so
    the backdrop follows whatever seam and orientation the caller chose.
    ``x_period`` is the axis width in the same units, used to duplicate
    segments that straddle the seam. Returns False if the data files are absent.
    """
    stars, lines = load()

    if lines is not None:
        x1 = np.asarray(ra_to_x(lines.ra1.values / 15.0), dtype=float)
        x2 = np.asarray(ra_to_x(lines.ra2.values / 15.0), dtype=float)
        y1, y2 = lines.dec1.values, lines.dec2.values
        x1 = _unwrap(x1, x2)
        segs = np.stack([np.column_stack([x1, y1]),
                         np.column_stack([x2, y2])], axis=1)
        # A segment pulled outside the axis by the unwrap still has to appear at
        # the opposite edge, so draw a shifted copy too and let clipping decide.
        off = segs.copy()
        off[:, :, 0] += np.where(segs[:, :, 0].min(axis=1) < 0,
                                 x_period, -x_period)[:, None]
        from matplotlib.collections import LineCollection
        for s in (segs, off):
            ax.add_collection(LineCollection(s, colors=FIG_COLOR, linewidths=0.7,
                                             alpha=0.55 * alpha, zorder=2))

    if stars is not None:
        s = stars[stars.vmag <= max_mag]
        ax.scatter(ra_to_x(s.ra_deg.values / 15.0), s.dec_deg.values,
                   s=star_sizes(s.vmag.values, max_mag), c=STAR_COLOR,
                   alpha=0.8 * alpha, linewidths=0, zorder=3)

    if labels and lines is not None:
        # Stellarium capitalises a couple of abbreviations its own way ("Cvn",
        # "Tra" against the IAU's "CVn", "TrA"), so match case-insensitively.
        by_lower = {k.lower(): v for k, v in NAMES.items()}
        y0, y1 = sorted(ax.get_ylim())
        x0, x1 = sorted(ax.get_xlim())
        for abbrev, g in lines.groupby("const"):
            # Mean of unit vectors, so constellations straddling the seam get a
            # centroid in the right place instead of halfway across the sky.
            ra = np.radians(np.concatenate([g.ra1.values, g.ra2.values]))
            dec = np.radians(np.concatenate([g.dec1.values, g.dec2.values]))
            v = np.array([np.cos(dec) * np.cos(ra), np.cos(dec) * np.sin(ra),
                          np.sin(dec)]).mean(axis=1)
            n = np.linalg.norm(v)
            if n < 1e-9:
                continue
            v /= n
            cra = np.degrees(np.arctan2(v[1], v[0])) % 360.0
            cdec = np.degrees(np.arcsin(np.clip(v[2], -1, 1)))
            cx = float(ra_to_x(cra / 15.0))
            # Matplotlib does not clip text by default, so a constellation whose
            # centre lies off the chart would otherwise print its name over the
            # axis labels and the caption.
            if not (y0 <= cdec <= y1 and x0 <= cx <= x1):
                continue
            # Near a vertical edge a centred label would be sliced through the
            # middle of a word, so anchor it inward instead.
            margin = 0.06 * (x1 - x0)
            ha = ("left" if cx < x0 + margin else
                  "right" if cx > x1 - margin else "center")
            ax.text(cx, cdec, by_lower.get(abbrev.lower(), abbrev),
                    color=NAME_COLOR, fontsize=label_size, ha=ha,
                    va="center", zorder=4, alpha=0.75 * alpha, clip_on=True)

    return stars is not None or lines is not None
