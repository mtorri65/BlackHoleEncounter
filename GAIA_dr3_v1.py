import numpy as np
import pandas as pd
import astropy.units as u
from astropy.time import Time
from astropy.coordinates import SkyCoord
from astroquery.gaia import Gaia

from pathlib import Path
import re

def resolve_bh_track_path(bh_track_filename: str,
                          project_root: Path | None = None) -> Path:
    """
    Resolve the full path to a BH RA/Dec Excel file given only its filename.

    Expected layout:
      <project_root>/simulations/<run_date>/<run_prefix>/<bh_track_filename>

    where:
      run_prefix is the filename without the trailing "__bh_radec__...xlsx"
      run_date is the first token of run_prefix (YYYYMMDD_HHMMSS)

    Falls back to recursive search under <project_root>/simulations.
    """
    bh_track_filename = bh_track_filename.strip()

    # project_root defaults to the directory containing this script
    if project_root is None:
        project_root = Path(__file__).resolve().parent

    sims_root = project_root / "simulations"

    # Extract run_prefix from filename
    # e.g. 20260121_114641__rp1p5__vinf10__inc20__toff2920__Om90__om180
    m = re.match(r"^(?P<prefix>.+?)__bh_radec__.+?\.xlsx$", bh_track_filename)
    if not m:
        # If pattern doesn't match, just search for the filename under simulations
        candidates = list(sims_root.rglob(bh_track_filename))
        if len(candidates) == 1:
            return candidates[0]
        if len(candidates) > 1:
            raise FileNotFoundError(
                f"Multiple matches for {bh_track_filename} under {sims_root}:\n"
                + "\n".join(str(c) for c in candidates[:10])
            )
        raise FileNotFoundError(f"Could not parse prefix and no file found for: {bh_track_filename}")

    run_prefix = m.group("prefix")

    # run_date is first token like 20260121_114641
    run_date = run_prefix.split("__", 1)[0]

    # Preferred expected path
    expected = sims_root / run_date / run_prefix / bh_track_filename
    if expected.exists():
        return expected

    # Fallback 1: maybe the prefix folder exists but file name differs slightly
    prefix_dir = sims_root / run_date / run_prefix
    if prefix_dir.exists():
        candidates = list(prefix_dir.glob("*bh_radec*.xlsx"))
        for c in candidates:
            if c.name == bh_track_filename:
                return c

    # Fallback 2: recursive search (slower but robust)
    candidates = list(sims_root.rglob(bh_track_filename))
    if len(candidates) == 1:
        return candidates[0]
    if len(candidates) > 1:
        raise FileNotFoundError(
            f"Multiple matches for {bh_track_filename} under {sims_root}:\n"
            + "\n".join(str(c) for c in candidates[:10])
        )

    raise FileNotFoundError(
        f"Could not find {bh_track_filename}.\n"
        f"Tried expected path: {expected}\n"
        f"and recursive search under: {sims_root}"
    )

# -----------------------
# USER SETTINGS
# -----------------------
SCHEDULE_CSV = "synthetic_rubin_schedule_REGENERATED_2026-05-04_to_2026-11-01_baseline3d_pairs30m.csv"

# Corridor around the BH track:
HALF_WIDTH_DEG = 0.08     # ~4.8 arcmin, same as we used in synthetic
PAD_ALONG_DEG  = 0.02     # extra padding along track endpoints (~1.2 arcmin)

# Detection rule:
K_SIGMA = 3.0             # require shift > K_SIGMA * sigma_ast
N_NIGHTS_REQUIRED = 5     # >= N distinct nights

# Gaia selection:
GAIA_MAG_LIMIT = 20.5     # tune; G<=20.5 gives manageable counts
ROW_LIMIT = 2_000_000     # safeguard; adjust as needed

# Astrometric noise model (mas per visit) vs Gaia G mag:
# (You can tune these constants to your preferred Rubin performance model.)
def sigma_ast_mas(Gmag: np.ndarray) -> np.ndarray:
    # ~10 mas at G=18, ~25 mas at G=20, ~63 mas at G=22
    return 10.0 * (10.0 ** (0.2 * (Gmag - 18.0)))

# -----------------------
# LOAD SCHEDULE
# -----------------------
sched = pd.read_csv(SCHEDULE_CSV)
sched["timestamp_utc"] = pd.to_datetime(sched["timestamp_utc"], utc=True)
sched = sched.sort_values("timestamp_utc").reset_index(drop=True)

import numpy as np
import pandas as pd

# --- Load BH track to get range_AU vs time ---
BH_TRACK_XLSX = "20260121_114641__rp1p5__vinf10__inc20__toff2920__Om90__om180__bh_radec__1p5__0p1.xlsx"
bh_path = resolve_bh_track_path(BH_TRACK_XLSX)
print(f"Using BH track: {bh_path}")
bh = pd.read_excel(bh_path)
bh["t"] = pd.to_datetime(bh["utc_iso"], utc=True)
bh = bh.sort_values("t").reset_index(drop=True)

tsec = bh["t"].astype("int64").to_numpy() / 1e9
range_au = bh["range_AU"].astype(float).to_numpy()

# interpolate range_AU to each schedule timestamp
x = sched["timestamp_utc"].astype("int64").to_numpy() / 1e9
sched["range_AU"] = np.interp(x, tsec, range_au)

# --- Compute thetaE_arcsec for M=0.1 Msun and Ds >> Dl ---
G = 6.67430e-11
c = 299792458.0
M_sun = 1.98847e30
AU = 1.495978707e11

M = 0.1 * M_sun
Dl_m = sched["range_AU"].to_numpy() * AU
thetaE_rad = np.sqrt((4 * G * M / c**2) / Dl_m)
sched["thetaE_arcsec"] = thetaE_rad * (180/np.pi) * 3600.0

# Lens track per visit
lens_ra = sched["ra_deg"].to_numpy()
lens_dec = sched["dec_deg"].to_numpy()
thetaE_arcsec = sched["thetaE_arcsec"].to_numpy()

# Unique nights
night = sched["timestamp_utc"].dt.floor("D")
unique_nights = np.sort(night.unique())
night_index = pd.Series(night).map({d:i for i,d in enumerate(unique_nights)}).to_numpy()
n_nights = len(unique_nights)

# Build a tangent-plane approximation around median track position (good for a few deg)
ra0 = float(np.median(lens_ra))
dec0 = float(np.median(lens_dec))
cosd0 = np.cos(np.deg2rad(dec0))

def radec_to_tan_deg(ra_deg, dec_deg):
    xi = (ra_deg - ra0) * cosd0
    eta = (dec_deg - dec0)
    return xi, eta

lens_xi, lens_eta = radec_to_tan_deg(lens_ra, lens_dec)

# Track polyline length in tangent plane (deg)
dx = np.diff(lens_xi); dy = np.diff(lens_eta)
seg_len = np.sqrt(dx*dx + dy*dy)
total_len = float(np.sum(seg_len))

# Corridor area (deg^2) approx
area_deg2 = (2 * HALF_WIDTH_DEG) * (total_len + 2*PAD_ALONG_DEG)

print(f"Track length (tan-plane): {total_len:.3f} deg")
print(f"Corridor area ~ {area_deg2:.3f} deg^2")

# -----------------------
# GAIA QUERY: pull real stars in a rectangular bounding box first (cheap),
# then filter to the corridor in Python (accurate).
# -----------------------
# Bounding box around track in RA/Dec:
ra_min = (lens_ra.min() - (HALF_WIDTH_DEG / cosd0) - PAD_ALONG_DEG) % 360
ra_max = (lens_ra.max() + (HALF_WIDTH_DEG / cosd0) + PAD_ALONG_DEG) % 360
dec_min = lens_dec.min() - HALF_WIDTH_DEG - PAD_ALONG_DEG
dec_max = lens_dec.max() + HALF_WIDTH_DEG + PAD_ALONG_DEG

# Handle RA wrap (if the box crosses 0/360)
wraps = ra_min > ra_max

Gaia.ROW_LIMIT = ROW_LIMIT

# Gaia DR3 table name in TAP
TABLE = "gaiadr3.gaia_source"

base_cols = """
source_id, ra, dec,
pmra, pmdec,
parallax,
phot_g_mean_mag, phot_bp_mean_mag, phot_rp_mean_mag
"""

mag_cut = f"phot_g_mean_mag <= {GAIA_MAG_LIMIT}"

if not wraps:
    adql = f"""
    SELECT {base_cols}
    FROM {TABLE}
    WHERE {mag_cut}
      AND ra BETWEEN {ra_min} AND {ra_max}
      AND dec BETWEEN {dec_min} AND {dec_max}
    """
else:
    # two pieces across 0°
    adql = f"""
    SELECT {base_cols}
    FROM {TABLE}
    WHERE {mag_cut}
      AND ((ra >= {ra_min}) OR (ra <= {ra_max}))
      AND dec BETWEEN {dec_min} AND {dec_max}
    """

print("Submitting Gaia ADQL query...")
job = Gaia.launch_job_async(adql)
gaia = job.get_results().to_pandas()
print(f"Fetched {len(gaia):,} Gaia stars (pre-corridor filter).")

# -----------------------
# FILTER STARS TO THE CORRIDOR AROUND THE TRACK
# Compute minimum distance from each star to the polyline in tangent plane.
# -----------------------
# Precompute segment endpoints
x0 = lens_xi[:-1]; y0 = lens_eta[:-1]
x1 = lens_xi[1:];  y1 = lens_eta[1:]
vx = x1 - x0; vy = y1 - y0
vv = vx*vx + vy*vy

# Stars tangent coords
sx, sy = radec_to_tan_deg(gaia["ra"].to_numpy(), gaia["dec"].to_numpy())

# Distance from point to polyline: min over segments
# Vectorized over stars but loop over segments (fast enough for typical sizes)
min_d = np.full(len(gaia), np.inf, dtype=float)
min_s = np.zeros(len(gaia), dtype=int)

for i in range(len(vx)):
    # projection t of star onto segment i
    px = sx - x0[i]; py = sy - y0[i]
    t = (px*vx[i] + py*vy[i]) / (vv[i] + 1e-30)
    t = np.clip(t, 0.0, 1.0)
    projx = x0[i] + t*vx[i]
    projy = y0[i] + t*vy[i]
    d = np.sqrt((sx - projx)**2 + (sy - projy)**2)  # deg
    better = d < min_d
    min_d[better] = d[better]
    min_s[better] = i

gaia["dist_to_track_deg"] = min_d
gaia = gaia[gaia["dist_to_track_deg"] <= HALF_WIDTH_DEG].copy()
print(f"Kept {len(gaia):,} Gaia stars in corridor (<= {HALF_WIDTH_DEG} deg).")

# -----------------------
# PROPAGATE GAIA POSITIONS TO EACH OBSERVATION TIME
# Gaia ra/dec are at reference epoch J2016.0 for DR3.
# We apply linear proper motion (ignoring parallax here; optional to add later).
# -----------------------
t_ref = Time("2016-01-01T00:00:00", scale="tcb")  # approx J2016; good enough for PM propagation
t_obs = Time(sched["timestamp_utc"].dt.strftime("%Y-%m-%dT%H:%M:%S").to_list(), scale="utc")

# Gaia values
ra = gaia["ra"].to_numpy()
dec = gaia["dec"].to_numpy()
pmra = gaia["pmra"].to_numpy()   # mas/yr (includes cos(dec))
pmdec = gaia["pmdec"].to_numpy() # mas/yr
Gmag = gaia["phot_g_mean_mag"].to_numpy()

# Years since ref
dt_yr = (t_obs.tcb.jyear - 2016.0).astype(float)  # shape (n_visits,)

# Propagate (deg): pmra is mas/yr along RA*cosdec in Gaia; convert to deg in RA
# dRA_deg = (pmra mas/yr) * dt / (cosdec * 3.6e6 mas/deg)
cosdec = np.cos(np.deg2rad(dec))
dra_deg = (pmra[:, None] * dt_yr[None, :]) / (cosdec[:, None] * 3.6e6)
ddec_deg = (pmdec[:, None] * dt_yr[None, :]) / (3.6e6)

ra_t = (ra[:, None] + dra_deg) % 360.0
dec_t = (dec[:, None] + ddec_deg)

# Convert stars + lens to tangent plane and separations (arcsec)
star_xi, star_eta = radec_to_tan_deg(ra_t, dec_t)  # deg arrays (n_star, n_visits)
lens_xi_v, lens_eta_v = radec_to_tan_deg(lens_ra, lens_dec)  # deg (n_visits,)

d_xi = (star_xi - lens_xi_v[None, :]) * 3600.0  # arcsec
d_eta = (star_eta - lens_eta_v[None, :]) * 3600.0
sep_arcsec = np.sqrt(d_xi*d_xi + d_eta*d_eta)

# Centroid shift magnitude (arcsec): delta = (u/(u^2+2))*thetaE
u_sep = sep_arcsec / thetaE_arcsec[None, :]
delta_arcsec = (u_sep / (u_sep*u_sep + 2.0)) * thetaE_arcsec[None, :]
delta_mas = delta_arcsec * 1000.0

# Per-visit detection: delta > K_SIGMA * sigma_ast(G)
sig = sigma_ast_mas(Gmag)
det_visit = delta_mas > (K_SIGMA * sig[:, None])

# Collapse to per-night detection (any visit that night counts)
det_night = np.zeros((len(gaia), n_nights), dtype=bool)
for j in range(n_nights):
    cols = np.where(night_index == j)[0]
    det_night[:, j] = np.any(det_visit[:, cols], axis=1)

epochs = det_night.sum(axis=1)

max_shift = delta_mas.max(axis=1)
max_snr = (delta_mas / sig[:, None]).max(axis=1)

passes = epochs >= N_NIGHTS_REQUIRED

# First/last detection nights
night_dates = pd.to_datetime(unique_nights, utc=True)
idx = np.arange(n_nights)

first_idx = np.where(det_night, idx[None, :], np.inf).min(axis=1)
last_idx  = np.where(det_night, idx[None, :], -np.inf).max(axis=1)

# initialize as NaT
first_dt = np.full(len(first_idx), np.datetime64("NaT"), dtype="datetime64[ns]")
last_dt  = np.full(len(last_idx),  np.datetime64("NaT"), dtype="datetime64[ns]")

# only assign for stars that have at least one detection night
mask = np.isfinite(first_idx)

first_dt[mask] = night_dates.to_numpy()[first_idx[mask].astype(int)]
last_dt[mask]  = night_dates.to_numpy()[last_idx[mask].astype(int)]


# Build catalog
cat = gaia.copy()
cat["epochs_detected"] = epochs
cat["first_detect_utc"] = pd.to_datetime(first_dt)
cat["last_detect_utc"] = pd.to_datetime(last_dt)
cat["sigma_ast_mas_per_visit"] = sig
cat["max_shift_mas"] = max_shift
cat["max_snr"] = max_snr
cat["passes_geN_nights"] = passes

# Add a list of detected nights for passing stars
det_list = []
for i in np.where(passes)[0]:
    nights = night_dates[det_night[i]]
    det_list.append(";".join([d.strftime("%Y-%m-%d") for d in nights]))
pass_nights = np.full(len(cat), "", dtype=object)
pass_nights[passes] = det_list
cat["detected_nights_utc"] = pass_nights

# Save
cat.to_csv("gaia_real_detection_catalog_all.csv", index=False)
inp = cat[cat["passes_geN_nights"]].to_csv("gaia_real_detection_catalog_pass_ge5nights_1.csv", index=False)

out = "gaia_real_detection_catalog_pass_ge5nights_EXCELSAFE.csv"

df = pd.read_csv(inp, dtype={"source_id": "string"})
df["source_id"] = "'" + df["source_id"]  # leading apostrophe forces Excel to treat as text
df.to_csv(out, index=False)
print("Wrote:", out)


print("Wrote:")
print("  gaia_real_detection_catalog_all.csv")
print("  gaia_real_detection_catalog_pass_ge5nights.csv")
print(cat["passes_geN_nights"].value_counts())
