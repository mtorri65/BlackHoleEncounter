import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

# ----------------
# FILES
# ----------------
SCHED_CSV = "synthetic_rubin_schedule_REGENERATED_2026-05-04_to_2026-11-01_baseline3d_pairs30m.csv"
PASS_CSV  = "gaia_real_detection_catalog_pass_ge5nights.csv"

# ----------------
# SETTINGS
# ----------------
HALF_WIDTH_DEG = 0.08
MAG_LIMIT = 12.1  # show only Gaia stars with G < 10
TARGET_ID = "6900394685110465408"  # mark this specific Gaia star
N_LABEL_BRIGHTEST = 5            # label the N brightest stars with their G mag

# ----------------
# LOAD DATA
# ----------------
sched = pd.read_csv(SCHED_CSV)
sched["timestamp_utc"] = pd.to_datetime(sched["timestamp_utc"], utc=True)
sched = sched.sort_values("timestamp_utc").reset_index(drop=True)

passcat = pd.read_csv(PASS_CSV, dtype={"source_id": "string"})
passcat["source_id"] = passcat["source_id"].str.replace(r"\.0$", "", regex=True)
passcat.loc[(passcat["phot_g_mean_mag"]-9.63).abs().idxmin(),
            ["source_id","phot_g_mean_mag","ra","dec","epochs_detected","max_snr"]]

# Filter to bright stars only
passcat = passcat[passcat["phot_g_mean_mag"] < MAG_LIMIT].copy()

ra = sched["ra_deg"].to_numpy()
dec = sched["dec_deg"].to_numpy()
t = sched["timestamp_utc"]

# ----------------
# KEY DATES
# ----------------
key_dates = [
    ("May 4", "2026-05-04"),
    ("Jun 1", "2026-06-01"),
    ("Jul 1", "2026-07-01"),
    ("Aug 1", "2026-08-01"),    
    ("Sep 1", "2026-09-01"),
    ("Oct 1", "2026-10-01"),
    ("Nov 1", "2026-11-01"),
]

key_idx = []
for label, d in key_dates:
    dt = pd.Timestamp(d, tz="UTC")
    i = np.argmin(np.abs((t - dt).to_numpy()))
    key_idx.append((label, i))

# ----------------
# TANGENT PLANE FOR CORRIDOR
# ----------------
ra0 = float(np.median(ra))
dec0 = float(np.median(dec))
cosd0 = np.cos(np.deg2rad(dec0))

def radec_to_tan(ra_deg, dec_deg):
    return (ra_deg - ra0) * cosd0, (dec_deg - dec0)

def tan_to_radec(xi, eta):
    return (ra0 + xi / cosd0) % 360.0, dec0 + eta

xi, eta = radec_to_tan(ra, dec)

dx = np.diff(xi)
dy = np.diff(eta)
seg_len = np.sqrt(dx*dx + dy*dy)
seg_len[seg_len == 0] = 1e-30

tx = dx / seg_len
ty = dy / seg_len
nx = -ty
ny = tx

nxi = np.zeros_like(xi)
neta = np.zeros_like(eta)
nxi[0], neta[0] = nx[0], ny[0]
nxi[-1], neta[-1] = nx[-1], ny[-1]

for i in range(1, len(xi)-1):
    vx = nx[i-1] + nx[i]
    vy = ny[i-1] + ny[i]
    norm = np.sqrt(vx*vx + vy*vy)
    if norm < 1e-30:
        vx, vy = nx[i], ny[i]
        norm = np.sqrt(vx*vx + vy*vy)
    nxi[i] = vx / norm
    neta[i] = vy / norm

xi_p = xi + HALF_WIDTH_DEG * nxi
eta_p = eta + HALF_WIDTH_DEG * neta
xi_m = xi - HALF_WIDTH_DEG * nxi
eta_m = eta - HALF_WIDTH_DEG * neta

ra_p, dec_p = tan_to_radec(xi_p, eta_p)
ra_m, dec_m = tan_to_radec(xi_m, eta_m)

# ----------------
# FIND TARGET STAR
# ----------------
print(passcat.columns)
print(passcat["phot_g_mean_mag"].describe())
print(passcat["source_id"].tolist())


target = passcat[passcat["source_id"] == TARGET_ID]

#target = passcat[passcat["source_id"] == TARGET_ID].copy()
if target.empty:
    print(f"⚠️ Target Gaia star {TARGET_ID} not found in the (G<{MAG_LIMIT}) pass catalog.")
else:
    row = target.iloc[0]
    print("✅ Target star found:")
    print(row[["source_id", "ra", "dec", "phot_g_mean_mag", "epochs_detected", "max_snr"]])

# ----------------
# PLOT
# ----------------
plt.figure(figsize=(14.5, 5.0))
bg_color = (0/255, 1/255, 32/255)

ax = plt.gca()
ax.set_facecolor(bg_color)
ax.patch.set_alpha(0.0)

plt.gcf().patch.set_facecolor(bg_color)
plt.gcf().patch.set_alpha(0.0)

# BH track
plt.plot(ra, dec, lw=1.8, color="yellow", label="BH track")

# Direction arrows
step = 6
for i in range(0, len(ra) - step, step):
    plt.arrow(
        ra[i],
        dec[i],
        ra[i + step] - ra[i],
        dec[i + step] - dec[i],
        head_width=0.015,
        head_length=0.03,
        fc="yellow",
        ec="yellow",
        alpha=0.5,
        length_includes_head=True,
        zorder=1,
    )

# Corridor
plt.plot(
    ra_p, dec_p,
    "--", lw=1.4, color="red",
    label=f"Detect corridor (±{HALF_WIDTH_DEG}°)"
)
plt.plot(
    ra_m, dec_m,
    "--", lw=1.4, color="red"
)

# Bright Gaia stars
plt.scatter(
    passcat["ra"],
    passcat["dec"],
    s=35,
    color="gold",
    edgecolor="k",
    alpha=0.9,
    label=f"Gaia stars (G < {MAG_LIMIT})",
    zorder=3
)

# Label N brightest stars
brightest = passcat.nsmallest(N_LABEL_BRIGHTEST, "phot_g_mean_mag")
for _, r in brightest.iterrows():
    plt.text(
        r["ra"] + 0.04,
        r["dec"] - 0.03,
        f"G={r['phot_g_mean_mag']:.2f}",
        fontsize=9,
        color="black",
        weight="bold",
        zorder=5,
    )

# Mark target star
if not target.empty:
    plt.scatter(
        target["ra"],
        target["dec"],
        s=220,
        marker="*",
        color="red",
        edgecolor="black",
        linewidth=1.2,
        zorder=6,
        label=f"Target Gaia {TARGET_ID}"
    )
    r = target.iloc[0]
    plt.text(
        r["ra"] + 0.06,
        r["dec"] - 0.05,
        f"G={r['phot_g_mean_mag']:.2f}",
        fontsize=10,
        weight="bold",
        color="red",
        zorder=7
    )

# Reference dates
for label, i in key_idx:
    plt.scatter(
        ra[i],
        dec[i],
        s=90,
        color="white",
        marker="x",
        linewidths=1.8,
        zorder=6
    )

    plt.text(
        ra[i] + 0.05,
        dec[i] + 0.02,
        label,
        fontsize=8,
        weight="bold",
        color="white",
        zorder=5
    )

plt.xlabel("Right Ascension (deg)", color="white")
plt.ylabel("Declination (deg)", color="white")
ax.tick_params(colors="white")
plt.title(
    "BH track + detect corridor + bright Gaia stars (May 4 – Nov 1, 2026)",
    color="white"
)
plt.grid(True, lw=0.5, alpha=0.35, color="white")
leg = plt.legend(framealpha=0.85)
for text in leg.get_texts():
    text.set_color("white")
leg.get_frame().set_facecolor(bg_color)
leg.get_frame().set_edgecolor("white")
plt.tight_layout()
plt.savefig(
    "overlay_plot.png",
    dpi=200,
    bbox_inches="tight",
    transparent=True
)
plt.show()
