import numpy as np
import matplotlib.pyplot as plt

# -----------------------------
# Black hole physical parameters
# -----------------------------
M_sun = 1.98847e30        # kg
M_bh = 0.1 * M_sun        # kg (your BH mass)
G = 6.67430e-11           # m^3/kg/s^2
c = 299792458             # m/s

# Schwarzschild radius
Rs = 2 * G * M_bh / c**2   # meters

# -----------------------------------
# "Safe distance" for spaghettification
# -----------------------------------
# Rough criterion: when differential tidal acceleration across 2 m human = 1g
human_height = 2.0        # m
g = 9.81                  # m/s^2

# Solve: 2GM * h / r^3 = g
r_safe = (2 * G * M_bh * human_height / g)**(1/3)

# -----------------------------------
# Domain: from Rs to just beyond the safe zone
# -----------------------------------
r_min = Rs
r_max = r_safe * 1.05   # 5% past safe distance so it is visible
r = np.linspace(r_min, r_max, 2000)

# -----------------------------------
# Compute tidal acceleration
# -----------------------------------
tidal_m_per_s2 = 2 * G * M_bh * human_height / r**3
tidal_g = tidal_m_per_s2 / g    # convert to g

# -----------------------------------
# Time dilation
# -----------------------------------
td = np.sqrt(1 - Rs / r)

# -----------------------------------
# Plot
# -----------------------------------
fig, ax1 = plt.subplots(figsize=(10,6))

# Tidal acceleration (left axis)
ax1.plot((r-Rs)/1000, tidal_g, color="red", label="Tidal acceleration (g)")
ax1.set_xlabel("Distance from event horizon (km)")
ax1.set_ylabel("Tidal acceleration (g)", color="red")
ax1.tick_params(axis='y', labelcolor="red")

# Mark safe distance
safe_dist_km = (r_safe - Rs) / 1000
ax1.axvline(safe_dist_km, color="purple", linestyle="--", label="Safe distance")
ax1.set_xlim(0, 2.0)

# Second y-axis for time dilation
ax2 = ax1.twinx()
ax2.plot((r-Rs)/1000, td, color="blue", label="Time dilation √(1−Rs/r)")
ax2.set_ylabel("Time dilation factor", color="blue")
ax2.tick_params(axis='y', labelcolor="blue")

# Labels and title
plt.title("Tidal Acceleration and Time Dilation vs. Distance from Event Horizon\n(0.1 Msun Black Hole)")
fig.tight_layout()

# Combine both legends
lines_1, labels_1 = ax1.get_legend_handles_labels()
lines_2, labels_2 = ax2.get_legend_handles_labels()
ax1.legend(lines_1 + lines_2, labels_1 + labels_2, loc="upper right")

plt.show()
