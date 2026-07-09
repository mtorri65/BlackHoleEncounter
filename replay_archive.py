#!/usr/bin/env python3
import time
import rebound
from pathlib import Path

# 1) Point this to the archive from the run you care about
run_dir = Path(r"simulations/20251112_224743/20251112_224743__rp10__vinf30")
archive = run_dir / "20251112_222344__rp10__vinf30__archive.bin"

print(f"Loading SimulationArchive from: {archive}")

# 2) Choose which snapshot to start from:
#    snapshot=-1  -> last snapshot in the archive (default)
#    snapshot=0   -> first snapshot
#    snapshot=42  -> 43rd snapshot, etc.
snapshot_index = -1

# This constructor reads directly from the SimulationArchive
sim = rebound.Simulation(str(archive), snapshot=snapshot_index)
print(f"Loaded snapshot {snapshot_index}, t = {sim.t:.3f} days")

# 3) Start the web visualizer on a different port
PORT = 1235
sim.start_server(port=PORT)
print(f"Visualizer running — open http://localhost:{PORT}  (press 'h' in viewer for help)")

# 4) Let the simulation run forward in time while you watch
try:
    while True:
        time.sleep(1.0)
except KeyboardInterrupt:
    print("\nReplay stopped.")
