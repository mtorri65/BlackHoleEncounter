# run_animation.py
from pathlib import Path
import numpy as np

# --- Backend: TkAgg ---
import matplotlib
matplotlib.use("TkAgg")  # must be before importing pyplot
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation


# ===================== CONFIG =====================
SNAP_DIR = Path("simulations/20251101_220642")  # <- adjust if needed
STEP = 10                # use every Nth snapshot
INITIAL_SPEED = 0.5      # 1.0 = sim-time; smaller = slower; bigger = faster
INITIAL_MIN_MS = 80      # clamp min ms between frames; increase to slow down
POINT_SIZE = 12          # scatter point size
# ==================================================


def extract_positions(d):
    x = np.asarray(d["x"]).ravel()
    y = np.asarray(d["y"]).ravel()
    z = np.asarray(d["z"]).ravel() if "z" in d else np.zeros_like(x)
    n = min(len(x), len(y), len(z))
    if n == 0:
        raise ValueError("Empty position arrays in snapshot.")
    return np.column_stack([x[:n], y[:n], z[:n]])


def load_frames_and_times(snap_dir: Path, step: int):
    frames = sorted(snap_dir.glob("snapshot_*.npz"))[::step]
    if not frames:
        raise FileNotFoundError(f"No snapshot_*.npz found in {snap_dir}")

    times = []
    for snap in frames:
        d = np.load(snap, allow_pickle=True)
        t_arr = np.array(d["t"])
        t_val = float(t_arr.ravel()[-1]) if t_arr.size else 0.0
        times.append(t_val)

    times = np.array(times, dtype=float)

    diffs = np.diff(times, prepend=times[0])
    if not np.all(np.isfinite(diffs)) or np.allclose(diffs, 0.0):
        diffs = np.ones_like(times)  # fall back if timestamps are weird
    return frames, times, diffs


def main():
    frames, times, dt_sim = load_frames_and_times(SNAP_DIR, STEP)

    fig, ax = plt.subplots()
    ax.set_aspect("equal", adjustable="box")
    ax.set_xlabel("x")
    ax.set_ylabel("y")

    d0 = np.load(frames[0], allow_pickle=True)
    pos0 = extract_positions(d0)
    sc = ax.scatter(pos0[:, 0], pos0[:, 1], s=POINT_SIZE)
    r0 = max(1.0, np.percentile(np.hypot(pos0[:, 0], pos0[:, 1]), 95))
    ax.set_xlim(-r0, r0)
    ax.set_ylim(-r0, r0)
    title = ax.set_title("")  # we’ll fill below

    # Playback state
    speed = INITIAL_SPEED            # affects scaling of sim Δt
    min_ms = INITIAL_MIN_MS          # hard clamp (dominates if Δt is tiny)
    paused = False

    def update_title(fi):
        try:
            fig.canvas.manager.set_window_title(
                f"Orbit animation — Space: pause | +/- or arrows: speed | [ ]: min clamp | "
                f"speed×{speed:.2f}, min={min_ms}ms, frame={fi+1}/{len(frames)}"
            )
        except Exception:
            pass
        title.set_text(f"{frames[fi].name}  |  t={times[fi]:.6g}  |  speed×{speed:.2f}, min={min_ms}ms")

    # initial interval
    dt_ms0 = max(min_ms, (dt_sim[0] / max(speed, 1e-9)) * 1000.0)

    def update(fi):
        nonlocal min_ms  # read by clamp below
        if paused:
            return sc, title

        d = np.load(frames[fi], allow_pickle=True)
        pos = extract_positions(d)
        sc.set_offsets(pos[:, :2])

        # compute the interval for THIS frame using current speed & clamp
        dt_ms = max(min_ms, (dt_sim[fi] / max(speed, 1e-9)) * 1000.0)
        ani.event_source.interval = int(dt_ms)

        update_title(fi)
        return sc, title

    ani = FuncAnimation(
        fig,
        update,
        frames=range(len(frames)),
        interval=int(dt_ms0),
        blit=False,
        repeat=True,
    )

    def on_key(event):
        nonlocal paused, speed, min_ms
        k = (event.key or "").lower()

        # Normalize several ways of pressing the same key across Tk/Win:
        is_plus  = k in {"+", "plus", "="}
        is_minus = k in {"-", "minus", "_"}
        is_up    = k in {"up"}
        is_down  = k in {"down"}
        is_openb = k in {"[", "braceleft"}
        is_closeb= k in {"]", "braceright"}

        if k == " ":
            paused = not paused

        # Speed controls (affect scaling of sim-time)
        elif is_plus or is_up:
            speed = min(64.0, speed * 1.25)   # cap to avoid overflow
        elif is_minus or is_down:
            speed = max(1e-3, speed / 1.25)   # don’t hit zero

        # Clamp controls (dominates when sim Δt is tiny)
        elif is_openb:
            min_ms = min(2000, int(min_ms * 1.25))  # slower (bigger clamp)
        elif is_closeb:
            min_ms = max(5, int(min_ms / 1.25))     # faster (smaller clamp)

        # Force an immediate interval recompute for the current frame
        cur = ani.frame_seq.__next__()  # advance generator to get an index
        # NOTE: next update() call will apply new speed/min_ms and adjust interval

        # Update title now to reflect new settings (frame index may not change yet)
        update_title(cur if isinstance(cur, int) else 0)

        # Also print to console for visibility
        print(f"[keys] speed×{speed:.3f}, min_ms={min_ms}ms")

    fig.canvas.mpl_connect("key_press_event", on_key)

    update_title(0)
    plt.show()


if __name__ == "__main__":
    main()
