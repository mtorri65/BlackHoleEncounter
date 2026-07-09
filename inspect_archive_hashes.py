import rebound
import sys
from pathlib import Path

def load_archive(path: str):
    if hasattr(rebound, "SimulationArchive"):
        return rebound.SimulationArchive(path)
    return rebound.Simulationarchive(path)

def resolve_archive_path(user_arg: str, simulations_dir: str = "simulations") -> Path:
    p = Path(user_arg)
    if p.is_file():
        return p.resolve()
    root = Path(simulations_dir)
    matches = list(root.rglob(p.name))
    if len(matches) != 1:
        raise RuntimeError(f"Need exactly one match for {p.name}, found {len(matches)}")
    return matches[0].resolve()

archive = resolve_archive_path(sys.argv[1])
sa = load_archive(str(archive))
sim = sa[0]

names = ["Sun","BH","Mercury","Venus","Earth","Moon","Mars","Jupiter","Saturn","Uranus","Neptune",
         "mercury","venus","earth","moon","mars","jupiter","saturn","uranus","neptune"]

print("Archive:", archive)
print("First snapshot t:", sim.t)
print("N particles:", sim.N)
print("\nHash lookups:")
for nm in names:
    try:
        _ = sim.particles[nm]
        print(f"  ✔ {nm}")
    except Exception:
        print(f"  ✘ {nm}")
