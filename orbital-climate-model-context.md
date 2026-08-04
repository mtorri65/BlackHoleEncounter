# Orbital Perturbation Climate Model — Conversation Export

Context document exported from a Claude chat (2026-07-26). Intended as background for a Claude Code session building a numerical energy balance model (EBM) of Earth's climate under a perturbed orbit.

---

## 1. Equilibrium temperature vs. distance from the Sun

Zeroth-order radiative equilibrium: absorbed power over cross-section πR² equals blackbody emission over surface 4πR²:

    L(1-A)/(4πd²) · πR² = 4πR²σT⁴

    T_eq = [L(1-A)/(16πσd²)]^(1/4)  ∝  d^(-1/2)

- Planet radius cancels; only stellar luminosity L, albedo A, and distance d matter.
- Earth (A ≈ 0.3): T_eq ≈ 255 K; actual mean surface T ≈ 288 K. The 33 K gap is the greenhouse effect.
- Actual mean temperatures vs. the curve: Mercury 440 K (0.39 AU), Venus 737 K (0.72 AU — runaway greenhouse, far above curve), Earth 288 K (1.0 AU), Mars 215 K (1.52 AU), Jupiter 165 K (5.2 AU — internal heat), Saturn 134 K (9.6 AU).

Relevant literature lineage:
- **Energy balance models:** Budyko (1969), Sellers (1969) — ice-albedo feedback, bifurcations, snowball hysteresis. Benchmark formulation: North, Cahalan & Coakley (1981).
- **Habitable zone theory:** Kasting, Whitmire & Reynolds (1993); updated by Kopparapu et al. (2013). Carbonate-silicate thermostat flattens the effective T(d) relation between runaway greenhouse (~0.95 AU) and CO₂ condensation (~1.7 AU) limits.
- **3-D GCM studies:** Wolf & Toon; Leconte group — clouds, heat transport, tidal locking.

## 2. The perturbed orbit scenario

Perturbation: perihelion moves in by 10%, aphelion out by 10%.

- Today: a = 1 AU, e ≈ 0.0167, r_p ≈ 0.983 AU, r_a ≈ 1.017 AU.
- Perturbed: r_p′ = 0.885 AU, r_a′ = 1.118 AU.
- Result: a′ ≈ 1.0017 AU (nearly unchanged), e′ ≈ 0.117 (~7× today; ~2× the Milankovitch maximum of ~0.058).
- Effectively a pure eccentricity injection.

### Annual-mean insolation (exact result — use as a unit test)

Time-averaging 1/r² over the orbit via angular momentum conservation (dθ/dt = h/r²) gives:

    ⟨S⟩ = L/(4πa²) · 1/√(1−e²)

Annual-mean insolation depends only on a and e, and increases with eccentricity (Jensen's inequality: ⟨1/r²⟩ > 1/⟨r⟩²).

- Eccentricity factor: 1/√(1−0.117²) ≈ +0.69%
- Semi-major axis factor: −0.33%
- Net: ~+0.35% insolation → ΔT ≈ +0.25 K global annual mean (T ∝ S^(1/4)). Essentially negligible.

### Seasonal structure (where the action is)

- Perihelion/aphelion flux ratio: (r_a/r_p)² goes from ~1.07 → ~1.60.
- Instantaneous-equilibrium swing T ∝ 1/√r ≈ ±12% (~±15 K), damped and lagged by ocean mixed-layer inertia (time constant ~1.6 yr vs. 1-yr forcing).
- Kepler II asymmetry: half-year imbalance Δt ≈ 4eP/π ≈ 54 days. Northern spring+summer ≈ 210 days; autumn+winter ≈ 156 days.
- Perihelion currently ~Jan 3 (longitude of perihelion ~283°): southern summers at perihelion (+28% flux), northern summers at aphelion (−20% flux).

### Key regional result: 65°N summer insolation

Daily-mean insolation computed via Kepler solver + declination + hour-angle integral: June peak at 65°N drops from ~480 to ~400 W/m² (−17%). Winter values near zero in both cases (perihelion boost wasted on polar night). ~400 W/m² sits at the glacial-inception end of the Milankovitch range (~390–550 W/m²) — sustained, this triggers northern hemisphere glaciation over 10–20 kyr.

## 3. Human consequences of a sudden orbital change

Timescale separation: atmosphere responds in months, ocean mixed layer in ~1–2 yr (most of new seasonal cycle within 5–8 yr), snowlines in decades, ice sheets in millennia, deep ocean in centuries.

- **Years 0–10:** Simultaneous cool-summer crop failure across all northern breadbaskets (Prairies, Ukraine, Siberia); global reserves ~3–4 months; food price crises. Southern hemisphere gets escalating perihelion-summer heat as ocean warms. Shift to short-season crops (barley, oats, potatoes, rye). Calendar reform forced (~366-day year, 54-day fast/slow asymmetry).
- **Years 10–20:** Northern grain belt migrates south several hundred km; Arctic sea ice expands, closing NW Passage/Northern Sea Route; perennial snowfields appear on Baffin Island / arctic Canada / Scandinavia (Laurentide nursery). Southern monsoons intensify; interior Australia approaches uninhabitable summers.
- **Years 20–50:** New normal established; snowfields thicken to firn/young glacier ice; ice-albedo feedback adds to the ~4–6 K direct orbital summer cooling at high northern latitudes; equatorward/coastward demographic drift of hundreds of millions; geopolitical weight shifts toward tropics and southern coasts.
- **Long term:** Apsidal precession relative to equinoxes (~21 kyr combined period) inverts the hemispheric configuration in ~10.5 kyr. Deliberate greenhouse forcing becomes a plausible lever to prop up northern summers.

## 4. Numerical model plan

Target: 1-D (latitude-resolved) seasonal EBM with orbital forcing, plus a sweep harness.

Core components:
1. **Kepler solver** — Newton iteration on E − e·sin E = M; true anomaly ν; r = a(1 − e·cos E).
2. **Insolation module** — solar longitude λ = ν + λ_p; declination δ = asin(sin ε · sin λ); hour angle H = acos(−tan φ · tan δ) (clamped); daily mean Q = (S₀/π)/r² · (H sin φ sin δ + cos φ cos δ sin H).
3. **EBM core** — North-style: C ∂T/∂t = Q·a(T) − (A + B·T) + D∇²T with ice-albedo step in a(T). Use an implicit/semi-implicit scheme for the diffusion term (stability).
4. **Spin-up** — run to a repeating seasonal cycle before applying perturbation; for "sudden change" experiments, switch orbital elements at t=0 and track transient.
5. **Sweep harness** — configs over (e, a, λ_p, D, C_land/C_ocean, albedo parameters); parallel runs; output NetCDF/parquet + plots.

Validation targets (build as tests from the start):
- Analytic annual mean: ⟨S⟩ = S₀/(4√(1−e²)) — exact, exercises the whole Kepler pipeline.
- Reproduce present-day ~288 K global mean with tuned greenhouse (A, B) terms.
- Recover North et al. (1981) EBM benchmarks.
- 65°N June-peak drop ~480 → ~400 W/m² for e = 0.117, λ_p = 283°.

## 5. Model selection for the coding work

Anthropic docs guidance (as of July 2026): start with Claude Opus 4.8 for complex agentic coding; Claude Fable 5 for highest available capability; default recommendation is the most intelligent generally available model with effort dialed to manage cost (cost-per-task often lower on stronger models). Some pages reference newer Opus 5 / Sonnet 5 releases — check https://platform.claude.com/docs/en/about-claude/models/overview for the current lineup.

Task mapping:
- **Physics core, discretization, and numerics review** (where subtle bugs hide — sign errors in hour angle, unstable explicit timesteps): strongest available model (Fable 5 / Opus 4.8), including a separate verification/review pass.
- **Sweep harness, config, plotting, CLI, parallelization:** Sonnet-class is sufficient and fast.
- Fits a worktree-parallelized multi-agent setup: physics core and experiment infrastructure are cleanly separable modules; use a planning/review agent on the strong model and iteration agents on Sonnet.
