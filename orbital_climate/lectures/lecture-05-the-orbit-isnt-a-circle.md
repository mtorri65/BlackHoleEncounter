# Lecture 5 — The orbit isn't a circle

*Act II begins. You'll need Lectures 1–4. There's more mathematics here than anywhere else in the course, and one equation that cannot be solved in closed form at all.*

---

## The lie we've been telling

Every lecture so far has quietly set `r = 1`. Earth, we assumed, stays exactly one astronomical unit from the Sun, all year, forever.

It doesn't. Earth's orbit is an ellipse — a very round one, but an ellipse. And the whole point of this project is worlds whose orbits are *not* round, because a passing black hole has rearranged them.

So the `1/r²` we've been ignoring has to come back. Doing that properly requires solving an equation that Kepler wrote down in 1609 and that nobody has solved in closed form since.

**Lead with the payoff, because it's the thing people get wrong.** By the end of this lecture you'll know that making an orbit *more* eccentric — more elongated, spending part of the year much further from the Sun — *increases* the average sunlight the planet receives.

Most people guess it decreases. Some guess it cancels out. Both are wrong, and the reason is Kepler's second law.

---

## Five numbers describe an orbit

We've met some of these. Here they are together, with what each one does to climate:

| Element | Symbol | What it controls |
|---|---|---|
| Semi-major axis | `a` | Overall distance → **the dominant control on temperature** |
| Eccentricity | `e` | How elongated → seasonal contrast, and (surprisingly) the annual mean |
| Obliquity | `ε` | Axial tilt → seasonal amplitude, and whether poles beat the equator (Lecture 2) |
| Longitude of perihelion | `λ_p` | **Which hemisphere's summer happens at closest approach** |
| Mean anomaly | `M` | Where the planet is right now |

The first four are the orbit's *shape and orientation* — fixed for a given world. Only `M` changes as time passes.

`λ_p` deserves a second look, because it's the one that sounds like bookkeeping and isn't. It decides whether your hemisphere's summer coincides with the close part of the orbit or the far part. Today `λ_p ≈ 283°`, which puts perihelion in early January — **southern** summer gets the close approach, northern summer gets the far one. Flip `λ_p` by 180° and you've swapped which hemisphere has mild seasons and which has harsh ones, without changing the orbit's shape at all.

---

## Three anomalies, and why you need all three

Here is the awkwardness at the heart of orbital mechanics.

**What advances uniformly is not what you can see.** A planet on an ellipse does not sweep out angle at a constant rate — it races through perihelion and dawdles at aphelion. But time, of course, ticks uniformly.

So we need three different angles, and each does a job the others can't:

**Mean anomaly `M`** — a fictitious angle that increases perfectly linearly with time. `M = 2π t / P`. It corresponds to nothing you could photograph; it's just "how far through the year we are." **This is what you step forward in a simulation**, because it is the one that tracks the clock.

**Eccentric anomaly `E`** — a geometric construction on the circle that circumscribes the ellipse. It has no physical meaning either. It exists because it's the bridge: it's the only one connecting `M` to real geometry through a manageable equation.

**True anomaly `ν`** — the actual, physical angle from perihelion to the planet, as seen from the Sun. **This is what the geometry needs.** It's what feeds the solar longitude `λ = ν + λ_p` from Lecture 2.

The workflow is therefore always the same:

$$
\text{time} \;\longrightarrow\; M \;\longrightarrow\; E \;\longrightarrow\; \nu \;\longrightarrow\; \text{geometry}
$$

You cannot skip the middle. `M` is the only one that tracks time; `ν` is the only one that describes reality; `E` is the only path between them.

---

## Kepler's equation

The bridge from `M` to `E` is:

$$
M = E - e \sin E
$$

Simple to write. **Impossible to invert in closed form.** There is no formula giving `E` in terms of `M` using elementary functions — this was proven, and four centuries of effort have produced series expansions and numerical methods but never a solution.

So we iterate. Newton's method: guess, measure the error, correct, repeat.

$$
E_{n+1} = E_n - \frac{E_n - e \sin E_n - M}{1 - e \cos E_n}
$$

It converges fast — a handful of iterations to machine precision. Our implementation reaches a residual of **1.8 × 10⁻¹⁵**, which is floating-point round-off. Kepler's equation is, in practice, solved.

Two implementation details worth knowing, because both are places people get burned:

- **Wrap `M` into [−π, π] before starting**, and add the removed whole turns back afterwards. Newton's method is well-behaved near zero and can wander badly if you hand it `M = 500 radians`.
- **The initial guess `E₀ = M + e·sin M`** is worth using. It's already close, and it keeps the iteration count low even at high eccentricity.

---

## Getting `ν`, and one trap

With `E` in hand, the true anomaly follows. You will find this formula in older references:

$$
\cos \nu = \frac{\cos E - e}{1 - e \cos E}
$$

**Don't use it.** Two failures, both silent:

1. **`arccos` only returns [0, π].** Half the orbit comes back with the wrong sign, and you won't notice until your seasons are mirrored.
2. **It loses precision through perihelion**, exactly where the planet moves fastest and you most need accuracy.

Use the half-angle form instead:

$$
\nu = 2 \arctan_2\!\left( \beta \sin(E/2),\; \cos(E/2) \right), \qquad \beta = \sqrt{\frac{1+e}{1-e}}
$$

`atan2` sees both components separately, so it gets the quadrant right everywhere and stays continuous through perihelion. Our round-trip test — `M → E → ν → M` for a thousand random inputs — closes to **1.8 × 10⁻¹⁵ radians**.

This is the same instinct as Lecture 2's clamp and Lecture 3's vanishing `(1−x²)`: **choose the formulation that has no bad cases**, rather than the one that needs guarding.

---

## The result you'll guess wrong

Now the payoff. Average the incoming flux over one complete orbit. The answer is:

$$
\langle S \rangle = \frac{S_0}{4 a^2 \sqrt{1 - e^2}}
$$

Stare at the `√(1−e²)` in the *denominator*. As `e` grows, that shrinks, so `⟨S⟩` **grows**. A more elongated orbit receives *more* average sunlight than a circular one of the same semi-major axis.

### Why, in one sentence

The flux goes as `1/r²`, and **the planet does not spend equal time at equal distances.**

Kepler's second law: it sweeps equal areas in equal times, so it moves fast near the Sun and slowly far away. You'd think that helps — less time spent in the bright region. But `1/r²` is *steeply* nonlinear. The flux gain from being close beats the time penalty of passing through quickly, and the arithmetic does not cancel.

More formally: the average of `1/r²` is not `1/⟨r⟩²`. For any convex function, averaging the function exceeds the function of the average — Jensen's inequality. The gap widens as the orbit gets more elongated.

### The size of it

For Earth's `e = 0.0167`, the factor `1/√(1−e²)` is **1.00014** — fourteen parts per hundred thousand. Utterly negligible.

For the perturbed orbit at the centre of this project, `e = 0.117`, it's **1.0069** — about +0.7%.

Also small! And that is itself the lesson: **eccentricity barely touches the annual mean.** What it does instead is redistribute sunlight violently *within* the year — which is where all the climate action is, and which the annual mean cannot see. We spent Lecture 3 learning that the global mean is blind to transport; here it's blind to eccentricity too.

---

## Milankovitch

Which brings us to why anyone cares.

**The Milankovitch hypothesis:** ice ages are paced by orbital variations — eccentricity, obliquity, and the precession of perihelion — acting *not* through the annual global mean, but through **high-latitude summer insolation**.

The mechanism is the one from Lecture 4, run in slow motion. Ice sheets grow when winter snow survives the summer. Winter snowfall varies little; what varies is whether summer is warm enough to melt it. So the controlling quantity is **summer sunlight at high northern latitudes** — conventionally **65 °N in June**, where the great northern ice sheets grew.

### The scenario driving this project

Move Earth's perihelion inward by 10% and its aphelion outward by 10%. The semi-major axis barely changes; the eccentricity goes from 0.0167 to about **0.117** — roughly **7× today's value, and about 2× the maximum Earth reaches in its natural Milankovitch cycles.**

The annual global mean, as we just computed, shifts by under a percent. Now look at what happens to June at 65 °N:

| | 65 °N June peak insolation |
|---|---|
| Today, `e = 0.0167` | **478.0 W/m²** |
| Perturbed, `e = 0.117` | **399.5 W/m²** |
| Change | **−16.4%** |

**A 0.7% change in the annual mean; a 16.4% collapse in the quantity that actually controls glaciation.**

For scale, the Milankovitch range over the ice-age cycles spans roughly 390–550 W/m² at this latitude. Our perturbed value of ~400 sits at the **glacial-inception end** of that range. Sustained, this is a configuration that grows ice sheets.

This single benchmark is the tightest validation target in the whole model. It exercises the Kepler solver, the anomaly conversions, the declination, the hour-angle clamp, and the daily-mean integral — and it has a known answer. If a change breaks any link in that chain, this number moves.

---

## What you now know

- An orbit is fixed by `a`, `e`, `ε`, `λ_p`; only the **mean anomaly `M`** changes with time.
- **`λ_p` decides which hemisphere gets summer at perihelion** — pure orientation, no change in shape.
- **Three anomalies, each necessary:** `M` is linear in time (step this), `E` bridges via Kepler's equation, `ν` is the real geometry.
- **Kepler's equation `M = E − e sin E` has no closed-form inverse.** Newton converges to machine precision (residual ~10⁻¹⁵).
- Get `ν` from the **half-angle `atan2` form**, never `arccos` — quadrants and perihelion precision.
- **Annual-mean insolation is `S₀/(4a²√(1−e²))`, and it *increases* with eccentricity** — because `1/r²` is convex and the planet doesn't spend equal time at equal distances.
- The effect on the annual mean is **tiny** (+0.7% at `e = 0.117`). Eccentricity's real work is *within* the year.
- **Milankovitch:** ice ages are paced by high-latitude summer insolation, not annual means.
- The benchmark: 65 °N June falls **478.0 → 399.5 W/m², −16.4%**, landing at the glacial-inception end of the natural range.

---

## Exercise 5 — Kepler and the annual mean

**Provide:** `orbital_climate/kepler.py`, `orbital_climate/insolation.py`.

**Part A.** Verify the Kepler solver residual `|E − e·sin E − M| < 10⁻¹¹` for `e ∈ {0, 0.1, 0.5, 0.9}`, sampling `M` across several full turns including negative values.

**Part B.** Round-trip `M → E → ν → M` for 1000 random `M`. Report the maximum error. Then deliberately swap in the `arccos` form for `ν` and re-run. Where do the failures appear, and can you explain their location from the shape of `arccos`?

**Part C.** Numerically integrate daily-mean insolation over a full year by **sampling mean anomaly uniformly**, and compare against `S₀/(4a²√(1−e²))`. You should match to a relative error around 10⁻⁶ or better.

**Part D — the important one.** Repeat Part C but sample the **true anomaly** uniformly instead. Compare again.

You will get an answer that is badly wrong, and increasingly wrong as `e` rises:

| e | Correct `⟨1/r²⟩` | Uniform-in-ν | Error |
|---|---|---|---|
| 0.0167 | 1.00014 | 1.00070 | +0.06% |
| 0.117 | 1.00692 | 1.03499 | **+2.8%** |
| 0.5 | 1.15470 | 2.00000 | **+73%** |
| 0.9 | 2.29416 | 38.91967 | **+1596%** |

Derive the closed form for the wrong answer — you should find `⟨1/r²⟩_ν = (1 + e²/2) / (a²(1−e²)²)` — and explain in one sentence why uniform steps in `ν` over-weight perihelion.

**Part E.** Reproduce the Milankovitch benchmark: peak daily-mean insolation at 65 °N for `e = 0.0167` and `e = 0.117`, both with `λ_p = 283°`. You should get 478.0 and 399.5 W/m².

Then set `λ_p = 103°` (perihelion moved half a year round) and repeat. What happened, and which hemisphere is now at risk of glaciation?

*Why Part D matters more than it looks: sampling `ν` uniformly **bypasses the Kepler solver entirely** — you'd never call it. So comparing that average against the analytic identity tests nothing; it's a tautology dressed as a validation. Sampling `M` forces the code through `M → E → ν`, which is exactly why that test is worth having. A test that cannot fail is not a test.*

---

**Next lecture:** we now have the full forcing — the right sunlight, at the right latitude, at the right moment. But we've been computing *equilibrium* temperatures, as though the planet responds instantly. It doesn't. Oceans take years to warm, land takes days, and a season only lasts a few months. When the forcing changes faster than the response, the answer depends on **timescales** rather than energy balance — and we'll find that the model as built has been quietly hiding an entire climate signal because of it.
