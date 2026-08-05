# Lecture 6 — Time constants

*The pivotal lecture. Lectures 1–5 assumed. This one is told as a failure — a fix that looks obvious, doesn't work, and is instructive precisely because of how it fails. Read it in order; skipping to the answer wastes it.*

---

## Act 1 — The symptom

Our model is finished, in the sense that it has everything: latitude, seasons, transport, ice, and now a real elliptical orbit.

So let's ask it a question we know the answer to. **What is the seasonal temperature range at 65 °N?**

| | Seasonal range at 65 °N |
|---|---|
| Our model | **~10 K** |
| Observed, continental interiors | **~40 K** |

We are wrong by a factor of four.

And this is not a peripheral quantity. Recall Lecture 5: glacial inception is controlled by whether **summer** is warm enough to melt the winter's snow. Seasonal extremes are the mechanism. A model that damps the seasonal cycle to a quarter of its real amplitude is a model that cannot speak about ice ages at all.

With everything else in place, this is now **the most load-bearing simplification remaining.**

Where does it come from? From one number we've barely examined. The heat capacity `C` appears in the model as a single value for the entire planet:

$$
C \frac{\partial T}{\partial t} = \dots
$$

and it is set to an ocean-like `1.05 × 10⁸ J/m²/°C`. Every square metre of Earth — Siberia, the Pacific, the Sahara — is being modelled as though it has the thermal inertia of 25 metres of seawater.

---

## The quantity that matters: τ = C/B

Before fixing anything, name the relevant quantity.

Strip the model down to a single point with no transport and no ice. It becomes:

$$
C \frac{dT}{dt} = \text{forcing} - B\,T
$$

which is the standard relaxation equation. Its **thermal time constant** is

$$
\tau = \frac{C}{B}
$$

τ is how long the system takes to respond. And the seasonal cycle is a forcing with a period of one year, so what governs the seasonal amplitude is the **ratio of τ to that year**:

- **τ ≪ 1 year** → the surface tracks the forcing almost exactly. Large seasonal swings.
- **τ ≫ 1 year** → the surface averages over the forcing and barely moves. Damped seasons.

For our model, `τ = 1.05 × 10⁸ / 2.09 ≈ **1.59 years**` — longer than the forcing period. The seasonal cycle is being averaged away by construction.

**Note carefully what τ does and does not do.** It sets how far the temperature swings *around* its mean. It has no influence at all on where that mean sits — the equilibrium is fixed by the radiation balance, exactly as in Lecture 3. Heat capacity is a *rate* parameter, not an energy one.

---

## Act 2 — The obvious fix

Land and ocean have wildly different thermal inertia. A metre of soil holds far less heat than fifty metres of ocean mixed layer. So make `C` depend on latitude, weighted by how much of that latitude band is land:

$$
C(\varphi) = f(\varphi)\, C_{\text{land}} + \left(1 - f(\varphi)\right) C_{\text{ocean}}
$$

with `f` the land fraction — Earth's real zonal profile, about 29% globally, most continental in the northern mid-latitudes.

This is genuinely a **two-line change**. The solver already stores `C` as a per-latitude array; nothing about the matrix structure has to change.

**Stop here and predict.** At 60 °N, which is 55% land, what do you expect τ to become? Commit to a number before reading on.

---

## Act 3 — It fails

| | τ at 60 °N |
|---|---|
| Uniform ocean-like `C` | 1.59 years |
| Area-blended `C` | **1.44 years** |

**A 9% change.**

Nine percent. We took a latitude that is more than half dry land, gave it a properly area-weighted heat capacity, and its response time moved by less than a tenth. The seasonal cycle is essentially as damped as it was. We have not fixed the problem; we have barely perturbed it.

If you predicted something dramatic, you are in good company — and the reason it fails is worth more than the fix.

---

## Act 4 — Why it fails

Here are the two heat capacities:

| Surface | C [J/m²/°C] | τ = C/B |
|---|---|---|
| Ocean, 50 m mixed layer | **2.1 × 10⁸** | 3.2 years |
| Land, ~1 m of soil | **1.2 × 10⁶** | **6.6 days** |

**A ratio of 175.**

Now compute the blend at 55% land:

$$
0.55 \times 1.2\times10^6 + 0.45 \times 2.1\times10^8 = 6.6\times10^5 + 9.45\times10^7
$$

The land term contributes **0.7%** of the total. The arithmetic mean of these two numbers is, to within a rounding error, *just the ocean value*. Weighting by area cannot rescue it: even if land covered 90% of the latitude, the ocean term would still dominate the sum.

So the blended model says: this latitude responds on a ~1.4-year timescale. But that describes **neither surface**. The real land responds in a week; the real ocean responds in three years. There is no physical surface anywhere on Earth with a 1.4-year time constant. We have computed an average that corresponds to nothing.

### The general lesson

This is worth lifting out of climate entirely:

> **Averaging a parameter is not the same as averaging the system's behaviour.**

You may only replace two components by their mean when the system's response is roughly *linear* in that parameter over the range involved. Here the response depends on τ in a strongly nonlinear way — fast surfaces track, slow surfaces damp — and the parameter spans **more than two orders of magnitude**.

The failure mode is universal, and it is at its worst exactly when the parameter range is widest. Whenever you are tempted to replace a heterogeneous population with its average, ask: *does anything in the real system actually have the average value?* Here, nothing does.

---

## Act 5 — The correct formulation

If one temperature cannot represent both surfaces, use two.

This is North & Coakley (1979). Each latitude carries a **land temperature** and an **ocean temperature**, each with its own heat capacity, exchanging energy with one another:

$$
C_l \frac{\partial T_l}{\partial t} = Q\,a(T_l) - (A + B T_l) + D\nabla^2 T_l + \nu\,(T_o - T_l)
$$

$$
C_o \frac{\partial T_o}{\partial t} = Q\,a(T_o) - (A + B T_o) + D\nabla^2 T_o - \nu\,(T_o - T_l)\cdot\frac{f}{1-f}
$$

Both surfaces see the same sunlight and radiate by the same law. What differs is `C` — and now the 175× contrast **survives**, because the two temperatures are never averaged together before the physics acts on them.

### Derive the f/(1−f) factor yourself

Before reading further, work out where that factor comes from. It takes ninety seconds and it is the kind of thing you remember for years.

The exchange term moves energy from ocean to land. In the land equation it appears as `ν(T_o − T_l)` — a flux **per square metre of land**. In the ocean equation it must appear as a flux **per square metre of ocean**.

But land and ocean occupy different areas. Energy conservation requires the *total* to balance:

$$
\underbrace{f \cdot \nu (T_o - T_l)}_{\text{gained by land}} = \underbrace{(1-f) \cdot X}_{\text{lost by ocean}}
\qquad\Longrightarrow\qquad
X = \nu (T_o - T_l)\,\frac{f}{1-f}
$$

That's it. **The factor is not a fudge — it's flux continuity.** Omit it and the model quietly creates or destroys energy at every latitude, in proportion to how far the land fraction is from a half.

### Implementation notes

Four things worth knowing, each of which bit during development:

- The **state becomes a 2N vector** `[T_land, T_ocean]`, and the implicit operator a **2N×2N block matrix**. It's still factorised once and reused — the cost is negligible (360×360 at typical resolution).
- **Ice-albedo is evaluated separately on each surface.** This is not a detail; see below.
- **Land fraction must be clamped** to `[10⁻³, 1−10⁻³]`. Antarctica is `f = 1`, which sends `f/(1−f)` to infinity. The clamp is the same "design the edge case away" move as Lecture 2's arccos.
- Blending the two temperatures for reporting must be **idempotent** — blending an already-blended field should be a no-op. Ours wasn't at first, and the resulting shape error was the only real bug in the implementation.

---

## Calibrating the coupling

`ν` controls how tightly the two surfaces are tied. Large `ν` drags land toward the ocean's mild cycle; small `ν` lets it run free. It has to be calibrated against observations — ~40 K seasonal range over land at 65 °N, ~8–9 K over ocean:

| ν [W/m²/°C] | Land range | Ocean range |
|---|---|---|
| 1.0 | 64.1 K | 7.3 K |
| 2.0 | 51.5 | 8.2 |
| 3.0 | 43.3 | 8.8 |
| **3.5** | **40.1** | **9.1** |
| 4.0 | 37.5 | 9.3 |
| 8.0 | 25.1 | 10.2 |

**ν = 3.5** matches both targets at once, and that is the value the model uses.

Note that no single `ν` was guaranteed to satisfy both constraints — the land and ocean ranges move in opposite directions as `ν` rises, so hitting both simultaneously is a genuine test rather than a fit. That it works is mild evidence the two-surface structure is right.

### And the global mean barely moves

| | Global mean T |
|---|---|
| Single surface | 288.14 K |
| Two surface (ν = 3.5) | **288.42 K** |

A shift of **0.28 K** after restructuring the entire thermal response of the planet.

This is **Lecture 3's lesson arriving from a new direction.** There, the global mean was untouched by transport. Here it is untouched by heat capacity. Both are rate parameters; equilibrium is set by radiation alone.

Worth recording: the report's author explicitly predicted that this change would require re-tuning the albedo to recover 288 K, and **was wrong** (report §7.1). The prediction was reasonable and the physics says otherwise. Note also that report §4.6 quotes 288.07 K for this comparison — that figure was computed before `ν` was calibrated, at `ν = 1.0`. With the calibrated `ν = 3.5` it is 288.42 K. Both are small; only one is current.

---

## Why this was worth five acts

Now the payoff. Take the pure eccentricity injection from Lecture 5 — `a` fixed at 1.0, `e` from 0.0167 to 0.117, the perturbation that drops 65 °N June insolation by 16.4% — and ask what happens to the **summer temperature**:

| Model | 65 °N summer peak change |
|---|---|
| Single surface | **−0.04 K** |
| Two-surface, **land** | **−7.53 K** |
| Two-surface, ocean | −0.02 K |

**The single-surface model reports essentially nothing.** A 16.4% collapse in summer sunlight — the canonical Milankovitch trigger — produces four hundredths of a degree.

The two-surface model, given the identical forcing, shows land summers cooling by **7.5 K**.

The reason is now obvious: with `τ = 6.6 days`, land tracks the insolation almost instantaneously and feels the full drop. With `τ = 3.2 years`, ocean averages it away to nothing. Blend them into one temperature and the ocean wins — which is precisely what the single-surface model was doing.

**For five lectures, the model has been silently unable to represent the mechanism this entire project exists to study.** Not wrong in a way that showed up as an error, or a bad fit, or an obviously silly number. It reported −0.04 K with the same composure it reports everything else.

That is what makes this the pivotal lecture. The bug was never in the code. It was in the **physical representation**, and no amount of testing the implementation would have found it.

---

## What you now know

- **τ = C/B** is the thermal time constant. Compared to the forcing period it sets seasonal amplitude — and nothing else.
- Heat capacity affects the **seasonal cycle but not the equilibrium mean**. It is a rate parameter.
- A single ocean-like `C` gives ~10 K seasonal range at 65 °N against ~40 K observed.
- **Area-blending the heat capacities fails** — only a 9% change in τ — because `C_ocean/C_land ≈ 175` and the arithmetic mean *is* the ocean term.
- **Averaging a parameter ≠ averaging the behaviour.** Worst when the parameter spans orders of magnitude, and the blended value may describe nothing that exists.
- **Two surfaces per latitude**, coupled by `ν`, preserve the contrast.
- The **f/(1−f)** factor is flux continuity, not a fudge — without it the model leaks energy.
- **ν = 3.5** hits both observational targets simultaneously: 40.1 K land, 9.1 K ocean.
- The global mean moves by 0.28 K — confirming Lecture 3 from a new angle.
- **Single-surface: −0.04 K. Two-surface land: −7.53 K.** The old model could not see the Milankovitch signal at all.

---

## Exercise 6 — Calibrate the coupling

**Provide:** the two-surface EBM with `ν` exposed as a free parameter, and the observational targets (~40 K land, ~8–9 K ocean at 65 °N). Don't look up the model's value of `ν` first.

**Part A.** Reproduce the failure. Implement the area-blended `C` and measure τ at 60 °N. Confirm it moves by roughly 9%. Then compute what fraction of the blended value comes from the land term.

**Part B.** Switch to the two-surface model. Sweep `ν ∈ {1.0, 2.0, 3.0, 3.5, 4.0, 6.0, 8.0}`.

**Part C.** For each `ν`, record the 65 °N seasonal range over land and over ocean separately. Plot both against `ν` on one figure.

**Part D.** Choose the `ν` that best matches both targets simultaneously, and justify your choice. Note that the two curves move in *opposite* directions — explain why, and why that makes satisfying both a meaningful test rather than a free fit.

**Part E.** Confirm the global mean temperature barely shifts relative to single-surface. Then explain why, referring back to Lecture 3. If your explanation invokes heat capacity at all, it's wrong.

**Extension.** With `ν` calibrated, rerun the pure eccentricity injection (`a = 1.0`, `e: 0.0167 → 0.117`) and measure the 65 °N summer peak change three ways: single-surface, two-surface land, two-surface ocean. You should find −0.04 K, −7.53 K, −0.02 K.

Then answer this: **the single-surface model was validated** — it conserved energy, reproduced 288 K, matched analytic solutions, and passed every test in the suite. How would you ever have discovered it was missing this signal? What kind of check finds a fault of this type?

*That last question is Lecture 8's subject, and it doesn't have a comfortable answer.*

---

**Next lecture:** Act III begins, and we turn from physics to the business of getting it into a computer that gives the right answer. The equation we've built is a PDE with a stiff diffusion term, and the obvious way to step it forward is unstable unless the timestep is punishingly small. We'll see why, meet the implicit scheme that fixes it, and find out why the model can take two-day steps when a naive method would need minutes — plus what to do when a nonlinear term won't fit into the matrix you factorised.
