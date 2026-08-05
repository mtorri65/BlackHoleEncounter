# Lecture 10 — The full pipeline

*Capstone. Everything assumed. This lecture assembles the whole thing and runs it, and its real subject is the decisions that made the science possible — none of which were about physics.*

---

## What we have

Nine lectures have produced a climate model. Separately, this project has an N-body simulation: a black hole passing through the solar system, integrated 672 times across a grid of encounter geometries, each run rearranging Earth's orbit differently.

Two independent programs. The pipeline is obvious in outline:

```
REBOUND sweep  ──►  Earth's post-flyby orbit  ──►  seasonal EBM  ──►  climate
  (672 runs)          (a, e, ε, λ_p, year)         (Lectures 1-7)
```

The middle arrow is where all the difficulty turned out to live.

---

## Problem 1 — The data doesn't contain what we need

The REBOUND runs write `a`, `e`, and `q` for every planet, before and after the encounter. That is the orbit's **size and shape**.

A seasonal climate model needs its **orientation**: the longitude of perihelion `λ_p` — which hemisphere gets summer at closest approach — and the obliquity `ε`. Lecture 5 established that `λ_p` decides whether an orbital change is a Milankovitch trigger or a curiosity.

**Neither is in the file.** Nor is either recoverable from `a`, `e`, `q`, which describe shape alone and carry no orientation information whatsoever.

The archives that would have contained full state vectors were disabled for this sweep. What remains is a per-run spreadsheet of positions and velocities — 220 MB each, 150 GB total.

---

## Problem 2 — The abstraction invites the wrong answer

Before solving that, a trap worth naming, because it is subtle and the tooling actively encourages falling into it.

The climate model has a sweep harness. You write:

```yaml
sweep:
  ecc:  [0.0167, 0.06, 0.117]
  a_au: [0.9, 1.0, 1.1]
```

and it runs every combination. Perfect for exploring a parameter space.

**Completely wrong here.** A sweep block computes a **Cartesian product**. Our parameters are **paired** — run 137 has *its own* `(a, e, ε, λ_p)`, and those four values belong together because one simulated encounter produced them.

Feed 662 values of `a` against 662 of `e` and you get **438,244 combinations**, of which 662 correspond to worlds that exist and 437,582 are arbitrary pairings of numbers from unrelated simulations.

So the bridge is a script, not configuration.

**The general form is worth carrying away.** The available abstraction fit the *shape* of the problem — many parameter values, many runs — while being wrong about its *structure*. That kind of mismatch produces confident, plentiful, meaningless output. Ask not only "can this tool express my problem?" but "does its notion of a run match mine?"

---

## The detective story

So: recover orientation from raw state vectors. Here is the chain of inference, because the reasoning is more transferable than the result.

**Step 1 — Test, don't assume.** The simulation's coordinate frame isn't documented. Rather than guess, measure something the frame determines. Compute the inclination of Earth's orbital plane relative to the frame's z-axis at t = 0:

$$
\text{inclination}(\hat z,\ \hat h) = \mathbf{23.457°}
$$

**Step 2 — Recognise it.** That is not an arbitrary number. It is Earth's **obliquity**, to three decimal places.

**Step 3 — Conclude the frame.** An orbit inclined to the z-axis by exactly the obliquity means the z-axis is the **celestial pole**, not the ecliptic pole. The simulation runs in the **equatorial J2000 frame** — and therefore *the frame's z-axis is Earth's spin axis.*

**Step 4 — Exploit the physics.** REBOUND treats Earth as a point mass. A point mass has no oblateness, so nothing exerts a torque on its spin axis. **The spin axis is fixed in inertial space** for the entire 308-year integration, while the black hole freely tilts the *orbital plane*.

**Step 5 — Read off both quantities.** With `ŝ = ẑ` fixed and `ĥ` the orbit normal:

$$
\varepsilon = \angle(\hat z, \hat h), \qquad
\hat e_{\text{eq}} = \hat h \times \hat s_\perp, \qquad
\lambda_p = \angle(\hat e_{\text{eq}} \to \vec{e})\ \text{about}\ \hat h
$$

where `ŝ⊥` is the spin axis projected into the orbital plane and `ē` is the eccentricity vector.

**Both missing quantities were recoverable from data that was never intended to record them** — because a single measured number identified the frame, and the frame carried the physics.

---

## Validation

A recovery this indirect demands proof. The test is available for free: **apply it at t = 0, where the answer is present-day Earth.**

| Quantity | Recovered | Expected |
|---|---|---|
| `a` | 0.99981 AU | 0.99981 (independent source) |
| `e` | 0.01717 | 0.01717 |
| `ε` | 23.457° | 23.44° |
| **`λ_p`** | **282.29°** | ~283° |

Then three cross-checks that could each have failed independently:

**1. Perihelion date.** `λ_p ≈ 282°` should place closest approach in early January. The simulation's minimum heliocentric distance falls on **1 January** — a semantic check, testing that the number *means* what we claim.

**2. Zero variance.** All 672 runs begin from the same epoch, so all must return identical pre-flyby values. Measured variance across 672 independent extractions: **~10⁻²²**. Not "small" — floating-point identical.

**3. Independent bound count.** 662 of 672 runs leave Earth bound, matching the count from the separate `a, e, q` files exactly.

Note the tiers from Lecture 8 at work: the t = 0 comparison is a **physical benchmark**, the perihelion date is **semantic**, and the zero variance is **internal consistency** — which alone would have proven nothing about correctness, but combined with the others closes off the plausible failure modes.

---

## The approximation, stated honestly

One assumption deserves scrutiny: **Earth's spin axis is not truly fixed.** It precesses on a 25,772-year cycle.

Over this simulation's span it was **~0.70° from `ẑ_J2000`** in 1873 and **~1.01° the other way** by 2181. The resulting drift in true `λ_p` is about **4.3°**.

Is that acceptable? Yes, for two independent reasons:

1. **Scale.** The flyby-induced `λ_p` spread is the full **0–360°**. A 4.3° systematic offset is ~1% of the signal.
2. **Unavoidability.** REBOUND models Earth as a point mass, so precession **cannot be represented** regardless of what frame we choose. This isn't an error the analysis introduced; it's a limit inherited from the simulation.

And since before-and-after both use the same fixed `ẑ`, **differences remain internally consistent** even if the absolute zero-point drifts.

*This is what an honest approximation statement looks like: what is assumed, how large the error is, why it's tolerable, and — critically — whether anything could have been done differently.*

---

## An aside that saved 124 GB

The extraction needed the first and last row of each 220 MB spreadsheet. The obvious approach — open each file with a spreadsheet library — takes 26 seconds per run: **4.8 hours** for the sweep.

But a modern spreadsheet is a ZIP archive of XML. Stream-decompress the sheet directly and keep only the head and tail bytes: **0.83 seconds per run**, a **30× speedup**, and the whole extraction finishes in ten minutes.

The same observation applied to storage. Those 150 GB of spreadsheets contained **17 columns of which only 7 were independent** — the rest were derived quantities and duplicate string copies of every number, all stored as XML text. Converted to a binary columnar format:

| | |
|---|---|
| Original | **150.7 GB** |
| Parquet, float32 | **26.4 GB** |
| Reduction | **5.7×** |

Verified across all 672 runs by re-deriving the orbital elements from the converted files and matching the originals to ~10⁻⁶ AU.

**Neither of these is science.** But 4.8 hours versus 10 minutes is the difference between an analysis you iterate on and one you run once and hope. Thinking about your data's format is not separate from doing the work.

---

## Results

| | |
|---|---|
| Pre-flyby baseline | **288.19 K** |
| Median outcome | **278.97 K (−9.2 K)** |
| Recognisably Earth-like (±10 K) | **138 runs (21%)** |
| Snowball | **250 runs (39%)** |
| Runs analysed | 644 of 672 |

The distribution is **bimodal**, with the forbidden gap from Lecture 4 — the bifurcation, observed from 644 independent directions.

And an **overlap zone** at `a = 1.073–1.122 AU` where *both* outcomes occur. Within that band, semi-major axis alone doesn't determine the answer; eccentricity, obliquity and `λ_p` decide which side of the cliff a world lands on.

Remember also Lecture 9's warnings before interpreting any of this: two-thirds of these runs sit outside the linear OLR's defensible range, and 16% describe worlds that would have run away entirely.

---

## The finding to end on

The results above are **equilibrium** states — where each world settles. That was the natural thing to compute.

Then we computed the **transient**: the year-by-year adjustment from Earth's pre-flyby climate to its new one. And found this:

> **198 runs — 31% of the sweep — pass through a peak hemispheric temperature asymmetry above 5 K that decays to under 1 K.** Median peak: 3.74 K. Median value at equilibrium: **0.00 K**. Typical year of the peak: **year 1**.

Nearly a third of these worlds spend their first years with one hemisphere dramatically out of balance with the other — and then that asymmetry **vanishes completely**, leaving no trace in the equilibrium state.

**The equilibrium sweep could not have discovered this.** Not because it was buggy, or under-resolved, or badly analysed. Because equilibrium states, by definition, contain no information about the paths taken to reach them. We ran 644 correct calculations that were structurally incapable of showing us this.

> **Choosing what to compute determines what you are capable of discovering.**

That is the last idea in the course, and it generalises past climate entirely. Every analysis embeds a choice about which quantities are worth producing, and that choice quietly bounds the space of possible findings. The equilibrium sweep wasn't wrong. It was **answering a different question than the one that had the interesting answer**, and nothing internal to it could have revealed that.

The only defence is to periodically ask: *what would this analysis be unable to tell me, no matter what the answer turned out to be?*

---

## What you now know

- Missing quantities can sometimes be **recovered from data never meant to record them** — if you can identify what the data implicitly encodes.
- **Test the frame, don't assume it.** One measured inclination (23.457°) identified the frame and unlocked both `ε` and `λ_p`.
- A **point-mass integration applies no torque**, so the spin axis is fixed — which is what made the recovery valid.
- **Validate indirect recovery against a known case.** t = 0 must return present-day Earth, plus independent semantic and consistency checks.
- State approximations with **size, justification, and whether an alternative existed** — here 4.3° of λ_p drift, ~1% of the signal, and unmodellable anyway.
- **A Cartesian-product sweep is wrong for paired parameters** — 438,244 combinations of which 662 are real. The abstraction fit the shape and not the structure.
- **Format choices are part of the work**: 30× faster extraction, 5.7× smaller storage, same science.
- Results: baseline 288.19 K, median −9.2 K, **21% Earth-like, 39% snowball**, bimodal with an overlap zone at `a = 1.073–1.122 AU`.
- **31% of runs pass through a >5 K hemispheric asymmetry that leaves no equilibrium trace** — and the equilibrium sweep was structurally incapable of finding it.
- **Choosing what to compute determines what you can discover.**

---

## Exercise 10 — Capstone

Choose **one**, written up as a short report with figures.

**(a) Recovery.** Given a synthetic orbit with known `ε` and `λ_p`, implement the frame-based recovery and demonstrate it returns the known values. Then perturb the assumed frame orientation by 1° and quantify the induced error in `λ_p`. Compare that to the 4.3° precession drift and comment on which dominates.

**(b) Regime map.** Using the sweep outputs, produce a two-dimensional map in `(a, e)` coloured by outcome class — temperate / snowball / outside validity band. Identify the overlap zone and determine what distinguishes temperate from snowball outcomes *within* it. Is the boundary sharp?

**(c) Transient hunt.** Find the runs with the largest peak hemispheric asymmetry. Determine what orbital configuration produces them and explain the mechanism. Verify the asymmetry decays to zero at equilibrium, and estimate how long a hypothetical observer would have to watch to detect it.

**(d) Break the model honestly.** Choose one limitation from the report's §9 and quantify it. For example: replace the linear OLR with a grey-body `σT⁴` form, re-run the sweep, and report how many outcome classifications change and in which direction.

Whichever you choose, the report should state clearly **what you would need to have computed differently** to answer a question your analysis cannot address. That section is the point of the capstone.

---

## Where this leaves you

You have built a climate model from the energy balance of a sphere, and you know exactly which parts of it to disbelieve.

More usefully, you have watched a real project **get things wrong and correct itself**: a heat-capacity average that was the wrong operation, a test that would have certified a broken solver, a correlation that was an artefact of how the data was generated, a comparison whose sign depended on the regime, and an entire class of finding that the original analysis was structurally unable to see.

None of those were failures of effort or care. They were failures of a kind that only show up when someone checks — and the habit of checking is the only thing this course was ever really about.

**Deliberately absent from this model:** the carbonate–silicate thermostat, clouds, ocean heat transport beyond diffusion, precession and obliquity cycling within a run, and any geometry beyond the zonal mean. Each is a natural extension, and each would change some of the numbers you have just spent ten lectures learning to distrust appropriately.
