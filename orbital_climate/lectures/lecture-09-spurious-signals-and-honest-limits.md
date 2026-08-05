# Lecture 9 — Spurious signals and honest limits

*Lectures 1–8 assumed. The model is built, validated, and honestly bounded. Now we use it on 644 worlds — and each of the three sections below is a different way of being wrong while doing everything correctly.*

---

## The situation

We have equilibrium climates for 644 simulated planets, each produced by a black-hole flyby that rearranged Earth's orbit. Every world has a semi-major axis, an eccentricity, an obliquity, a longitude of perihelion, and a resulting temperature.

This is a dataset. The temptation is to treat it like one — correlate everything against everything and report what comes out.

That temptation is the subject of Part 1.

---

## Part 1 — A correlation that isn't there

Correlate equilibrium temperature against each orbital element:

| Predictor | Correlation with T |
|---|---|
| Annual-mean insolation `S_mean` | **+0.966** |
| Semi-major axis `a` | −0.790 |
| **Longitude of perihelion `λ_p`** | **+0.593** |
| Eccentricity `e` | +0.110 |
| Obliquity `ε` | +0.067 |

Look at `λ_p`. **+0.593** across 644 samples — statistically overwhelming, stronger than eccentricity and obliquity combined. Reported as-is, this says the orientation of a planet's orbit substantially controls its temperature.

**But we know that's false before running anything.**

`λ_p` says where perihelion sits relative to the seasons. From Lecture 5, annual-mean insolation is `S₀/(4a²√(1−e²))` — a function of `a` and `e` **only**. `λ_p` does not appear. It cannot appear. Rotating an orbit in its own plane changes *when* sunlight arrives, not *how much* arrives per year.

So a +0.593 correlation with the annual-mean temperature is not a discovery. It is a warning that something is wrong with the analysis.

### The confounder

Check whether the predictors are independent of each other:

$$
\text{corr}(\lambda_p,\ a) = -0.417
$$

They are not. In this sweep `λ_p` and `a` are **42% correlated** — not because of any physics, but because of how the data was generated. A black-hole flyby that strongly changes a planet's orbital orientation also tends to change its size; the encounter geometry couples them.

So `λ_p` is acting as a **proxy for `a`**, and `a` genuinely does control temperature (−0.790). The `λ_p` correlation is inherited, not earned.

### Controlling for it

Restrict to a narrow band of annual-mean insolation — hold the real cause roughly fixed and see what `λ_p` does on its own:

| Within `S_mean ∈ [330, 360]` W/m², n = 68 | |
|---|---|
| `corr(λ_p, T)` | **−0.220** |
| `corr(S_mean, T)` | **+0.998** |

The `λ_p` correlation collapses from +0.593 to −0.220 — **and changes sign.** Meanwhile `S_mean` explains essentially everything: +0.998.

### What makes this example unusually good

In most fields, confounding is diagnosed statistically and argued about afterwards. Here **the physics tells you the answer in advance.**

We could predict, from Lecture 5 alone, that `λ_p` must have near-zero direct effect on the annual mean. The statistics then either agree — confirming the analysis was done correctly — or disagree, indicating a mistake. The theory validated the statistics, not the other way round.

That is the strongest position you can occupy when analysing model output, and it is available far more often than people use it. **Before correlating anything, write down what the physics requires.** Anything the data says beyond that needs explaining.

### And the residual is real

Don't discard the −0.220. It is small, it is the right size, and it is **physically expected**: `λ_p` controls whether a hemisphere's summer falls at perihelion or aphelion, which changes seasonal extremes, which changes where ice forms, which — through the ice-albedo feedback of Lecture 4 — slightly changes the annual mean.

The naive analysis found a large effect that doesn't exist. Careful analysis finds a small effect that does. **Controlling for a confounder does not make an effect vanish; it reveals its actual size.**

---

## Part 2 — When a better model gives a "worse" answer

Now something subtler, and more uncomfortable.

Lecture 6 established that the two-surface land/ocean model is unambiguously more physical than the single-surface one. It resolves a 175× contrast in thermal inertia that the single-surface version averages away, and it recovers the Milankovitch signal the old model missed entirely (−7.53 K vs −0.04 K).

So: apply both models across all 220 habitable-band runs and compare summer cooling.

| | Median 65 °N summer change | Runs cooling > 4 K |
|---|---|---|
| Single surface | **−6.31 K** | 138 (63%) |
| Two-surface land | **−3.37 K** | 104 (47%) |

**The better model shows *less* cooling.** Roughly half as much, and in a third fewer runs.

If you expected the improvement to make the Milankovitch signal *stronger* everywhere — as Lecture 6 might have led you to — this looks like a contradiction.

### Why

In this sweep, semi-major axis ranges from **0.63 to 2.8 AU**. That is an enormous spread: annual-mean insolation varies by a factor of twenty, and **annual-mean changes dominate everything else.**

Now recall the time constants. The single-surface model has `τ ≈ 1.6 years` — longer than a year, so its seasonal cycle is heavily damped and its **summer temperature is essentially its annual mean.** When the annual mean drops 20 K, its "summer" drops 20 K too.

Land has `τ ≈ 6.6 days`. It tracks the seasonal cycle properly, so its summer sits **well above** the annual mean. When the annual mean falls, land summer falls too — but from a higher starting point, and buffered by the large seasonal amplitude.

So in a sweep dominated by orbital *size* changes, the single-surface model **over**-predicts summer cooling, because it conflates summer with the annual mean.

### Isolating the mechanism

Restrict to the 72 runs where `a` stayed within 2% of 1.0 AU — orbital **shape** changed, size didn't, annual mean roughly preserved:

| | Median 65 °N summer change |
|---|---|
| Single surface | **+0.48 K** (warming) |
| Two-surface land | **−0.98 K** (cooling) |

**A sign flip.** With the annual mean held fixed, the single-surface model reports summer *warming* while land correctly reports *cooling* — exactly the Milankovitch mechanism, and exactly the Lecture 6 result at smaller amplitude.

### The accurate statement

Write this out in your own words; it's the hardest idea in the lecture:

> **The two-surface model is correct in both regimes. But the *direction* of the correction it applies depends on what changed.**
>
> For orbital **shape** changes it reveals a summer signal the simpler model gets backwards. For orbital **size** changes it moderates an overestimate.

A more complete model is not uniformly "more" or "less" of anything. It is *different*, in a way that depends on the regime you apply it to. Anyone who summarised this as "the two-surface model shows stronger glacial signals" would be wrong across most of this dataset — and anyone who summarised it as "shows weaker signals" would be wrong in exactly the regime the project cares about.

**Report the regime, or don't report the comparison.**

---

## Part 3 — Numbers the model isn't entitled to

Finally, the limit that was flagged in Lecture 1 and has been waiting ever since.

The linear OLR `A + B·T` was fitted to satellite observations clustered near 288 K. Across our 644 worlds:

| | |
|---|---|
| Runs in the defensible 250–300 K band | **210 (33%)** |
| Runs outside it | **434 (67%)** |

Two-thirds of the results are extrapolations.

At the hot extreme the model reports **559.7 K**. At that temperature:

| | |
|---|---|
| `A + B·T` gives | **802 W/m²** |
| Blackbody physics demands | **5565 W/m²** |
| Factor wrong | **6.9** |

Lecture 1 promised this exact number, in week one, before any of the machinery existed:

> *"Some will come out at 560 K. At 560 K, `A + B·T` gives roughly 800 W/m² of outgoing radiation. Actual blackbody physics demands around 5500. The model will be wrong by a factor of seven. And it will not warn you."*

It didn't warn us. It returned `559.7 K` — three significant figures, a decimal place, formatted identically to the numbers that are correct.

### What to do about it

Three responses, in increasing order of honesty:

**1. Report the number.** Wrong, and indefensible once you know.

**2. Report it as ordinal, not quantitative.** "This world is hotter than that one" survives extrapolation even when "this world is 559.7 K" doesn't. Rankings are more robust than values, because monotonicity is preserved even where calibration fails.

**3. Detect the invalid regime and refuse.** Real moist atmospheres cannot radiate more than the Simpson–Nakajima ceiling (~300 W/m²); above it there is **no equilibrium at all** — the planet runs away. Compute the absorbed flux, compare it to that ceiling, and flag those runs as unmodellable rather than reporting a temperature.

The model now does the third. **102 runs (16%) are flagged as runaway** — worlds for which the honest output is not a number but "this planet has no stable climate."

And note the distinction that made option 3 available: a *cap* on OLR would have manufactured a stable state that doesn't exist. A *flag* says the model has left its domain. **Refusing to answer is a legitimate output.**

### A note on inherited numbers

The syllabus for this lecture quotes 33% in the 250–300 K band, which is what the report originally said and what I've reproduced above. Later analysis refined the window to **230–300 K** — because the linear form's failure below 230 K is qualitative, not merely quantitative — giving 35%, and separately identified the 16% runaway population. Both figures are defensible; they answer slightly different questions.

Mention this only to make the habit explicit one last time: **when you inherit a number, find out what question it answered.**

---

## What you now know

- A **+0.593 correlation** can be entirely spurious. `λ_p` inherited it from `a`, with which it is 42% collinear by construction of the dataset.
- **Control for the confounder** and it collapses to −0.220 while `S_mean` explains +0.998.
- **Physics can predict what the statistics must show.** Use theory to validate the analysis, not the reverse. Write the prediction down first.
- A residual after controlling is often **real and informative** — here, the genuine small seasonal/ice-albedo effect of `λ_p`.
- **Regime dependence:** the two-surface model shows *less* summer cooling across the full sweep (−3.37 vs −6.31 K) but *more* — with a **sign flip** — when restricted to shape-only orbit changes (−0.98 vs +0.48 K).
- A more complete model is **not uniformly "more" or "less"** of anything. Report the regime.
- **Only 33% of runs** sit in the linear OLR's defensible band. At 559.7 K the model is wrong by a factor of **6.9** and says nothing about it.
- Outside the calibration range, report results as **ordinal, not quantitative** — or detect the invalid regime and **refuse to answer**.
- **Flag, don't cap.** Capping invents a state that doesn't exist.

---

## Exercise 9 — Find the confounder yourself

**Provide:** `<STAMP>_climate.csv` and `<STAMP>_earth_elements.csv`. No guidance as to which predictors matter.

**Part A.** Compute pairwise correlations between equilibrium temperature and each of `a`, `e`, `ε`, `λ_p`. Rank them.

**Part B.** *Before any further analysis*, write down which you believe are causal and why. Commit to it in writing — the exercise depends on you being on record.

**Part C.** Compute the correlation matrix **among the predictors themselves**. Identify collinear pairs. Which of your Part B answers now look unsafe?

**Part D.** Construct the derived predictor `S₀/(4a²√(1−e²))` and correlate it with temperature. How does it compare to the best raw predictor, and why should you have expected that?

**Part E.** Partial out annual-mean insolation and recompute the `λ_p` correlation.

**Part F — the graded item.** *From physics alone*, predict what Part E should give — before looking at your answer. Then compare.

Your prediction should be reasoned, not guessed: state what `λ_p` does and does not affect, conclude what its direct correlation with annual-mean temperature must be, and predict the sign and rough size of any residual. **A correct prediction that precedes the calculation is worth more than a correct calculation.**

**Extension.** Restrict to runs with `|a − 1.0| < 2%` and reproduce the sign flip in summer-peak change between the single- and two-surface models. Then explain, in your own words, why "the two-surface model shows stronger glacial signals" is a false summary of this dataset — and why "shows weaker signals" is also false.

---

**Next lecture:** the capstone. We assemble everything — the N-body flyby, the recovery of orbital elements from a simulation that never recorded them, the climate model, and the honest limits — into the full pipeline, and run it end to end. You'll see how a frame convention discovered by accident made the whole thing possible, why 30 seconds of thought about file formats saved 124 GB, and what 644 rearranged worlds actually look like when you finally have the right to interpret them.
