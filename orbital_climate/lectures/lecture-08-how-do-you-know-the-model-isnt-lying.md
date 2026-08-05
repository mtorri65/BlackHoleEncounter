# Lecture 8 — How do you know the model isn't lying to you?

*Lectures 1–7 assumed. No new physics and no new numerics. This lecture is about the question that decides whether any of the previous seven were worth doing, and it almost never gets taught.*

---

## The uncomfortable position

You have a model. It runs. It produces smooth, plausible numbers with decimal places.

**How do you know any of them are right?**

Not "does the code run without crashing" — that's nothing. Not "does it look reasonable" — Lecture 6 destroyed that standard, when a model that passed every test in its suite turned out to be structurally incapable of representing the one signal the project existed to measure. It reported **−0.04 K** with exactly the composure it reports correct answers.

The honest position is this: **you cannot validate a model against the truth, because if you knew the truth you wouldn't need the model.** What you can do is construct situations where you *do* know the answer, and demand the model reproduce them.

Not all such demands are equally strong. Here is a hierarchy.

---

## Tier 1 — Analytic limits

**The strongest test available.** Disable mechanisms until the equations have a closed-form solution, then require the code to reproduce it exactly.

In this model: set `D = 0` and suppress the ice feedback. Every latitude is now radiatively independent, and the answer is known:

$$
T(\varphi) = \frac{Q(\varphi)\,a_0 - A}{B}
$$

Measured agreement: **0.0002 °C.**

Look at what that single test exercises. To produce that number the code must correctly execute the Kepler solver, the anomaly conversions, the declination, the hour-angle clamp, the daily-mean integral, the grid construction, the time-stepping scheme, and the radiation terms. **Any error in that entire chain moves the result.** And the expected value came from algebra, not from the code.

That last property is what makes it Tier 1: **the answer has an independent origin.**

Analytic limits are always available if you look. Turn off a term. Take a parameter to zero or infinity. Impose a symmetry. The model becomes less interesting and more checkable, which is exactly the trade you want.

---

## Tier 2 — Conservation laws and identities

Quantities that *must* hold, by structure, regardless of parameters:

| Test | Result |
|---|---|
| Diffusion annihilates a constant field | 1.8 × 10⁻¹² |
| Area-weighted operator sum vanishes | 10⁻¹⁵–10⁻¹⁷ relative |
| Global budget closes: `⟨absorbed⟩ = A + B⟨T⟩` | to 0.05 W/m² |
| Round-trip `M → E → ν → M` | 1.8 × 10⁻¹⁵ rad |
| Annual mean `= S₀/(4a²√(1−e²))` | 5.7 × 10⁻⁷ relative |

These are nearly as strong as Tier 1 and often easier to write. They share the crucial property: **the expected value is 0, or 1, or a closed-form expression** — something no amount of coding error can influence.

A round-trip test is particularly cheap and particularly good. `M → E → ν → M` must return what it started with; you need know nothing about orbital mechanics to assert that.

But note the limitation, because it's the theme of this lecture: a round-trip test verifies that your forward and inverse transforms are **mutually consistent**. Both could be wrong in compensating ways. Consistency is not correctness.

---

## Tier 3 — Physical benchmarks

Comparisons against the real world:

| Benchmark | Model | Reality |
|---|---|---|
| Present-day global mean | 288.15 K | 288 K |
| Zonal land fraction, integrated | 0.290 | 0.29 |
| Perihelion date | Jan 1 | early Jan |
| 65 °N seasonal range, land | 40.1 K | ~40 K |

Weaker, for two reasons. **First, tuning.** The 288.15 K is not a prediction — `coalbedo_a0` was *adjusted* until it came out. A tuned quantity cannot validate the model that was tuned to produce it; it only confirms the tuning converged.

**Second, tolerance.** "Within a few percent of observations" accommodates a great deal of wrongness. Real Earth's global mean is uncertain at the few-tenths level and varies by year.

They aren't worthless — the land-fraction check would catch a genuine data error, and the 40.1 K land range was *not* tuned to two targets simultaneously — but treat them as sanity checks, not proofs.

---

## Tier 4 — Internal consistency

The weakest tier, and the easiest to over-trust:

- Pre-flyby orbital elements show **zero variance across all 672 runs** — as they must, since every run starts from the same epoch.
- The bound/unbound count, **662/672**, matches an independently computed source.

These catch **pipeline** bugs: mis-indexing, dropped rows, a file read twice, results silently overwritten. In a project processing 672 simulations they are genuinely valuable.

**But they cannot catch a physics error at all.** If the model's physics is wrong, it will be wrong with perfect consistency across all 672 runs, and every Tier 4 test will pass. Consistency measures reproducibility, not truth.

---

## The punchline: a test that would have certified a broken solver

Now the story this lecture exists for. It's from this project's own record.

An early test asserted the eccentric anomaly for `e = 0.2, M = 1.0`:

```python
E = solve_eccentric_anomaly(1.0, 0.2)
assert E == pytest.approx(1.1934205, abs=1e-6)
```

Textbook shape. Analytic-flavoured. Tight tolerance.

**The correct value is 1.1853242.**

| | |
|---|---|
| Asserted | 1.1934205 |
| Correct | **1.1853242** |
| Discrepancy | 0.0081 rad = **0.46°** |

Where did the wrong number come from? It was *written down as though known* — plausible, right magnitude, seven digits of false precision.

Now consider what that test does. With `abs=1e-6`, it **rejects the correct solver.** A developer running it sees a failure and concludes the solver is broken. The natural response is to go and "fix" working code until it matches.

And what would a passing solver look like? Solving backwards: a solver that produced 1.1934205 would be one using **`e = 0.20806` instead of 0.2** — a 4% error in eccentricity, silently applied to every orbit the model ever integrates.

**The test would not have caught a bug. It would have created one.**

The general principle:

> **A test is code. Code has bugs. A test whose expected value came from the thing it tests is a tautology wearing a lab coat.**

The fix is not cleverness, it's provenance. Every expected value must have an **independent origin**: algebra you did yourself, a published table, a different implementation, a closed-form identity. In this case the correct value was recovered by writing a four-line Newton iteration in a separate process and comparing.

Ask of every assertion in your suite: **where did this number come from?** If the answer is "I ran the code and pasted the output," you have written a change detector, not a test. Change detectors are useful — they catch regressions — but they certify nothing about correctness, and they *lock in* whatever was wrong on the day you wrote them.

---

## An experiment: what does a good suite actually catch?

The model has 61 tests, written across four tiers by someone trying to be careful. Let's measure them. Three bugs were injected, one at a time, and the full suite run against each.

**Bug 1 — sign error in the diffusion flux divergence.** `L[i,i] = +(wl+wr)/dx²` instead of negative. Heat now flows *up* the gradient.

> **Caught immediately.** `test_diffusion_annihilates_constant` fails. A Tier 2 conservation test, and the fastest in the suite.

**Bug 3 — off-by-one in the seasonal index.** Advance the orbital position *before* stepping instead of after, so temperature is computed against tomorrow's insolation.

> **Caught.** `test_local_radiative_equilibrium_no_diffusion` fails — the Tier 1 analytic limit. A subtle error, and the strongest test found it.

**Bug 2 — `cos φ` weighting applied to the already-equal-area grid.** Exactly the Lecture 7 bug that the grid design was supposed to make impossible.

> **All 61 tests pass.**

Read that again. A wrong global-mean weighting, undetected by a suite specifically built to check this model.

### Why it escaped — and this is the real lesson

Investigating revealed something neither I nor the test suite had noticed: **the model computes its global mean in two different places.**

- `EBM.global_mean()` — used for spin-up convergence and a recorded diagnostic
- `T_annual.mean()`, written out longhand in `experiment.py` — the one that produces **every reported result**

The bug went into the first. Nothing that any test asserts on flows through it.

Inject the *same conceptual bug* into the reporting path instead, and **three tests fail immediately**, including the 288 K benchmark — the reported mean jumps to 292.07 K.

So the identical error is caught or missed depending on **which of two duplicated code paths it lands in.**

> **Test coverage is a property of which code is executed, not which concepts you believed you checked.**

This is why coverage tooling exists, and why "we test the global mean" was a false statement about this suite. The duplication itself is the underlying fault: two implementations of one idea, only one of them tested, free to diverge silently.

---

## A taxonomy of real bugs

From this project's actual defect record, ordered by how hard they were to notice:

**Loud — crashes.** A regex assumed XML attribute order (`Id` before `Target`); the writer emitted the reverse. Nothing parsed. *Easy: the failure is total.*

**Medium — wrong shape.** A blend applied to already-blended data threw a dimension mismatch. *Easy: arrays have shapes, and shapes are checked.*

**Quiet — plausible wrong numbers.** A `--limit` flag truncated the dataset *before* statistics were computed, so a request for 8 runs reported **"654 outside validity band."**

That third class is the dangerous one. It didn't crash. It printed a well-formatted, correctly-typed, entirely wrong number, in a summary designed to be skimmed. It was caught only because 654 was implausible *for a request of 8* — a coincidence of noticing.

**The bugs that hurt you are the ones that return plausible output.** Everything else announces itself.

---

## What you now know

- **You cannot validate against truth.** You construct situations with known answers and demand agreement.
- **Tier 1, analytic limits** — disable mechanisms until closed form exists. Strongest, because the expected value has independent origin. Ours matches to 0.0002 °C.
- **Tier 2, conservation and identities** — expected values of 0, 1, or a formula. Nearly as strong.
- **Tier 3, physical benchmarks** — weakened by tuning and loose tolerances. Sanity checks.
- **Tier 4, internal consistency** — catches pipeline bugs; **cannot catch physics errors**, which are consistent too.
- **A test whose expected value came from the code under test is a tautology.** The 1.1934205 assertion would have rejected a correct solver and certified one with a 4% eccentricity error.
- **Every expected value needs provenance.** "I pasted the output" makes a change detector, not a test.
- The same bug is caught or missed depending on **which code path it lands in**. Coverage is about executed lines, not intended concepts.
- **Duplicated logic is a testing hazard**: two implementations, one tested, free to diverge.
- The dangerous bugs **return plausible numbers**.

---

## Exercise 8 — Write the tests, then break the code

**Provide:** the EBM package, plus a copy with three deliberately injected bugs.

**Part A.** *Before looking at the buggy version*, write four tests — one Tier 1 analytic limit, one Tier 2 conservation check, one round-trip, one Tier 3 physical benchmark. For each, write one sentence stating **where your expected value came from**. If you cannot answer that, the test isn't finished.

**Part B.** Run your suite against the buggy version. How many bugs does it catch?

**Part C.** For every bug you missed, write the test that would have caught it. Then ask why you didn't write it the first time — was the concept unconsidered, or considered but untested?

**Part D.** Take an existing test from `orbital_climate/tests/` and determine whether its expected value was derived independently or is a tautology. Some are genuinely one, some genuinely the other.

**Part E — extension.** Find a quantity in the codebase computed in more than one place. (One is documented in this lecture; there may be others.) Write a test asserting the two agree. Then consider: would that test have caught the `cos φ` bug in *either* location?

*Assessment note: a low catch rate in Part B is a legitimate and interesting outcome, not a failure. Most people catch the sign error and miss the weighting bug — the same result the model's own 61-test suite produced. Discovering that a careful suite misses a real defect is the point of the exercise, not an embarrassment.*

---

**Next lecture:** we now have a validated model and honest limits on it. So let's use it — and immediately run into trouble. We'll find a correlation of +0.59 between two variables that have no causal connection whatsoever, discover that a conclusion drawn from a controlled experiment reverses sign when applied to the full dataset, and see the model report confident temperatures for planets that physically cannot exist. Each is a different way of being wrong while doing everything right.
