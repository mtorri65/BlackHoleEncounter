# Lecture 7 — Turning a PDE into a program

*Act III begins. Lectures 1–6 assumed. This lecture is about numerical methods, and its argument is that the good choices here are not optimisations — they are the difference between a program that can be trusted and one that merely runs.*

---

## The whole equation, at last

Here is everything, in one line:

$$
C \frac{\partial T}{\partial t} = \underbrace{Q(x,t)}_{\text{L2, L5}} \cdot \underbrace{a(x,T)}_{\text{L3, L4}} - \underbrace{(A + B\,T)}_{\text{L1}} + \underbrace{D \frac{\partial}{\partial x}\left[(1-x^2)\frac{\partial T}{\partial x}\right]}_{\text{L3}}
$$

with `x = sin φ`, and `C` either one heat capacity or two (Lecture 6).

**Before reading on, name where each term came from.** If any of them is unfamiliar, that's the lecture to revisit — from here we stop adding physics and start worrying about arithmetic.

It's a partial differential equation: continuous in latitude, continuous in time, with no analytic solution once the ice feedback makes `a` depend on `T`. To learn anything we must discretise it. And **the choices made in that discretisation determine which bugs are possible.**

---

## Grid design: making a bug category impossible

Every global mean in this model is an average over the sphere. The area of a latitude band is proportional to `cos φ`, so on a uniform-in-latitude grid, every average must be written:

```python
T_global = np.sum(T * np.cos(phi)) / np.sum(np.cos(phi))
```

That is not hard. But it appears in the equilibrium diagnostic, the energy-budget check, the transient reduction, the sweep aggregation, the plotting code — **a dozen places**, each an opportunity to write `np.mean(T)` in a hurry. And such a bug is quiet: the number stays plausible, just wrong by a few percent, biased toward whichever latitudes are over-weighted.

So don't use that grid. Use cells that are **equal in `x = sin φ`**:

$$
x_i = -1 + \left(i + \tfrac{1}{2}\right)\Delta x, \qquad \Delta x = \frac{2}{n}
$$

Because `d(area) ∝ cos φ dφ = dx`, **equal steps in `x` are equal areas on the sphere**. Verified: the spread in cell area across the grid is `2 × 10⁻¹⁶` — machine precision.

The consequence:

```python
T_global = np.mean(T)     # already area-weighted
```

**The `cos φ` weighting bug cannot be written**, because there is no weighting step to get wrong. This is the third instance of a pattern this course keeps returning to — Lecture 2's clamp, Lecture 3's vanishing `(1−x²)`, and now this. *Arrange the structure so the error has nowhere to live.*

---

## Conservation as a structural property

The diffusion term is a **flux divergence**: energy leaving one cell enters its neighbour. Discretise it so that this is literally true of the arithmetic, by evaluating the weight `w = 1 − x²` at cell **interfaces** and differencing the fluxes:

$$
(L T)_i = \frac{w_{i+1/2}\,(T_{i+1} - T_i) - w_{i-1/2}\,(T_i - T_{i-1})}{\Delta x^2}
$$

Each interface flux appears **exactly twice** — added to one cell, subtracted from its neighbour. Sum over all cells and every interior term cancels algebraically. What's left is the two boundary fluxes, at `x = ±1` — where `w = 1 − x² = 0`.

So the total is **identically zero**, for any temperature field whatsoever. Measured:

| Field | `Σ (L·T)` | relative to `Σ |L·T|` |
|---|---|---|
| Smooth, n=180 | −9.1 × 10⁻¹¹ | 2.6 × 10⁻¹⁵ |
| Random, n=180 | +1.4 × 10⁻⁹ | **2.3 × 10⁻¹⁷** |

The absolute residual looks larger for the random field only because the tendencies themselves are larger. **Relative to the magnitudes being summed, the cancellation is at machine precision** — this is floating-point round-off, not physics.

And the operator annihilates a constant field to `1.8 × 10⁻¹²`: a planet at uniform temperature has no heat transport, as it must not.

**This is worth being precise about.** Energy conservation here is not a property we hoped for, tested, and were relieved to find. It is a **consequence of how the operator was built** — provable from the structure, before running anything. The test in the suite exists to catch someone breaking it later, not to discover whether it holds.

That distinction matters. A model that conserves energy *approximately*, because the errors happen to be small, will eventually meet a case where they aren't.

---

## Stiffness, and why the obvious method fails

Now time-stepping. The obvious approach — **explicit** Euler — evaluates everything at the current temperature and steps forward:

$$
T^{n+1} = T^n + \frac{\Delta t}{C}\left[ Q\,a(T^n) - (A + B T^n) + D\,L\,T^n \right]
$$

Simple, no linear algebra. And unusable here.

The problem is **stiffness**: the equation contains processes with wildly different natural timescales, and an explicit method is limited by the *fastest* one no matter which you care about.

The fast process is diffusion between adjacent cells. Its rate scales as `1/Δx²`, so halving the cell size *quadruples* it. Measured, by bisecting for the largest `Δt` that survives three simulated years:

| Resolution | Δx | Max stable explicit Δt |
|---|---|---|
| n = 45 | 0.0444 | **2.205 days** |
| n = 90 | 0.0222 | **0.529 days** |
| n = 180 | 0.0111 | **0.130 days** |

Each doubling of resolution cuts the allowed step by a factor of **4.17, then 4.05** — the `Δt ∝ Δx² ∝ 1/n²` scaling, confirmed.

At the resolution we actually use, an explicit scheme would need steps of about **three hours**, forever, in every run of every sweep — not because three-hour resolution tells us anything, but because anything longer explodes.

**This is the distinction to hold on to.** A timestep chosen for **accuracy** is one you picked because you want to resolve something. A timestep forced by **stability** is one the method imposed on you regardless of what you need. The second is a tax.

---

## The implicit scheme

Treat the troublesome terms at the **new** time instead:

$$
\left(\frac{C}{\Delta t} + B\right) T^{n+1} - D\,L\,T^{n+1} = \frac{C}{\Delta t}T^n + Q\,a(T^n) - A
$$

The left side is now a linear system. It has to be solved each step — but the matrix

$$
M_{\text{imp}} = \text{diag}\!\left(\frac{C}{\Delta t} + B\right) - D\,L
$$

**depends on nothing that changes.** So it's LU-factorised **once**, at construction, and every subsequent step is a cheap back-substitution. The per-step cost is comparable to explicit, and the stability limit is gone entirely.

### Why not implicit everything?

Because `a(x, T)` is the ice-albedo switch from Lecture 4, and it is **nonlinear and discontinuous**. Putting it inside the matrix would mean re-factorising every step — and worse, iterating to convergence against a step function.

So we split. This is **IMEX** — implicit-explicit:

- **Implicit:** diffusion and the linear OLR. Stiff, linear, constant.
- **Explicit:** `Q·a(T)`. Nonlinear, but *slow* — the seasonal cycle and the ice line move on timescales of months and years, not hours.

The rule is simple: **implicit for what's stiff, explicit for what's awkward.** It works because in this system those two sets don't overlap.

### The payoff, measured

| Δt | Global mean | 65 °N seasonal range |
|---|---|---|
| 1 day | 288.151 K | 11.54 K |
| 5 days | 288.131 | 11.45 |
| 10 days | 288.091 | 11.48 |
| **30 days** | **287.973** | **10.76** |

A **thirtyfold** change in timestep moves the global mean by **0.18 K**. Nothing blows up; nothing oscillates. The differences are *discretisation error* — they shrink as `Δt` shrinks, exactly as accuracy errors should.

At n = 180: explicit needs **2810 steps per year**; IMEX at `Δt = 30 days` needs **12**. A factor of **231**.

The default is 2 days — chosen for accuracy, with plenty of margin, which is a choice we are free to make precisely because stability isn't making it for us.

---

## When the split stops working

One honest complication, because it arose in this project after the scheme was built.

Lecture 1's linear OLR turned out to be unusable far from 288 K, and was replaced (optionally) by a nonlinear form. But `OLR` was in the **implicit** half — the constant matrix. A nonlinear `OLR(T)` would force re-factorisation every step, throwing away the scheme's main advantage.

The fix keeps the linear term in the matrix and puts the *correction* in the explicit source:

$$
\text{rhs} = \frac{C}{\Delta t}T^n + Q\,a(T^n) - \text{OLR}(T^n) + B\,T^n
$$

Watch what happens at convergence. When `T^{n+1} = T^n = T`, the `(C/Δt)T` terms cancel, and so do the two `B·T` terms — one from the matrix, one from the source. What remains is:

$$
D\,L\,T + Q\,a(T) - \text{OLR}(T) = 0
$$

**the true balance, with no trace of `B`.** The linear term has become a pure *preconditioner*: it stabilises the iteration and contributes nothing to the answer.

Two things worth noticing. First, setting `OLR(T) = A + B·T` collapses the source back to `Q·a(T) − A` — the original scheme exactly, which is why adding this option left all 54 existing tests passing untouched. Second, this is a general and reusable trick: **keep a linear approximation implicit for stability, and carry the difference explicitly.** You'll meet it again in any stiff nonlinear problem.

---

## What you now know

- An **equal-area grid in `x = sin φ`** makes global means plain cell means — so the `cos φ` weighting bug *cannot be written*.
- The **flux-divergence** form with `w` at interfaces conserves energy **identically**: interior fluxes cancel algebraically, boundary fluxes vanish because `w(±1) = 0`.
- Conservation is **structural, not empirical** — provable before running. The test guards it; it doesn't establish it.
- **Stiff** means the equation contains a process much faster than the one you care about. Diffusion here is stiff: its rate goes as `1/Δx²`.
- Explicit stability scales as **Δt ∝ Δx² ∝ 1/n²** — measured 2.205 / 0.529 / 0.130 days at n = 45 / 90 / 180.
- **IMEX:** implicit for stiff-and-linear, explicit for nonlinear-and-slow. The matrix is constant, so it's factorised once.
- The result: **Δt becomes an accuracy knob, not a stability limit.** A 30× change in Δt moves the answer 0.18 K; the scheme is 231× cheaper than explicit at n = 180.
- A nonlinear term can be handled by keeping a **linear preconditioner implicit** and carrying the difference explicitly — the linear part cancels exactly at convergence.

---

## Exercise 7 — Break the stepper

**Provide:** the EBM with the time-stepping scheme swappable between explicit and IMEX.

**Part A.** With IMEX, run identical scenarios at `Δt = 1, 5, 10, 30` days. Compare equilibrium global means and 65 °N seasonal ranges. Are the differences consistent with discretisation error, and how would you tell that from incipient instability?

**Part B.** Switch to fully explicit. Find the largest `Δt` that remains stable at `n_lat = 180`. Bisection is the efficient approach; define "stable" precisely before you start, since a run can look fine for a year and then diverge.

**Part C.** Repeat at `n_lat = 45` and `n_lat = 90`. Deduce the scaling of the stability limit with resolution, and check it against `Δx²`. Then predict the limit at `n_lat = 360` and verify.

**Part D.** Verify conservation directly: build a random temperature field, apply the diffusion operator, and compute the area-weighted sum. Report it **both** absolutely and relative to `Σ|L·T|`. Explain why the relative figure is the meaningful one.

**Part E.** Verify the operator annihilates a constant field, and explain physically why it must.

**Part F — extension.** Deliberately break conservation: replace the interface weights `w_{i±1/2}` with cell-centre values `w_i`. The operator will still look reasonable. Measure the conservation residual now, and run a long integration to see how the global mean drifts. Roughly how many simulated years before the error exceeds the 0.28 K signal from Lecture 6?

*Part F is the point of the exercise. A discretisation can be consistent — converging to the right PDE as `Δx → 0` — and still fail to conserve at finite resolution. Consistency and conservation are different properties, and only one of them is visible in a short run.*

---

**Next lecture:** we have a model that is physically motivated, numerically sound, and conserves energy to machine precision. It is also, as Lecture 6 demonstrated, entirely capable of being **confidently and quietly wrong** — the single-surface model passed every test in the suite while missing the signal the project existed to measure. So: how do you test a program when you don't know the right answer? What kinds of error do tests catch, and which kinds can they never catch?
