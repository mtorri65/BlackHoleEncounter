# Lecture 3 — Why don't the poles freeze and the tropics boil?

*You'll need Lectures 1 and 2. Some familiarity with the idea of a derivative helps for one section, but nothing is computed by hand.*

---

## Assembling what we have

Two lectures in, we own both halves of a real climate model.

From Lecture 1: absorbed sunlight balances emitted infrared, and we write the outgoing side as the empirical line `OLR = A + B·T`.

From Lecture 2: we can compute the sunlight arriving at any latitude on any day.

So let's just connect them. At each latitude, set what arrives equal to what leaves and solve for temperature:

$$
T(\varphi) = \frac{Q(\varphi)\,a(\varphi) - A}{B}
$$

Every latitude minds its own business. No latitude knows any other exists. This is called **local radiative equilibrium**, and it is the most natural thing to try.

It is also catastrophically wrong, and the way it fails tells us exactly what's missing.

---

## First, the absorbed fraction

One detail we glossed over. In Lecture 1 albedo was a single number, α ≈ 0.30. But the poles are bright with ice and take sunlight at a punishing slant, while the tropics are darker and lit from overhead. Absorption varies with latitude.

The model writes the **coalbedo** — the fraction *absorbed*, so 1 − α — as

$$
a(x) = a_0 + a_2 P_2(x), \qquad P_2(x) = \tfrac{1}{2}(3x^2 - 1)
$$

with `a₀ = 0.676`, `a₂ = −0.200`, and `x = sin(latitude)` as always. `P₂` is the second Legendre polynomial; for our purposes it's simply a smooth shape that is low at the equator and high at the poles.

| Latitude | coalbedo *a* | albedo α |
|---|---|---|
| Equator | 0.776 | 0.224 |
| Pole | 0.479 | 0.521 |

The poles reflect over half of what reaches them; the tropics reflect under a quarter.

**And here's the tidy part.** `P₂` averages to exactly zero over the globe, so the *global mean* coalbedo is precisely `a₀ = 0.676`. If you did Exercise 1 Part D, that number is the one you solved for. It was never arbitrary — it's the global-mean absorbed fraction, and now it has latitude structure wrapped around it.

---

## Running the broken model

Now compute `T(φ) = (Q·a − A)/B` at every latitude, with ice effects switched off so we see the geometry cleanly.

| | Local radiative equilibrium | Observed Earth |
|---|---|---|
| Equator | **+57 °C** | ~+26 °C |
| Pole | **−57 °C** | ~−25 °C |
| Equator-to-pole difference | **113 K** | ~40–45 K |

The equator is hotter than the hottest place ever recorded on Earth. The pole is colder than almost anywhere outside Antarctica's interior. And the gradient between them is roughly **two and a half times** what the real planet exhibits.

*(A note on the syllabus, if you're following it: it describes this gradient as "roughly double" the observed one. Computed properly it comes out at about 2.6×. Worth checking numbers rather than inheriting them — including mine.)*

The model isn't broken in a random direction. It is **too extreme everywhere at once**: too hot where there's surplus energy, too cold where there's deficit. That is the signature of a specific missing ingredient.

Nothing is moving energy from where there's too much to where there's too little.

Which, of course, is what winds and ocean currents do all day.

---

## The transport term

We add one term to the balance:

$$
D \frac{\partial}{\partial x}\left[ (1 - x^2) \frac{\partial T}{\partial x} \right]
$$

This is a **diffusion** term. Its logic is the simplest possible: energy flows down the temperature gradient, at a rate proportional to how steep that gradient is. Steeper gradient, faster flow. Where the gradient is flat, nothing moves.

The full model, at last:

$$
C \frac{\partial T}{\partial t} = Q(x,t)\,a(x,T) - (A + B\,T) + D \frac{\partial}{\partial x}\left[(1-x^2)\frac{\partial T}{\partial x}\right]
$$

That equation is the entire remaining course. Everything from here refines a piece of it.

### Be honest about what this is

You should feel some discomfort, and I want to name it rather than let it fester.

Real poleward heat transport is done by **baroclinic eddies** — the mid-latitude weather systems on any forecast map — plus the **Hadley circulation** in the tropics and **ocean gyres** like the Gulf Stream. These are large-scale, organised, three-dimensional flows. Some carry heat as latent energy in water vapour.

**None of that is diffusion.** Diffusion describes the random jostling of molecules. Weather systems are neither random nor molecular.

So this term is a **parameterisation**, not a derivation. We are asserting that the *net* effect of enormously complicated fluid dynamics can be approximated by "flow proportional to gradient," and then choosing `D` to make the answer come out right.

`D = 0.58 W/m²/°C` is not measured. It is **fitted** — the value that reproduces the observed equator-to-pole gradient.

This is the **second-largest approximation in the model**, after the linear OLR from Lecture 1. Recall the discipline we set there: know where your approximations were fitted, and refuse to trust results from outside that range. The same caution applies here. A planet whose circulation differs radically from Earth's has no business being modelled with Earth's `D`, and several of the worlds later in this course are exactly that.

### The factor that vanishes

Look at `(1 − x²)`. Since `x = sin φ`, that is exactly `cos²φ` — and it goes to **zero at both poles**.

That's not decoration. It's a statement that the *area available for heat to flow through* shrinks to nothing at a pole, because a pole is a point.

The consequence is quietly excellent: **the no-flux boundary condition comes for free.** We never have to tell the model "don't let energy leak out of the north pole." The geometry makes leakage impossible, because the term it would flow through is zero there.

Compare this to the arccos clamp in Lecture 2. Same idea in a different costume: **arrange the mathematics so the awkward case cannot arise**, rather than detecting and handling it. Lecture 7 returns to why this matters enormously once we discretise onto a grid.

---

## Turning the dial

Now vary `D` and watch. Ice is still disabled, so this is transport alone.

| D | Equator | Pole | Equator − pole | **Global mean** |
|---|---|---|---|---|
| 0.00 | +57.0 °C | −56.5 °C | 113.5 K | **15.902 °C** |
| 0.10 | +47.6 | −42.2 | 89.8 | **15.902** |
| 0.30 | +37.7 | −24.8 | 62.5 | **15.902** |
| **0.58** | **+31.1** | **−12.6** | **43.7** | **15.902** |
| 1.50 | +23.5 | +1.5 | 22.0 | **15.902** |
| 5.00 | +18.5 | +10.9 | 7.6 | **15.902** |

Two things to see.

**First, the calibration works.** At `D = 0.58` the gradient is **43.7 K**, landing squarely in the observed 40–45 K range. The equator sits at +31 °C and the pole at −13 °C. Not perfect, but recognisably Earth.

At `D = 5.0` the planet is nearly isothermal — pole and equator within 8 K, a world of endless mild weather. At `D = 0` we have the broken model we started with.

**Second — and this is the real lesson — look at the last column.**

---

## The column that doesn't move

The global mean temperature is **15.902 °C for every single value of D.** Not approximately. Identically, to every digit the model produces.

Most people find this genuinely surprising, and expect a well-mixed planet to be warmer or cooler than a stratified one. It is neither.

The reason is worth stating carefully:

> **Diffusion redistributes energy. It cannot create or destroy it.**

The transport term moves heat between latitudes. Summed over the whole globe, everything it takes from one place it delivers to another, so its total contribution is **exactly zero**. And if it contributes nothing globally, it cannot shift the global energy balance — which is what sets the global mean temperature.

The global mean is fixed by radiation alone: absorbed in, emitted out, exactly as in Lecture 1. **All of Lecture 2's geometry and all of this lecture's transport rearrange that temperature across the planet without changing its average.**

This isn't a numerical accident either. Recall from Lecture 1 that the planet's radius cancelled out. This is the same kind of fact: a structural property of the equations, provable rather than measured. In Lecture 8 we'll turn it into a *test* — because a property that must hold exactly is the best possible check that your code is correct.

---

## What you now know

- **Local radiative equilibrium** — each latitude balancing alone — gives a 113 K equator-to-pole gradient against an observed ~43 K. Too hot in the tropics, too cold at the poles.
- The **coalbedo** varies with latitude as `a₀ + a₂P₂(x)`: 0.776 absorbed at the equator, 0.479 at the poles. Its global mean is exactly `a₀ = 0.676`.
- Heat transport enters as a **diffusion term** `D ∂/∂x[(1−x²) ∂T/∂x]`.
- This is a **parameterisation, not a derivation.** Real transport is baroclinic eddies, the Hadley cell and ocean gyres — none of which is diffusion. `D = 0.58` is fitted to the observed gradient. It is the model's second-largest approximation.
- `(1 − x²) = cos²φ` **vanishes at the poles**, so the no-flux boundary condition is automatic — the same "design the edge case away" move as Lecture 2's clamp.
- Raising `D` flattens the planet; lowering it sharpens the contrast.
- **The global mean is completely unaffected by D.** Diffusion moves energy; it cannot make or destroy any.

---

## Exercise 3 — The transport dial

**Provide:** `orbital_climate/ebm.py`, with the ice-albedo feedback disabled (set `T_ice_degC` far below any achievable temperature).

**Part A.** Run to equilibrium with `D = 0`. Verify each latitude matches the analytic local-equilibrium formula `T = (Q·a − A)/B` to within 0.05 °C. If it doesn't, your solver has a bug — find it before continuing, because everything downstream inherits it.

**Part B.** Sweep `D` over `{0, 0.1, 0.3, 0.58, 1.5, 5.0}` and overlay the six equilibrium temperature profiles on one figure, temperature against latitude. Use a single-hue sequential colour scale ordered by `D` — this is a magnitude, so let lightness carry it.

**Part C.** Plot equator-minus-pole temperature difference against `D`. Is the relationship linear? What happens as `D` grows large, and what physical state is the planet approaching?

**Part D — the point of the exercise.** For every `D`, compute the global mean temperature. Confirm it is unchanged to within 0.01 K.

Then explain *why*, in two or three sentences, without using the word "diffusion". If your explanation would work equally well for any energy-conserving transport process, you've understood it.

**Part E — extension.** Real Earth's transport is not down-gradient everywhere; the Hadley circulation moves energy *poleward across the tropics* by a mechanism that a diffusion law represents poorly. Read off your `D = 0.58` profile at 20° and 30° latitude, compare to observed zonal-mean temperatures, and comment on where this model should be least trusted.

*Part D is the one that matters. Students who expected a stirred planet to be warmer have learned something structural about what transport can and cannot do.*

---

**Next lecture:** we disabled ice throughout this lecture, and that was hiding something. Ice is bright, so it reflects sunlight, which cools the surface, which makes more ice — a loop that feeds itself. Feedbacks like that don't merely adjust an answer; they can produce *several* valid answers for identical conditions, and abrupt jumps between them. We'll find that this model has a temperature range it simply cannot occupy, and that the planet you get can depend on the planet you started with.
