# Lecture 4 — Ice, feedback, and tipping points

*Lectures 1–3 assumed. No new mathematics; the new ingredient is a single conditional, and the consequences are out of all proportion to it.*

---

## One line of code

Every lecture so far has held albedo fixed. Bright poles, dark tropics, and that structure never changed no matter what the temperature did.

But ice is white. If a place gets cold enough to freeze, it becomes far more reflective — and reflectivity is precisely what controls how much sunlight the planet keeps.

So we make the coalbedo depend on temperature:

$$
a(x, T) = \begin{cases}
a_{\text{ice}} = 0.38 & \text{if } T < T_{\text{ice}} = -10\ ^\circ\text{C} \\[4pt]
a_0 + a_2 P_2(x) & \text{otherwise}
\end{cases}
$$

That's it. That's the whole change. Where it's cold enough, the ground absorbs 38% of the sunlight reaching it instead of 48–78%.

**And with that one conditional, the model becomes nonlinear** — the temperature now feeds back into the very thing that determines the temperature. Everything remarkable in this lecture follows from that.

---

## The loop

Trace it slowly, because the whole lecture is contained in these five arrows:

```
   cooling  →  more ice  →  higher albedo  →  less sunlight absorbed  →  more cooling
       ↑                                                                      │
       └──────────────────────────────────────────────────────────────────────┘
```

The output feeds the input **with the same sign**. Cooling causes more cooling. That is a **positive feedback**, and the name is a trap worth defusing now: "positive" here means *self-reinforcing*, not *beneficial*. This particular feedback is capable of killing a biosphere.

Run it the other way and it works just as hard: warming melts ice, which darkens the surface, which absorbs more, which warms further.

Compare this to the feedback we've had all along without naming it. In Lecture 1, a warmer planet radiated more (`σT⁴`, or our linear stand-in), which cooled it back down — output opposing input, a **negative** feedback. That's what makes planetary temperature stable at all. It's the stiff spring.

**We now have both springs in the same model**, pulling against each other. When the negative feedback wins, the planet is stable. This lecture is about what happens when it loses.

---

## Dimming the Sun

Here's the experiment. Take present-day Earth and slowly turn down the Sun — multiply `S₀` by a factor `s`, and walk `s` down in 1% steps. After each step, let the model settle, then step again. Each run starts from where the last one finished, exactly as a real planet would experience a gradual change.

| s | Global mean T | Ice edge |
|---|---|---|
| 1.00 | 288.13 K | 70.8° |
| 0.98 | 285.43 | 64.2° |
| 0.96 | 282.55 | 61.4° |
| 0.94 | 279.47 | 56.4° |
| 0.92 | 276.06 | 52.1° |
| 0.90 | 272.17 | 46.2° |
| 0.89 | 269.71 | 42.7° |
| **0.88** | **266.67** | **39.3°** |
| **0.87** | **229.71** | **ice everywhere** |
| 0.86 | 229.09 | ice everywhere |

Read the first rows and nothing alarms you. Dim the Sun by 1%, lose a couple of degrees, ice creeps a few degrees equatorward. Perfectly proportionate. You could extrapolate confidently.

You would be extrapolating straight off a cliff.

Between `s = 0.88` and `s = 0.87` — one more 1% step, the same size as all the others — the planet drops **37 K** and freezes to the equator. The ice edge doesn't march from 39° to 30°. It goes to zero. Every latitude, frozen.

**Notice what the ice edge was doing just before.** From 70.8° down to 39.3°, each step moved it a little further than the one before: 6.6°, then 3.4°, then 3.5°... The acceleration was visible in the data, if you knew to look. That's the feedback gaining on the restoring force. By 39° the ice has reached latitudes that receive real sunlight, so each new band of ice removes far more absorbed energy than a polar band did — and the loop closes on itself.

---

## Nobody wrote this

Now the part that matters most in this course.

**Go and look at the code.** There is no snowball branch in it. There is no `if planet_should_freeze:`. Nobody wrote a rule saying "when the ice line passes 39 degrees, collapse to a frozen state." No threshold of that kind appears anywhere.

What is in the code is the four-line albedo conditional at the top of this lecture, plus the diffusion and radiation from Lectures 1–3. The cliff is not a feature of the program. **It is a consequence of the equations.**

This distinction — between behaviour that is **imposed** and behaviour that is **emergent** — is arguably the single most important idea in this course, and this is the cleanest demonstration of it you will get.

An imposed behaviour tells you only that someone anticipated it. An emergent one is a genuine prediction: the model is telling you something you did not put in, which means it *can* be wrong, which means it is worth testing.

The same signature shows up in the black-hole sweep this course is built around. Across 644 simulated worlds, the equilibrium temperatures come out **bimodal** — a cluster of temperate planets, a cluster of frozen ones, and a **forbidden gap around 235–260 K containing almost nothing**. Nobody sorted those planets into bins. They were run independently, from orbital elements produced by a completely separate N-body simulation, and they landed in two groups because the underlying equations admit two stable states and very little in between.

You are looking at the same cliff, seen from 644 different directions.

---

## Going back up

If dimming the Sun by 12% freezes the planet, then brightening it by 12% should thaw it. Symmetry.

Let's check. Start from a fully glaciated world and ramp `s` back **up**:

| s | Global mean T | Ice edge |
|---|---|---|
| 0.87 | 229.71 K | ice everywhere |
| 0.95 | 234.66 | ice everywhere |
| 1.00 | 237.75 | ice everywhere |
| 1.11 | 244.56 | ice everywhere |
| 1.21 | 250.74 | ice everywhere |
| **1.27** | **254.46** | **ice everywhere** |
| **1.29** | **321.87** | **ice-free** |

It does not thaw at 0.88. It does not thaw at 1.00. **It stays frozen until the Sun is 29% brighter than today**, then leaps 67 K in a single step to an ice-free hothouse.

Freezing needs `s ≈ 0.87`. Thawing needs `s ≈ 1.29`. **The path down and the path back are not the same path.** That's **hysteresis**, and the gap between the two thresholds is enormous — 41 percentage points of solar output.

The physical reason is straightforward once stated: a frozen planet is *white*. It reflects most of what arrives, so you must supply drastically more sunlight to warm it past freezing than you needed to remove to freeze it in the first place. The system defends whichever state it's in.

### What this means for Earth right now

Look again at the row for `s = 1.00`.

At **exactly today's sunlight**, this model has **two** stable answers. A temperate Earth at 288 K with ice edge near 71°, and a frozen Earth at 238 K with ice to the equator. Both are genuine equilibria. Both persist indefinitely. The physics does not prefer one.

**Which one you get depends only on which one you started from.** Here is the same present-day configuration reached from four different initial states:

| Starting state | Result |
|---|---|
| Uniform +15 °C | 288.13 K, ice edge 70.8° |
| Earth-like profile | 288.13 K, ice edge 70.8° |
| Uniform −30 °C | **237.75 K, fully glaciated** |
| Uniform −60 °C | **237.75 K, fully glaciated** |

Same equations, same forcing, same solar constant. Two destinations.

This is not a quirk of our implementation — it's the classic Budyko–Sellers result from 1969, and it's part of why the geological evidence for Snowball Earth episodes was resisted for so long. The problem was never getting into one. It was getting out.

---

## An honest correction

I want to show you something from this project's own development record, because it's a better lesson than any I could invent.

Report §7.1 asked whether the model's equilibria depend on initial conditions. It tested 12 configurations straddling the bifurcation, running each from two different starting states, and found **zero** cases where the outcome differed. It concluded the results were robust to initial conditions.

You have just seen, on this page, that the model is *massively* bistable at present-day forcing.

**Both statements are true, and here's the resolution.** The two initial conditions in that test were "uniform 15 °C" and "Earth's equilibrium profile." Both are *warm* states. Both sit on the warm branch. The test compared two points on the same side of the divide and, unsurprisingly, found they landed in the same place.

The bistability is real. The test simply wasn't sensitive to it. To find the cold branch you have to *start* on the cold branch — as the −30 °C and −60 °C rows above do.

Nothing about that earlier conclusion was dishonest, and the arithmetic was correct. It was **under-tested in a way that wasn't obvious until you knew what to look for**. That is the ordinary condition of computational science, and the reason results get re-examined rather than filed.

Ask of any negative result: *what would this test have failed to detect?*

---

## What we should not claim

The bifurcation's **existence** is robust. It follows from a positive feedback overwhelming a negative one, it appears across the whole Budyko–Sellers literature, and it does not depend on our particular numbers.

The bifurcation's **location** is not. It sits where it does because we drew a hard step at `T_ice = −10 °C` and dropped the coalbedo to exactly `0.38`. Real sea ice doesn't appear discontinuously at a threshold temperature; its reflectivity depends on thickness, snow cover, melt ponds, and the angle of the Sun. Substitute a gradual transition and the cliff softens and shifts.

So: **"this model has a tipping point" is a defensible claim. "Earth tips at 88% of present sunlight" is not.** The first is physics; the second is our parameterisation.

Report §9.2 records this limitation explicitly, and it's worth internalising the general form: *know which digits of your result are physics and which are your choices.*

---

## What you now know

- Making coalbedo temperature-dependent — one conditional — makes the model **nonlinear**, and everything here follows from that.
- **Ice-albedo is a positive feedback**: cooling → ice → higher albedo → less absorption → more cooling. "Positive" means self-reinforcing, not good.
- Gradual dimming produces gradual cooling until it doesn't. Between `s = 0.88` and `0.87` the planet drops **37 K** and freezes completely.
- **Nobody programmed that collapse.** It is emergent, not imposed — the distinction that separates a prediction from a restatement.
- The same signature appears as a **bimodal outcome distribution with a forbidden gap** across 644 independent simulated worlds.
- **Hysteresis:** freezing needs `s ≈ 0.87`, thawing needs `s ≈ 1.29`. A frozen planet defends its state.
- **At present-day sunlight the model is bistable** — 288 K or 238 K, decided entirely by where you started.
- A test that finds no initial-condition sensitivity may simply have probed one branch twice.
- The bifurcation's **existence** is robust physics; its **location** is a parameterisation choice.

---

## Exercise 4 — Find the cliff

**Provide:** the full single-surface EBM with ice feedback enabled.

**Part A.** Scale `S₀` by a factor `s`. Sweep `s` downward from 1.0 in steps of 0.01, recording the equilibrium global mean temperature at each. **Important:** start each run from the previous run's final state, not from a fixed guess — you are modelling a planet experiencing gradual change, and that requires memory.

**Part B.** Identify the critical `s` at which the planet collapses to a snowball. How large is the temperature jump across that single 1% step, compared to the steps either side of it?

**Part C.** Now sweep *upward* from a fully glaciated state. Does it de-glaciate at the same `s`? Report both thresholds.

**Part D.** Plot ice-edge latitude against `s` for both sweep directions on one figure, with arrows showing the direction of travel. You should get the classic two-branch loop.

**Part E.** At `s = 1.00`, run the model from four starting states: uniform +15 °C, uniform 0 °C, uniform −30 °C, uniform −60 °C. How many distinct equilibria do you find? Locate the approximate starting temperature that separates the two outcomes.

**Discussion.** If Earth had ever fallen into a snowball state, what could possibly have got it out? Our model says it would need a 29% brighter Sun — which never happened. So either the snowball episodes in the geological record didn't occur, or something is missing from this model.

*Something is missing. Volcanoes continue to emit CO₂ regardless of surface temperature, but the weathering reactions that normally remove it require liquid water — which a frozen planet has none of. So CO₂ accumulates for millions of years until the greenhouse effect overwhelms the ice. This is the **carbonate–silicate thermostat**, it operates on ~10⁶-year timescales, and it is listed in report §9.6 as absent from our model. Nothing in this course will ever thaw a snowball, and now you know why.*

---

**Next lecture:** Act I is complete — we have a planet with latitude, seasons, transport, and ice. Every bit of it has assumed a circular orbit. Real orbits are ellipses, which means the distance to the Sun varies through the year, which means the `1/r²` we've been quietly setting to 1 starts to matter. We'll need to solve an equation that has no closed-form solution and has been bothering people since Kepler wrote it down in 1609 — and then discover that eccentricity does something to the annual-mean sunlight that most people guess backwards.
