# Lecture 1 — Why is Earth 288 Kelvin?

*No prior climate science assumed. You'll need to be comfortable with algebra and with the idea of a fourth power. That's genuinely it.*

---

## A question you can't answer yet

Earth's average surface temperature is about 15 °C — 288 Kelvin, if we use the scale that starts at absolute zero. We'll use Kelvin throughout, because the physics of radiation only works in Kelvin.

Here's the question: **why that number?**

Not "why is it comfortable" — that's backwards, we evolved to find it comfortable. The question is what physical process picks 288 rather than 200 or 400. Something must be setting it.

By the end of this lecture you'll be able to calculate that number from scratch. You'll get it *wrong* on your first attempt, by about 33 degrees, and that failure will turn out to be the most interesting thing in the lecture.

---

## The only idea in this course

A planet sits in space. Sunlight arrives and warms it. The planet, being warm, glows — not visibly, but in infrared — and that glow carries energy away.

If more arrives than leaves, the planet warms up. But a warmer planet glows *more brightly*, so it loses energy faster. If more leaves than arrives, the planet cools, glows more dimly, and loses energy more slowly.

Either way it settles at the temperature where the two exactly balance:

> **Energy in = Energy out.**

That's it. That's the whole subject. Everything in the next nine lectures is a refinement of that one sentence — adding latitude, seasons, ice, oceans, orbits — but the sentence never changes, and every time you get lost, come back to it.

Let's make both sides concrete.

---

## Energy out: how brightly does a warm object glow?

Every object above absolute zero radiates energy. A hot object radiates more. The relationship isn't gentle — it's a fourth power:

$$ \text{energy radiated per square metre} = \sigma T^4 $$

where σ (sigma) is a constant of nature, the Stefan–Boltzmann constant:

$$ \sigma = 5.67 \times 10^{-8} \ \text{W m}^{-2}\text{K}^{-4} $$

The units say what it does: give it a temperature in Kelvin, get back watts per square metre.

**Sit with the fourth power for a moment.** Double an object's temperature and it radiates sixteen times as much. This is the reason planetary temperatures are stable rather than runaway. If a planet warms slightly, its energy loss shoots up steeply and drags it back. The fourth power is a very stiff spring.

For Earth at 288 K:

$$ \sigma T^4 = 5.67 \times 10^{-8} \times 288^4 \approx 390 \ \text{W/m}^2 $$

Hold onto 390. It comes back.

---

## Energy in: the awkward factor of 4

At Earth's distance, the Sun delivers about **1361 W/m²**. This number is called the *solar constant*, S₀. It's measured by satellites and it's the amount of power crossing one square metre held face-on to the Sun.

Now, the trap. You might think each square metre of Earth receives 1361 W. It doesn't, and the reason is geometry.

**The catching side.** From the Sun's point of view, Earth is a flat circular target of radius R. The sunlight it intercepts is whatever passes through a disc of area:

$$ \pi R^2 $$

So total power intercepted = $S_0 \times \pi R^2$.

**The radiating side.** But Earth doesn't radiate from a disc. It radiates from its entire spherical surface — day side, night side, poles, everywhere — with area:

$$ 4\pi R^2 $$

So when we ask "how much solar power per square metre of *planet surface*", we divide the intercepted power by the whole sphere:

$$ \frac{S_0 \pi R^2}{4 \pi R^2} = \frac{S_0}{4} = \frac{1361}{4} \approx 340 \ \text{W/m}^2 $$

**A warning, because nearly everyone makes this mistake:** the 4 is *not* because half the planet is in darkness. That would be a factor of 2. The 4 comes from the ratio of a sphere's surface area to its shadow's area — it accounts for night *and* for the fact that even in daylight, sunlight strikes most of the surface at a slant rather than head-on. A square metre at high latitude in the afternoon catches far less than a square metre with the Sun directly overhead.

Notice also that R cancelled. **A planet's size doesn't affect its temperature at all.** A bigger planet catches more sunlight, but it has proportionally more surface to radiate from. This surprises most people, and it's worth pausing on — it's the first time in this course that an obviously-relevant variable turns out to be irrelevant. It won't be the last.

---

## Not all of it sticks

Some sunlight bounces straight back to space without ever warming anything — off clouds, off ice, off bright desert. The fraction reflected is called the **albedo**, α.

For Earth, α ≈ 0.30. Nearly a third of incoming sunlight is simply returned, and clouds do most of that work.

So the energy actually absorbed, per square metre of planet:

$$ \text{absorbed} = \frac{S_0}{4}(1 - \alpha) = 340 \times 0.70 \approx 238 \ \text{W/m}^2 $$

---

## Your first calculation

Now set the two sides equal:

$$ \frac{S_0}{4}(1-\alpha) = \sigma T^4 $$

and solve for T:

$$ T = \left[ \frac{S_0(1-\alpha)}{4\sigma} \right]^{1/4} $$

Put the numbers in:

$$ T = \left[ \frac{238}{5.67 \times 10^{-8}} \right]^{1/4} = \left[ 4.20 \times 10^{9} \right]^{1/4} \approx 255 \ \text{K} $$

**255 K. That's −18 °C.**

You just calculated a planet's temperature from three numbers and one idea. That's a real achievement, and you should be pleased with it.

You should also notice that it's wrong.

---

## The 33-degree hole

Observed: 288 K. Calculated: 255 K. We are short by **33 K**.

This is not a rounding error or a sloppy constant. It's a third of the way to freezing the planet solid. Something structural is missing.

And here's the thing worth appreciating: our calculation isn't wrong so much as *incomplete in a specific, locatable way*. Remember we said Earth radiates 390 W/m² at 288 K, but only absorbs 238 W/m²? Those numbers don't balance, and they're both correct. The surface really does radiate 390. The planet really does only absorb 238.

The resolution is that **the surface is not what space sees.**

---

## What's actually going on

Certain gases in the atmosphere — water vapour, carbon dioxide, methane — are nearly transparent to incoming sunlight but strongly absorb outgoing infrared. So the infrared radiated by the ground doesn't escape directly. It's absorbed a few metres up, re-radiated in all directions, absorbed again, and so on, stumbling upward through the atmosphere.

Radiation only escapes freely from high up, where the air is thin enough that there's nothing left above to absorb it. And high up, the air is *cold*.

So space doesn't see a 288 K surface. It sees an effective emitting layer several kilometres up, at roughly 255 K — which is exactly the number we calculated. **Our calculation was right about the wrong surface.**

The atmosphere cools with height at about 6.5 K per kilometre. If the emission level sits around 5 km, the ground beneath it must be roughly 33 K warmer than the emission level. That's our missing 33 K, and it's why the greenhouse effect is better described as *raising the altitude of escape* than as "trapping heat."

**A common misconception to discard now:** the greenhouse effect is not a blanket that stops energy leaving. In steady state, exactly as much energy leaves as arrives — it must, or the planet wouldn't be in steady state. What changes is the *temperature at which* the escaping radiation is emitted, and therefore how warm the surface underneath has to be.

---

## The decision that shapes this entire course

We now face a choice. To model climate properly, we need to know how much infrared escapes to space at a given surface temperature. Call this the **outgoing longwave radiation**, or OLR.

We could calculate it honestly. That means tracking hundreds of thousands of individual absorption lines for each gas, at every altitude, at every wavelength. It's called line-by-line radiative transfer, it's a career, and it's how serious climate models do it.

We're not going to do that. Instead we're going to cheat, deliberately and openly.

Look at the observations: measure the surface temperature at many places on Earth, measure the outgoing infrared above each, and plot one against the other. Over the range of temperatures Earth actually exhibits, the relationship comes out **very nearly a straight line**. So we write:

$$ \text{OLR} = A + B \cdot T $$

with, for Earth,

$$ A = 203.3 \ \text{W/m}^2, \qquad B = 2.09 \ \text{W/m}^2/°\text{C} $$

(Note T here is in **degrees Celsius**, not Kelvin — this parameterisation is conventionally written that way, and mixing them up is a classic first-week bug.)

Check that it works. At 15 °C:

$$ 203.3 + 2.09 \times 15 = 234.7 \ \text{W/m}^2 $$

Compare to the 238 W/m² we said Earth absorbs. Close — within a couple of watts, which is about the accuracy of the constants themselves.

**What have we actually done here?** We've replaced all of atmospheric radiative physics — every molecule, every spectral line, the entire vertical structure — with two numbers fitted to observations. In exchange we get a model that runs in seconds instead of hours and that you can fully understand.

Notice the two-way trade. Compare `σT⁴`, which is a law of physics valid at any temperature, against `A + B·T`, which is a line fitted through a cloud of data points clustered around 288 K. The first is true everywhere. The second is true *near where it was fitted*, and makes no promises elsewhere.

---

## A promise about how this will fail

I want to flag something now, in week one, so that when it happens in Lecture 9 you recognise it rather than being blindsided.

Later in this course we'll use this model to explore planets on wildly different orbits — some much closer to their star, some much further. Some will come out at 560 K.

At 560 K, `A + B·T` gives roughly 800 W/m² of outgoing radiation. Actual blackbody physics at that temperature demands around 5500 W/m².

The model will be wrong by a factor of seven. And it will not warn you. It will return `559.7 K` with a decimal place, looking every bit as authoritative as the numbers that are correct.

**This is the most important thing in the lecture, and it isn't about climate at all.** A model gives you numbers whether or not it's entitled to. Knowing where your approximations were fitted — and refusing to trust results from outside that range — isn't pedantry. It's the difference between doing science and generating plausible digits.

Every model you ever build will have a version of this line drawn around it. Your job is to know where the line is.

---

## What you now know

- **Energy in = energy out.** Everything else is refinement.
- A warm object radiates as **σT⁴**. The fourth power is what makes planetary temperature stable.
- The **factor of 4** comes from sphere-versus-disc geometry — not from night.
- Planetary radius **cancels**. Size doesn't set temperature.
- **Albedo** (α ≈ 0.30) removes a third of the incoming energy before it can warm anything.
- These give **255 K** for Earth, against an observed **288 K**.
- The **33 K gap** exists because space sees a cold, high emission layer rather than the warm surface.
- We replace all of that physics with **OLR = A + B·T**, a straight line fitted near 288 K.
- That line will produce confident, precise, completely fictitious numbers far from 288 K.

---

## Exercise 1 — A zero-dimensional planet

You need nothing but a calculator or a few lines of Python.

**Part A.** Compute the effective temperature for:

| Planet | S₀ (W/m²) | albedo α | observed surface T |
|---|---|---|---|
| Venus | 2601 | 0.77 | 737 K |
| Earth | 1361 | 0.30 | 288 K |
| Mars | 586 | 0.25 | 210 K |

**Part B.** For each, compute the greenhouse gap (observed minus calculated). Rank the planets by gap size.

**Part C.** Venus has a very high albedo — it reflects 77% of incoming sunlight, more than any other planet in the solar system. Yet it's by far the hottest. Explain how both facts are true at once.

**Part D.** Using `OLR = A + B·T` with T in °C, find the value of absorbed radiation that yields exactly 15.0 °C. Then work backwards: what albedo would Earth need for the simple `S₀(1−α)/4` calculation to deliver that value? Compare your answer to the 0.30 we used, and to the value 0.676 that appears in the model you'll meet in Lecture 3 as a *co*-albedo (the fraction absorbed, i.e. 1 − α).

**Part E — the one that matters.** You found Venus's greenhouse gap in Part B. It's roughly 500 K. Our linear OLR was fitted to observations spanning maybe 250–300 K on Earth. Write two or three sentences on what you would expect if someone used Earth's `A` and `B` values to model Venus — and on how they would know something had gone wrong.

*Part E has no single right answer. It's asking you to reason about the limits of a tool, which is the skill this entire course is built around.*

---

**Next lecture:** we've treated Earth as a single point with one temperature. But the tropics are hot, the poles are frozen, and there are seasons. None of that is in what we've built. Next time we add geometry — and discover, before doing any physics at all, that the sunniest place on Earth is the summer pole.
