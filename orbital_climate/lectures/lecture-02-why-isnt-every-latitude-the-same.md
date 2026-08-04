# Lecture 2 — Why isn't every latitude the same?

*You'll need trigonometry — sines, cosines, and the inverse functions. No physics beyond Lecture 1. In fact there's no new physics in this lecture at all, which is the point.*

---

## What we broke last time

Last lecture we computed Earth's temperature and got a single number. One planet, one temperature, 255 K.

But you know that's false. The Sahara and Antarctica are not the same temperature. July and January are not the same temperature. Our model has no way to express any of that, because we collapsed the entire planet into a point.

So this lecture adds **geometry**. And here is the striking thing: we will add no physics whatsoever. No new constants, no new mechanisms, nothing about atmospheres. Just the shape of the Earth, the tilt of its axis, and where the sunlight lands.

That alone will produce seasons. It will produce polar night. And it will produce a result that most people flatly refuse to believe the first time they see it.

---

## Kill this idea first

Before we start, let's destroy a misconception that would otherwise quietly wreck your understanding:

> **"Seasons happen because Earth is closer to the Sun in summer."**

No. This is wrong, and it's wrong in a way you can check yourself.

Earth is actually *closest* to the Sun in early January — northern winter. If distance drove seasons, January would be the northern hemisphere's summer. It isn't.

More decisively: **everything in this lecture assumes a perfectly circular orbit.** The distance to the Sun will be constant. And we will still get seasons — strong ones. If distance were the cause, our circular-orbit model would have no seasons at all.

Seasons are caused by **axial tilt**. That's the whole story, and by the end of this lecture you'll be able to derive it.

*(Distance isn't entirely irrelevant — it's a modest modulation, and in Lecture 5 we'll put the eccentricity back in and see exactly how much it matters. But it is not the cause.)*

---

## Three angles

We need three quantities. Take them slowly; each is simple, and the combination is where the interest lies.

### 1. Where Earth is in its orbit — the solar longitude λ

As Earth travels around the Sun, we need a way to say "how far around". Call that angle the **solar longitude**, λ.

$$ \lambda = \nu + \lambda_p $$

Here ν (nu) is the **true anomaly** — the angle measured from perihelion, the closest point to the Sun — and λ_p is the **longitude of perihelion**, which says where perihelion sits relative to the seasons.

For now, with a circular orbit, ν just increases steadily: 360° over one year, about one degree per day. The machinery for computing ν on a genuinely elliptical orbit is Lecture 4's business, and it's more interesting than you'd expect.

Set λ = 0 at the northern **vernal equinox** — around 20 March. Then λ = 90° is the June solstice, 180° the September equinox, 270° the December solstice.

### 2. How high the Sun rides — the declination δ

Earth's axis is tilted by an angle ε (epsilon), the **obliquity**, currently **23.44°**. Crucially, that axis points in a *fixed* direction in space while Earth orbits. It does not swing around to follow the Sun.

That fixity is what creates seasons. For half the year the northern end leans sunward; for the other half it leans away.

The **declination** δ is the latitude at which the Sun is directly overhead at noon:

$$ \delta = \arcsin(\sin \varepsilon \cdot \sin \lambda) $$

Check the extremes, because they should feel familiar:

| λ | δ | What it is |
|---|---|---|
| 0° | 0° | March equinox — Sun overhead at the equator |
| 90° | **+23.44°** | June solstice — Sun overhead at the Tropic of Cancer |
| 180° | 0° | September equinox |
| 270° | **−23.44°** | December solstice — Tropic of Capricorn |

**The Tropics are not an arbitrary line drawn by geographers.** They are exactly ±ε — the furthest the overhead Sun ever wanders. If Earth's tilt changed, the Tropics would move. In some of the orbits we'll study later, they move a *lot*.

Note that δ can never exceed ε. The obliquity is a hard ceiling on the declination, which is worth remembering when we get to planets tilted at 74°.

### 3. How long the Sun is up — the hour angle H₀

Now fix a latitude φ (phi) and ask: on this day, what fraction of the 24 hours is the Sun above the horizon?

The answer comes from spherical trigonometry, and I'll state it rather than derive it:

$$ H_0 = \arccos(-\tan \varphi \cdot \tan \delta) $$

H₀ is the **sunset hour angle** — half the length of daylight, expressed as an angle. H₀ = π means 24-hour daylight; H₀ = 0 means the Sun never rises; H₀ = π/2 means a twelve-hour day.

Test it at the equator, φ = 0. Then tan φ = 0, so H₀ = arccos(0) = π/2. Twelve hours of daylight, every single day of the year, regardless of season. That is correct, and it is why equatorial regions have no meaningful seasons.

---

## The clamp — the most important paragraph in this lecture

Look hard at that arccos.

Its argument is `−tan φ · tan δ`. Near the poles, tan φ becomes enormous. Push φ to 80° in northern summer and that product sails past 1.

**And arccos of anything greater than 1 is undefined.** Feed it to a calculator and you get an error. Feed it to Python and you get `nan`, which then silently contaminates everything downstream.

So what does it *mean* physically? It means there is no sunset. The Sun circles the sky without ever touching the horizon. You've hit the midnight Sun, and the equation is telling you that the question "when does the Sun set?" has no answer.

The naive fix is to write branching logic:

```python
if -tan(phi)*tan(delta) > 1:
    H0 = 0.0          # polar night
elif -tan(phi)*tan(delta) < -1:
    H0 = pi           # polar day
else:
    H0 = arccos(-tan(phi)*tan(delta))
```

Don't. Do this instead:

```python
H0 = arccos(clip(-tan(phi)*tan(delta), -1.0, 1.0))
```

**Clamping the argument to [−1, 1] gives exactly the right answer on both branches.** When the product exceeds +1, the clamp returns arccos(1) = 0 — polar night. When it drops below −1, arccos(−1) = π — polar day. The two "special cases" are simply the two saturated ends of one continuous formula.

This is worth naming as a habit of mind, because you'll use it for the rest of your career:

> **Design the edge case out of existence rather than testing for it.**

Branching logic is where bugs live. Every `if` is a place where one path gets exercised and the other rots quietly for two years until the day it runs. A clamp has no branches to get wrong. It's also faster, and it works on whole arrays at once without any special handling.

The polar day/night tests in this model's test suite exist precisely to lock this behaviour down — not because anyone doubts arccos, but because someone will eventually be tempted to "fix" the clamp.

---

## Putting it together: the daily-mean insolation

We want the energy arriving per square metre of horizontal ground, averaged over a full 24-hour rotation. Combining the pieces:

$$ Q = \frac{S_0/\pi}{r^2}\left( H_0 \sin\varphi \sin\delta + \cos\varphi \cos\delta \sin H_0 \right) $$

The bracket does two jobs at once, and it's worth seeing them separately:

- **`H₀ sin φ sin δ`** — how *long* the Sun is up. It's positive when your hemisphere is tilted sunward and negative when it isn't.
- **`cos φ cos δ sin H₀`** — how *high* the Sun climbs. Light striking at a slant is spread over more ground.

Seasonal contrast is these two effects *reinforcing* each other. Summer gives you both a longer day and a higher Sun. Winter takes both away. That compounding is why seasons are as strong as they are.

The `r²` is the distance to the Sun. For this lecture r = 1 always, since we assumed a circular orbit. It returns in Lecture 5.

### One check you can do by hand

Put yourself at the equator at equinox: φ = 0, δ = 0. Then sin φ = 0, cos φ = 1, cos δ = 1, and H₀ = π/2, so sin H₀ = 1. Everything collapses:

$$ Q = \frac{S_0}{\pi}\left(0 + 1 \cdot 1 \cdot 1\right) = \frac{S_0}{\pi} = \frac{1361}{\pi} \approx 433 \ \text{W/m}^2 $$

Clean, exact, and easy to verify in code. **Any implementation that fails this test is broken**, and it's the first thing to check when your plot looks strange.

Notice this is *larger* than the 340 W/m² from Lecture 1. No contradiction: 340 was the average over the whole planet including the night side and the slanted high latitudes. This is one specific well-favoured place.

---

## The result nobody believes

Now predict, before computing. On the June solstice, where on Earth does the most sunlight land per square metre per day?

Most people say the Tropic of Cancer, since that's where the Sun is directly overhead.

Run the numbers and you get this:

| Location, June solstice | Daily-mean insolation |
|---|---|
| **North Pole** | **541 W/m²** |
| Equator | 398 W/m² |

**The North Pole. By a margin of 36%.**

It is the single sunniest place on the planet, and it beats the equator decisively.

The reason is in the two terms we separated. The pole loses badly on Sun *height* — the Sun crawls around at 23.44° elevation, never higher. But it wins overwhelmingly on *duration*: **24 hours of continuous daylight**, against the equator's unvarying twelve. Length beats angle.

And this is not a curiosity. It is the mechanism behind ice ages. Whether ice sheets grow or melt depends on whether high-latitude snow survives the summer, which depends on exactly this quantity. In Lecture 6 we'll watch 65 °N summer insolation fall from 480 to 400 W/m² under a perturbed orbit — and that 17% drop is enough to start a glaciation.

If polar summer insolation feels like an obscure number, it isn't. It's the number.

---

## Turning the tilt knob

The clearest way to see that tilt *is* the seasons is to change it.

**Set ε = 0.** The declination becomes `asin(0 · sin λ) = 0` — permanently. The Sun sits over the equator all year, every day is twelve hours everywhere, and the insolation pattern is **completely uniform in time**. No seasons anywhere on the planet, ever, despite Earth orbiting exactly as before. Latitude still matters enormously — the poles are dim and the tropics bright — but nothing changes through the year.

**Set ε = 90°**, tipping the axis fully onto its side. Now the poles alternate between pointing straight at the Sun and straight away. Seasons become violent, and the annual-mean insolation at the pole *exceeds* the equator's:

| Obliquity | Pole (annual mean) | Equator (annual mean) |
|---|---|---|
| 23.44° (today) | 172 W/m² | 416 W/m² |
| **54°** | **350** | **350** ← crossover |
| 74° | 416 | 299 |
| 90° | 433 | 276 |

Somewhere near **54°** the ranking flips, and the poles become the warmest part of the planet in the annual mean. That's not a hypothetical: the black-hole flyby simulations this course is built around contain Earth-like planets left at **74° obliquity**, comfortably past the crossover. Those worlds have frozen equators and temperate poles.

Uranus, incidentally, sits at 98°.

---

## What you now know

- **Seasons come from axial tilt, not distance.** A circular orbit still has full seasons.
- **Solar longitude** λ says where in the orbit we are; **declination** δ = asin(sin ε sin λ) says which latitude has the Sun overhead.
- δ can never exceed ε — the Tropics sit at exactly ±ε.
- The **sunset hour angle** H₀ = arccos(−tan φ tan δ) gives day length.
- **Polar day and polar night are not special cases.** They're the saturated branches of one formula, and clamping the arccos argument to [−1, 1] handles both with no branching.
- Daily-mean insolation combines **day length** and **Sun height**, which reinforce each other seasonally.
- At the equator at equinox, Q = **S₀/π ≈ 433 W/m²** exactly — your correctness check.
- **The summer pole is the sunniest place on Earth**, beating the equator by 36%.
- Past ~54° obliquity, the poles beat the equator in the *annual* mean too.

---

## Exercise 2 — The insolation surface

**Provide:** `orbital_climate/insolation.py`, circular orbit (`ecc = 0.0`).

**Part A.** Produce a filled contour plot of Q(φ, day) over a full year, with latitude from −90° to +90° on one axis and day of year on the other. Use a sequential colour scale — this is a magnitude, so light-to-dark of a single hue, not a rainbow.

**Part B.** Mark the latitude and day of the global maximum on your plot. Where is it? Was that what you expected before plotting?

**Part C.** Verify numerically that the equator-at-equinox value equals S₀/π. Get it exact, not approximate — if it isn't, something in your pipeline is wrong and you want to know now rather than in Lecture 6.

**Part D.** Set ε = 0 and re-plot. Then set ε = 90° and re-plot. Describe in one sentence each what happened to the seasonal structure, and explain both using the declination formula alone.

**Part E — the extension.** At what obliquity does the *annual-mean* insolation at the pole first exceed that at the equator? Find it numerically to within a degree.

Then consider: the simulation sweep behind this course contains runs where Earth was left at 74° obliquity by a passing black hole. Using your answer, describe qualitatively what that planet's climate zones look like. Which latitudes freeze?

*Part E is the one to spend time on. It's the first point in the course where you'll compute something about a world that doesn't exist, using tools validated on one that does.*

---

## A note on what we've built

Stop and notice something. We have not written down a single new physical law this lecture. Everything came from three angles and the geometry of a sphere.

Yet we now have seasons, polar night, tropics, and a quantitative prediction about ice ages. **Geometry is doing an enormous amount of work here** — and it will keep doing so. When we get to orbital perturbations in Lecture 5, the changes that matter most will again be geometric: not "how much energy arrives" but "when and where it arrives."

---

**Next lecture:** we have insolation everywhere on the planet, and from Lecture 1 we know how to turn absorbed energy into a temperature. Time to connect them. When we do, the answer comes out badly wrong — a tropics far hotter than any place on Earth, poles far colder, and an equator-to-pole gradient more than twice the real one. The missing ingredient is that energy does not stay where it lands: winds and ocean currents carry it polewards, and we'll need some way to represent that. We'll also discover something that surprises almost everyone about what all that heat transport does to the planet's *average* temperature.
