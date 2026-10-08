# Flat Panel Speakers: Physics, Topology & Tools

*From violin bass bars to finite element analysis — a self-contained lesson, and the reasoning behind Panel Speaker Designer.*

---

## 01 — Foundation: The Problem Every Acoustic Surface Shares

Whether you're looking at a violin top plate, a speaker cone, or a flat panel on a wall, the fundamental engineering problem is the same: a transducer applies force to a surface, and that surface must move in a controlled way to radiate sound. The question is always *how does energy travel across the surface*, and *how do you shape that travel*?

### The Violin as a Teaching Model

A violin has two internal structures that solve this problem. The **soundpost** — a tight spruce dowel wedged between the top and back plates — sits under the treble (E-string) foot of the bridge and acts as a near-rigid coupling point. The **bass bar** — a spruce strip glued along the inside of the top — runs under the bass (G-string) foot.

Together they create a deliberate asymmetry. The soundpost side is stiff; the bass bar side is free to flex, but in a controlled way. The bar runs parallel to the wood grain, and spruce is dramatically stiffer along the grain than across it — so the bar efficiently channels vibration longitudinally down the plate rather than letting it dissipate locally. The whole lower bout vibrates as a coherent unit instead of just the patch under the bridge foot.

> **Core Principle:** A stiffener doesn't just reinforce — it shapes the *path* that vibrational energy takes. Where energy goes determines what resonates, and what resonates determines what you hear.

### Conventional Speakers: Suppressing Breakup

A pistonic speaker cone tries to be a perfect rigid piston — all parts moving together. At higher frequencies the cone "breaks up," meaning different regions flex out of phase and introduce coloration. Cone geometry, radial ribs, and dust cap design all push that breakup point as high as possible. The goal is to suppress modal behavior. Radial ribs, like a bass bar, stiffen along a specific axis to make the surface behave as a more unified whole.

---

## 02 — Distributed Mode Loudspeakers (DML)

Flat panel speakers — developed commercially under the NXT/DML banner in the late 1990s — invert the pistonic philosophy entirely. Instead of suppressing resonant breakup modes, they *deliberately excite as many as possible*, densely enough that they statistically average into something approaching flat response. The panel is never trying to move as a piston. It's a chaotic but managed resonant system.

### Why Exciter Placement Matters So Much

Every resonant mode of a panel has a characteristic shape — regions that move a lot (antinodes) and lines where movement is zero (node lines). If you place your exciter on a node line of a given mode, that mode doesn't get driven at all. It's absent from the output.

Center placement on a rectangular panel is nearly the worst possible choice. A rectangle's modes fall into four symmetry families: symmetric or antisymmetric about the vertical centre line, and symmetric or antisymmetric about the horizontal one. Every mode that is antisymmetric about a centre line has a node line along it. So anywhere on one centre line misses about half the modes, and the centre — on both lines — misses about three quarters. **Almost any informed placement beats center placement.** In Panel Speaker Designer's heat map this shows up as a dark cross along both centre lines.

```
BAD: center placement          BETTER: offset placement       OPTIMAL: scored placement

┌─────────────┐                ┌─────────────┐                ┌─────────────┐
│      │      │                │      │      │                │  ╌╌╌│╌╌╌╌╌ │
│      │      │                │      │      │                │╌╌╌╌╌│╌╌╌╌╌╌│
│──────●──────│                │──────│──────│                │  ╌╌╌│╌╌ ●  │
│      │      │                │      │  ●   │                │╌╌╌╌╌│╌╌╌╌╌╌│
│      │      │                │      │      │                │  ╌╌╌│╌╌╌╌╌ │
└─────────────┘                └─────────────┘                └─────────────┘
Center sits on both            Offset avoids the              Scored to avoid
centre-line node lines         centre lines                   many modes' node lines
```

### The Role of Panel Material

Bending wave speed in a panel depends on stiffness, density, and frequency. The dispersion relation — how wave speed varies with frequency — determines where modes fall. For an **isotropic** material (same properties in all directions), this is a single clean equation. For **orthotropic** materials (wood, plywood, many composites), bending stiffness differs along each axis, and the mode structure becomes directionally asymmetric.

This matters more than it sounds. Balsa is roughly 30× stiffer along the grain than across it, so a balsa panel's modes look nothing like an acrylic panel's of the same size. Panel Speaker Designer takes stiffness along the width (`E_x`), along the height (`E_y`) and in shear (`G`, which governs twisting modes); for wood, put the grain along the width.

---

## 03 — Level One: What You Can Calculate Analytically

For a rectangular panel, the modes can be calculated without a mesh. This is the right place to start, and it's where this project started.

### The Simple Model (and Why It Isn't Enough)

The first version of the software used the textbook shortcut: treat each free-edge mode as `cos(mπx/Lx)·cos(nπy/Ly)` and give it the frequency of a simply-supported plate.

```python
# Simplified sketch of the original model. Easy to write, but see below.

import numpy as np

def mode_frequency(m, n, Lx, Ly, h, E, rho, nu):
    D = (E * h**3) / (12 * (1 - nu**2))          # bending stiffness
    kx, ky = m * np.pi / Lx, n * np.pi / Ly
    omega = np.sqrt(D / (rho * h)) * (kx**2 + ky**2)
    return omega / (2 * np.pi)                    # Hz

def score_position(x, y, Lx, Ly, modes):
    # Higher score = position moves more in more modes
    return sum(abs(np.cos(m * np.pi * x / Lx) * np.cos(n * np.pi * y / Ly))
               for m, n in modes)
```

A brute-force grid search over positions gives a heat map, and that heat map already clears the "beat centre placement" bar. But the shortcut gets free edges wrong. A free edge is not where a cosine peaks; real free-edge mode shapes bend differently near the edges, and their frequencies come from different wavenumbers (a free-free beam's first bending mode uses about 1.506π/L, not 2π/L). On a 300 × 200 × 3 mm acrylic panel, the shortcut puts the first two modes at 87 and 107 Hz; the correct values are 52 and 57 Hz. For balsa, where the shortcut also ignores the grain, it is off by 4–10×.

### What the App Does Instead

**Rayleigh-Ritz.** The panel's bending and kinetic energies are written out in full (including the grain-dependent stiffnesses) and minimised over a family of smooth trial shapes (Legendre polynomials). The result is the true free-edge modes, matching Leissa's published values for a free square plate to within 0.5%. Simply-supported edges have an exact closed form, which the app uses directly.

**The exciter isn't a point.** An exciter drives the panel around its voice-coil ring (typically 20–30 mm across), so its coupling to a mode is the mode's motion averaged around that ring. Modes smaller than the coil partly cancel themselves out, which a point model can't see.

**Edges are excluded.** Free edges move more than anywhere else on a panel, and they're impractical mounting points, so the search ignores a margin of 10% of the shorter side.

---

## 04 — Scoring Placement: Response Flatness

Counting how much a position moves in each mode was the first scoring rule. It has a blind spot: every mode counts the same, so a position that drives ten crowded modes around 2 kHz but leaves a hole at 300 Hz can still score well. What a listener hears is the *response* — how strongly the panel vibrates at each frequency.

### From Modes to a Response

Drive the panel with a steady force at one point and each mode responds like a damped resonator: strongly near its own frequency, weakly elsewhere, with the peak width set by the material's damping (the loss factor η). The panel's total vibration at any frequency is the sum of all those contributions, each weighted by how well the exciter couples to that mode. Averaged over the panel, the cross terms between modes cancel, which leaves a clean sum the app can evaluate for every position on the heat map at once.

### Raggedness

That response always slopes downward overall — that's the panel's physics, not something placement can fix. What placement *does* control is the peaks and dips on top of the slope. So the app scores each position by **raggedness**: how far, in dB RMS, the response wanders from its own one-octave-smoothed trend. Lower is flatter.

On the 300 × 200 × 3 mm acrylic panel, the best position measures ±2.9 dB; the centre measures ±4.3 dB. Click anywhere on the heat map to see that position's response next to the best one.

### Modal Density and Gaps

Below the first few hundred hertz, any panel has only a handful of modes, and no exciter position can fill the gaps between them. The app's *modes per ⅓ octave* chart shows where those gaps are. The default acrylic panel has two empty bands, at 63 and 80 Hz, between its first pair of modes and the next. Bigger, thinner or less stiff panels push the gaps lower.

### What This Still Can't Tell You

The response here is how much the panel *vibrates*, not how loud it sounds at a listening position. Turning vibration into sound pressure needs radiation efficiency, which depends on frequency, panel size and baffling. This tier tells you how evenly a position drives the panel — very useful for placement and design comparisons — but measured output is the final word.

---

## 05 — Topology as Acoustic Engineering

Once you move beyond a plain rectangle, the interesting design space opens up. Two tools are available: **stiffeners** (adding material) and **cutouts** (removing it). Both modify the path bending waves take across the panel — but they work differently.

### Stiffeners

A stiffener — a rib, bar, or bonded strip — raises bending stiffness along its axis. Bending waves traveling parallel to the stiffener propagate faster in that region. Waves crossing the stiffener encounter a stiffness discontinuity that reflects and scatters energy. The net effect is to redistribute mode frequencies and, crucially, to engineer anisotropy into an otherwise isotropic panel — exactly what the violin's bass bar does to spruce.

In the app, on a 300 × 200 × 3 mm panel:

- A spruce bar (5 × 8 mm) along the width of an **acrylic** panel raises the first bending mode along the bar from 56 Hz to 100 Hz. The twisting mode barely moves (52 → 53 Hz), because a thin bar adds little resistance to twist.
- Three spruce ribs across the grain of a **balsa** panel raise its weak cross-grain mode from 64 Hz to 155 Hz — the bass-bar principle used to correct a material's weak direction. The panel then has fewer modes below 5 kHz (92 instead of 108), a trade-off worth watching in the modal-density chart.

### Cutouts

A cutout doesn't just remove mass — it *severs a wave propagation path*. A slot forces bending waves to detour around the cut ends, increasing effective path length. Longer paths mean lower effective wave speed in that direction, shifting modes downward and creating directional asymmetry. A panel with parallel slots becomes strongly anisotropic: waves traveling along the slots propagate freely; waves crossing them must detour significantly.

On the same acrylic panel:

| Panel | 1st mode | 2nd mode | Modes below 5 kHz |
|---|---|---|---|
| Plain | 52 Hz | 56 Hz | 115 |
| 40 mm round hole | 51 Hz | 55 Hz | 117 |
| One slot across the middle | 46 Hz | 51 Hz | 120 |
| Three parallel slots | 42 Hz | 50 Hz | 126 |

A round hole barely matters. Slots lower the modes that bend across them and pack in more modes — more modal density is exactly what a DML panel wants.

> **The slotted panel interpreted:** A rectangular panel with internal slots creates a series of coupled resonating fingers — like a marimba where the bars are still joined at both ends. Each finger has its own resonant character, but couples energy with its neighbors. The panel preferentially radiates energy along the slot axis, biasing the modal density directionally. This is engineered anisotropy from topology alone — no exotic materials required.

### Directional Bias

Both stiffeners and cutouts create directional bias — the panel radiates differently depending on the axis of the structure. This isn't necessarily a problem; it can be a design handle. Orienting that bias toward the listening area, or using it to suppress room modes in a particular axis, are real strategies. The key is knowing which direction you're biasing toward, which requires either calculation or measurement.

### A Note on Tuning Forks

The same principle applies at small scale: a tuning fork's resonant frequency can be trimmed by removing material from the tips (lowers frequency) or the base (raises it). Precision forks are adjusted this way after casting. Cutouts in a panel are doing the same thing — just distributed across a 2D surface rather than at a single point.

---

## 06 — When You Need Finite Element Analysis

Once the outline isn't a rectangle, or there are cutouts or stiffeners, there are no closed-form mode shapes. You need numerical methods — finite element analysis (FEA).

### What FEA Actually Does

FEA divides the panel into many small elements, each described by simple local equations. It assembles these into a global system and solves for the mode shapes and frequencies of the whole structure at once — an eigenvalue problem.

### How Panel Speaker Designer Does It

The whole pipeline runs inside the app, in Rust:

```
Outline (built-in shape or your Bézier curves) + holes + slots + rib lines
        │
Mesh: constrained Delaunay triangulation, refined to a quality mesh
        │
Elements: DKT thin-plate triangles (with grain-dependent stiffness),
          plus beam elements along each rib
        │
Eigen-solve: sparse factorisation + block Lanczos → modes up to your max frequency
        │
Same scoring as before: heat map, response, node lines
```

The mesh is sized from the frequency range — about five elements per half-wavelength at the top of the band — which gives about 1–2% frequency error there and much less for the low modes. The solver was checked against the analytic rectangle solver (both isotropic acrylic and grain-heavy balsa), the exact simply-supported solution, and Leissa's values for a free circular plate, including the matched pairs of modes a circle has.

Plain rectangles still use the faster analytic solver. Everything else takes about a second per solve, which is why the app has a **Solve** button rather than recalculating on every keystroke.

### Stiffeners in FEA

A rib is modelled as a beam running along mesh edges, bending along its length and twisting, with its mass added to the panel. Because it's glued to one face, a rib bends together with a strip of the panel as a T-shaped section — much stiffer than the bar alone. The app uses the standard engineering approximation for how wide that strip is (the rib width plus ten panel thicknesses either side). Real composite action depends on the panel's in-plane stiffness, which a thin-plate model doesn't include, so treat rib results as good for comparing designs and confirm with measurement.

### Going Further in Python

If you want to explore beyond the app — anisotropic composites with arbitrary layups, thick plates, or coupled acoustics — the Python ecosystem has the tools: Gmsh for meshing, and FEniCSx or SfePy for general-purpose FEA. They're far more flexible and far more work.

### What FEA Still Can't Tell You

FEA gives you structural modes — how the panel vibrates. It doesn't directly give you acoustic output (SPL vs. frequency at a listening point). For that you need acoustic modelling — for a flat panel in a baffle, the Rayleigh integral is a tractable next step. In practice: use structural analysis to design the panel, then *measure* the acoustic result. Simulation guides the design; measurement validates it.

---

## 07 — Measurement and Calibration

Simulation is only as good as its material constants, and those are the weak link. Foam and balsa stiffness can vary by ±50% from one sheet to the next, so a precise model fed a textbook value gives a precise wrong answer. Two measurements close the loop.

### Tap Test: Calibrate the Material

Hang the bare panel from two threads (that's a free edge), tap it lightly, and record the sound with a phone or measurement mic. The spectrum shows sharp peaks at the panel's lowest modes. Their *pattern* is set by the shape; their *absolute frequencies* are set by stiffness and density. Weigh the panel for density, and the measured peaks pin down the stiffness of your actual sheet.

### Frequency Response: Validate the Design

Once a panel is built, measure its frequency response with a measurement mic and REW (Room EQ Wizard, free). Where measurement and prediction agree, trust the model; where they don't, the gap points to what's missing — damping, mounting, the exciter's own mass.

Neither is in the app yet (roadmap phase C, below), but the tap test is worth doing on any sheet before you trust the numbers.

---

## 08 — Roadmap: The Sensible Order of Operations

| Phase | Goal | Status in the app |
|---|---|---|
| 1 — Analytical placement | Beat centre placement on a rectangle | **Done** |
| A — Accurate modes | True free-edge modes, wood grain, exciter size | **Done** |
| B — Response scoring | Score placement by response flatness; show modal density and gaps | **Done** |
| D — Geometry (FEA) | Any outline, cutouts, stiffeners, custom Bézier shapes | **Done** |
| C — Calibration | Fit stiffness to tap-test peaks; overlay REW measurements | Waiting on a physical panel |
| E — Radiation | Estimate sound pressure from a baffled panel (Rayleigh integral) | Optional |

> Phase 1 is an afternoon. FEA is a project. Measurement is where you find out what you missed. Start with a plain rectangle — it already tells you something real.
