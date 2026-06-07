# Flat Panel Speakers: Physics, Topology & Tools

*From violin bass bars to finite element analysis — a self-contained lesson.*

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

Center placement on a rectangular panel is nearly the worst possible choice: it sits on the symmetry axes, which means it can't excite any antisymmetric mode — roughly half of all available modes. **Almost any informed placement beats center placement.** This is the low bar your first software only needs to clear.

```
BAD: center placement          BETTER: offset placement       OPTIMAL: scored placement

┌─────────────┐                ┌─────────────┐                ┌─────────────┐
│      │      │                │      │      │                │  ╌╌╌│╌╌╌╌╌ │
│      │      │                │      │      │                │╌╌╌╌╌│╌╌╌╌╌╌│
│──────●──────│                │──────│──────│                │  ╌╌╌│╌╌ ●  │
│      │      │                │      │  ●   │                │╌╌╌╌╌│╌╌╌╌╌╌│
│      │      │                │      │      │                │  ╌╌╌│╌╌╌╌╌ │
└─────────────┘                └─────────────┘                └─────────────┘
Center sits on node            Offset avoids (1,1)            Scored to avoid
lines of (1,1) mode            node lines                     multiple modes
```

### The Role of Panel Material

Bending wave speed in a panel depends on stiffness, density, and frequency. The dispersion relation — how wave speed varies with frequency — determines where modes fall. For an **isotropic** material (same properties in all directions), this is a single clean equation. For **anisotropic** materials (wood, some composites), bending wave speed differs along each axis, and the mode structure becomes directionally asymmetric.

---

## 03 — Level One Software: What You Can Calculate Analytically

For a rectangular, isotropic, free-edge panel, the mode shapes and frequencies are known analytically. You don't need simulation — you need arithmetic. This is the right place to start.

### What the Code Does

Given panel dimensions (width, height, thickness) and material constants (Young's modulus, density, Poisson's ratio), you can calculate: the frequency of each mode; the node line pattern of each mode; and a score for any candidate exciter position based on how close it sits to node lines across all modes.

**Key libraries:** NumPy (array math, mode frequencies), SciPy (eigenvalue solvers for boundary variations), Matplotlib (node line and score visualisation).

```python
# Simplified plate mode frequency — isotropic free plate
# f_mn depends on mode indices m, n and material/geometry constants

import numpy as np

def mode_frequency(m, n, Lx, Ly, h, E, rho, nu):
    # Bending stiffness D
    D = (E * h**3) / (12 * (1 - nu**2))
    # Approximate frequency for free rectangular plate
    kx = (m * np.pi) / Lx
    ky = (n * np.pi) / Ly
    omega = np.sqrt(D / rho / h) * (kx**2 + ky**2)
    return omega / (2 * np.pi)  # Hz

def score_position(x, y, Lx, Ly, modes):
    # Score how well a position couples to each mode
    # Higher score = drives more modes = better placement
    score = 0
    for m, n in modes:
        amplitude = np.abs(
            np.cos(m * np.pi * x / Lx) *
            np.cos(n * np.pi * y / Ly)
        )
        score += amplitude
    return score
```

A brute-force grid search over candidate positions, scored this way, gives you a heat map of good and bad placements. That's a genuinely useful result from maybe two hours of code.

### What It Can't Tell You

This tier doesn't predict actual SPL response, radiation efficiency, or how the panel sounds at a specific listening position. It tells you where *not* to put the exciter, which is already valuable. Actual acoustic output you measure — simulation at this level is for placement strategy, not final prediction.

---

## 04 — Advanced Physics: Topology as Acoustic Engineering

Once you move beyond a simple homogeneous rectangle, the interesting design space opens up. Two tools are available: **stiffeners** (adding material) and **cutouts** (removing it). Both modify the path bending waves take across the panel — but they work differently.

### Stiffeners

A stiffener — a rib, bar, or bonded strip — raises local bending stiffness along its axis. Bending waves traveling parallel to the stiffener propagate faster in that region. Waves crossing the stiffener encounter a stiffness discontinuity that reflects and scatters energy. The net effect is to redistribute mode frequencies and, crucially, to engineer anisotropy into an otherwise isotropic panel — exactly what the violin's bass bar does to spruce.

### Cutouts

A cutout doesn't just remove mass — it *severs a wave propagation path*. A slot forces bending waves to detour around the cut ends, increasing effective path length. Longer paths mean lower effective wave speed in that direction, shifting modes downward and creating directional asymmetry. A panel with parallel slots becomes strongly anisotropic: waves traveling along the slots propagate freely; waves crossing them must detour significantly.

> **The slotted panel interpreted:** A rectangular panel with internal slots creates a series of coupled resonating fingers — like a marimba where the bars are still joined at both ends. Each finger has its own resonant character, but couples energy with its neighbors. The panel preferentially radiates energy along the slot axis, biasing the modal density directionally. This is engineered anisotropy from topology alone — no exotic materials required.

### Directional Bias

Both stiffeners and cutouts create directional bias — the panel radiates differently depending on the axis of the structure. This isn't necessarily a problem; it can be a design handle. Orienting that bias toward the listening area, or using it to suppress room modes in a particular axis, are real strategies. The key is knowing which direction you're biasing toward, which requires either calculation or measurement.

### A Note on Tuning Forks

The same principle applies at small scale: a tuning fork's resonant frequency can be trimmed by removing material from the tips (lowers frequency) or the base (raises it). Precision forks are adjusted this way after casting. Cutouts in a panel are doing the same thing — just distributed across a 2D surface rather than at a single point.

---

## 05 — Advanced Software: When You Need Finite Element Analysis

Once the geometry is non-rectangular, the material is anisotropic, or stiffeners and cutouts are involved, the analytical solutions break down. The mode shapes no longer have closed-form expressions and you need numerical methods — finite element analysis (FEA).

### What FEA Actually Does

FEA divides your panel into many small elements, each described by simple local equations. It assembles these into a global system and solves for the mode shapes and frequencies of the whole structure simultaneously — an eigenvalue problem. For a thin plate, the relevant FEA formulation uses shell or plate elements specifically designed for bending behavior.

### The Full Pipeline

```
Geometry Input
(SVG / CAD / STEP / DXF)
        │
Geometry Processing
(svgpathtools / CadQuery / pythonOCC → boundary representation)
        │
Mesh Generation
(Gmsh Python API → triangulated mesh of panel surface)
        │
FEA Eigenvalue Solve
(FEniCSx or SfePy → mode frequencies + shapes)
        │
Analysis & Output
(node line visualisation, exciter scoring, modal density plots)
```

### Choosing a Library

| Library | Best For | Learning Curve |
|---|---|---|
| `FEniCSx` | Research-grade, arbitrary physics, well-documented plate formulations | Steep — requires comfort with weak form PDEs |
| `SfePy` | More approachable, good shell/plate element support, Python-native | Moderate |
| `Gmsh` | Meshing only — the keystone between geometry and FEA solver | Low for basic use |
| `Calculix` | Full solver, Abaqus-like capability, open source | High — traditional solver workflow |

### Stiffeners in FEA

Adding a bass-bar-style stiffener in FEA is modeled as a line of elements with different (higher) stiffness properties — or as beam elements superimposed on the plate mesh. This is where the perturbation approach breaks down and FEA earns its keep: the stiffness discontinuity along the bar creates mode shape changes that can't be approximated analytically.

### What FEA Still Can't Tell You

FEA gives you structural modes — how the panel vibrates. It doesn't directly give you acoustic output (SPL vs. frequency at a listening point). For that you need acoustic FEA or boundary element methods, which are significantly more complex. In practice: use structural FEA to optimize the panel, then *measure* the acoustic result. Simulation guides the design; measurement validates it.

---

## 06 — Roadmap: The Sensible Order of Operations

| Phase | Goal | Tools | Output |
|---|---|---|---|
| 1 — Analytical | Beat center placement on a rectangle | NumPy, Matplotlib | Position heat map, mode frequency list |
| 2 — Modal density | Check spectral coverage, find gaps | SciPy, NumPy | Modes-per-octave plot, coverage score |
| 3 — Simple FEA | Handle non-rectangular shapes | Gmsh + SfePy | Mode shapes for arbitrary geometry |
| 4 — Full FEA | Model stiffeners, cutouts, anisotropy | Gmsh + FEniCSx | Topology-aware mode optimisation |
| 5 — Measure | Validate against reality | REW, microphone, exciter | Waterfall plots, actual frequency response |

> Phase 1 is an afternoon. Phase 3 is a weekend. Phase 4 is a project. Phase 5 is where you find out what you missed. Start at 1 — it already tells you something real.
