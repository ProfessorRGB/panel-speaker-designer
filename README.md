# Panel Speaker Designer

A native desktop app for designing flat-panel (DML) speakers: work out where to mount the exciter, and see how the panel's shape, material, cutouts and stiffeners change its behaviour.

Built with [Tauri 2](https://tauri.app) — Rust backend, TypeScript/Canvas frontend. Builds natively on macOS, Windows, and Linux.

---

## What it does

Distributed Mode Loudspeaker (DML) panels work by exciting as many resonant modes as possible across a flat panel. Where you mount the exciter determines which modes get driven — place it on a node line of a given mode and that mode contributes nothing to the output.

The app calculates the panel's vibration modes, predicts how the panel responds across frequency when driven from every possible position, and scores each position by how flat that response is. The result is a heat map of good (warm) and bad (cool) positions with the best one marked, plus charts of the predicted response and of how the modes are spread across frequency.

Panels can be rectangles, rounded rectangles, ellipses, regular polygons or outlines you draw yourself with Bézier curves, with holes, slots and glued-on stiffening ribs.

---

## Using it

1. **Set up the panel** in the sidebar: shape, size and thickness, material (presets or your own values), and any cutouts or ribs.
2. **Press Solve** (or Enter in any field). Edits fade the heat map until you solve again. The **Auto** checkbox solves after each change instead — handy for plain rectangles, which solve instantly; other shapes take about a second.
3. **Read the results:**
   - The heat map shows each position's score; the crosshair marks the best one. The dimmed band around the edges is excluded from the search.
   - Hover to read the score at any point. **Click** a point to plot its response against the best position's.
   - **Mode overlay** draws any mode's node lines on the panel.
4. **Custom outlines:** choose *Custom (Bézier)*. The outline starts as a copy of the previous shape. Drag points and handles, double-click the outline to add a point, select a point and press Delete to remove it, Option-click a point to switch between a corner and a smooth curve, and ⌘Z to undo. Press *Done editing*, then Solve.

---

## Physics

The panel is a thin (Kirchhoff) plate of orthotropic material: stiffness can differ along the width (`E_x`, e.g. along the grain) and height (`E_y`), with shear modulus `G` controlling twisting. Isotropic materials are the special case `E_x = E_y`, `G = E/2(1+ν)`. Thin-plate theory needs the panel to be much larger than it is thick, so width and height must be at least 10× the thickness.

**Rectangles (analytic):** free-edge modes come from a Rayleigh-Ritz solution using Legendre polynomials. This gives the true free-edge modes rather than the common `cos·cos` approximation, which misplaces free-edge frequencies by up to ~2× (and much more for wood). Results match Leissa's published values for a free square plate to within 0.5%. Simply-supported rectangles use the exact `sin·sin` solution.

**Everything else (FEA):** the outline, holes, slots and rib lines are meshed into triangles and solved with DKT thin-plate elements, beam elements for ribs, and a sparse eigen-solver, all in Rust. The mesh is sized from the frequency range (five elements per half-wavelength), giving about 1–2% frequency error at the top of the band. Verified against the rectangle solver (acrylic and balsa), the exact simply-supported solution, and Leissa's values for a free circular plate.

**Stiffeners:** a rib glued to one face bends together with a strip of panel as a T-section. The app uses a standard effective-width approximation (rib width plus 10 panel thicknesses each side), so rib results are approximate; they're best used to compare designs.

**Coupling:** an exciter drives the panel around its voice-coil ring, so its coupling to each mode is the mode shape averaged around that ring.

**Response and placement score:** from the couplings, the damped modal response gives the panel's average vibration level across frequency (in 1/12-octave bands, with the loss factor η setting damping). Placement is scored by *raggedness*: how far the response wanders, in dB RMS, from its own octave-smoothed trend. The overall slope is set by the panel; peaks and dips are what exciter position controls. The older score (sum of `|coupling|`) is still available.

**Edge margin:** free edges move far more than anywhere practical to mount an exciter, so the search excludes positions within 10% of the shorter side from the outer edge, and within the exciter's radius plus 3 mm of a cutout. The heat map's colour scale is set by the searchable area for the same reason.

### Limits

- The response is a vibration measure, not sound pressure at a listening position.
- Material constants are the weak link: foam and balsa vary a lot between sheets. Measure if you can (see the lesson).
- Rib stiffness is an approximation; ribs are straight only.
- Designs aren't saved between sessions.
- Large problems (foam panels, 1 mm sheets) take 2–6 s per solve.

---

## Features

- **Panel shapes** — rectangle (exact solver), plus rounded rectangle, ellipse/circle, regular polygon and custom Bézier outlines drawn on the panel view (FEA)
- **7 material presets** — XPS foam, EPS foam, balsa, birch plywood, acrylic, aluminium, carbon fibre — with fields that auto-fill and remain editable; balsa and plywood are orthotropic (grain along the width)
- **Cutouts** — round holes and slots at any angle; slots cut the bending path across them, lowering those modes and adding modal density
- **Stiffeners** — straight ribs glued to one face (spruce, carbon fibre, aluminium or the panel's own material), like a violin's bass bar
- **Exciter size** — voice-coil diameter, which averages out modes smaller than the coil
- **Damping** — loss factor per material, which sets how sharp the resonance peaks are
- **Both boundary conditions** — free edges (realistic for DML) and simply supported (outer edge; cutout edges are always free)
- **Heat map** — placement score on a proportional canvas, clipped to the panel's real outline
- **Response chart** — predicted response at the best position; click anywhere on the panel to compare that position
- **Modal density** — modes per ⅓ octave, with empty bands flagged
- **Mode node-line overlay** — select any mode to see its node lines drawn on the panel
- **Solve on demand** — Solve button and Enter, with optional auto-solve

---

## Requirements

- [Rust](https://rustup.rs) (stable)
- Node.js 18+
- macOS 11+, Windows 10+, or a modern Linux desktop

---

## Development

```bash
# Install frontend dependencies
npm install

# Run in dev mode — hot-reloads HTML/CSS/TS, recompiles Rust on change
npm run tauri dev

# Run the solver tests (release mode: the FEA tests are slow in debug)
cd src-tauri && cargo test --release --lib

# Build a release .app bundle
npm run tauri build -- --bundles app
```

The first dev build after a fresh clone compiles all dependencies with optimisation, which takes a couple of minutes.

---

## Roadmap

| Phase | Status | Scope |
|---|---|---|
| 1 — Analytical | **Done** | Rectangular panel, placement heat map |
| A — Accurate modes | **Done** | Rayleigh-Ritz free-edge modes, orthotropic materials, exciter footprint |
| B — Response scoring | **Done** | Damped modal response; score placement by flatness; modes-per-octave plot |
| D1 — FEA engine | **Done** | Native Rust FEA (DKT elements, sparse eigen-solver); rounded-rectangle, ellipse/circle and polygon panels |
| D2 — Cutouts | **Done** | Round holes and rounded-end slots at any angle |
| D3 — Stiffeners | **Done** | Straight ribs glued to one face |
| D4 — Custom outlines | **Done** | Bézier outline editor on the panel view |
| C — Calibration | Planned | Fit stiffness to tap-test frequencies; overlay REW measurements |
| E — Radiation | Optional | Baffled-panel SPL estimate via the Rayleigh integral |

The original plan (separate SfePy and FEniCSx phases via a Python sidecar) was replaced by one native Rust FEA solver. Calibration (C) is waiting on a physical panel to measure.

---

## Reference

The physics and design rationale behind this tool, including how to calibrate against a real panel, are written up as a self-contained lesson in [`docs/panel-speaker-lesson.md`](docs/panel-speaker-lesson.md).
