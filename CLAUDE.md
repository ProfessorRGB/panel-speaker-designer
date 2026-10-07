# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Commands

```bash
# Run in dev mode (starts Vite + compiles Rust + opens native window)
npm run tauri dev

# Type-check TypeScript only (no emit)
npx tsc --noEmit

# Check Rust without building the full binary
cd src-tauri && cargo check

# Build release bundle (.app)
npm run tauri build
```

Hot-reload applies to frontend changes (HTML/CSS/TS) automatically. Rust changes trigger a recompile and window reload via Tauri's file watcher.

## Architecture

This is a **Tauri 2.0** app: a Rust backend exposed as IPC commands, a Vite/TypeScript frontend rendered in WKWebView.

### Data flow

```
User input (HTML inputs) → getParams() → invoke("compute_heatmap") → Rust
    → CalculationResult (grid[], modes[], optimal_x/y) → render() on <canvas>
```

Calculation runs when the user presses **Solve** (or Enter in any sidebar field). Edits only mark the result stale (`scheduleCalculate` → `markStale`, faded heat map); with the **Auto** checkbox on, edits solve after a 400 ms pause. Live solving on every keystroke used to send half-typed values (e.g. a 4 mm width while typing "450") that the FEA solver choked on. `editVersion` keeps edits made during a solve marked stale. Resize re-renders without recalculating.

### Rust backend

`src-tauri/src/lib.rs` — Tauri commands, scoring and caching:
- `compute_heatmap(params) -> CalculationResult` — heat map, mode list, optimal position
- `mode_shape(params, index, n) -> Vec<f64>` — one mode sampled on an n×n grid, for the node-line overlay
- `response_at(params, x, y) -> ResponseCurve` — response (dB per band) and raggedness at a clicked position
- Both run off the main thread (`#[tauri::command(async)]`) and share the last eigen-solve through a static cache keyed on `SolveKey`
- **Coupling** to each mode = its shape averaged around the exciter's voice-coil ring (`exciter_d`; 0 = point drive)
- **Placement score** (`params.score`): `"flatness"` (default) = raggedness of the damped response, lower is better; `"coupling"` = sum of |coupling| over modes ≤ `freq_max` (the v0.1 score). `grid` is normalised so 1 = best; `grid_raw` keeps the score's units
- The solve is sized for `1.3 × freq_max` (`RESPONSE_HEADROOM`) so the response near `freq_max` includes the tails of modes just above it
- Grid cells are classified in `region`: 0 = search, 1 = edge margin, 2 = outside the panel (scores are NaN → JSON `null`). The margin is 10% of the shorter bounding-box side from the outer edge (or the exciter radius if larger), and exciter radius + 3 mm from cutout edges. Search and colour range use region 0 only, because free edges move far more than anywhere an exciter can go
- `outline` (normalised rings) and `solver` (`"analytic"` | `"fea"`) are returned for drawing and status
- `grid_n` is clamped to 100 on both the JS and Rust sides to prevent UI lock-up
- Guards against pathological input: `Plate::validate` rejects panels smaller than 10× their thickness (thin-plate theory), FEA solves only up to the requested frequency even when the mesh resolves more, and the Lanczos basis is capped (error instead of unbounded growth)
- Command fns must remain non-`pub` — Tauri 2's `#[tauri::command]` macro generates duplicate names when applied to `pub fn`

`src-tauri/src/model.rs` — one interface over both solvers. `ModelKey` picks the analytic solver for plain rectangles and FEA for every other shape or any panel with cutouts or stiffeners; `Model::solve` validates cutouts and stiffeners first; `section_properties` turns a rib into beam properties (T-section with an effective panel flange of rib width + 20 × panel thickness, minus the flange's own bending, which the plate already carries — an approximation); `Model::grid_eval` / `point_eval` sample mode shapes (NaN outside the panel). The model is cached on its key in `lib.rs`.

`src-tauri/src/geometry.rs` — outlines as closed paths of line and cubic Bézier segments (`Path`), built-in shapes (`Panel`: rectangle, rounded rectangle, ellipse, regular polygon) and cutouts (`Cutout`: round hole, rounded-end slot with angle) generated as paths, `Shape::validate` (cutouts inside the outline, ≥ 1 mm from edges and each other), small holes flattened to ≥ 24 edges, straight `Stiffener` ribs with `Shape::validate_stiffeners` (inside the material ≥ 1 mm from all edges, no crossings), adaptive flattening, and `Outline` inside/distance queries (holes supported by crossing parity). Coordinates in metres, origin top-left, y down.

`src-tauri/src/mesh.rs` — constrained Delaunay triangulation with refinement (`spade`): max area from target element size, 25° minimum angle, holes excluded, capped at 40k nodes.

`src-tauri/src/fem.rs` — FEA for arbitrary shapes:
- DKT plate element (nodal w, ∂w/∂x, ∂w/∂y), built constructively from its definition; orthotropic bending stiffness
- Lumped mass plus rotary inertia; simply-supported edges via a stiff spring on w at outline nodes
- Shift-invert block Lanczos (block 8, full reorthogonalisation) on a `faer` sparse Cholesky of K + σM; block size handles the repeated frequencies of symmetric shapes
- Element size: 5 elements per half-wavelength at 1.3 × `freq_max`, node count capped (`MAX_FEA_NODES`); the resolvable limit is reported as `truncated_above`
- Hole edges are always free; simply-supported applies to the outer outline only
- Stiffeners (`Beam`): Euler-Bernoulli bending along the rib plus St Venant torsion on the plate's nodal unknowns, rib mass lumped on w. The mesher embeds rib lines as constraint edges; with lines present it classifies triangles by centroid, because spade's inside/outside parity is flipped by open constraint lines
- Tests: rectangle vs the Ritz solver (isotropic and balsa), simply supported vs exact, free circle vs Leissa (including double modes); in `model.rs`, a tiny hole barely changes modes, a slotted panel converges under mesh refinement, a negligible rib changes nothing, and a ribbed strip matches free-free beam theory

`src-tauri/src/plate.rs` — analytic modal solver for rectangles (pure physics, with unit tests):
- Orthotropic material: `E_x`, `E_y`, `G`, `ν_xy`, principal axes along the edges (x = width = grain)
- **Free edges**: Rayleigh-Ritz with orthonormal Legendre polynomials, split into four symmetry blocks, solved with `nalgebra::SymmetricEigen`. Rigid-body modes are filtered out. Basis size per axis is chosen from `freq_max` (≈1.6 per half-wave + 10, capped at 48); modes above the resolvable limit are dropped and reported via `truncated_above`
- **Simply supported**: closed-form orthotropic `sin·sin`
- Mode `(m, n)` labels on free plates are approximate: median count of sign changes along several lines
- Tests check against Leissa's free-square values, convergence, symmetry and normalisation: `cd src-tauri && cargo test --release --lib`

`src-tauri/src/response.rs` — damped response and flatness:
- Spatially averaged mean-square velocity `Σ c_k² ω²/|ω_k²−ω²+iηω_k²|²` (mode orthogonality removes cross terms). A vibration proxy, not radiated SPL
- Evaluated as averages over 1/12-octave bands (Lorentzian band integral with a tail correction), so results don't depend on peak sharpness vs band spacing
- Heat map: `C² · T` as one matrix product (cells × modes times modes × bands)
- **Raggedness**: RMS dB deviation from the response's own one-octave moving average (window narrowed symmetrically at the ends)

### Frontend — `src/main.ts`

All state is module-level. Key globals: `lastResult`, `selectedModeIdx`, `selectedShape`, `requestId` (drops stale async results).

- `getParams()` reads inputs and converts units (mm→m, MPa→Pa)
- "Isotropic" checkbox derives `E_y = E_x` and `G = E/2(1+ν)`; presets with `ey`/`g` set are orthotropic
- `render()` draws the heat map on `<canvas id="heatmap">` using a 6-stop colormap (dark blue → red), dimming the excluded edge margin
- `drawNodeLines()` traces the zero contour of `selectedShape` (marching squares) in yellow
- Shape selector shows corner-radius / sides fields as needed. Cutouts are edited as cards (`renderCutouts`, `CUTOUT_FIELDS`) in mm and sent to Rust in metres. Stiffeners likewise (`renderStiffeners`, `RIB_MATERIALS`; "Same as panel" uses the panel's directional modulus along the rib) and are drawn on the heat map from local state. `render()` clips the heat map to the returned `outline` (Path2D, even-odd) and fills cells straddling curved edges from neighbours (`fillOutside`)
- Clicking the heat map sets `probe`; its response is fetched with `response_at` and re-fetched on every recalculation
- `src/charts.ts`: `ResponseChart` (log-frequency line chart, crosshair tooltip) and `ModeDensityChart` (modes per ⅓ octave, empty bands marked "0"). Canvases sit in `.chart-plot` containers because their backing store is sized at device pixel ratio
- Canvas is sized to preserve the panel's physical aspect ratio within the available container

### Roadmap

Phase 1, Phase A (accurate free-edge modes, orthotropy, exciter footprint), Phase B (response-flatness scoring, response and modal-density charts) Phase D1 (FEA engine, built-in shapes), D2 (holes and slots) and D3 (straight stiffeners) are done. Next in D:, D4 custom Bézier outlines (path model already supports cubic segments). C — calibrate stiffness from tap-test frequencies and overlay REW measurements; D — one FEA solver for arbitrary geometry, cutouts and stiffeners (preferably native Rust, avoiding a Python sidecar); E (optional) — Rayleigh-integral SPL estimate. See README.
