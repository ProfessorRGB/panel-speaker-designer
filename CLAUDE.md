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

Calculation is triggered on every input change, debounced 250 ms. Resize re-renders without recalculating.

### Rust backend

`src-tauri/src/lib.rs` — Tauri commands, scoring and caching:
- `compute_heatmap(params) -> CalculationResult` — heat map, mode list, optimal position
- `mode_shape(params, index, n) -> Vec<f64>` — one mode sampled on an n×n grid, for the node-line overlay
- Both run off the main thread (`#[tauri::command(async)]`) and share the last eigen-solve through a static cache keyed on `SolveKey`
- **Placement score**: sum over modes of |mode shape averaged around the exciter's voice-coil ring| (`exciter_d`; 0 = point drive)
- **Optimal search** and the heat map's colour range use only the interior: a 10% edge margin (widened to keep the exciter on the panel) is excluded, because free edges move far more than anywhere an exciter can go
- `grid_n` is clamped to 100 on both the JS and Rust sides to prevent UI lock-up
- Command fns must remain non-`pub` — Tauri 2's `#[tauri::command]` macro generates duplicate names when applied to `pub fn`

`src-tauri/src/plate.rs` — modal solver (pure physics, with unit tests):
- Orthotropic material: `E_x`, `E_y`, `G`, `ν_xy`, principal axes along the edges (x = width = grain)
- **Free edges**: Rayleigh-Ritz with orthonormal Legendre polynomials, split into four symmetry blocks, solved with `nalgebra::SymmetricEigen`. Rigid-body modes are filtered out. Basis size per axis is chosen from `freq_max` (≈1.6 per half-wave + 10, capped at 48); modes above the resolvable limit are dropped and reported via `truncated_above`
- **Simply supported**: closed-form orthotropic `sin·sin`
- Mode `(m, n)` labels on free plates are approximate: median count of sign changes along several lines
- Tests check against Leissa's free-square values, convergence, symmetry and normalisation: `cd src-tauri && cargo test --release --lib`

### Frontend — `src/main.ts`

All state is module-level. Key globals: `lastResult`, `selectedModeIdx`, `selectedShape`, `requestId` (drops stale async results).

- `getParams()` reads inputs and converts units (mm→m, MPa→Pa)
- "Isotropic" checkbox derives `E_y = E_x` and `G = E/2(1+ν)`; presets with `ey`/`g` set are orthotropic
- `render()` draws the heat map on `<canvas id="heatmap">` using a 6-stop colormap (dark blue → red), dimming the excluded edge margin
- `drawNodeLines()` traces the zero contour of `selectedShape` (marching squares) in yellow
- Canvas is sized to preserve the panel's physical aspect ratio within the available container

### Roadmap

Phase 1 and Phase A (accurate free-edge modes, orthotropy, exciter footprint) are done. Next: B — score placement by damped modal response flatness (modes-per-octave falls out); C — calibrate stiffness from tap-test frequencies and overlay REW measurements; D — one FEA solver for arbitrary geometry, cutouts and stiffeners (preferably native Rust, avoiding a Python sidecar); E (optional) — Rayleigh-integral SPL estimate. See README.
