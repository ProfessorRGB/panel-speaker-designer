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

### Rust backend — `src-tauri/src/lib.rs`

Single Tauri command: `compute_heatmap(params: PanelParams) -> CalculationResult`

Key physics:
- **Bending stiffness**: `D = Eh³ / 12(1−ν²)`
- **Mode frequency**: `f_mn = (1/2π) · √(D/ρh) · (kx² + ky²)` where `kx = mπ/Lx`
- **Mode shape amplitude** — free BC: `cos(mπx/Lx)·cos(nπy/Ly)` (cosine approximation, standard for DML); simply-supported: `sin·sin`
- **Placement score**: sum of `|amplitude|` over all modes below `freq_max`
- **Optimal search**: excludes a 10% edge margin on all sides — without this, corners always win because `cos(mπ·0/L) = 1` for every mode
- Free BC skips modes where `m+n < 2` (rigid-body modes)
- `grid_n` is clamped to 100 on both the JS and Rust sides to prevent UI lock-up

`fn compute_heatmap` must remain non-`pub` — Tauri 2's `#[tauri::command]` macro generates duplicate names when applied to `pub fn`.

### Frontend — `src/main.ts`

All state is module-level. Key globals: `lastResult`, `selectedModeIdx`.

- `getParams()` reads inputs and converts units (mm→m, MPa→Pa)
- `render()` draws the heat map on `<canvas id="heatmap">` using a 6-stop colormap (dark blue → red)
- `drawNodeLines()` overlays mode node lines in yellow when a mode is selected from the toolbar dropdown
- Canvas is sized to preserve the panel's physical aspect ratio within the available container

### Planned future phases (from the reference document)

| Phase | Scope |
|---|---|
| 2 | Modal density / modes-per-octave analysis |
| 3 | Non-rectangular geometry via Gmsh mesh + SfePy FEA |
| 4 | Stiffeners, cutouts, anisotropic materials via FEniCSx |

FEA will be added as a Python sidecar process invoked via Tauri's sidecar API, keeping the same `invoke` pattern in the frontend.
