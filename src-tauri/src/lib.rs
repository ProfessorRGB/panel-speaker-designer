mod plate;

use plate::{Boundary, GridBasis, Plate, SolveKey, Solution};
use serde::{Deserialize, Serialize};
use std::f64::consts::PI;
use std::sync::{Arc, Mutex};

#[derive(Deserialize, Clone)]
pub struct PanelParams {
    pub lx: f64,          // panel width  [m]
    pub ly: f64,          // panel height [m]
    pub h: f64,           // thickness    [m]
    pub ex: f64,          // Young's modulus along x (width / grain) [Pa]
    pub ey: f64,          // Young's modulus along y [Pa]
    pub g: f64,           // in-plane shear modulus [Pa]
    pub nu: f64,          // major Poisson's ratio ν_xy
    pub rho: f64,         // density [kg/m³]
    pub boundary: String, // "free" | "simply_supported"
    pub freq_max: f64,    // upper frequency limit [Hz]
    pub grid_n: usize,    // heat-map grid resolution (N×N)
    pub exciter_d: f64,   // exciter voice-coil diameter [m]; 0 = point drive
}

#[derive(Serialize, Clone)]
pub struct ModeInfo {
    pub m: u32,
    pub n: u32,
    pub freq: f64,
}

#[derive(Serialize)]
pub struct CalculationResult {
    // Row-major NxN, normalised so the interior search region spans [0, 1];
    // cells inside the edge margin can exceed 1.
    pub grid: Vec<f64>,
    pub grid_n: usize,
    pub modes: Vec<ModeInfo>,
    pub mode_count: usize,
    pub optimal_x: f64,      // normalised [0, 1]
    pub optimal_y: f64,
    pub optimal_score_raw: f64,
    pub margin_x: f64,       // edge margin excluded from the search, normalised
    pub margin_y: f64,
    // Set when the solver could not resolve modes all the way to freq_max;
    // modes above this frequency are omitted.
    pub truncated_above: Option<f64>,
}

impl PanelParams {
    fn key(&self) -> Result<SolveKey, String> {
        let plate = Plate {
            lx: self.lx,
            ly: self.ly,
            h: self.h,
            ex: self.ex,
            ey: self.ey,
            g: self.g,
            nu_xy: self.nu,
            rho: self.rho,
        };
        plate.validate()?;
        if !self.freq_max.is_finite() || self.freq_max <= 0.0 {
            return Err("Max frequency must be positive".into());
        }
        let boundary = match self.boundary.as_str() {
            "free" => Boundary::Free,
            "simply_supported" => Boundary::SimplySupported,
            other => return Err(format!("Unknown boundary condition: {other}")),
        };
        Ok(SolveKey::for_freq(plate, boundary, self.freq_max))
    }
}

// The eigen-solve depends only on the panel, not on the grid or exciter, so
// keep the last one for re-scoring and mode overlays.
static CACHE: Mutex<Option<Arc<Solution>>> = Mutex::new(None);

fn solve_cached(key: SolveKey) -> Arc<Solution> {
    let mut cache = CACHE.lock().unwrap_or_else(|e| e.into_inner());
    if let Some(sol) = cache.as_ref().filter(|s| s.key == key) {
        return sol.clone();
    }
    let sol = Arc::new(plate::solve(key));
    *cache = Some(sol.clone());
    sol
}

fn modes_in_range(sol: &Solution, freq_max: f64) -> &[plate::Mode] {
    let count = sol.modes.partition_point(|m| m.freq <= freq_max);
    &sol.modes[..count]
}

// Bilinear sample taps for the exciter's voice-coil ring around each grid
// cell: (cell index, weight) pairs that together average over the ring.
const RING_POINTS: usize = 8;

fn ring_stencil(n: usize, rx: f64, ry: f64) -> Vec<Vec<(usize, f64)>> {
    // Below a quarter cell the ring is indistinguishable from a point.
    if rx * (n as f64) < 0.25 && ry * (n as f64) < 0.25 {
        return (0..n * n).map(|i| vec![(i, 1.0)]).collect();
    }
    let max = (n - 1) as f64;
    let w = 1.0 / RING_POINTS as f64;
    let mut stencil = Vec::with_capacity(n * n);
    for row in 0..n {
        for col in 0..n {
            let mut taps = Vec::with_capacity(RING_POINTS * 4);
            for k in 0..RING_POINTS {
                let theta = 2.0 * PI * k as f64 / RING_POINTS as f64;
                // Grid coordinates (cell centres at integers), clamped to the panel.
                let gx = (col as f64 + rx * n as f64 * theta.cos()).clamp(0.0, max);
                let gy = (row as f64 + ry * n as f64 * theta.sin()).clamp(0.0, max);
                let (x0, y0) = (gx.floor() as usize, gy.floor() as usize);
                let (x1, y1) = ((x0 + 1).min(n - 1), (y0 + 1).min(n - 1));
                let (fx, fy) = (gx - x0 as f64, gy - y0 as f64);
                taps.push((y0 * n + x0, w * (1.0 - fx) * (1.0 - fy)));
                taps.push((y0 * n + x1, w * fx * (1.0 - fy)));
                taps.push((y1 * n + x0, w * (1.0 - fx) * fy));
                taps.push((y1 * n + x1, w * fx * fy));
            }
            stencil.push(taps);
        }
    }
    stencil
}

#[tauri::command(async)]
fn compute_heatmap(params: PanelParams) -> Result<CalculationResult, String> {
    let sol = solve_cached(params.key()?);
    let modes = modes_in_range(&sol, params.freq_max);
    let mode_count = modes.len();
    let n = params.grid_n.clamp(4, 100);
    let truncated_above = (sol.freq_limit < params.freq_max).then_some(sol.freq_limit);

    // Optimal search excludes a 10% edge margin (widened if needed to keep the
    // whole exciter on the panel). Free edges always move a lot, and edges and
    // corners are impractical mounting locations.
    const EDGE_MARGIN: f64 = 0.10;
    let radius = params.exciter_d.max(0.0) / 2.0;
    let margin_x = EDGE_MARGIN.max(radius / params.lx);
    let margin_y = EDGE_MARGIN.max(radius / params.ly);

    let mode_infos = modes.iter().map(|m| ModeInfo { m: m.m, n: m.n, freq: m.freq }).collect();

    if mode_count == 0 {
        return Ok(CalculationResult {
            grid: vec![0.0; n * n],
            grid_n: n,
            modes: mode_infos,
            mode_count: 0,
            optimal_x: 0.5,
            optimal_y: 0.5,
            optimal_score_raw: 0.0,
            margin_x,
            margin_y,
            truncated_above,
        });
    }

    // --- Score every grid cell ---
    // Coupling to a mode is its displacement averaged over the voice-coil ring;
    // the score is the sum of |coupling| over all modes.
    let stencil = ring_stencil(n, radius / params.lx, radius / params.ly);
    let basis = GridBasis::new(&sol, n);
    let mut grid = vec![0.0f64; n * n];
    for mode in modes {
        let shape = basis.eval(&sol, mode);
        for (cell, taps) in grid.iter_mut().zip(&stencil) {
            let coupling: f64 = taps.iter().map(|&(i, w)| w * shape[i]).sum();
            *cell += coupling.abs();
        }
    }

    // The colour range comes from the interior only: free edges move far more
    // than anywhere an exciter can go, and would otherwise flatten the interior.
    let mut max_score = f64::NEG_INFINITY;
    let mut min_score = f64::INFINITY;
    let mut opt_x = 0.5f64;
    let mut opt_y = 0.5f64;
    let mut opt_score = f64::NEG_INFINITY;

    for row in 0..n {
        let norm_y = (row as f64 + 0.5) / n as f64;
        for col in 0..n {
            let norm_x = (col as f64 + 0.5) / n as f64;
            let score = grid[row * n + col];

            let interior = norm_x >= margin_x && norm_x <= 1.0 - margin_x
                        && norm_y >= margin_y && norm_y <= 1.0 - margin_y;
            if !interior {
                continue;
            }

            max_score = max_score.max(score);
            min_score = min_score.min(score);
            if score > opt_score {
                opt_score = score;
                opt_x = norm_x;
                opt_y = norm_y;
            }
        }
    }

    // Normalise to [0, 1]
    let range = max_score - min_score;
    if range.is_finite() && range > 1e-12 {
        for v in &mut grid {
            *v = (*v - min_score) / range;
        }
    }

    Ok(CalculationResult {
        grid,
        grid_n: n,
        modes: mode_infos,
        mode_count,
        optimal_x: opt_x,
        optimal_y: opt_y,
        optimal_score_raw: opt_score,
        margin_x,
        margin_y,
        truncated_above,
    })
}

/// Shape of mode `index` (in the order `compute_heatmap` returned) sampled on
/// an n×n grid, row-major. Used to draw node lines.
#[tauri::command(async)]
fn mode_shape(params: PanelParams, index: usize, n: usize) -> Result<Vec<f64>, String> {
    let sol = solve_cached(params.key()?);
    let mode = modes_in_range(&sol, params.freq_max)
        .get(index)
        .ok_or("Mode index out of range")?;
    Ok(plate::shape_grid(&sol, mode, n.clamp(4, 200)))
}

#[cfg_attr(mobile, tauri::mobile_entry_point)]
pub fn run() {
    tauri::Builder::default()
        .plugin(tauri_plugin_opener::init())
        .invoke_handler(tauri::generate_handler![compute_heatmap, mode_shape])
        .run(tauri::generate_context!())
        .expect("error while running tauri application");
}
