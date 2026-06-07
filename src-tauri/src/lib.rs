use serde::{Deserialize, Serialize};
use std::f64::consts::PI;

#[derive(Deserialize, Clone)]
pub struct PanelParams {
    pub lx: f64,          // panel width  [m]
    pub ly: f64,          // panel height [m]
    pub h: f64,           // thickness    [m]
    pub e: f64,           // Young's modulus [Pa]
    pub rho: f64,         // density [kg/m³]
    pub nu: f64,          // Poisson's ratio
    pub boundary: String, // "free" | "simply_supported"
    pub freq_max: f64,    // upper frequency limit [Hz]
    pub grid_n: usize,    // heat-map grid resolution (N×N)
}

#[derive(Serialize, Clone)]
pub struct ModeInfo {
    pub m: u32,
    pub n: u32,
    pub freq: f64,
}

#[derive(Serialize)]
pub struct CalculationResult {
    pub grid: Vec<f64>,      // row-major NxN, normalised to [0, 1]
    pub grid_n: usize,
    pub modes: Vec<ModeInfo>,
    pub mode_count: usize,
    pub optimal_x: f64,      // normalised [0, 1]
    pub optimal_y: f64,
    pub optimal_score_raw: f64,
}

fn bending_stiffness(e: f64, h: f64, nu: f64) -> f64 {
    e * h.powi(3) / (12.0 * (1.0 - nu * nu))
}

// Approximate plate mode frequency using the thin-plate dispersion relation.
// Exact for simply-supported; a well-accepted approximation for free edges
// used throughout the DML literature.
fn mode_frequency(m: u32, n: u32, lx: f64, ly: f64, h: f64, e: f64, rho: f64, nu: f64) -> f64 {
    let d = bending_stiffness(e, h, nu);
    let kx = m as f64 * PI / lx;
    let ky = n as f64 * PI / ly;
    let omega = (d / (rho * h)).sqrt() * (kx * kx + ky * ky);
    omega / (2.0 * PI)
}

// Signed coupling amplitude between an exciter at (x, y) and mode (m, n).
// Free BC  : cosine mode shapes — W_mn = φ_m(x)·φ_n(y),  φ_0=1, φ_k=cos(kπx/L)
// SS BC    : sine mode shapes   — W_mn = sin(mπx/Lx)·sin(nπy/Ly)
fn mode_amplitude(m: u32, n: u32, x: f64, y: f64, lx: f64, ly: f64, boundary: &str) -> f64 {
    if boundary == "simply_supported" {
        (m as f64 * PI * x / lx).sin() * (n as f64 * PI * y / ly).sin()
    } else {
        let phi_x = if m == 0 { 1.0 } else { (m as f64 * PI * x / lx).cos() };
        let phi_y = if n == 0 { 1.0 } else { (n as f64 * PI * y / ly).cos() };
        phi_x * phi_y
    }
}

#[tauri::command]
fn compute_heatmap(params: PanelParams) -> CalculationResult {
    // --- Collect valid modes ---
    let mut modes: Vec<ModeInfo> = Vec::new();

    for m in 0u32..=24 {
        for n in 0u32..=24 {
            // Simply-supported: modes start at (1,1)
            if params.boundary == "simply_supported" && (m == 0 || n == 0) {
                continue;
            }
            // Free: skip rigid-body modes — (0,0), (1,0), (0,1)
            if params.boundary == "free" && m + n < 2 {
                continue;
            }

            let freq = mode_frequency(
                m, n, params.lx, params.ly, params.h,
                params.e, params.rho, params.nu,
            );

            if freq > 0.0 && freq <= params.freq_max {
                modes.push(ModeInfo { m, n, freq });
            }
        }
    }

    modes.sort_by(|a, b| a.freq.partial_cmp(&b.freq).unwrap_or(std::cmp::Ordering::Equal));

    let mode_count = modes.len();

    let n = params.grid_n.clamp(4, 100);

    if mode_count == 0 {
        return CalculationResult {
            grid: vec![0.0; n * n],
            grid_n: n,
            modes,
            mode_count: 0,
            optimal_x: 0.5,
            optimal_y: 0.5,
            optimal_score_raw: 0.0,
        };
    }

    // --- Score every grid cell ---
    // Optimal search excludes a 10% edge margin — corners always win under the
    // cosine approximation (cos(mπ·0/L)=1 for all m), but are physically
    // impractical mounting locations.
    const EDGE_MARGIN: f64 = 0.10;

    let mut grid = vec![0.0f64; n * n];
    let mut max_score = f64::NEG_INFINITY;
    let mut min_score = f64::INFINITY;
    let mut opt_x = 0.5f64;
    let mut opt_y = 0.5f64;
    let mut opt_score = f64::NEG_INFINITY;

    for row in 0..n {
        let norm_y = (row as f64 + 0.5) / n as f64;
        let y = norm_y * params.ly;
        for col in 0..n {
            let norm_x = (col as f64 + 0.5) / n as f64;
            let x = norm_x * params.lx;

            let score: f64 = modes.iter().map(|mode| {
                mode_amplitude(mode.m, mode.n, x, y, params.lx, params.ly, &params.boundary)
                    .abs()
            }).sum();

            grid[row * n + col] = score;

            if score > max_score { max_score = score; }
            if score < min_score { min_score = score; }

            let interior = norm_x >= EDGE_MARGIN && norm_x <= 1.0 - EDGE_MARGIN
                        && norm_y >= EDGE_MARGIN && norm_y <= 1.0 - EDGE_MARGIN;

            if interior && score > opt_score {
                opt_score = score;
                opt_x = norm_x;
                opt_y = norm_y;
            }
        }
    }

    // Normalise to [0, 1]
    let range = max_score - min_score;
    if range > 1e-12 {
        for v in &mut grid {
            *v = (*v - min_score) / range;
        }
    }

    CalculationResult {
        grid,
        grid_n: n,
        modes,
        mode_count,
        optimal_x: opt_x,
        optimal_y: opt_y,
        optimal_score_raw: opt_score,
    }
}

#[cfg_attr(mobile, tauri::mobile_entry_point)]
pub fn run() {
    tauri::Builder::default()
        .plugin(tauri_plugin_opener::init())
        .invoke_handler(tauri::generate_handler![compute_heatmap])
        .run(tauri::generate_context!())
        .expect("error while running tauri application");
}
