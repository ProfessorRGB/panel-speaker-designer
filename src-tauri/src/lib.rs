mod fem;
mod geometry;
mod mesh;
mod model;
mod plate;
mod response;

use geometry::{Panel, Pt};
use model::{Model, ModelKey};
use nalgebra::DMatrix;
use plate::{Boundary, Plate};
use response::Bands;
use serde::{Deserialize, Serialize};
use std::f64::consts::PI;
use std::sync::{Arc, Mutex};

#[derive(Deserialize, Clone)]
pub struct PanelParams {
    pub shape: String,    // "rectangle" | "rounded_rectangle" | "ellipse" | "polygon"
    pub corner_r: f64,    // rounded-rectangle corner radius [m]
    pub sides: usize,     // polygon side count
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
    pub eta: f64,         // damping loss factor
    pub score: String,    // "flatness" | "coupling"
}

#[derive(Serialize, Clone)]
pub struct ModeInfo {
    // Nodal-line counts; absent for shapes where they aren't meaningful.
    pub m: Option<u32>,
    pub n: Option<u32>,
    pub freq: f64,
}

// Grid cell classes for the heat map.
const REGION_SEARCH: u8 = 0;   // candidate exciter position
const REGION_MARGIN: u8 = 1;   // on the panel, too close to an edge
const REGION_OUTSIDE: u8 = 2;  // not panel material

#[derive(Serialize)]
pub struct CalculationResult {
    // Row-major NxN over the bounding box, normalised so the search region
    // spans [0, 1] with 1 = best. Null (NaN) outside the panel.
    pub grid: Vec<f64>,
    // The same cells in the score's own units: raggedness in dB (lower is
    // better) or summed coupling (higher is better).
    pub grid_raw: Vec<f64>,
    pub region: Vec<u8>,     // per cell: 0 search, 1 edge margin, 2 outside
    pub grid_n: usize,
    // Outline and hole rings in normalised coordinates, for drawing.
    pub outline: Vec<Vec<Pt>>,
    pub modes: Vec<ModeInfo>,
    pub mode_count: usize,
    pub optimal_x: f64,      // normalised [0, 1]
    pub optimal_y: f64,
    pub optimal_score_raw: f64,
    pub bands: Vec<f64>,         // response band centres [Hz]
    pub response_opt: Vec<f64>,  // response at the optimal position [dB]
    pub raggedness_opt: f64,     // its raggedness [dB]
    pub solver: &'static str,    // "analytic" | "fea"
    // Set when the solver could not resolve modes all the way to freq_max;
    // modes above this frequency are omitted.
    pub truncated_above: Option<f64>,
}

impl PanelParams {
    fn key(&self) -> Result<ModelKey, String> {
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
        let panel = match self.shape.as_str() {
            "rectangle" => Panel::Rectangle,
            // A zero radius is just a rectangle; use the faster exact solver.
            "rounded_rectangle" if self.corner_r <= 0.0 => Panel::Rectangle,
            "rounded_rectangle" => Panel::RoundedRectangle { radius: self.corner_r },
            "ellipse" => Panel::Ellipse,
            "polygon" => Panel::Polygon { sides: self.sides.clamp(3, 64) },
            other => return Err(format!("Unknown shape: {other}")),
        };
        // Resolve modes somewhat above freq_max, so the response near the top
        // of the band includes the tails of the modes just beyond it.
        Ok(ModelKey::for_freq(plate, boundary, panel, self.freq_max * RESPONSE_HEADROOM))
    }

    /// Minimum distance from any edge for a candidate exciter position.
    /// Free edges always move a lot, and edges are impractical mounting
    /// locations; the margin also keeps the whole exciter on the panel.
    fn edge_margin(&self) -> f64 {
        const EDGE_MARGIN: f64 = 0.10;
        (EDGE_MARGIN * self.lx.min(self.ly)).max(self.exciter_d.max(0.0) / 2.0)
    }
}

const RESPONSE_HEADROOM: f64 = 1.3;

/// Everything needed to evaluate the response at any exciter position.
struct ResponseModel {
    mode_count: usize,
    bands: Bands,
    transfer: DMatrix<f64>,  // modes × bands
}

impl ResponseModel {
    fn new(model: &Model, params: &PanelParams) -> ResponseModel {
        let mode_count = model.modes.partition_point(|m| m.freq <= params.freq_max * RESPONSE_HEADROOM);
        // From just below the first mode to freq_max, or lower if the solver
        // ran out of resolution.
        let f_start = model.modes.first().map_or(params.freq_max, |m| m.freq * 2f64.powf(-1.0 / 24.0));
        let f_end = params.freq_max.min(model.freq_limit / RESPONSE_HEADROOM);
        let bands = Bands::new(f_start, f_end);
        let freqs: Vec<f64> = model.modes[..mode_count].iter().map(|m| m.freq).collect();
        let transfer = response::transfer(&freqs, params.eta, &bands);
        ResponseModel { mode_count, bands, transfer }
    }

    /// Response (dB) for an exciter with the given modal couplings.
    fn curve(&self, couplings: &[f64]) -> Vec<f64> {
        let power: Vec<f64> = (0..self.bands.len())
            .map(|b| couplings.iter().enumerate().map(|(k, c)| c * c * self.transfer[(k, b)]).sum())
            .collect();
        response::to_db(&power)
    }

    /// Response at a normalised panel position, evaluating mode shapes
    /// directly around the voice-coil ring rather than interpolating the grid.
    fn curve_at(&self, model: &Model, params: &PanelParams, x: f64, y: f64) -> Vec<f64> {
        let radius = params.exciter_d.max(0.0) / 2.0;
        let points: Vec<Pt> = if radius > 0.0 {
            (0..RING_POINTS).map(|k| {
                let theta = 2.0 * PI * k as f64 / RING_POINTS as f64;
                [
                    (x + radius * theta.cos() / params.lx).clamp(0.0, 1.0),
                    (y + radius * theta.sin() / params.ly).clamp(0.0, 1.0),
                ]
            }).collect()
        } else {
            vec![[x, y]]
        };
        let eval = model.point_eval(points);
        let couplings: Vec<f64> = (0..self.mode_count)
            .map(|k| finite_mean(eval.eval(k).iter().map(|&v| (v, 1.0))))
            .collect();
        self.curve(&couplings)
    }
}

/// Weighted mean over the finite values only (ring points that fall off the
/// panel are skipped). NaN if none are finite.
fn finite_mean(values: impl Iterator<Item = (f64, f64)>) -> f64 {
    let (mut sum, mut weight) = (0.0, 0.0);
    for (v, w) in values {
        if v.is_finite() {
            sum += w * v;
            weight += w;
        }
    }
    if weight > 0.0 { sum / weight } else { f64::NAN }
}

// The eigen-solve depends only on the panel, not on the grid or exciter, so
// keep the last one for re-scoring and mode overlays.
static CACHE: Mutex<Option<Arc<Model>>> = Mutex::new(None);

fn solve_cached(key: ModelKey) -> Result<Arc<Model>, String> {
    let mut cache = CACHE.lock().unwrap_or_else(|e| e.into_inner());
    if let Some(model) = cache.as_ref().filter(|m| m.key == key) {
        return Ok(model.clone());
    }
    let model = Arc::new(Model::solve(key)?);
    *cache = Some(model.clone());
    Ok(model)
}

fn modes_in_range(model: &Model, freq_max: f64) -> usize {
    model.modes.partition_point(|m| m.freq <= freq_max)
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
                // Grid coordinates (cell centres at integers), clamped to the box.
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
    let model = solve_cached(params.key()?)?;
    let mode_count = modes_in_range(&model, params.freq_max);
    let n = params.grid_n.clamp(4, 100);
    let truncated_above = (model.freq_limit < params.freq_max).then_some(model.freq_limit);
    let flatness = params.score != "coupling";
    let solver = if model.key.panel == Panel::Rectangle { "analytic" } else { "fea" };

    let modes = model.modes[..mode_count].iter()
        .map(|m| ModeInfo { m: m.label.map(|l| l.0), n: m.label.map(|l| l.1), freq: m.freq })
        .collect();
    let outline = model.outline.rings().iter()
        .map(|ring| ring.iter().map(|p| [p[0] / params.lx, p[1] / params.ly]).collect())
        .collect();

    // Classify grid cells by position relative to the panel edges.
    let margin = params.edge_margin();
    let region: Vec<u8> = (0..n * n).map(|i| {
        let p = [((i % n) as f64 + 0.5) / n as f64 * params.lx, ((i / n) as f64 + 0.5) / n as f64 * params.ly];
        if !model.outline.contains(p) {
            REGION_OUTSIDE
        } else if model.outline.distance_to_edge(p) < margin {
            REGION_MARGIN
        } else {
            REGION_SEARCH
        }
    }).collect();

    let empty = |grid: Vec<f64>| CalculationResult {
        grid: grid.clone(),
        grid_raw: grid,
        region: region.clone(),
        grid_n: n,
        outline: Vec::new(),
        modes: Vec::new(),
        mode_count: 0,
        optimal_x: 0.5,
        optimal_y: 0.5,
        optimal_score_raw: 0.0,
        bands: Vec::new(),
        response_opt: Vec::new(),
        raggedness_opt: 0.0,
        solver,
        truncated_above,
    };
    if mode_count == 0 {
        let grid = region.iter().map(|&r| if r == REGION_OUTSIDE { f64::NAN } else { 0.0 }).collect();
        return Ok(CalculationResult { outline, modes, ..empty(grid) });
    }

    // --- Couplings: each mode's shape averaged around the voice-coil ring ---
    let response_model = ResponseModel::new(&model, &params);
    let radius = params.exciter_d.max(0.0) / 2.0;
    let stencil = ring_stencil(n, radius / params.lx, radius / params.ly);
    let eval = model.grid_eval(n);
    let mut couplings = DMatrix::<f64>::zeros(n * n, response_model.mode_count);
    for k in 0..response_model.mode_count {
        let shape = eval.eval(k);
        for (cell, taps) in stencil.iter().enumerate() {
            couplings[(cell, k)] = if region[cell] == REGION_OUTSIDE {
                0.0
            } else {
                // Ring taps off the panel are skipped.
                let c = finite_mean(taps.iter().map(|&(i, w)| (shape[i], w)));
                if c.is_finite() { c } else { 0.0 }
            };
        }
    }

    // --- Score every grid cell ---
    // Flatness: raggedness of the damped response (lower is better).
    // Coupling: sum of |coupling| over modes up to freq_max (higher is better).
    let grid_raw: Vec<f64> = if flatness {
        let power = couplings.map(|c| c * c) * &response_model.transfer;  // cells × bands
        (0..n * n)
            .map(|cell| {
                if region[cell] == REGION_OUTSIDE {
                    return f64::NAN;
                }
                let row: Vec<f64> = power.row(cell).iter().copied().collect();
                response::raggedness(&response::to_db(&row))
            })
            .collect()
    } else {
        (0..n * n)
            .map(|cell| {
                if region[cell] == REGION_OUTSIDE {
                    return f64::NAN;
                }
                (0..mode_count).map(|k| couplings[(cell, k)].abs()).sum()
            })
            .collect()
    };
    // Higher is better for the search and the colour scale.
    let goodness = |raw: f64| if flatness { -raw } else { raw };

    // The colour range comes from the search region only: free edges move far
    // more than anywhere an exciter can go, and would otherwise flatten it.
    let mut max_score = f64::NEG_INFINITY;
    let mut min_score = f64::INFINITY;
    let mut opt = None;
    let mut opt_score = f64::NEG_INFINITY;
    for cell in 0..n * n {
        if region[cell] != REGION_SEARCH {
            continue;
        }
        let score = goodness(grid_raw[cell]);
        max_score = max_score.max(score);
        min_score = min_score.min(score);
        if score > opt_score {
            opt_score = score;
            opt = Some(cell);
        }
    }
    // A shape too small for the margin has no search region; fall back to
    // the cell furthest from the edges.
    let opt = opt.unwrap_or_else(|| {
        (0..n * n)
            .filter(|&c| region[c] != REGION_OUTSIDE)
            .max_by(|&a, &b| {
                let p = |c: usize| [((c % n) as f64 + 0.5) / n as f64 * params.lx, ((c / n) as f64 + 0.5) / n as f64 * params.ly];
                model.outline.distance_to_edge(p(a)).total_cmp(&model.outline.distance_to_edge(p(b)))
            })
            .unwrap_or(n * n / 2 + n / 2)
    });
    let (opt_x, opt_y) = (((opt % n) as f64 + 0.5) / n as f64, ((opt / n) as f64 + 0.5) / n as f64);

    // Normalise to [0, 1]
    let range = max_score - min_score;
    let grid = grid_raw.iter()
        .map(|&raw| {
            if !raw.is_finite() {
                f64::NAN
            } else if range.is_finite() && range > 1e-12 {
                (goodness(raw) - min_score) / range
            } else {
                0.5
            }
        })
        .collect();

    let response_opt = response_model.curve_at(&model, &params, opt_x, opt_y);
    let raggedness_opt = response::raggedness(&response_opt);

    Ok(CalculationResult {
        grid,
        optimal_score_raw: grid_raw[opt],
        grid_raw,
        region,
        grid_n: n,
        outline,
        modes,
        mode_count,
        optimal_x: opt_x,
        optimal_y: opt_y,
        bands: response_model.bands.centres.clone(),
        response_opt,
        raggedness_opt,
        solver,
        truncated_above,
    })
}

/// Shape of mode `index` (in the order `compute_heatmap` returned) sampled on
/// an n×n grid, row-major; null outside the panel. Used to draw node lines.
#[tauri::command(async)]
fn mode_shape(params: PanelParams, index: usize, n: usize) -> Result<Vec<f64>, String> {
    let model = solve_cached(params.key()?)?;
    if index >= modes_in_range(&model, params.freq_max) {
        return Err("Mode index out of range".into());
    }
    Ok(model.grid_eval(n.clamp(4, 200)).eval(index))
}

#[derive(Serialize)]
pub struct ResponseCurve {
    pub db: Vec<f64>,       // per band, same bands as compute_heatmap
    pub raggedness: f64,    // [dB]
}

/// Response with the exciter at a normalised panel position (x, y in 0..1).
#[tauri::command(async)]
fn response_at(params: PanelParams, x: f64, y: f64) -> Result<ResponseCurve, String> {
    let model = solve_cached(params.key()?)?;
    let p = [x.clamp(0.0, 1.0) * params.lx, y.clamp(0.0, 1.0) * params.ly];
    if !model.outline.contains(p) {
        return Err("That position isn't on the panel".into());
    }
    let response_model = ResponseModel::new(&model, &params);
    let db = response_model.curve_at(&model, &params, x.clamp(0.0, 1.0), y.clamp(0.0, 1.0));
    let raggedness = response::raggedness(&db);
    Ok(ResponseCurve { db, raggedness })
}

#[cfg_attr(mobile, tauri::mobile_entry_point)]
pub fn run() {
    tauri::Builder::default()
        .plugin(tauri_plugin_opener::init())
        .invoke_handler(tauri::generate_handler![compute_heatmap, mode_shape, response_at])
        .run(tauri::generate_context!())
        .expect("error while running tauri application");
}
