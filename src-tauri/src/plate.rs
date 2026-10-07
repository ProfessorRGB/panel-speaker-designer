// Plate modal solver.
//
// Free edges   : Rayleigh-Ritz with a tensor-product basis of orthonormal
//                Legendre polynomials. Free edges impose no geometric
//                constraints, so any complete basis is admissible, and the
//                Ritz method gives true free-edge mode shapes rather than the
//                cos·cos approximation.
// Simply supp. : closed-form sin·sin solution (exact for orthotropic plates
//                with no bending-twisting coupling).
//
// Both support orthotropic materials with principal axes along the panel
// edges (x = width, e.g. the grain direction of wood).

use nalgebra::{DMatrix, SymmetricEigen};
use std::f64::consts::PI;

#[derive(Clone, Copy, Debug, PartialEq)]
pub struct Plate {
    pub lx: f64,    // width  [m]
    pub ly: f64,    // height [m]
    pub h: f64,     // thickness [m]
    pub ex: f64,    // Young's modulus along x [Pa]
    pub ey: f64,    // Young's modulus along y [Pa]
    pub g: f64,     // in-plane shear modulus [Pa]
    pub nu_xy: f64, // major Poisson's ratio
    pub rho: f64,   // density [kg/m³]
}

#[derive(Clone, Copy, Debug, PartialEq)]
pub enum Boundary {
    Free,
    SimplySupported,
}

pub struct Rigidities {
    pub d11: f64,
    pub d22: f64,
    pub d12: f64,
    pub d66: f64,
}

impl Plate {
    pub fn validate(&self) -> Result<(), String> {
        let vals = [self.lx, self.ly, self.h, self.ex, self.ey, self.g, self.rho];
        if vals.iter().any(|v| !v.is_finite() || *v <= 0.0) {
            return Err("All dimensions and material constants must be positive".into());
        }
        // Positive-definite stiffness requires ν_xy·ν_yx < 1
        if !self.nu_xy.is_finite() || self.nu_xy < 0.0 || self.nu_xy * self.nu_xy * self.ey / self.ex >= 1.0 {
            return Err("Poisson's ratio is out of range for these moduli".into());
        }
        Ok(())
    }

    pub fn rigidities(&self) -> Rigidities {
        let nu_yx = self.nu_xy * self.ey / self.ex;
        let denom = 12.0 * (1.0 - self.nu_xy * nu_yx);
        let h3 = self.h.powi(3);
        let d22 = self.ey * h3 / denom;
        Rigidities {
            d11: self.ex * h3 / denom,
            d22,
            d12: self.nu_xy * d22,
            d66: self.g * h3 / 12.0,
        }
    }

    fn mass_per_area(&self) -> f64 {
        self.rho * self.h
    }

    // Upper bound on the number of half-waves along x and y for any mode at
    // or below `freq`: ω²ρh = D11·kx⁴ + 2(D12+2D66)·kx²ky² + D22·ky⁴ ≥ D11·kx⁴.
    fn max_half_waves(&self, freq: f64) -> (f64, f64) {
        let d = self.rigidities();
        let omega = 2.0 * PI * freq;
        let k4 = omega * omega * self.mass_per_area();
        (
            self.lx / PI * (k4 / d.d11).powf(0.25),
            self.ly / PI * (k4 / d.d22).powf(0.25),
        )
    }

    // Inverse of `max_half_waves`: frequency at which `m` half-waves fit along x.
    fn freq_for_half_waves_x(&self, m: f64) -> f64 {
        let k = m * PI / self.lx;
        (self.rigidities().d11 * k.powi(4) / self.mass_per_area()).sqrt() / (2.0 * PI)
    }

    fn freq_for_half_waves_y(&self, n: f64) -> f64 {
        let k = n * PI / self.ly;
        (self.rigidities().d22 * k.powi(4) / self.mass_per_area()).sqrt() / (2.0 * PI)
    }
}

// ── Basis sizing ─────────────────────────────────────────────────────────────

// Polynomial degree needed to resolve `m` half-waves to well under 1% error
// in frequency (verified by the convergence test below).
const BASIS_PER_HALF_WAVE: f64 = 1.6;
const BASIS_OFFSET: f64 = 10.0;
const BASIS_MIN: usize = 12;
const BASIS_MAX: usize = 48; // per axis — four parity blocks of at most 24×24 = 576 dof
const SS_MAX: usize = 400;

/// Discretisation used for a solve. Two solves with the same plate, boundary
/// and size produce identical modes, so this doubles as a cache key.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct SolveKey {
    pub plate: Plate,
    pub boundary: Boundary,
    pub nx: usize,
    pub ny: usize,
}

impl SolveKey {
    pub fn for_freq(plate: Plate, boundary: Boundary, freq_max: f64) -> SolveKey {
        let (mx, my) = plate.max_half_waves(freq_max);
        let (nx, ny) = match boundary {
            Boundary::Free => {
                let size = |m: f64| ((BASIS_PER_HALF_WAVE * m + BASIS_OFFSET).ceil() as usize)
                    .clamp(BASIS_MIN, BASIS_MAX);
                (size(mx), size(my))
            }
            Boundary::SimplySupported => (
                (mx.ceil() as usize + 1).min(SS_MAX),
                (my.ceil() as usize + 1).min(SS_MAX),
            ),
        };
        SolveKey { plate, boundary, nx, ny }
    }

    /// Highest frequency this discretisation resolves reliably.
    fn freq_limit(&self) -> f64 {
        match self.boundary {
            Boundary::Free => {
                let half_waves = |n: usize| (n as f64 - BASIS_OFFSET) / BASIS_PER_HALF_WAVE;
                self.plate.freq_for_half_waves_x(half_waves(self.nx))
                    .min(self.plate.freq_for_half_waves_y(half_waves(self.ny)))
            }
            Boundary::SimplySupported => self.plate.freq_for_half_waves_x(self.nx as f64)
                .min(self.plate.freq_for_half_waves_y(self.ny as f64)),
        }
    }
}

// ── Modes ────────────────────────────────────────────────────────────────────

pub enum Shape {
    Sine { m: u32, n: u32 },
    // Coefficients over basis functions φ_i(ξ)·ψ_j(η) with i ≡ px, j ≡ py (mod 2),
    // stored row-major as [i_sub * ny_sub + j_sub].
    Legendre { px: usize, py: usize, coeffs: Vec<f64> },
}

pub struct Mode {
    pub freq: f64,
    // Nodal-line counts along x and y. Exact for simply-supported plates; for
    // free plates an approximate label, since many modes mix several patterns.
    pub m: u32,
    pub n: u32,
    pub shape: Shape,
}

pub struct Solution {
    pub key: SolveKey,
    pub modes: Vec<Mode>,  // elastic modes, ascending frequency, all ≤ freq_limit
    pub freq_limit: f64,
}

pub fn solve(key: SolveKey) -> Solution {
    let freq_limit = key.freq_limit();
    let mut modes = match key.boundary {
        Boundary::Free => solve_free(&key),
        Boundary::SimplySupported => solve_simply_supported(&key),
    };
    modes.retain(|m| m.freq <= freq_limit);
    modes.sort_by(|a, b| a.freq.total_cmp(&b.freq));
    if key.boundary == Boundary::Free {
        let labeller = Labeller::new(&key);
        for mode in &mut modes {
            (mode.m, mode.n) = labeller.label(mode);
        }
    }
    Solution { key, modes, freq_limit }
}

fn solve_simply_supported(key: &SolveKey) -> Vec<Mode> {
    let p = &key.plate;
    let d = p.rigidities();
    let mut modes = Vec::new();
    for m in 1..=key.nx as u32 {
        for n in 1..=key.ny as u32 {
            let kx = m as f64 * PI / p.lx;
            let ky = n as f64 * PI / p.ly;
            let stiff = d.d11 * kx.powi(4)
                + 2.0 * (d.d12 + 2.0 * d.d66) * kx * kx * ky * ky
                + d.d22 * ky.powi(4);
            let freq = (stiff / p.mass_per_area()).sqrt() / (2.0 * PI);
            modes.push(Mode { freq, m, n, shape: Shape::Sine { m, n } });
        }
    }
    modes
}

// ── Legendre basis ───────────────────────────────────────────────────────────

/// Orthonormal Legendre polynomials φ_k = √((2k+1)/2)·P_k on [-1, 1] and their
/// first two derivatives, evaluated at each point: tables[k][point].
struct LegendreTable {
    val: Vec<Vec<f64>>,
    d1: Vec<Vec<f64>>,
    d2: Vec<Vec<f64>>,
}

fn legendre_table(count: usize, xs: &[f64]) -> LegendreTable {
    let np = xs.len();
    let mut val = vec![vec![0.0; np]; count];
    let mut d1 = vec![vec![0.0; np]; count];
    let mut d2 = vec![vec![0.0; np]; count];
    for (p, &x) in xs.iter().enumerate() {
        // Unnormalised P_k via Bonnet's recurrence; derivatives via
        // P'_{k+1} = P'_{k-1} + (2k+1)·P_k (and the same for P'').
        let (mut p0, mut p1) = (1.0, x);
        let (mut dp0, mut dp1) = (0.0, 1.0);
        let (mut ddp0, mut ddp1) = (0.0, 0.0);
        for k in 0..count {
            let (pk, dpk, ddpk) = if k == 0 { (p0, dp0, ddp0) } else { (p1, dp1, ddp1) };
            let s = ((2 * k + 1) as f64 / 2.0).sqrt();
            val[k][p] = s * pk;
            d1[k][p] = s * dpk;
            d2[k][p] = s * ddpk;
            if k >= 1 {
                let kf = k as f64;
                let p2 = ((2.0 * kf + 1.0) * x * p1 - kf * p0) / (kf + 1.0);
                let dp2 = dp0 + (2.0 * kf + 1.0) * p1;
                let ddp2 = ddp0 + (2.0 * kf + 1.0) * dp1;
                (p0, p1, dp0, dp1, ddp0, ddp1) = (p1, p2, dp1, dp2, ddp1, ddp2);
            }
        }
    }
    LegendreTable { val, d1, d2 }
}

/// Gauss-Legendre nodes and weights on [-1, 1].
fn gauss_legendre(q: usize) -> (Vec<f64>, Vec<f64>) {
    let mut nodes = vec![0.0; q];
    let mut weights = vec![0.0; q];
    for i in 0..q {
        let mut x = (PI * (i as f64 + 0.75) / (q as f64 + 0.5)).cos();
        let mut dp = 1.0;
        for _ in 0..100 {
            let (mut p0, mut p1) = (1.0, x);
            for k in 1..q {
                let kf = k as f64;
                let p2 = ((2.0 * kf + 1.0) * x * p1 - kf * p0) / (kf + 1.0);
                p0 = p1;
                p1 = p2;
            }
            dp = q as f64 * (x * p1 - p0) / (x * x - 1.0);
            let dx = p1 / dp;
            x -= dx;
            if dx.abs() < 1e-15 {
                break;
            }
        }
        nodes[i] = x;
        weights[i] = 2.0 / ((1.0 - x * x) * dp * dp);
    }
    (nodes, weights)
}

/// 1-D integrals ∫ φ_i^(a) φ_k^(b) dξ needed for plate bending energy.
struct Integrals {
    e11: Vec<Vec<f64>>, // ∫ φ_i' φ_k'
    e22: Vec<Vec<f64>>, // ∫ φ_i'' φ_k''
    e20: Vec<Vec<f64>>, // ∫ φ_i'' φ_k
}

fn integrals(count: usize) -> Integrals {
    // Integrands have degree < 2·count, so `count + 1` Gauss points are exact.
    let (xs, ws) = gauss_legendre(count + 1);
    let t = legendre_table(count, &xs);
    let integrate = |a: &Vec<Vec<f64>>, b: &Vec<Vec<f64>>| -> Vec<Vec<f64>> {
        (0..count)
            .map(|i| (0..count)
                .map(|k| ws.iter().enumerate().map(|(g, w)| w * a[i][g] * b[k][g]).sum())
                .collect())
            .collect()
    };
    Integrals {
        e11: integrate(&t.d1, &t.d1),
        e22: integrate(&t.d2, &t.d2),
        e20: integrate(&t.d2, &t.val),
    }
}

fn parity_indices(count: usize, parity: usize) -> Vec<usize> {
    (parity..count).step_by(2).collect()
}

fn solve_free(key: &SolveKey) -> Vec<Mode> {
    let p = &key.plate;
    let d = p.rigidities();
    let ix = integrals(key.nx);
    let iy = integrals(key.ny);

    // Map ξ = 2x/Lx − 1, η = 2y/Ly − 1. The Jacobian Lx·Ly/4 appears in both
    // stiffness and mass, so it cancels from ω².
    let ax2 = (2.0 / p.lx).powi(2);
    let ay2 = (2.0 / p.ly).powi(2);
    let c11 = d.d11 * ax2 * ax2;
    let c22 = d.d22 * ay2 * ay2;
    let c12 = d.d12 * ax2 * ay2;
    let c66 = 4.0 * d.d66 * ax2 * ay2;

    // Legendre polynomials have definite parity, and the rectangle is symmetric
    // about both centre lines, so the problem splits into four independent blocks.
    let mut modes = Vec::new();
    for px in 0..2 {
        for py in 0..2 {
            let xi = parity_indices(key.nx, px);
            let yi = parity_indices(key.ny, py);
            let size = xi.len() * yi.len();
            let k = DMatrix::from_fn(size, size, |r, c| {
                let (i, j) = (xi[r / yi.len()], yi[r % yi.len()]);
                let (kk, l) = (xi[c / yi.len()], yi[c % yi.len()]);
                // The basis is orthonormal, so ∫φ_iφ_k = δ_ik.
                let dx = (i == kk) as u8 as f64;
                let dy = (j == l) as u8 as f64;
                c11 * ix.e22[i][kk] * dy
                    + c22 * dx * iy.e22[j][l]
                    + c12 * (ix.e20[i][kk] * iy.e20[l][j] + ix.e20[kk][i] * iy.e20[j][l])
                    + c66 * ix.e11[i][kk] * iy.e11[j][l]
            });
            let eig = SymmetricEigen::new(k);
            for (idx, &lambda) in eig.eigenvalues.iter().enumerate() {
                let omega2 = lambda.max(0.0) / p.mass_per_area();
                let freq = omega2.sqrt() / (2.0 * PI);
                let mut coeffs: Vec<f64> = eig.eigenvectors.column(idx).iter().copied().collect();
                // Deterministic sign: largest coefficient positive.
                let big = coeffs.iter().copied().fold(0.0, |a: f64, v| if v.abs() > a.abs() { v } else { a });
                if big < 0.0 {
                    coeffs.iter_mut().for_each(|v| *v = -*v);
                }
                modes.push(Mode { freq, m: 0, n: 0, shape: Shape::Legendre { px, py, coeffs } });
            }
        }
    }

    // Drop the three rigid-body modes (translation and two rotations), which
    // have zero strain energy. Elastic modes sit far above this threshold.
    let d_min = d.d11.min(d.d22).min(d.d66);
    let f_ref = (d_min / p.mass_per_area()).sqrt() * (PI / p.lx.max(p.ly)).powi(2) / (2.0 * PI);
    modes.retain(|m| m.freq > 0.05 * f_ref);
    modes
}

// ── Evaluating mode shapes ───────────────────────────────────────────────────

/// Mode shape sampled on an n×n grid of cell centres, row-major with rows
/// along y. Shapes are normalised to RMS 0.5 over the panel (the RMS of a
/// unit-amplitude cos·cos or sin·sin mode), so modes are comparable.
pub fn shape_grid(sol: &Solution, mode: &Mode, n: usize) -> Vec<f64> {
    let centres: Vec<f64> = (0..n).map(|i| (i as f64 + 0.5) / n as f64).collect();
    match &mode.shape {
        Shape::Sine { m, n: nn } => {
            let sx: Vec<f64> = centres.iter().map(|x| (*m as f64 * PI * x).sin()).collect();
            let sy: Vec<f64> = centres.iter().map(|y| (*nn as f64 * PI * y).sin()).collect();
            sy.iter().flat_map(|y| sx.iter().map(move |x| x * y)).collect()
        }
        Shape::Legendre { px, py, coeffs } => {
            let xs: Vec<f64> = centres.iter().map(|t| 2.0 * t - 1.0).collect();
            let tx = legendre_table(sol.key.nx, &xs).val;
            let ty = legendre_table(sol.key.ny, &xs).val;
            separable_eval(coeffs, &parity_indices(sol.key.nx, *px), &parity_indices(sol.key.ny, *py), &tx, &ty)
        }
    }
}

/// Precomputed basis tables for evaluating many modes on the same grid.
pub struct GridBasis {
    n: usize,
    tx: Vec<Vec<f64>>,
    ty: Vec<Vec<f64>>,
}

impl GridBasis {
    pub fn new(sol: &Solution, n: usize) -> GridBasis {
        let xs: Vec<f64> = (0..n).map(|i| 2.0 * (i as f64 + 0.5) / n as f64 - 1.0).collect();
        let (tx, ty) = match sol.key.boundary {
            Boundary::Free => (legendre_table(sol.key.nx, &xs).val, legendre_table(sol.key.ny, &xs).val),
            Boundary::SimplySupported => (Vec::new(), Vec::new()),
        };
        GridBasis { n, tx, ty }
    }

    pub fn eval(&self, sol: &Solution, mode: &Mode) -> Vec<f64> {
        match &mode.shape {
            Shape::Sine { .. } => shape_grid(sol, mode, self.n),
            Shape::Legendre { px, py, coeffs } => separable_eval(
                coeffs,
                &parity_indices(sol.key.nx, *px),
                &parity_indices(sol.key.ny, *py),
                &self.tx,
                &self.ty,
            ),
        }
    }
}

// w[row][col] = Σ_i Σ_j c_ij · φ_i(x_col) · ψ_j(y_row), computed as two
// matrix products so cost is O(n·Nx·Ny + n²·Ny) rather than O(n²·Nx·Ny).
fn separable_eval(
    coeffs: &[f64],
    xi: &[usize],
    yi: &[usize],
    tx: &[Vec<f64>],
    ty: &[Vec<f64>],
) -> Vec<f64> {
    let nxp = tx.first().map_or(0, |r| r.len());
    let nyp = ty.first().map_or(0, |r| r.len());
    // partial[j][col] = Σ_i c_ij φ_i(x_col)
    let mut partial = vec![vec![0.0; nxp]; yi.len()];
    for (a, &i) in xi.iter().enumerate() {
        for (b, row) in partial.iter_mut().enumerate() {
            let c = coeffs[a * yi.len() + b];
            for (col, v) in row.iter_mut().enumerate() {
                *v += c * tx[i][col];
            }
        }
    }
    let mut out = vec![0.0; nxp * nyp];
    for (b, &j) in yi.iter().enumerate() {
        for row in 0..nyp {
            let w = ty[j][row];
            let dst = &mut out[row * nxp..(row + 1) * nxp];
            for (v, p) in dst.iter_mut().zip(&partial[b]) {
                *v += w * p;
            }
        }
    }
    out
}

// ── Mode labels ──────────────────────────────────────────────────────────────

// Counts sign changes along several lines parallel to each axis and takes the
// median, since lines close to a node line or a free edge can over- or
// under-count. The lines avoid the centre so they don't coincide with
// symmetry node lines.
const LABEL_LINES: [f64; 7] = [0.137, 0.233, 0.389, 0.461, 0.617, 0.743, 0.853];
const LABEL_SAMPLES: usize = 240;

struct Labeller {
    nx: usize,
    ny: usize,
    along_x: Vec<Vec<f64>>,  // φ_i at samples along a line
    along_y: Vec<Vec<f64>>,  // ψ_j at samples along a line
    at_lines_x: Vec<Vec<f64>>, // φ_i at the fixed x of each vertical line
    at_lines_y: Vec<Vec<f64>>, // ψ_j at the fixed y of each horizontal line
}

impl Labeller {
    fn new(key: &SolveKey) -> Labeller {
        let samples: Vec<f64> = (0..LABEL_SAMPLES)
            .map(|i| 2.0 * (i as f64 + 0.5) / LABEL_SAMPLES as f64 - 1.0)
            .collect();
        let lines: Vec<f64> = LABEL_LINES.iter().map(|t| 2.0 * t - 1.0).collect();
        Labeller {
            nx: key.nx,
            ny: key.ny,
            along_x: legendre_table(key.nx, &samples).val,
            along_y: legendre_table(key.ny, &samples).val,
            at_lines_x: legendre_table(key.nx, &lines).val,
            at_lines_y: legendre_table(key.ny, &lines).val,
        }
    }

    fn label(&self, mode: &Mode) -> (u32, u32) {
        let Shape::Legendre { px, py, coeffs } = &mode.shape else {
            return (mode.m, mode.n);
        };
        let xi = parity_indices(self.nx, *px);
        let yi = parity_indices(self.ny, *py);
        let c = |a: usize, b: usize| coeffs[a * yi.len() + b];
        let mut m_counts = Vec::with_capacity(LABEL_LINES.len());
        let mut n_counts = Vec::with_capacity(LABEL_LINES.len());
        for line in 0..LABEL_LINES.len() {
            // Horizontal line at y = y_line: w(x) = Σ_i φ_i(x) · (Σ_j c_ij ψ_j(y_line))
            let weights: Vec<f64> = (0..xi.len())
                .map(|a| yi.iter().enumerate().map(|(b, &j)| c(a, b) * self.at_lines_y[j][line]).sum())
                .collect();
            let row: Vec<f64> = (0..LABEL_SAMPLES)
                .map(|s| xi.iter().zip(&weights).map(|(&i, w)| w * self.along_x[i][s]).sum())
                .collect();
            m_counts.push(sign_changes(&row));

            let weights: Vec<f64> = (0..yi.len())
                .map(|b| xi.iter().enumerate().map(|(a, &i)| c(a, b) * self.at_lines_x[i][line]).sum())
                .collect();
            let col: Vec<f64> = (0..LABEL_SAMPLES)
                .map(|s| yi.iter().zip(&weights).map(|(&j, w)| w * self.along_y[j][s]).sum())
                .collect();
            n_counts.push(sign_changes(&col));
        }
        (median(m_counts), median(n_counts))
    }
}

fn median(mut v: Vec<u32>) -> u32 {
    v.sort_unstable();
    v[v.len() / 2]
}

fn sign_changes(values: &[f64]) -> u32 {
    let peak = values.iter().fold(0.0f64, |a, v| a.max(v.abs()));
    let threshold = 0.02 * peak;
    let mut last = 0.0f64;
    let mut count = 0;
    for &v in values {
        if v.abs() < threshold {
            continue;
        }
        if last != 0.0 && v.signum() != last {
            count += 1;
        }
        last = v.signum();
    }
    count
}

// ── Tests ────────────────────────────────────────────────────────────────────

#[cfg(test)]
mod tests {
    use super::*;

    fn isotropic(lx: f64, ly: f64, nu: f64) -> Plate {
        let e = 70e9;
        Plate { lx, ly, h: 0.002, ex: e, ey: e, g: e / (2.0 * (1.0 + nu)), nu_xy: nu, rho: 2700.0 }
    }

    // Non-dimensional frequency λ = ω·a²·√(ρh/D).
    fn lambda(p: &Plate, freq: f64) -> f64 {
        let d = p.rigidities().d11;
        2.0 * PI * freq * p.lx * p.lx * (p.mass_per_area() / d).sqrt()
    }

    fn solve_up_to(p: Plate, boundary: Boundary, count: usize) -> Solution {
        // Pick a cut-off comfortably above the first `count` modes.
        let mut f = 100.0;
        loop {
            let sol = solve(SolveKey::for_freq(p, boundary, f));
            if sol.modes.len() >= count {
                return sol;
            }
            f *= 2.0;
        }
    }

    #[test]
    fn free_square_matches_leissa() {
        // Leissa, "Vibration of Plates" (NASA SP-160), free square plate, ν = 0.3.
        let p = isotropic(0.3, 0.3, 0.3);
        let sol = solve_up_to(p, Boundary::Free, 6);
        let expected = [13.468, 19.596, 24.270, 34.801, 34.801];
        for (mode, want) in sol.modes.iter().zip(expected) {
            let got = lambda(&p, mode.freq);
            assert!((got - want).abs() / want < 0.005, "λ = {got:.3}, expected {want}");
        }
    }

    #[test]
    fn free_plate_excludes_rigid_body_modes() {
        let p = isotropic(0.3, 0.3, 0.3);
        let sol = solve_up_to(p, Boundary::Free, 3);
        assert!(lambda(&p, sol.modes[0].freq) > 5.0);
        // The first mode of a free square is the (1,1) twisting mode.
        assert_eq!((sol.modes[0].m, sol.modes[0].n), (1, 1));
    }

    #[test]
    fn simply_supported_fundamental_is_exact() {
        let p = isotropic(0.3, 0.3, 0.3);
        let sol = solve_up_to(p, Boundary::SimplySupported, 1);
        let got = lambda(&p, sol.modes[0].freq);
        assert!((got - 2.0 * PI * PI).abs() < 1e-9);
    }

    #[test]
    fn rotating_an_orthotropic_plate_preserves_frequencies() {
        let a = Plate { lx: 0.3, ly: 0.2, h: 0.003, ex: 10e9, ey: 2e9, g: 0.8e9, nu_xy: 0.3, rho: 600.0 };
        let b = Plate { lx: a.ly, ly: a.lx, ex: a.ey, ey: a.ex, nu_xy: a.nu_xy * a.ey / a.ex, ..a };
        let sa = solve(SolveKey::for_freq(a, Boundary::Free, 3000.0));
        let sb = solve(SolveKey::for_freq(b, Boundary::Free, 3000.0));
        for (ma, mb) in sa.modes.iter().zip(&sb.modes).take(12) {
            assert!((ma.freq - mb.freq).abs() / ma.freq < 1e-6, "{} vs {}", ma.freq, mb.freq);
            assert_eq!((ma.m, ma.n), (mb.n, mb.m));
        }
    }

    #[test]
    fn free_modes_are_converged_up_to_the_reported_limit() {
        let p = Plate { lx: 0.3, ly: 0.2, h: 0.003, ex: 3.2e9, ey: 3.2e9, g: 1.17e9, nu_xy: 0.37, rho: 1190.0 };
        let coarse = solve(SolveKey::for_freq(p, Boundary::Free, 5000.0));
        let fine = solve(SolveKey { nx: coarse.key.nx + 12, ny: coarse.key.ny + 12, ..coarse.key });
        for (c, f) in coarse.modes.iter().zip(&fine.modes) {
            assert!((c.freq - f.freq) / f.freq < 0.01, "{} vs {}", c.freq, f.freq);
        }
    }

    #[test]
    fn shapes_are_normalised() {
        let p = isotropic(0.3, 0.2, 0.3);
        let sol = solve_up_to(p, Boundary::Free, 10);
        let n = 200;
        for mode in &sol.modes {
            let w = shape_grid(&sol, mode, n);
            let rms = (w.iter().map(|v| v * v).sum::<f64>() / w.len() as f64).sqrt();
            assert!((rms - 0.5).abs() < 0.01, "rms {rms}");
        }
    }
}
