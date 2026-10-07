// Damped modal response and placement flatness.
//
// With the exciter at x0 driving with unit force, the panel's spatially
// averaged mean-square velocity is
//
//     <v²>(ω) ∝ Σ_k c_k² · ω² / |ω_k² − ω² + iηω_k²|²
//
// where c_k is the exciter's coupling to mode k (its shape averaged around the
// voice-coil ring). Cross terms vanish in the spatial average because modes
// are orthogonal, and every mode shape has the same normalisation, so modal
// mass is a common factor. This tracks how evenly the panel is driven across
// frequency; it is not the radiated sound pressure (that needs radiation
// efficiency — see roadmap phase E).

use nalgebra::DMatrix;
use std::f64::consts::PI;

pub const BANDS_PER_OCTAVE: usize = 12;

// The flatness reference is the response smoothed over one octave (±½ octave).
const SMOOTH_HALF_WIDTH: usize = BANDS_PER_OCTAVE / 2;

/// Fractional-octave bands, contiguous, from `f_start` up to at most `f_end`.
pub struct Bands {
    pub centres: Vec<f64>,
    lo: Vec<f64>,
    hi: Vec<f64>,
}

impl Bands {
    pub fn new(f_start: f64, f_end: f64) -> Bands {
        let step = 2f64.powf(1.0 / BANDS_PER_OCTAVE as f64);
        let (mut centres, mut lo, mut hi) = (Vec::new(), Vec::new(), Vec::new());
        let mut a = f_start;
        while a * step <= f_end {
            lo.push(a);
            hi.push(a * step);
            centres.push(a * step.sqrt());
            a *= step;
        }
        Bands { centres, lo, hi }
    }

    pub fn len(&self) -> usize {
        self.centres.len()
    }
}

/// Band-averaged velocity response of each mode for unit coupling:
/// a (modes × bands) matrix T, so the response at a point is c² · T.
///
/// Averaging over the band (rather than sampling at its centre) makes the
/// result independent of how sharp the peaks are relative to the band width.
/// Each mode is integrated as a Lorentzian, which is exact near resonance;
/// the factor ω_c²/(ω_k + ω_c)² restores the correct far-off-resonance tails
/// (ω²/ω_k⁴ well below, 1/ω² well above).
pub fn transfer(mode_freqs: &[f64], eta: f64, bands: &Bands) -> DMatrix<f64> {
    let eta = eta.max(1e-4);
    DMatrix::from_fn(mode_freqs.len(), bands.len(), |k, b| {
        let wk = 2.0 * PI * mode_freqs[k];
        let (wa, wb, wc) = (2.0 * PI * bands.lo[b], 2.0 * PI * bands.hi[b], 2.0 * PI * bands.centres[b]);
        let gamma = eta * wk / 2.0;
        let lorentz = (((wb - wk) / gamma).atan() - ((wa - wk) / gamma).atan()) / (gamma * (wb - wa));
        lorentz * wc * wc / ((wk + wc) * (wk + wc))
    })
}

pub fn to_db(power: &[f64]) -> Vec<f64> {
    power.iter().map(|p| 10.0 * p.max(1e-300).log10()).collect()
}

/// RMS deviation (dB) of a response from its own one-octave-smoothed trend.
/// Placement barely changes the overall trend but strongly changes the peaks
/// and dips, so this isolates what the exciter position controls.
pub fn raggedness(db: &[f64]) -> f64 {
    let n = db.len();
    if n < 3 {
        return 0.0;
    }
    let sum_sq: f64 = (0..n).map(|i| {
        // Symmetric window, narrowed at the ends so a straight slope has no
        // deviation anywhere.
        let half = SMOOTH_HALF_WIDTH.min(i).min(n - 1 - i);
        let window = &db[i - half..=i + half];
        let trend = window.iter().sum::<f64>() / window.len() as f64;
        (db[i] - trend).powi(2)
    }).sum();
    (sum_sq / n as f64).sqrt()
}

#[cfg(test)]
mod tests {
    use super::*;

    // Exact band average of ω²/|ω_k² − ω² + iηω_k²|² by numerical integration.
    fn exact_band_average(fk: f64, eta: f64, lo: f64, hi: f64) -> f64 {
        let wk = 2.0 * PI * fk;
        let steps = 200_000;
        let (wa, wb) = (2.0 * PI * lo, 2.0 * PI * hi);
        let dw = (wb - wa) / steps as f64;
        (0..steps).map(|i| {
            let w = wa + (i as f64 + 0.5) * dw;
            let re = wk * wk - w * w;
            let im = eta * wk * wk;
            w * w / (re * re + im * im)
        }).sum::<f64>() / steps as f64
    }

    #[test]
    fn band_transfer_matches_exact_integral() {
        let bands = Bands::new(50.0, 2000.0);
        let mode = [300.0];
        for eta in [0.01, 0.05] {
            let t = transfer(&mode, eta, &bands);
            for b in 0..bands.len() {
                let exact = exact_band_average(mode[0], eta, bands.lo[b], bands.hi[b]);
                let err_db = 10.0 * (t[(0, b)] / exact).log10();
                assert!(err_db.abs() < 0.5, "band {} Hz, η {eta}: {err_db:.2} dB", bands.centres[b]);
            }
        }
    }

    #[test]
    fn bands_are_contiguous_and_within_range() {
        let bands = Bands::new(100.0, 1000.0);
        assert_eq!(bands.len(), 39); // 3.32 octaves × 12, whole bands only
        assert!(*bands.hi.last().unwrap() <= 1000.0);
        for b in 1..bands.len() {
            assert!((bands.lo[b] - bands.hi[b - 1]).abs() < 1e-9);
        }
    }

    #[test]
    fn raggedness_ignores_trend_but_sees_dips() {
        let sloped: Vec<f64> = (0..48).map(|i| -0.25 * i as f64).collect();
        assert!(raggedness(&sloped) < 0.2);
        let mut dipped = sloped.clone();
        dipped[24] -= 12.0;
        assert!(raggedness(&dipped) > 1.5);
    }
}
