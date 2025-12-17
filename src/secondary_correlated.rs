use serde::{Serialize, Deserialize};
use crate::reaction_product::TabulatedProbability;
use rand::Rng;

/// Interpolation type for tabular distributions
#[derive(Clone, Copy, Debug, PartialEq, Eq, Serialize, Deserialize)]
pub enum Interpolation {
    Histogram,
    LinLin,
}

/// Tabular distribution for outgoing energy or angle
#[derive(Clone, Debug, Serialize, Deserialize)]
pub struct Tabular {
    pub x: Vec<f64>,
    pub p: Vec<f64>,
    pub c: Vec<f64>,
    pub interpolation: Interpolation,
    pub n_discrete: usize, // Number of discrete lines (for outgoing energy)
}

impl Tabular {
    /// Sample from the tabular distribution (ACE/OpenMC style)
    pub fn sample<R: Rng>(&self, rng: &mut R) -> f64 {
        let r = rng.gen::<f64>();
        let n = self.x.len();
        let mut k = 0;
        let mut c_k = self.c[0];
        let mut end = n - 2;
        // Discrete portion
        for j in 0..self.n_discrete {
            k = j;
            c_k = self.c[k];
            if r < c_k {
                end = j;
                break;
            }
        }
        // Continuous portion
        let mut c_k1 = 0.0;
        for j in self.n_discrete..end {
            k = j;
            c_k1 = self.c[k + 1];
            if r < c_k1 {
                break;
            }
            k = j + 1;
            c_k = c_k1;
        }
        let x_k = self.x[k];
        let p_k = self.p[k];
        if self.interpolation == Interpolation::Histogram {
            if p_k > 0.0 && k >= self.n_discrete {
                x_k + (r - c_k) / p_k
            } else {
                x_k
            }
        } else {
            // LinLin
            let x_k1 = self.x[k + 1];
            let p_k1 = self.p[k + 1];
            let frac = (p_k1 - p_k) / (x_k1 - x_k);
            if frac == 0.0 {
                x_k + (r - c_k) / p_k
            } else {
                x_k + ((p_k * p_k + 2.0 * frac * (r - c_k)).max(0.0).sqrt() - p_k) / frac
            }
        }
    }
}

/// Correlated angle-energy distribution (OpenMC style)
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct CorrelatedAngleEnergy {
    pub energy: Vec<f64>, // incident energy grid
    pub distributions: Vec<CorrTable>,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct CorrTable {
    pub interpolation: Interpolation,
    pub n_discrete: usize,
    pub e_out: Vec<f64>,
    pub p: Vec<f64>,
    pub c: Vec<f64>,
    pub angle: Vec<Tabular>, // angular distributions for each outgoing energy
}

impl CorrelatedAngleEnergy {
    /// Sample outgoing energy and angle given incident energy
    pub fn sample<R: Rng>(&self, e_in: f64, rng: &mut R) -> (f64, f64) {
        let n_energy = self.energy.len();
        if n_energy < 2 {
            // Not enough data, fallback
            return (e_in, 2.0 * rng.gen::<f64>() - 1.0);
        }
        let (i, r) = if e_in < self.energy[0] {
            (0, 0.0)
        } else if e_in > self.energy[n_energy - 1] {
            (n_energy - 2, 1.0)
        } else {
            let i = lower_bound_index(&self.energy, e_in);
            let r = (e_in - self.energy[i]) / (self.energy[i + 1] - self.energy[i]);
            (i, r)
        };
        let l = if r > rng.gen::<f64>() { i + 1 } else { i };
        let dist = &self.distributions[l];
        let n_energy_out = dist.e_out.len();
        let n_discrete = dist.n_discrete;
        let e_out = {
            let tab = Tabular {
                x: dist.e_out.clone(),
                p: dist.p.clone(),
                c: dist.c.clone(),
                interpolation: dist.interpolation,
                n_discrete: dist.n_discrete,
            };
            tab.sample(rng)
        };
        // Find outgoing energy bin
        let mut k = 0;
        for j in 0..n_energy_out - 1 {
            if e_out >= dist.e_out[j] && e_out <= dist.e_out[j + 1] {
                k = j;
                break;
            }
        }
        // Interpolate mu between angle[k] and angle[k+1] if needed
        let mu = if k + 1 < dist.angle.len() {
            let mu0 = dist.angle[k].sample(rng);
            let mu1 = dist.angle[k + 1].sample(rng);
            let e0 = dist.e_out[k];
            let e1 = dist.e_out[k + 1];
            let f = if e1 > e0 { (e_out - e0) / (e1 - e0) } else { 0.0 };
            mu0 * (1.0 - f) + mu1 * f
        } else {
            dist.angle[k].sample(rng)
        };
        (e_out, mu.max(-1.0).min(1.0))
    }
}


/// Find the lower bound index for interpolation
fn lower_bound_index(grid: &[f64], value: f64) -> usize {
    for i in 0..grid.len() - 1 {
        if value < grid[i + 1] {
            return i;
        }
    }
    grid.len() - 2
}
