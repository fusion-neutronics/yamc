
use rand::Rng;
use serde::{Serialize, Deserialize};

#[derive(Debug, Clone, Serialize, Deserialize)]
pub enum Interpolation {
    Histogram,
    LinLin,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct KalbachMann {
    pub energy: Vec<f64>, // incident energy grid
    pub distributions: Vec<KMTable>,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct KMTable {
    pub interpolation: Interpolation,
    pub n_discrete: usize,
    pub e_out: Vec<f64>,
    pub p: Vec<f64>,
    pub c: Vec<f64>,
    pub r: Vec<f64>,
    pub a: Vec<f64>,
}

impl KalbachMann {
    pub fn sample<R: Rng>(&self, e_in: f64, rng: &mut R) -> (f64, f64) {
        let n_energy = self.energy.len();
        if n_energy < 2 {
            return (e_in, 2.0 * rng.gen::<f64>() - 1.0);
        }
        let (i, r_interp) = if e_in < self.energy[0] {
            (0, 0.0)
        } else if e_in > self.energy[n_energy - 1] {
            (n_energy - 2, 1.0)
        } else {
            let i = lower_bound_index(&self.energy, e_in);
            let r = (e_in - self.energy[i]) / (self.energy[i + 1] - self.energy[i]);
            (i, r)
        };
        let l = if r_interp > rng.gen::<f64>() { i + 1 } else { i };
        let dist = &self.distributions[l];
        let n_energy_out = dist.e_out.len();
        let n_discrete = dist.n_discrete;
        // Outgoing energy sampling (OpenMC style)
        let e_out = {
            let r1 = rng.gen::<f64>();
            let mut k = 0;
            let mut c_k = dist.c[0];
            let mut end = n_energy_out - 2;
            // Discrete portion
            for j in 0..n_discrete {
                k = j;
                c_k = dist.c[k];
                if r1 < c_k {
                    end = j;
                    break;
                }
            }
            // Continuous portion
            let mut c_k1 = 0.0;
            for j in n_discrete..end {
                k = j;
                c_k1 = dist.c[k + 1];
                if r1 < c_k1 {
                    break;
                }
                k = j + 1;
                c_k = c_k1;
            }
            let e_l_k = dist.e_out[k];
            let p_l_k = dist.p[k];
            match dist.interpolation {
                Interpolation::Histogram => {
                    if p_l_k > 0.0 && k >= n_discrete {
                        e_l_k + (r1 - c_k) / p_l_k
                    } else {
                        e_l_k
                    }
                }
                Interpolation::LinLin => {
                    let e_l_k1 = dist.e_out[k + 1];
                    let p_l_k1 = dist.p[k + 1];
                    let frac = (p_l_k1 - p_l_k) / (e_l_k1 - e_l_k);
                    if frac == 0.0 {
                        e_l_k + (r1 - c_k) / p_l_k
                    } else {
                        e_l_k + ((p_l_k * p_l_k + 2.0 * frac * (r1 - c_k)).max(0.0).sqrt() - p_l_k) / frac
                    }
                }
            }
        };
        // Find outgoing energy bin
        let mut k = 0;
        for j in 0..n_energy_out - 1 {
            if e_out >= dist.e_out[j] && e_out <= dist.e_out[j + 1] {
                k = j;
                break;
            }
        }
        // Interpolate r and a for e_out
        let (r_val, a_val) = if k + 1 < dist.e_out.len() {
            let e0 = dist.e_out[k];
            let e1 = dist.e_out[k + 1];
            let f = if e1 > e0 { (e_out - e0) / (e1 - e0) } else { 0.0 };
            let r_interp = dist.r[k] * (1.0 - f) + dist.r[k + 1] * f;
            let a_interp = dist.a[k] * (1.0 - f) + dist.a[k + 1] * f;
            (r_interp, a_interp)
        } else {
            (dist.r[k], dist.a[k])
        };
        // Mixture: with probability r, sample precompound; else compound
        let mu = if rng.gen::<f64>() > r_val {
            // Compound: exp(a*mu) on [-1,1]
            if a_val.abs() < 1e-6 {
                2.0 * rng.gen::<f64>() - 1.0
            } else {
                let r = rng.gen::<f64>();
                if a_val > 0.0 {
                    ((1.0 - r) * (-a_val).exp() + r * a_val.exp()).ln() / a_val
                } else {
                    -((1.0 - r) * a_val.exp() + r * (-a_val).exp()).ln() / a_val
                }
            }
        } else {
            // Precompound: sample as in OpenMC
            let T = (2.0 * rng.gen::<f64>() - 1.0) * a_val.sinh();
            (T + (T * T + 1.0).sqrt()).ln() / a_val
        };
        (e_out, mu.max(-1.0).min(1.0))
    }
}

fn lower_bound_index(grid: &[f64], value: f64) -> usize {
    for i in 0..grid.len() - 1 {
        if value < grid[i + 1] {
            return i;
        }
    }
    grid.len() - 2
}

// Tests should be rewritten to use the new KalbachMann struct and its sample method.
