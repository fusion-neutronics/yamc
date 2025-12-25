use serde::{Deserialize, Serialize};
use crate::reaction_product::ReactionProduct;

/// Represents a single reaction channel (identified by ENDF/MT number) for a
/// specific nuclide at a given temperature.
///
/// The reaction may either provide its own truncated energy grid or rely on
/// the parent nuclide's top‑level temperature grid (offset by `threshold_idx`).
/// `cross_section` values correspond 1‑to‑1 with the reaction's effective
/// energy grid (either its own `energy` or a slice of the parent grid).
#[derive(Debug, Serialize, Deserialize, Clone)]
pub struct Reaction {
    /// Cross section values in barns for the reaction energy grid.
    pub cross_section: Vec<f64>,
    /// Index into the parent (top‑level) energy grid where this reaction becomes active.
    pub threshold_idx: usize,
    /// Interpolation flags (currently informational; retained from source data).
    pub interpolation: Vec<i32>,
    /// Reaction‑specific energy grid (may be empty until synthesized).
    #[serde(skip, default)]
    pub energy: Vec<f64>, // Reaction-specific energy grid
    /// ENDF/MT reaction identifier.
    pub mt_number: i32, // The MT number for this reaction
    /// Q-value of the reaction in eV (energy released/absorbed in the reaction)
    #[serde(default)]
    pub q_value: f64,
    /// Products emitted by this reaction (e.g., neutrons, photons, fragments)
    #[serde(default)]
    pub products: Vec<ReactionProduct>,
    /// Whether scattering is in center-of-mass frame (requires CM to LAB conversion)
    #[serde(default)]
    pub scatter_in_cm: bool,
}

impl Reaction {
    /// Returns the cross section value for a given neutron energy using linear interpolation.
    /// If the energy is below the grid, returns the first value.
    /// If above, returns the last value.
    /// Otherwise, performs linear interpolation between grid points.
    #[inline]
    pub fn cross_section_at(&self, energy: f64) -> Option<f64> {
        if self.energy.is_empty() || self.cross_section.is_empty() {
            return None;
        }

        let n = self.energy.len();

        // Handle boundary cases
        if energy <= self.energy[0] {
            return Some(self.cross_section[0]);
        }
        if energy >= self.energy[n - 1] {
            return Some(self.cross_section[n - 1]);
        }

        // Binary search for the interval
        match self.energy.binary_search_by(|e| e.partial_cmp(&energy).unwrap()) {
            Ok(idx) => Some(self.cross_section[idx]),
            Err(idx) => {
                // idx is the insertion point, so energy is between [idx-1] and [idx]
                let i = idx - 1;
                let e0 = self.energy[i];
                let e1 = self.energy[idx];
                let xs0 = self.cross_section[i];
                let xs1 = self.cross_section[idx];

                // Linear interpolation
                let t = (energy - e0) / (e1 - e0);
                Some(xs0 + t * (xs1 - xs0))
            }
        }
    }
}
#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_cross_section_at() {
        let reaction = Reaction {
            cross_section: vec![1.0, 2.0, 3.0, 4.0],
            threshold_idx: 0,
            interpolation: vec![],
            energy: vec![0.5, 1.0, 2.0, 5.0],
            mt_number: 102,
            q_value: 0.0,
            products: vec![],
            scatter_in_cm: false,
        };

        // Below grid
        assert_eq!(reaction.cross_section_at(0.1), Some(1.0));
        // Exact match
        assert_eq!(reaction.cross_section_at(1.0), Some(2.0));
        // Between grid points (linear interpolation: 2.0 + 0.5*(3.0-2.0) = 2.5)
        assert_eq!(reaction.cross_section_at(1.5), Some(2.5));
        // Above grid
        assert_eq!(reaction.cross_section_at(10.0), Some(4.0));
    }
}
