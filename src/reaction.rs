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
}

impl Reaction {
    /// Returns the cross section value for a given neutron energy.
    /// Uses the interpolation type specified in the reaction data.
    pub fn cross_section_at(&self, energy: f64) -> Option<f64> {
        if self.energy.is_empty() || self.cross_section.is_empty() {
            return None;
        }

        // Only allow known interpolation types
        let interp_type = self.interpolation.get(0).copied().unwrap_or(2);
        match interp_type {
            2 => { // lin-lin interpolation
                Some(crate::utilities::interpolate_linear(&self.energy, &self.cross_section, energy))
            }
            // Extend here for other interpolation types (e.g., log-log)
            // 5 => Some(crate::utilities::interpolate_log_log(&self.energy, &self.cross_section, energy)),
            _ => {
                panic!("Unknown interpolation type {} in Reaction::cross_section_at. Please implement support for this type.", interp_type);
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
        };

        // Below grid
        assert_eq!(reaction.cross_section_at(0.1), Some(1.0));
        // Exact match
        assert_eq!(reaction.cross_section_at(1.0), Some(2.0));
        // Between grid points
        assert_eq!(reaction.cross_section_at(1.5), Some(2.0));
        // Above grid
        assert_eq!(reaction.cross_section_at(10.0), Some(4.0));
    }
}
