#[cfg(test)]
mod tests {
    use super::*;
    use crate::cell::Cell;
    use crate::tally::Filter;

    fn dummy_cell(cell_id: i32) -> Cell {
        // Minimal cell with just an ID for filter testing
            use crate::region::{Region, HalfspaceType};
            use std::sync::Arc;
            use crate::surface::{Surface, SurfaceKind, BoundaryType};
            let dummy_surface = Arc::new(Surface {
                surface_id: Some(1),
                kind: SurfaceKind::Plane { a: 1.0, b: 0.0, c: 0.0, d: 0.0 },
                boundary_type: BoundaryType::default(),
            });
            let region = Region::new_from_halfspace(HalfspaceType::Above(dummy_surface));
            Cell {
                cell_id: Some(cell_id as u32),
                name: None,
                region,
                material: None,
            }
    }

    #[test]
    fn test_score_event_direct_mt() {
        use std::sync::atomic::Ordering;
        let mut tally = Tally::new();
        tally.set_score(101); // Absorption
        tally.initialize_batches(1); // Start batch
        let cell = dummy_cell(1);
        assert!(tally.score_event(101, &cell, Some(42), 0));
        assert_eq!(tally.batch_data.lock().unwrap()[0].load(Ordering::Relaxed), 1, "Should increment for direct MT match");
    }

    #[test]
    fn test_score_event_mt4_inelastic_constituent() {
        use std::sync::atomic::Ordering;
        let mut tally = Tally::new();
        tally.set_score(4); // MT 4 (inelastic)
        tally.initialize_batches(1); // Start batch
        let cell = dummy_cell(1);
        assert!(tally.score_event(53, &cell, Some(42), 0)); // MT 53 is inelastic constituent
        assert_eq!(tally.batch_data.lock().unwrap()[0].load(Ordering::Relaxed), 1, "Should increment for MT 4 when constituent MT occurs");
    }

    #[test]
    fn test_score_event_inelastic_constituent_only() {
        use std::sync::atomic::Ordering;
        let mut tally = Tally::new();
        tally.set_score(53); // MT 53
        tally.initialize_batches(1); // Start batch
        let cell = dummy_cell(1);
        assert!(tally.score_event(53, &cell, Some(42), 0));
        assert_eq!(tally.batch_data.lock().unwrap()[0].load(Ordering::Relaxed), 1, "Should increment for direct constituent MT match");
        assert!(!tally.score_event(4, &cell, Some(42), 0));
        assert_eq!(tally.batch_data.lock().unwrap()[0].load(Ordering::Relaxed), 1, "Should not increment for MT 4 if tally is for constituent");
    }

    #[test]
    fn test_score_event_with_cell_filter() {
        use std::sync::atomic::Ordering;
        let mut tally = Tally::new();
        tally.set_score(101);
        tally.initialize_batches(1);
        // Add a cell filter that matches cell_id 1
        tally.filters.push(Filter::Cell(crate::tally::CellFilter { cell_id: 1 }));
        let cell = dummy_cell(1);
        assert!(tally.score_event(101, &cell, Some(42), 0));
        assert_eq!(tally.batch_data.lock().unwrap()[0].load(Ordering::Relaxed), 1, "Should increment when cell filter matches");
        let cell2 = dummy_cell(2);
        assert!(!tally.score_event(101, &cell2, Some(42), 0));
        assert_eq!(tally.batch_data.lock().unwrap()[0].load(Ordering::Relaxed), 1, "Should not increment when cell filter does not match");
    }
}
use std::fmt;
use std::collections::HashMap;
use std::sync::{Arc, Mutex};
use std::sync::atomic::{AtomicU32, AtomicU64, Ordering};
use crate::filters::{CellFilter, MaterialFilter};

/// Unified filter enum for tallies
#[derive(Debug, Clone, PartialEq)]
pub enum Filter {
    Cell(CellFilter),
    Material(MaterialFilter),
}

impl Filter {
    /// Get the type name of this filter for validation
    pub fn type_name(&self) -> &'static str {
        match self {
            Filter::Cell(_) => "CellFilter",
            Filter::Material(_) => "MaterialFilter",
        }
    }
}

/// Unified tally structure serving as both input specification and results container
/// Note: Cannot derive Clone because AtomicU64 is not cloneable (use Arc<Tally> for sharing)
#[derive(Debug)]
pub struct Tally {
    // Input specification fields
    pub id: Option<u32>,           // Optional tally ID
    pub name: Option<String>,      // Optional tally name
    pub score: i32,               // MT number to score (e.g., 101 for absorption)
    pub filters: Vec<Filter>, // Filters for this tally

    // Results fields (populated during simulation)
    pub units: String,          // Units (e.g., "particles", "collisions", "reactions")
    pub batch_data: Mutex<Arc<Vec<AtomicU64>>>,   // Atomic counts, Mutex only for init (lock-free during scoring!)
    // Statistics (calculated from batch_data, use Atomic types for thread-safety)
    pub mean: AtomicU64,              // Mean per source particle (stored as f64 bits)
    pub std_dev: AtomicU64,           // Standard deviation per source particle (stored as f64 bits)
    pub rel_error: AtomicU64,         // Relative error (coefficient of variation) (stored as f64 bits)
    pub n_batches: AtomicU32,         // Number of batches
    pub particles_per_batch: AtomicU32, // Particles per batch for normalization
}

impl Tally {
    /// Score a reaction event for this tally, including MT 4/inelastic constituent logic
    pub fn score_event(&self, reaction_mt: i32, cell: &crate::cell::Cell, material_id: Option<u32>, batch_index: usize) -> bool {
    let mut should_score = self.score == reaction_mt;
        // If this is an inelastic constituent (MT 50-91), also score for MT 4
        if (50..92).contains(&reaction_mt) && self.score == 4 {
            should_score = true;
        }
        // If this is a constituent absorption (explicit list and continuum), also score for MT 101
        let is_absorption_constituent =
            self.score == 101 && (
                matches!(reaction_mt,
                    102 | 103 | 104 | 105 | 106 | 107 | 108 | 109 | 111 | 112 | 113 | 114 |
                    115 | 116 | 117 | 155 | 182 | 191 | 192 | 193 | 197
                ) ||
                (600..650).contains(&reaction_mt) ||
                (650..700).contains(&reaction_mt) ||
                (700..750).contains(&reaction_mt) ||
                (750..800).contains(&reaction_mt) ||
                (800..850).contains(&reaction_mt)
            );
        if is_absorption_constituent && self.score == 101 {
            should_score = true;
        }
        if should_score {
            let passes_filters = if self.filters.is_empty() {
                true
            } else {
                self.filters.iter().all(|filter| match filter {
                    crate::tally::Filter::Cell(cell_filter) => {
                        cell.cell_id.map_or(false, |id| cell_filter.matches(id))
                    },
                    crate::tally::Filter::Material(material_filter) => {
                        material_filter.matches(material_id)
                    }
                })
            };
            if passes_filters {
                // Clone the Arc (cheap - just increments ref count) so we don't hold the Mutex lock
                let batch_data = self.batch_data.lock().unwrap().clone();
                if let Some(atomic) = batch_data.get(batch_index) {
                    atomic.fetch_add(1, Ordering::Relaxed);
                }
                return true;
            }
        }
        false
    }
    /// Create a new tally specification
    pub fn new() -> Self {
        Self {
            id: None,
            name: None,
            score: 0,
            filters: Vec::new(),
            units: String::new(),
            batch_data: Mutex::new(Arc::new(Vec::new())),
            mean: AtomicU64::new(0.0_f64.to_bits()),
            std_dev: AtomicU64::new(0.0_f64.to_bits()),
            rel_error: AtomicU64::new(0.0_f64.to_bits()),
            n_batches: AtomicU32::new(0),
            particles_per_batch: AtomicU32::new(0),
        }
    }

    /// Initialize batch_data for a given number of batches
    pub fn initialize_batches(&self, n_batches: usize) {
        let new_batch_data: Vec<AtomicU64> = (0..n_batches)
            .map(|_| AtomicU64::new(0))
            .collect();
        *self.batch_data.lock().unwrap() = Arc::new(new_batch_data);
    }

    /// Clone the tally specification (without batch_data) for creating a new working copy
    pub fn clone_spec(&self) -> Self {
        Self {
            id: self.id,
            name: self.name.clone(),
            score: self.score,
            filters: self.filters.clone(),
            units: self.units.clone(),
            batch_data: Mutex::new(Arc::new(Vec::new())),  // Empty - will be initialized separately
            mean: AtomicU64::new(0.0_f64.to_bits()),
            std_dev: AtomicU64::new(0.0_f64.to_bits()),
            rel_error: AtomicU64::new(0.0_f64.to_bits()),
            n_batches: AtomicU32::new(0),
            particles_per_batch: AtomicU32::new(0),
        }
    }

    /// Create a new tally with name and units (for simulation results)
    pub fn with_name_and_units(name: &str, units: &str) -> Self {
        Self {
            id: None,
            name: Some(name.to_string()),
            score: 0,
            filters: Vec::new(),
            units: units.to_string(),
            batch_data: Mutex::new(Arc::new(Vec::new())),
            mean: AtomicU64::new(0.0_f64.to_bits()),
            std_dev: AtomicU64::new(0.0_f64.to_bits()),
            rel_error: AtomicU64::new(0.0_f64.to_bits()),
            n_batches: AtomicU32::new(0),
            particles_per_batch: AtomicU32::new(0),
        }
    }
    
    /// Set the MT number to score
    pub fn set_score(&mut self, mt_number: i32) {
        self.score = mt_number;
    }
    
    /// Validate that the tally configuration is valid
    pub fn validate(&self) -> Result<(), String> {
        // Check for duplicate filter types
        let mut filter_type_counts: HashMap<&str, u32> = HashMap::new();
        
        for filter in &self.filters {
            let type_name = filter.type_name();
            *filter_type_counts.entry(type_name).or_insert(0) += 1;
        }
        
        for (filter_type, count) in filter_type_counts {
            if count > 1 {
                return Err(format!(
                    "Tally '{}' has {} filters of type {}. Multiple filters of the same type are not allowed as they create impossible conditions (particles cannot be in multiple cells or materials simultaneously).", 
                    self.display_name(),
                    count,
                    filter_type
                ));
            }
        }
        
        Ok(())
    }

    /// Add a batch result to the tally (stores value in existing batch slot)
    pub fn add_batch(&self, batch_index: usize, count: u32, particles_per_batch: u32) {
        let batch_data = self.batch_data.lock().unwrap();
        if let Some(atomic) = batch_data.get(batch_index) {
            atomic.store(count as u64, Ordering::Relaxed);
        }
    }

    /// Update statistics from current batch data
    pub fn update_statistics(&self, particles_per_batch: u32) {
        let batch_data = self.batch_data.lock().unwrap();
        if batch_data.is_empty() || particles_per_batch == 0 {
            self.mean.store(0.0_f64.to_bits(), Ordering::Relaxed);
            self.std_dev.store(0.0_f64.to_bits(), Ordering::Relaxed);
            self.rel_error.store(0.0_f64.to_bits(), Ordering::Relaxed);
            self.n_batches.store(0, Ordering::Relaxed);
            self.particles_per_batch.store(particles_per_batch, Ordering::Relaxed);
            return;
        }

        let n = batch_data.len() as f64;
        let particles_per_batch_f64 = particles_per_batch as f64;

        // Normalize each batch by particles per batch to get per-source-particle values
        let normalized_data: Vec<f64> = batch_data
            .iter()
            .map(|atomic_count| atomic_count.load(Ordering::Relaxed) as f64 / particles_per_batch_f64)
            .collect();

        let mean_val = normalized_data.iter().sum::<f64>() / n;
        let variance = normalized_data
            .iter()
            .map(|x| (x - mean_val).powi(2))
            .sum::<f64>() / (n - 1.0).max(1.0);
        let std_dev_val = variance.sqrt();
        let rel_error_val = if mean_val > 0.0 { std_dev_val / mean_val } else { 0.0 };

        self.mean.store(mean_val.to_bits(), Ordering::Relaxed);
        self.std_dev.store(std_dev_val.to_bits(), Ordering::Relaxed);
        self.rel_error.store(rel_error_val.to_bits(), Ordering::Relaxed);
        self.n_batches.store(n as u32, Ordering::Relaxed);
        self.particles_per_batch.store(particles_per_batch, Ordering::Relaxed);
    }

    /// Get the total count across all batches
    pub fn total_count(&self) -> u64 {
        self.batch_data.lock().unwrap().iter().map(|a| a.load(Ordering::Relaxed)).sum()
    }

    /// Get the display name of the tally
    pub fn display_name(&self) -> String {
        self.name.clone().unwrap_or_else(|| "Unnamed Tally".to_string())
    }

    /// Helper to get mean as f64
    pub fn get_mean(&self) -> f64 {
        f64::from_bits(self.mean.load(Ordering::Relaxed))
    }

    /// Helper to get std_dev as f64
    pub fn get_std_dev(&self) -> f64 {
        f64::from_bits(self.std_dev.load(Ordering::Relaxed))
    }

    /// Helper to get rel_error as f64
    pub fn get_rel_error(&self) -> f64 {
        f64::from_bits(self.rel_error.load(Ordering::Relaxed))
    }
}

impl fmt::Display for Tally {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        writeln!(f, "Tally name: {}", self.display_name())?;
        writeln!(f, "  Mean: {:.6} per particle", self.get_mean())?;
        writeln!(f, "    Std Dev: {:.6} per particle", self.get_std_dev())?;
        writeln!(f, "    Rel Error: {:.4} ({:.2}%)", self.get_rel_error(), self.get_rel_error() * 100.0)?;
        writeln!(f, "    Batches: {}", self.n_batches.load(Ordering::Relaxed))?;
        writeln!(f, "    Particles per batch: {}", self.particles_per_batch.load(Ordering::Relaxed))?;
        writeln!(f, "  Total {}: {}", self.units.to_lowercase(), self.total_count())?;
        // Print batch_data values (not the AtomicU64 debug format)
        let batch_values: Vec<u64> = self.batch_data.lock().unwrap().iter().map(|a| a.load(Ordering::Relaxed)).collect();
        write!(f, "  Batch data: {:?}", batch_values)
    }
}

/// Convenience type alias for leakage tallies (backwards compatibility)
pub type LeakageTally = Tally;
pub type CountTally = Tally;

/// Initialize tallies from user tally specifications
pub fn create_tallies_from_specs(tally_specs: &[Tally]) -> Vec<Tally> {
    let mut tallies = Vec::new();
    
    // Always add leakage tally as the first tally
    tallies.push(Tally::with_name_and_units("Leakage", "particles"));
    
    // Add user-specified tallies
    for (i, spec) in tally_specs.iter().enumerate() {
        let name = spec.name.clone().unwrap_or_else(|| {
            format!("Tally {}", i + 1)
        });
        
        let units = "events";
        
        let mut tally = Tally::with_name_and_units(&name, units);
        tally.id = spec.id;
        tally.score = spec.score;
        tally.filters = spec.filters.clone();
        
        // Validate the tally configuration
        if let Err(err) = tally.validate() {
            panic!("Invalid tally configuration: {}", err);
        }
        
        tallies.push(tally);
    }
    
    tallies
}