use std::collections::HashMap;
use std::fmt;
use std::sync::atomic::{AtomicU32, AtomicU64, Ordering};
use std::sync::{Arc, Mutex};
// ...existing code...
use crate::filter::Filter;

/// Special score value for flux (track-length estimator)

/// Represents a score - either an MT number or a named score like "flux"
#[derive(Debug, Clone, PartialEq)]
pub enum Score {
    MT(i32),
    Flux,
}

impl Score {
    /// Convert to i32 representation
    pub fn to_i32(&self) -> i32 {
        match self {
            Score::MT(mt) => *mt,
            Score::Flux => panic!("Direct integer value for flux score is not supported. Use 'flux' string only."),
        }
    }
    
    /// Parse from string
    pub fn from_str(s: &str) -> Result<Self, String> {
        match s {
            "flux" => Ok(Score::Flux),
            _ => Err(format!("Unknown score: '{}'", s)),
        }
    }
}

/// Unified tally structure serving as both input specification and results container
/// Note: Cannot derive Clone because AtomicU64 is not cloneable (use Arc<Tally> for sharing)
#[derive(Debug)]
pub struct Tally {
    // Input specification fields
    pub id: Option<u32>,      // Optional tally ID
    pub name: Option<String>, // Optional tally name
    pub scores: Vec<Score>,   // Score enum (MT or Flux)
    pub filters: Vec<Filter>, // Filters for this tally

    // Results fields (populated during simulation)
    pub units: String, // Units (e.g., "particles", "collisions", "reactions")
    // Mutex ONLY for initialization; scoring is lock-free via Arc clones
    // Length = num_scores * num_energy_bins (or just num_scores if no energy filter)
    pub batch_data: Vec<Mutex<Arc<Vec<AtomicU64>>>>, // One per score*energy_bin, Mutex only for init
    // Statistics (calculated from batch_data, use Atomic types for thread-safety)
    pub mean: Vec<AtomicU64>,      // Mean per score*energy_bin (stored as f64 bits)
    pub std_dev: Vec<AtomicU64>,   // Std dev per score*energy_bin
    pub rel_error: Vec<AtomicU64>, // Rel error per score*energy_bin
    pub n_batches: AtomicU32,      // Number of batches
    pub particles_per_batch: AtomicU32, // Particles per batch for normalization
}

impl Tally {
    /// Get the number of energy bins from the energy filter (1 if no energy filter)
    pub fn num_energy_bins(&self) -> usize {
        self.filters
            .iter()
            .find_map(|f| {
                if let Filter::Energy(ef) = f {
                    Some(ef.num_bins())
                } else {
                    None
                }
            })
            .unwrap_or(1)
    }

    /// Get the flat index for a score and energy bin
    /// Returns None if out of range
    pub fn get_bin_index(&self, score_index: usize, energy_bin: usize) -> Option<usize> {
        let num_energy_bins = self.num_energy_bins();
        if score_index >= self.scores.len() || energy_bin >= num_energy_bins {
            return None;
        }
        Some(score_index * num_energy_bins + energy_bin)
    }

    /// Get the energy filter if present
    pub fn get_energy_filter(&self) -> Option<&crate::filters::EnergyFilter> {
        self.filters.iter().find_map(|f| {
            if let Filter::Energy(ef) = f {
                Some(ef)
            } else {
                None
            }
        })
    }
    /// Score a reaction event for this tally, including MT 4/inelastic constituent logic
    pub fn score_event(
        &self,
        reaction_mt: i32,
        cell: &crate::cell::Cell,
        material_id: Option<u32>,
        energy: f64,
        batch_index: usize,
    ) -> bool {
        let mut triggered = false;
        
        // Get the energy bin if there's an energy filter
        let energy_bin = if let Some(energy_filter) = self.get_energy_filter() {
            match energy_filter.get_bin(energy) {
                Some(bin) => bin,
                None => return false, // Energy outside filter range, don't score
            }
        } else {
            0 // No energy filter, use bin 0
        };
        
        // For each score, check if it should be incremented
        for (score_idx, score) in self.scores.iter().enumerate() {
            let mut should_score = false;
            match score {
                Score::MT(mt) => {
                    if *mt == reaction_mt {
                        should_score = true;
                    }
                    // If this is an inelastic constituent (MT 50-91), also score for MT 4
                    if (50..92).contains(&reaction_mt) && *mt == 4 {
                        should_score = true;
                    }
                    // If this is a constituent absorption (explicit list and continuum), also score for MT 101
                    let is_absorption_constituent = *mt == 101
                        && (matches!(
                            reaction_mt,
                            102 | 103
                                | 104
                                | 105
                                | 106
                                | 107
                                | 108
                                | 109
                                | 111
                                | 112
                                | 113
                                | 114
                                | 115
                                | 116
                                | 117
                                | 155
                                | 182
                                | 191
                                | 192
                                | 193
                                | 197
                        ) || (600..650).contains(&reaction_mt)
                            || (650..700).contains(&reaction_mt)
                            || (700..750).contains(&reaction_mt)
                            || (750..800).contains(&reaction_mt)
                            || (800..850).contains(&reaction_mt));
                    if is_absorption_constituent {
                        should_score = true;
                    }
                }
                Score::Flux => {
                    // No event-based scoring for flux
                }
            }
            if should_score {
                let passes_filters = if self.filters.is_empty() {
                    true
                } else {
                    self.filters.iter().all(|filter| match filter {
                        Filter::Cell(cell_filter) => {
                            cell.cell_id.map_or(false, |id| cell_filter.matches(id))
                        }
                        Filter::Material(material_filter) => material_filter.matches(material_id),
                        Filter::Energy(_) => true, // Energy already checked above
                    })
                };
                if passes_filters {
                    // Get the flat bin index combining score and energy bin
                    if let Some(bin_idx) = self.get_bin_index(score_idx, energy_bin) {
                        // Lock-free: clone the Arc to get shared access to atomics
                        let batch_arc = self.batch_data[bin_idx].lock().unwrap().clone();
                        if let Some(atomic) = batch_arc.get(batch_index) {
                            atomic.fetch_add(1, Ordering::Relaxed);
                            triggered = true;
                        }
                    }
                }
            }
        }
        triggered
    }

    /// Score a track-length contribution for flux tallies
    /// This is called every time a particle moves through a cell
    pub fn score_track_length(
        &self,
        track_length: f64,
        cell: &crate::cell::Cell,
        material_id: Option<u32>,
        energy: f64,
        batch_index: usize,
    ) {
        // Get the energy bin if there's an energy filter
        let energy_bin = if let Some(energy_filter) = self.get_energy_filter() {
            match energy_filter.get_bin(energy) {
                Some(bin) => bin,
                None => return, // Energy outside filter range, don't score
            }
        } else {
            0 // No energy filter, use bin 0
        };
        
        // For each score, check if it's a flux score
        for (score_idx, score) in self.scores.iter().enumerate() {
            if let Score::Flux = score {
                let passes_filters = if self.filters.is_empty() {
                    true
                } else {
                    self.filters.iter().all(|filter| match filter {
                        Filter::Cell(cell_filter) => {
                            cell.cell_id.map_or(false, |id| cell_filter.matches(id))
                        }
                        Filter::Material(material_filter) => material_filter.matches(material_id),
                        Filter::Energy(_) => true, // Energy already checked above
                    })
                };
                if passes_filters {
                    // Get the flat bin index combining score and energy bin
                    if let Some(bin_idx) = self.get_bin_index(score_idx, energy_bin) {
                        // Lock-free: clone the Arc to get shared access to atomics
                        let batch_arc = self.batch_data[bin_idx].lock().unwrap().clone();
                        if let Some(atomic) = batch_arc.get(batch_index) {
                            // Atomically add to floating point accumulator using f64 bits
                            // This is lock-free but requires a compare-exchange loop
                            let mut current = atomic.load(Ordering::Relaxed);
                            loop {
                                let current_f64 = f64::from_bits(current);
                                let new_f64 = current_f64 + track_length;
                                let new_bits = new_f64.to_bits();
                                match atomic.compare_exchange_weak(
                                    current,
                                    new_bits,
                                    Ordering::Relaxed,
                                    Ordering::Relaxed,
                                ) {
                                    Ok(_) => break,
                                    Err(x) => current = x,
                                }
                            }
                        }
                    }
                }
            }
        }
    }

    /// Create a new tally specification
    pub fn new() -> Self {
        Self {
            id: None,
            name: None,
            scores: Vec::new(),
            filters: Vec::new(),
            units: String::new(),
            batch_data: Vec::new(),
            mean: Vec::new(),
            std_dev: Vec::new(),
            rel_error: Vec::new(),
            n_batches: AtomicU32::new(0),
            particles_per_batch: AtomicU32::new(0),
        }
    }

    /// Initialize batch_data for a given number of batches (mutable version)
    pub fn initialize_batches(&mut self, n_batches: usize) {
        let num_bins = self.scores.len() * self.num_energy_bins();
        
        self.batch_data = (0..num_bins)
            .map(|_| {
                Mutex::new(Arc::new(
                    (0..n_batches).map(|_| AtomicU64::new(0)).collect(),
                ))
            })
            .collect();
        self.mean = (0..num_bins)
            .map(|_| AtomicU64::new(0.0_f64.to_bits()))
            .collect();
        self.std_dev = (0..num_bins)
            .map(|_| AtomicU64::new(0.0_f64.to_bits()))
            .collect();
        self.rel_error = (0..num_bins)
            .map(|_| AtomicU64::new(0.0_f64.to_bits()))
            .collect();
        self.n_batches.store(n_batches as u32, Ordering::Relaxed);
    }

    /// Initialize through shared reference (for use with Arc<Tally>)
    /// This uses Mutex for initialization, but scoring remains lock-free
    pub fn initialize_batches_shared(&self, n_batches: usize) {
        let num_bins = self.scores.len() * self.num_energy_bins();
        
        // Initialize each score*energy_bin's batch data
        for (i, batch_mutex) in self.batch_data.iter().enumerate() {
            if i < num_bins {
                let new_arc = Arc::new((0..n_batches).map(|_| AtomicU64::new(0)).collect());
                *batch_mutex.lock().unwrap() = new_arc;
            }
        }
        // Initialize statistics
        let num_bins = self.scores.len() * self.num_energy_bins();
        for i in 0..num_bins {
            if let Some(m) = self.mean.get(i) {
                m.store(0.0_f64.to_bits(), Ordering::Relaxed);
            }
            if let Some(s) = self.std_dev.get(i) {
                s.store(0.0_f64.to_bits(), Ordering::Relaxed);
            }
            if let Some(r) = self.rel_error.get(i) {
                r.store(0.0_f64.to_bits(), Ordering::Relaxed);
            }
        }
        self.n_batches.store(n_batches as u32, Ordering::Relaxed);
    }

    /// Clone the tally specification (without batch_data) for creating a new working copy
    pub fn clone_spec(&self) -> Self {
        let num_bins = self.scores.len() * self.num_energy_bins();
        
        Self {
            id: self.id,
            name: self.name.clone(),
            scores: self.scores.clone(),
            filters: self.filters.clone(),
            units: self.units.clone(),
            batch_data: Vec::new(), // Empty - will be initialized separately
            mean: (0..num_bins)
                .map(|_| AtomicU64::new(0.0_f64.to_bits()))
                .collect(),
            std_dev: (0..num_bins)
                .map(|_| AtomicU64::new(0.0_f64.to_bits()))
                .collect(),
            rel_error: (0..num_bins)
                .map(|_| AtomicU64::new(0.0_f64.to_bits()))
                .collect(),
            n_batches: AtomicU32::new(0),
            particles_per_batch: AtomicU32::new(0),
        }
    }

    /// Create a new tally with name and units (for simulation results)
    pub fn with_name_and_units(name: &str, units: &str) -> Self {
        Self {
            id: None,
            name: Some(name.to_string()),
            scores: Vec::new(),
            filters: Vec::new(),
            units: units.to_string(),
            batch_data: Vec::new(),
            mean: Vec::new(),
            std_dev: Vec::new(),
            rel_error: Vec::new(),
            n_batches: AtomicU32::new(0),
            particles_per_batch: AtomicU32::new(0),
        }
    }

    /// Set scores from a mix of integers and score names
    pub fn set_scores_mixed(&mut self, scores: Vec<Score>) {
        self.scores = scores;
        let num_bins = self.scores.len() * self.num_energy_bins();
        
        // Re-initialize storage
        self.batch_data = (0..num_bins)
            .map(|_| Mutex::new(Arc::new(Vec::new())))
            .collect();
        self.mean = (0..num_bins)
            .map(|_| AtomicU64::new(0))
            .collect();
        self.std_dev = (0..num_bins)
            .map(|_| AtomicU64::new(0))
            .collect();
        self.rel_error = (0..num_bins)
            .map(|_| AtomicU64::new(0))
            .collect();
    }

    /// Set the MT number to score
    pub fn set_scores(&mut self, mt_numbers: Vec<i32>) {
        self.scores = mt_numbers.into_iter().map(Score::MT).collect();
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
    pub fn add_batch(
        &self,
        score_index: usize,
        batch_index: usize,
        count: u32,
        _particles_per_batch: u32,
    ) {
        if let Some(batch_mutex) = self.batch_data.get(score_index) {
            let batch_arc = batch_mutex.lock().unwrap();
            if let Some(atomic) = batch_arc.get(batch_index) {
                atomic.store(count as u64, Ordering::Relaxed);
            }
        }
    }

    /// Update statistics from current batch data
    pub fn update_statistics(&self, particles_per_batch: u32) {
        for (i, batch_mutex) in self.batch_data.iter().enumerate() {
            let batch_arc = batch_mutex.lock().unwrap();
            if batch_arc.is_empty() || particles_per_batch == 0 {
                if let Some(m) = self.mean.get(i) {
                    m.store(0.0_f64.to_bits(), Ordering::Relaxed);
                }
                if let Some(s) = self.std_dev.get(i) {
                    s.store(0.0_f64.to_bits(), Ordering::Relaxed);
                }
                if let Some(r) = self.rel_error.get(i) {
                    r.store(0.0_f64.to_bits(), Ordering::Relaxed);
                }
                self.n_batches.store(0, Ordering::Relaxed);
                self.particles_per_batch
                    .store(particles_per_batch, Ordering::Relaxed);
                continue;
            }

            let n = batch_arc.len() as f64;
            let particles_per_batch_f64 = particles_per_batch as f64;

            // Check if this is a flux score (needs special handling)
            let is_flux_score = matches!(self.scores.get(i), Some(Score::Flux));

            // Normalize each batch by particles per batch to get per-source-particle values
            let normalized_data: Vec<f64> = batch_arc
                .iter()
                .map(|atomic_count| {
                    // For flux scores, value is stored as f64 bits
                    // For reaction scores, value is stored as integer count
                    let count = if is_flux_score {
                        f64::from_bits(atomic_count.load(Ordering::Relaxed))
                    } else {
                        atomic_count.load(Ordering::Relaxed) as f64
                    };
                    count / particles_per_batch_f64
                })
                .collect();

            let mean_val = normalized_data.iter().sum::<f64>() / n;
            let variance = normalized_data
                .iter()
                .map(|x| (x - mean_val).powi(2))
                .sum::<f64>()
                / (n - 1.0).max(1.0);
            let std_dev_val = variance.sqrt();
            let rel_error_val = if mean_val > 0.0 {
                std_dev_val / mean_val
            } else {
                0.0
            };

            if let Some(m) = self.mean.get(i) {
                m.store(mean_val.to_bits(), Ordering::Relaxed);
            }
            if let Some(s) = self.std_dev.get(i) {
                s.store(std_dev_val.to_bits(), Ordering::Relaxed);
            }
            if let Some(r) = self.rel_error.get(i) {
                r.store(rel_error_val.to_bits(), Ordering::Relaxed);
            }
            self.n_batches.store(n as u32, Ordering::Relaxed);
            self.particles_per_batch
                .store(particles_per_batch, Ordering::Relaxed);
        }
    }

    /// Get the total count across all batches
    pub fn total_count(&self) -> Vec<u64> {
        self.batch_data
            .iter()
            .map(|batch_mutex| {
                batch_mutex
                    .lock()
                    .unwrap()
                    .iter()
                    .map(|a| a.load(Ordering::Relaxed))
                    .sum()
            })
            .collect()
    }

    /// Get the display name of the tally
    pub fn display_name(&self) -> String {
        self.name
            .clone()
            .unwrap_or_else(|| "Unnamed Tally".to_string())
    }

    /// Helper to get mean as f64 for each score
    pub fn get_mean(&self) -> Vec<f64> {
        self.mean
            .iter()
            .map(|m| f64::from_bits(m.load(Ordering::Relaxed)))
            .collect()
    }

    /// Helper to get std_dev as f64 for each score
    pub fn get_std_dev(&self) -> Vec<f64> {
        self.std_dev
            .iter()
            .map(|s| f64::from_bits(s.load(Ordering::Relaxed)))
            .collect()
    }

    /// Helper to get rel_error as f64 for each score
    pub fn get_rel_error(&self) -> Vec<f64> {
        self.rel_error
            .iter()
            .map(|r| f64::from_bits(r.load(Ordering::Relaxed)))
            .collect()
    }

    /// Get batch data for a specific score index
    /// Returns a vector of normalized scores (per particle) for each batch
    pub fn get_batch_data(&self, score_index: usize) -> Vec<f64> {
        if score_index >= self.batch_data.len() {
            return vec![];
        }

        let batch_arc = self.batch_data[score_index].lock().unwrap().clone();
        let particles_per_batch = self.particles_per_batch.load(Ordering::Relaxed) as f64;

        batch_arc
            .iter()
            .map(|atomic_count| {
                let count = atomic_count.load(Ordering::Relaxed) as f64;
                if particles_per_batch > 0.0 {
                    count / particles_per_batch
                } else {
                    count
                }
            })
            .collect()
    }
}

impl fmt::Display for Tally {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        writeln!(f, "Tally name: {}", self.display_name())?;
        let means = self.get_mean();
        let std_devs = self.get_std_dev();
        let rel_errors = self.get_rel_error();
        let total_counts = self.total_count();

        for (i, score) in self.scores.iter().enumerate() {
            let score_name = match score {
                Score::MT(mt) => format!("MT {}", mt),
                Score::Flux => "Flux".to_string(),
            };
            writeln!(f, "  Score {}:", score_name)?;
            let mean = means.get(i).copied().unwrap_or(0.0);
            let std_dev = std_devs.get(i).copied().unwrap_or(0.0);
            let rel_error = rel_errors.get(i).copied().unwrap_or(0.0);
            let total = total_counts.get(i).copied().unwrap_or(0);

            writeln!(f, "    Mean: {:.6} per particle", mean)?;
            writeln!(f, "    Std Dev: {:.6} per particle", std_dev)?;
            writeln!(
                f,
                "    Rel Error: {:.4} ({:.2}%)",
                rel_error,
                rel_error * 100.0
            )?;
            writeln!(f, "    Total {}: {}", self.units.to_lowercase(), total)?;
            if let Some(batch_mutex) = self.batch_data.get(i) {
                let batch_arc = batch_mutex.lock().unwrap();
                let batch_values: Vec<u64> = batch_arc
                    .iter()
                    .map(|a| a.load(Ordering::Relaxed))
                    .collect();
                writeln!(f, "    Batch data: {:?}", batch_values)?;
            }
        }
        writeln!(f, "    Batches: {}", self.n_batches.load(Ordering::Relaxed))?;
        writeln!(
            f,
            "    Particles per batch: {}",
            self.particles_per_batch.load(Ordering::Relaxed)
        )
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
        let name = spec
            .name
            .clone()
            .unwrap_or_else(|| format!("Tally {}", i + 1));

        let units = "events";

        let mut tally = Tally::with_name_and_units(&name, units);
        tally.id = spec.id;
        tally.scores = spec.scores.clone();
        tally.filters = spec.filters.clone();

        // Validate the tally configuration
        if let Err(err) = tally.validate() {
            panic!("Invalid tally configuration: {}", err);
        }

        tallies.push(tally);
    }

    tallies
}

#[cfg(test)]
mod tests {
    #[test]
    fn test_duplicate_scores_in_tally() {
        use std::sync::atomic::Ordering;
        let mut tally = Tally::new();
        // Add the same score twice
        tally.scores = vec![Score::MT(101), Score::MT(101)];
        tally.initialize_batches(1);
        let cell = dummy_cell(1);
        // Trigger the score event
        assert!(tally.score_event(101, &cell, Some(42), 1e6, 0));
        // Should only increment once
        let batch_arc = tally.batch_data[0].lock().unwrap();
        assert_eq!(
            batch_arc[0].load(Ordering::Relaxed),
            1,
            "Should increment only once for duplicate scores"
        );
        drop(batch_arc);
        tally.update_statistics(1);
        let means = tally.get_mean();
        assert_eq!(means.len(), 2, "Should have two means for two scores");
        assert_eq!(means[0], 1.0, "First mean should be 1.0");
        assert_eq!(means[1], 1.0, "Second mean should be 1.0");
        assert_eq!(means[0], means[1], "Means should be equal");
    }
    use super::*;
    use crate::cell::Cell;

    fn dummy_cell(cell_id: i32) -> Cell {
        // Minimal cell with just an ID for filter testing
        use crate::region::{HalfspaceType, Region};
        use crate::surface::{BoundaryType, Surface, SurfaceKind};
        use std::sync::Arc;
        let dummy_surface = Arc::new(Surface {
            surface_id: Some(1),
            kind: SurfaceKind::Plane {
                a: 1.0,
                b: 0.0,
                c: 0.0,
                d: 0.0,
            },
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
        tally.set_scores(vec![101]); // Absorption
        tally.initialize_batches(1); // Start batch
        let cell = dummy_cell(1);
        assert!(tally.score_event(101, &cell, Some(42), 1e6, 0));
        let batch_arc = tally.batch_data[0].lock().unwrap();
        assert_eq!(
            batch_arc[0].load(Ordering::Relaxed),
            1,
            "Should increment for direct MT match"
        );
    }

    #[test]
    fn test_score_event_mt4_inelastic_constituent() {
        use std::sync::atomic::Ordering;
        let mut tally = Tally::new();
        tally.set_scores(vec![4]); // MT 4 (inelastic)
        tally.initialize_batches(1); // Start batch
        let cell = dummy_cell(1);
        assert!(tally.score_event(53, &cell, Some(42), 1e6, 0)); // MT 53 is inelastic constituent
        let batch_arc = tally.batch_data[0].lock().unwrap();
        assert_eq!(
            batch_arc[0].load(Ordering::Relaxed),
            1,
            "Should increment for MT 4 when constituent MT occurs"
        );
    }

    #[test]
    fn test_score_event_inelastic_constituent_only() {
        use std::sync::atomic::Ordering;
        let mut tally = Tally::new();
        tally.set_scores(vec![53]); // MT 53
        tally.initialize_batches(1); // Start batch
        let cell = dummy_cell(1);
        assert!(tally.score_event(53, &cell, Some(42), 1e6, 0));
        {
            let batch_arc = tally.batch_data[0].lock().unwrap();
            assert_eq!(
                batch_arc[0].load(Ordering::Relaxed),
                1,
                "Should increment for direct constituent MT match"
            );
        }
        assert!(!tally.score_event(4, &cell, Some(42), 1e6, 0));
        let batch_arc = tally.batch_data[0].lock().unwrap();
        assert_eq!(
            batch_arc[0].load(Ordering::Relaxed),
            1,
            "Should not increment for MT 4 if tally is for constituent"
        );
    }

    #[test]
    fn test_score_event_with_cell_filter() {
        use std::sync::atomic::Ordering;
        let mut tally = Tally::new();
        tally.set_scores(vec![101]);
        tally.initialize_batches(1);
        // Add a cell filter that matches cell_id 1
        tally
            .filters
            .push(Filter::Cell(crate::filters::CellFilter { cell_id: 1 }));
        let cell = dummy_cell(1);
        assert!(tally.score_event(101, &cell, Some(42), 1e6, 0));
        {
            let batch_arc = tally.batch_data[0].lock().unwrap();
            assert_eq!(
                batch_arc[0].load(Ordering::Relaxed),
                1,
                "Should increment when cell filter matches"
            );
        }
        let cell2 = dummy_cell(2);
        assert!(!tally.score_event(101, &cell2, Some(42), 1e6, 0));
        let batch_arc = tally.batch_data[0].lock().unwrap();
        assert_eq!(
            batch_arc[0].load(Ordering::Relaxed),
            1,
            "Should not increment when cell filter does not match"
        );
    }
}
