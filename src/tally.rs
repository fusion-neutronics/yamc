use std::fmt;
use crate::filters::{CellFilter, MaterialFilter};

/// Unified filter enum for tallies
#[derive(Debug, Clone, PartialEq)]
pub enum Filter {
    Cell(CellFilter),
    Material(MaterialFilter),
}

/// Unified tally structure serving as both input specification and results container
#[derive(Debug, Clone)]
pub struct Tally {
    // Input specification fields
    pub id: Option<u32>,           // Optional tally ID
    pub name: Option<String>,      // Optional tally name  
    pub score: i32,               // MT number to score (e.g., 101 for absorption)
    pub filters: Vec<Filter>, // Filters for this tally
    
    // Results fields (populated during simulation)
    pub units: String,          // Units (e.g., "particles", "collisions", "reactions")  
    pub batch_data: Vec<u32>,   // Raw counts per batch
    // Statistics (calculated from batch_data)
    pub mean: f64,              // Mean per source particle
    pub std_dev: f64,           // Standard deviation per source particle
    pub rel_error: f64,         // Relative error (coefficient of variation)
    pub n_batches: u32,         // Number of batches
    pub particles_per_batch: u32, // Particles per batch for normalization
}

impl Tally {
    /// Create a new tally specification
    pub fn new() -> Self {
        Self {
            id: None,
            name: None,
            score: 0,
            filters: Vec::new(),
            units: String::new(),
            batch_data: Vec::new(),
            mean: 0.0,
            std_dev: 0.0,
            rel_error: 0.0,
            n_batches: 0,
            particles_per_batch: 0,
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
            batch_data: Vec::new(),
            mean: 0.0,
            std_dev: 0.0,
            rel_error: 0.0,
            n_batches: 0,
            particles_per_batch: 0,
        }
    }
    
    /// Set the MT number to score
    pub fn set_score(&mut self, mt_number: i32) {
        self.score = mt_number;
    }

    /// Add a batch result to the tally
    pub fn add_batch(&mut self, count: u32, particles_per_batch: u32) {
        self.batch_data.push(count);
        self.update_statistics(particles_per_batch);
    }

    /// Update statistics from current batch data
    fn update_statistics(&mut self, particles_per_batch: u32) {
        if self.batch_data.is_empty() || particles_per_batch == 0 {
            self.mean = 0.0;
            self.std_dev = 0.0;
            self.rel_error = 0.0;
            self.n_batches = 0;
            self.particles_per_batch = particles_per_batch;
            return;
        }

        let n = self.batch_data.len() as f64;
        let particles_per_batch_f64 = particles_per_batch as f64;
        
        // Normalize each batch by particles per batch to get per-source-particle values
        let normalized_data: Vec<f64> = self.batch_data
            .iter()
            .map(|&count| count as f64 / particles_per_batch_f64)
            .collect();
        
        self.mean = normalized_data.iter().sum::<f64>() / n;
        let variance = normalized_data
            .iter()
            .map(|x| (x - self.mean).powi(2))
            .sum::<f64>() / (n - 1.0).max(1.0);
        self.std_dev = variance.sqrt();
        self.rel_error = if self.mean > 0.0 { self.std_dev / self.mean } else { 0.0 };
        self.n_batches = n as u32;
        self.particles_per_batch = particles_per_batch;
    }

    /// Get the total count across all batches
    pub fn total_count(&self) -> u32 {
        self.batch_data.iter().sum()
    }
    
    /// Get the display name of the tally
    pub fn display_name(&self) -> String {
        self.name.clone().unwrap_or_else(|| "Unnamed Tally".to_string())
    }
}

impl fmt::Display for Tally {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        writeln!(f, "Tally name: {}", self.display_name())?;
        writeln!(f, "  Mean: {:.6} per particle", self.mean)?;
        writeln!(f, "    Std Dev: {:.6} per particle", self.std_dev)?;
        writeln!(f, "    Rel Error: {:.4} ({:.2}%)", self.rel_error, self.rel_error * 100.0)?;
        writeln!(f, "    Batches: {}", self.n_batches)?;
        writeln!(f, "    Particles per batch: {}", self.particles_per_batch)?;
        writeln!(f, "  Total {}: {}", self.units.to_lowercase(), self.total_count())?;
        write!(f, "  Batch data: {:?}", self.batch_data)
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
        
        tallies.push(tally);
    }
    
    tallies
}