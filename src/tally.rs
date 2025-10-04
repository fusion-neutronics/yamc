use std::fmt;

/// Generic tally for counting events (leakage, collisions, absorptions, etc.)
#[derive(Debug, Clone)]
pub struct CountTally {
    pub name: String,           // Name of the tally (e.g., "Leakage", "Absorption", "Flux")
    pub units: String,          // Units (e.g., "particles", "collisions", "reactions")  
    pub batch_data: Vec<u32>,   // Raw counts per batch
    // Statistics (calculated from batch_data)
    pub mean: f64,              // Mean per source particle
    pub std_dev: f64,           // Standard deviation per source particle
    pub rel_error: f64,         // Relative error (coefficient of variation)
    pub n_batches: u32,         // Number of batches
    pub particles_per_batch: u32, // Particles per batch for normalization
}

impl CountTally {
    /// Create a new count tally with a name and units
    pub fn new(name: &str, units: &str) -> Self {
        Self {
            name: name.to_string(),
            units: units.to_string(),
            batch_data: Vec::new(),
            mean: 0.0,
            std_dev: 0.0,
            rel_error: 0.0,
            n_batches: 0,
            particles_per_batch: 0,
        }
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
}

impl fmt::Display for CountTally {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        writeln!(f, "Tally name: {}", self.name)?;
        writeln!(f, "  Mean: {:.6} per particle", self.mean)?;
        writeln!(f, "    Std Dev: {:.6} per particle", self.std_dev)?;
        writeln!(f, "    Rel Error: {:.4} ({:.2}%)", self.rel_error, self.rel_error * 100.0)?;
        writeln!(f, "    Batches: {}", self.n_batches)?;
        writeln!(f, "    Particles per batch: {}", self.particles_per_batch)?;
        writeln!(f, "  Total {}: {}", self.units.to_lowercase(), self.total_count())?;
        write!(f, "  Batch data: {:?}", self.batch_data)
    }
}

/// Convenience type alias for leakage tallies
pub type LeakageTally = CountTally;

/// User-configurable tally specification
#[derive(Debug, Clone)]
pub struct Tally {
    pub id: Option<u32>,           // Optional tally ID
    pub name: Option<String>,      // Optional tally name
    pub score: Vec<i32>,          // MT numbers to score (e.g., [101] for absorption)
    pub filters: Vec<String>,     // Future: energy, material, cell filters
}

impl Tally {
    /// Create a new tally specification
    pub fn new() -> Self {
        Self {
            id: None,
            name: None,
            score: Vec::new(),
            filters: Vec::new(),
        }
    }
    
    /// Set the MT numbers to score
    pub fn set_score(&mut self, mt_numbers: Vec<i32>) {
        self.score = mt_numbers;
    }
    
    /// Add an MT number to score
    pub fn add_score(&mut self, mt_number: i32) {
        self.score.push(mt_number);
    }
}

/// Initialize user tallies from tally specifications
pub fn create_tallies_from_specs(tally_specs: &[Tally]) -> Vec<CountTally> {
    let mut tallies = Vec::new();
    
    // Always add leakage tally as the first tally
    tallies.push(CountTally::new("Leakage", "particles"));
    
    // Add user-specified tallies
    for (i, spec) in tally_specs.iter().enumerate() {
        let name = spec.name.clone().unwrap_or_else(|| {
            format!("Tally {}", i + 1)
        });
        
        let units = "events";
        
        tallies.push(CountTally::new(&name, units));
    }
    
    tallies
}