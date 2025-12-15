use crate::source::IndependentSource;

#[derive(Debug, Clone)]
pub struct Settings {
    pub particles: usize,
    pub batches: usize,
    pub source: IndependentSource,
    /// Random number seed for reproducibility. If None, a random seed is used.
    /// Default: Some(1) for reproducible results
    pub seed: Option<u64>,
}

impl Settings {
    /// Create a new Settings with default seed (1)
    pub fn new(particles: usize, batches: usize, source: IndependentSource) -> Self {
        Self {
            particles,
            batches,
            source,
            seed: Some(1),
        }
    }

    /// Get the seed to use, defaulting to 1 if None
    pub fn get_seed(&self) -> u64 {
        self.seed.unwrap_or(1)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::source::IndependentSource;

    #[test]
    fn test_settings_construction() {
        let src = IndependentSource {
            space: [0.0, 0.0, 0.0],
            angle: crate::distribution_multi::AngularDistribution::new_monodirectional(0.0, 1.0, 0.0),
            energy: 1e5,
        };
        let settings = Settings {
            particles: 100,
            batches: 10,
            source: src.clone(),
            seed: Some(1),
        };
        assert_eq!(settings.particles, 100);
        assert_eq!(settings.batches, 10);
        assert_eq!(settings.source.space, [0.0, 0.0, 0.0]);
        assert_eq!(settings.seed, Some(1));
        assert_eq!(settings.get_seed(), 1);
    }

    #[test]
    fn test_settings_source_assignment() {
        let src = IndependentSource {
            space: [1.0, 1.0, 1.0],
            angle: crate::distribution_multi::AngularDistribution::new_monodirectional(0.0, 0.0, 1.0),
            energy: 1e6,
        };
        let mut settings = Settings {
            particles: 50,
            batches: 5,
            source: src.clone(),
            seed: Some(42),
        };
        settings.source = src.clone();
        assert_eq!(settings.source.space, [1.0, 1.0, 1.0]);
        assert_eq!(settings.get_seed(), 42);
    }
}
