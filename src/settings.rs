use crate::source::IndependentSource;

#[derive(Debug, Clone)]
pub struct Settings {
    pub particles: usize,
    pub batches: usize,
    pub source: IndependentSource,
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::source::IndependentSource;

    #[test]
    fn test_settings_construction() {
        use crate::stats::Monodirectional;
        let src = IndependentSource {
            space: [0.0, 0.0, 0.0],
            angle: Box::new(Monodirectional::new(0.0, 1.0, 0.0)),
            energy: 1e5,
        };
        let settings = Settings {
            particles: 100,
            batches: 10,
            source: src.clone(),
        };
        assert_eq!(settings.particles, 100);
        assert_eq!(settings.batches, 10);
        assert_eq!(settings.source.space, [0.0, 0.0, 0.0]);
    }

    #[test]
    fn test_settings_source_assignment() {
        use crate::stats::Monodirectional;
        let src = IndependentSource {
            space: [1.0, 1.0, 1.0],
            angle: Box::new(Monodirectional::new(0.0, 0.0, 1.0)),
            energy: 1e6,
        };
        let mut settings = Settings {
            particles: 50,
            batches: 5,
            source: src.clone(),
        };
        settings.source = src.clone();
        assert_eq!(settings.source.space, [1.0, 1.0, 1.0]);
    }
}
