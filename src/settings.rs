use crate::source::Source;

#[derive(Debug, Clone)]
pub struct Settings {
    pub particles: usize,
    pub batches: usize,
    pub source: Source,
}


#[cfg(test)]
mod tests {
    use super::*;
    use crate::source::Source;

    #[test]
    fn test_settings_construction() {
        let src = Source {
            position: [0.0, 0.0, 0.0],
            direction: [0.0, 1.0, 0.0],
            energy: 1e5,
        };
        let settings = Settings {
            particles: 100,
            batches: 10,
            source: src.clone(),
        };
        assert_eq!(settings.particles, 100);
        assert_eq!(settings.batches, 10);
        assert_eq!(settings.source.position, [0.0, 0.0, 0.0]);
    }

    #[test]
    fn test_settings_source_assignment() {
        let src = Source {
            position: [1.0, 1.0, 1.0],
            direction: [0.0, 0.0, 1.0],
            energy: 1e6,
        };
        let mut settings = Settings {
            particles: 50,
            batches: 5,
            source: src.clone(),
        };
        settings.source = src.clone();
        assert_eq!(settings.source.position, [1.0, 1.0, 1.0]);
    }
}