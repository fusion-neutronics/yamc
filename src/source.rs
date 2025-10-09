#[derive(Debug, Clone)]
pub struct IndependentSource {
    pub space: [f64; 3],
    pub angle: [f64; 3],
    pub energy: f64,
}

impl IndependentSource {
    pub fn sample(&self) -> crate::particle::Particle {
        crate::particle::Particle::new(self.space, self.angle, self.energy)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_source_construction() {
        let s = IndependentSource {
            space: [1.0, 2.0, 3.0],
            angle: [0.0, 0.0, 1.0],
            energy: 2e6,
        };
        let p = s.sample();
        assert_eq!(p.position, [1.0, 2.0, 3.0]);
        assert_eq!(p.direction, [0.0, 0.0, 1.0]);
        assert_eq!(p.energy, 2e6);
        assert!(p.alive);
    }
}
