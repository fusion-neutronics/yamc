#[derive(Debug, Clone)]
pub struct Source {
    pub position: [f64; 3],
    pub direction: [f64; 3],
    pub energy: f64,
}

impl Source {
    pub fn sample(&self) -> crate::particle::Particle {
        crate::particle::Particle::new(self.position, self.direction, self.energy)
    }
}


#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_source_construction() {
        let s = Source {
            position: [1.0, 2.0, 3.0],
            direction: [0.0, 0.0, 1.0],
            energy: 2e6,
        };
    let p = s.sample();
    assert_eq!(p.position, [1.0, 2.0, 3.0]);
    assert_eq!(p.direction, [0.0, 0.0, 1.0]);
    assert_eq!(p.energy, 2e6);
    assert!(p.alive);
    }
}