#[derive(Debug, Clone)]
pub struct Particle {
    pub position: [f64; 3],
    pub direction: [f64; 3],
    pub energy: f64,
    pub alive: bool,
}

impl Particle {
    pub fn new(position: [f64; 3], direction: [f64; 3], energy: f64) -> Self {
        Self {
            position,
            direction,
            energy,
            alive: true,
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_particle_construction() {
        let p = Particle::new([0.0, 1.0, 2.0], [1.0, 0.0, 0.0], 1e6);
        assert_eq!(p.position, [0.0, 1.0, 2.0]);
        assert_eq!(p.direction, [1.0, 0.0, 0.0]);
        assert_eq!(p.energy, 1e6);
        assert!(p.alive);
    }
}
