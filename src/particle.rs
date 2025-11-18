#[derive(Debug, Clone)]
pub struct Particle {
    pub position: [f64; 3],
    pub direction: [f64; 3],
    pub energy: f64,
    pub alive: bool,
    pub id: usize, // Particle ID for tracking
    pub current_cell_index: Option<usize>, // Cached cell index for performance
}

impl Particle {
    pub fn new(position: [f64; 3], direction: [f64; 3], energy: f64) -> Self {
        Self {
            position,
            direction,
            energy,
            alive: true,
            id: 0, // Default ID
            current_cell_index: None, // Will be set on first transport step
        }
    }

    /// Move the particle along its current direction by the specified distance
    pub fn move_by(&mut self, distance: f64) {
        for i in 0..3 {
            self.position[i] += self.direction[i] * distance;
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

    #[test]
    fn test_particle_move_by() {
        let mut p = Particle::new([0.0, 0.0, 0.0], [1.0, 0.0, 0.0], 1e6);

        // Move 2 units in x direction
        p.move_by(2.0);
        assert_eq!(p.position, [2.0, 0.0, 0.0]);

        // Move another 1.5 units
        p.move_by(1.5);
        assert_eq!(p.position, [3.5, 0.0, 0.0]);

        // Test with diagonal direction (normalized)
        let sqrt_2_inv = 1.0 / 2.0_f64.sqrt();
        let mut p2 = Particle::new([0.0, 0.0, 0.0], [sqrt_2_inv, sqrt_2_inv, 0.0], 1e6);
        p2.move_by(2.0_f64.sqrt()); // Move sqrt(2) units to get 1 unit in both x and y

        // Allow for floating point precision
        assert!((p2.position[0] - 1.0).abs() < 1e-10);
        assert!((p2.position[1] - 1.0).abs() < 1e-10);
        assert_eq!(p2.position[2], 0.0);
    }
}
