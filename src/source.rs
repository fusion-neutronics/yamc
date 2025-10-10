#[derive(Debug, Clone)]
pub struct IndependentSource {
    pub space: [f64; 3],
    pub angle: Box<dyn crate::stats::AngularDistribution>,
    pub energy: f64,
}

impl IndependentSource {
    pub fn sample(&self) -> crate::particle::Particle {
        let sampled_angle = self.angle.sample();
        crate::particle::Particle::new(self.space, sampled_angle, self.energy)
    }
    
    pub fn new() -> Self {
        use crate::stats::Isotropic;
        Self {
            space: [0.0, 0.0, 0.0],
            angle: Box::new(Isotropic),
            energy: 14.06e6,
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_source_construction() {
        use crate::stats::{Monodirectional, Isotropic};
        
        let mut s = IndependentSource::new();
        s.space = [1.0, 2.0, 3.0];
        s.angle = Box::new(Monodirectional::new(0.0, 0.0, 1.0));
        s.energy = 2e6;
        
        let p = s.sample();
        assert_eq!(p.position, [1.0, 2.0, 3.0]);
        assert_eq!(p.direction, [0.0, 0.0, 1.0]);
        assert_eq!(p.energy, 2e6);
        assert!(p.alive);
    }
    
    #[test]
    fn test_default_source() {
        let s = IndependentSource::new();
        let p = s.sample();
        
        // Check that we get valid values
        assert_eq!(p.position, [0.0, 0.0, 0.0]);
        assert_eq!(p.energy, 14.06e6);
        assert!(p.alive);
        
        // Check that direction is normalized (isotropic sampling)
        let mag = (p.direction[0] * p.direction[0] + p.direction[1] * p.direction[1] + p.direction[2] * p.direction[2]).sqrt();
        assert!((mag - 1.0).abs() < 1e-10);
    }

    #[test]
    fn test_source_space_modification() {
        let mut s = IndependentSource::new();
        
        // Test different space values
        let test_positions = [
            [1.0, 2.0, 3.0],
            [-5.0, 0.0, 10.0],
            [0.5, -0.5, 0.0],
        ];
        
        for &position in &test_positions {
            s.space = position;
            let p = s.sample();
            assert_eq!(p.position, position);
        }
    }

    #[test]
    fn test_source_energy_modification() {
        let mut s = IndependentSource::new();
        
        // Test different energy values
        let test_energies = [1e6, 2.5e6, 14.1e6, 20e6];
        
        for &energy in &test_energies {
            s.energy = energy;
            let p = s.sample();
            assert_eq!(p.energy, energy);
        }
    }

    #[test]
    fn test_source_angle_switching() {
        use crate::stats::{Monodirectional, Isotropic};
        
        let mut s = IndependentSource::new();
        
        // Start with monodirectional
        s.angle = Box::new(Monodirectional::new(1.0, 0.0, 0.0));
        let p1 = s.sample();
        assert_eq!(p1.direction, [1.0, 0.0, 0.0]);
        
        // Switch to different monodirectional
        s.angle = Box::new(Monodirectional::new(0.0, 1.0, 0.0));
        let p2 = s.sample();
        assert_eq!(p2.direction, [0.0, 1.0, 0.0]);
        
        // Switch to isotropic
        s.angle = Box::new(Isotropic);
        let p3 = s.sample();
        let p4 = s.sample();
        
        // Both should be normalized
        let p3_mag = (p3.direction[0] * p3.direction[0] + p3.direction[1] * p3.direction[1] + p3.direction[2] * p3.direction[2]).sqrt();
        let p4_mag = (p4.direction[0] * p4.direction[0] + p4.direction[1] * p4.direction[1] + p4.direction[2] * p4.direction[2]).sqrt();
        assert!((p3_mag - 1.0).abs() < 1e-10);
        assert!((p4_mag - 1.0).abs() < 1e-10);
        
        // Very unlikely to be identical with isotropic sampling
        assert_ne!(p3.direction, p4.direction);
    }

    #[test]
    fn test_source_consistency() {
        use crate::stats::Monodirectional;
        
        let mut s = IndependentSource::new();
        s.space = [1.0, 2.0, 3.0];
        s.energy = 5e6;
        s.angle = Box::new(Monodirectional::new(0.0, 0.0, 1.0));
        
        // Multiple samples should be consistent for monodirectional
        for _ in 0..10 {
            let p = s.sample();
            assert_eq!(p.position, [1.0, 2.0, 3.0]);
            assert_eq!(p.energy, 5e6);
            assert_eq!(p.direction, [0.0, 0.0, 1.0]);
            assert!(p.alive);
        }
    }

    #[test]
    fn test_source_isotropic_variation() {
        let s = IndependentSource::new(); // Uses isotropic by default
        
        // Sample many particles and check for variation
        let mut directions = Vec::new();
        for _ in 0..100 {
            let p = s.sample();
            
            // All should have same position and energy (deterministic)
            assert_eq!(p.position, [0.0, 0.0, 0.0]);
            assert_eq!(p.energy, 14.06e6);
            assert!(p.alive);
            
            // Direction should be normalized
            let mag = (p.direction[0] * p.direction[0] + p.direction[1] * p.direction[1] + p.direction[2] * p.direction[2]).sqrt();
            assert!((mag - 1.0).abs() < 1e-10);
            
            directions.push(p.direction);
        }
        
        // Should have variation in directions (isotropic)
        let first_direction = directions[0];
        let all_same = directions.iter().all(|&d| d == first_direction);
        assert!(!all_same, "Isotropic source should produce varying directions");
    }
}
