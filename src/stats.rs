use rand::Rng;

/// Trait for statistical distributions that can be sampled for 3D direction vectors
pub trait AngularDistribution: std::fmt::Debug + Send + Sync {
    fn sample(&self) -> [f64; 3];
    fn clone_box(&self) -> Box<dyn AngularDistribution>;
    fn type_name(&self) -> &'static str;
    fn as_any(&self) -> &dyn std::any::Any;
}

impl Clone for Box<dyn AngularDistribution> {
    fn clone(&self) -> Self {
        self.clone_box()
    }
}

/// Isotropic distribution (uniform on sphere)
#[derive(Debug, Clone)]
pub struct Isotropic;

impl AngularDistribution for Isotropic {
    fn sample(&self) -> [f64; 3] {
        // Sample isotropic direction using OpenMC's method
        let mut rng = rand::thread_rng();
        // Generate uniform random numbers
        let xi1: f64 = rng.gen();
        let xi2: f64 = rng.gen();
        
        // Convert to spherical coordinates
        let mu = 2.0 * xi1 - 1.0; // cosine of polar angle
        let phi = 2.0 * std::f64::consts::PI * xi2; // azimuthal angle
        
        // Convert to Cartesian coordinates
        let sin_theta = (1.0 - mu * mu).sqrt();
        [
            sin_theta * phi.cos(),
            sin_theta * phi.sin(),
            mu,
        ]
    }
    
    fn clone_box(&self) -> Box<dyn AngularDistribution> {
        Box::new(self.clone())
    }

    fn type_name(&self) -> &'static str {
        "Isotropic"
    }

    fn as_any(&self) -> &dyn std::any::Any {
        self
    }
}

/// Monodirectional distribution (fixed direction)
#[derive(Debug, Clone)]
pub struct Monodirectional {
    pub reference_uvw: [f64; 3],
}

impl Monodirectional {
    pub fn new(u: f64, v: f64, w: f64) -> Self {
        // Normalize the direction vector
        let mag = (u * u + v * v + w * w).sqrt();
        if mag == 0.0 {
            panic!("Direction vector cannot be zero");
        }
        Self {
            reference_uvw: [u / mag, v / mag, w / mag],
        }
    }
}

impl AngularDistribution for Monodirectional {
    fn sample(&self) -> [f64; 3] {
        self.reference_uvw
    }
    
    fn clone_box(&self) -> Box<dyn AngularDistribution> {
        Box::new(self.clone())
    }

    fn type_name(&self) -> &'static str {
        "Monodirectional"
    }

    fn as_any(&self) -> &dyn std::any::Any {
        self
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_monodirectional_distribution() {
        let mono = Monodirectional::new(0.0, 0.0, 1.0);
        let sample = mono.sample();
        assert_eq!(sample, [0.0, 0.0, 1.0]);
    }

    #[test]
    fn test_monodirectional_different_directions() {
        let test_directions = [
            [1.0, 0.0, 0.0],
            [0.0, 1.0, 0.0],
            [0.0, 0.0, 1.0],
            [-1.0, 0.0, 0.0],
            [0.5773502691896257, 0.5773502691896257, 0.5773502691896257], // normalized (1,1,1)
        ];

        for &direction in &test_directions {
            let mono = Monodirectional::new(direction[0], direction[1], direction[2]);
            let sample = mono.sample();
            assert_eq!(sample, direction);
        }
    }

    #[test]
    fn test_monodirectional_consistency() {
        let mono = Monodirectional::new(1.0, 0.0, 0.0);
        
        // Should return the same direction consistently
        for _ in 0..100 {
            let sample = mono.sample();
            assert_eq!(sample, [1.0, 0.0, 0.0]);
        }
    }

    #[test]
    fn test_isotropic_distribution() {
        let iso = Isotropic;
        let sample = iso.sample();
        
        // Check that the vector is normalized
        let mag = (sample[0] * sample[0] + sample[1] * sample[1] + sample[2] * sample[2]).sqrt();
        assert!((mag - 1.0).abs() < 1e-10);
    }

    #[test]
    fn test_isotropic_randomness() {
        let iso = Isotropic;
        
        // Sample many directions and check they're all normalized and vary
        let mut samples = Vec::new();
        for _ in 0..1000 {
            let sample = iso.sample();
            
            // Check normalization
            let mag = (sample[0] * sample[0] + sample[1] * sample[1] + sample[2] * sample[2]).sqrt();
            assert!((mag - 1.0).abs() < 1e-10);
            
            samples.push(sample);
        }
        
        // Check that we have variation (not all the same)
        let first_sample = samples[0];
        let all_same = samples.iter().all(|&s| s == first_sample);
        assert!(!all_same, "Isotropic samples should vary");
    }

    #[test]
    fn test_trait_object_safety() {
        // Test that we can use trait objects (Box<dyn AngularDistribution>)
        let iso: Box<dyn AngularDistribution> = Box::new(Isotropic);
        let mono: Box<dyn AngularDistribution> = Box::new(Monodirectional::new(0.0, 0.0, 1.0));
        
        // Should be able to call sample through trait objects
        let iso_sample = iso.sample();
        let mono_sample = mono.sample();
        
        // Isotropic should be normalized
        let mag = (iso_sample[0] * iso_sample[0] + iso_sample[1] * iso_sample[1] + iso_sample[2] * iso_sample[2]).sqrt();
        assert!((mag - 1.0).abs() < 1e-10);
        
        // Monodirectional should match expected
        assert_eq!(mono_sample, [0.0, 0.0, 1.0]);
    }

    #[test]
    fn test_send_sync_bounds() {
        // Test that our distributions implement Send + Sync for threading
        fn assert_send<T: Send>() {}
        fn assert_sync<T: Sync>() {}
        
        assert_send::<Isotropic>();
        assert_sync::<Isotropic>();
        assert_send::<Monodirectional>();
        assert_sync::<Monodirectional>();
        
        // Test with trait objects too
        assert_send::<Box<dyn AngularDistribution>>();
        assert_sync::<Box<dyn AngularDistribution>>();
    }
}