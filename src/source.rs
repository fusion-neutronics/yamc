#[derive(Debug, Clone)]
pub struct Source {
    pub position: [f64; 3],
    pub direction: [f64; 3],
    pub energy: f64,
}

impl Source {
    pub fn sample(&self) -> ([f64; 3], [f64; 3], f64) {
        // For now, just return the fixed values
        (self.position, self.direction, self.energy)
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
        let (pos, dir, en) = s.sample();
        assert_eq!(pos, [1.0, 2.0, 3.0]);
        assert_eq!(dir, [0.0, 0.0, 1.0]);
        assert_eq!(en, 2e6);
    }
}