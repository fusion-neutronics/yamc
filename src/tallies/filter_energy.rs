/// Energy filter for tallies - filters events based on particle energy
/// Energy bins are defined by bin edges [E0, E1, E2, ..., En]
/// This creates bins: [E0, E1), [E1, E2), ..., [En-1, En]
#[derive(Debug, Clone, PartialEq)]
pub struct EnergyFilter {
    /// Energy bin boundaries in eV, must be in ascending order
    /// For n boundaries, creates n-1 bins
    pub bins: Vec<f64>,
}

impl EnergyFilter {
    /// Create a new EnergyFilter with the given bin boundaries
    ///
    /// # Arguments
    /// * `bins` - Energy bin boundaries in eV, must be in ascending order
    ///
    /// # Returns
    /// A new EnergyFilter
    ///
    /// # Panics
    /// Panics if bins are not in ascending order or if there are fewer than 2 bins
    pub fn new(bins: Vec<f64>) -> Self {
        if bins.len() < 2 {
            panic!("EnergyFilter requires at least 2 bin boundaries (to create at least 1 bin)");
        }
        
        // Verify bins are in ascending order
        for i in 1..bins.len() {
            if bins[i] <= bins[i - 1] {
                panic!("Energy bins must be in strictly ascending order");
            }
        }
        
        Self { bins }
    }
    
    /// Get the bin index for a given energy
    ///
    /// # Arguments
    /// * `energy` - The particle energy in eV
    ///
    /// # Returns
    /// `Some(bin_index)` if the energy falls within the filter range, `None` otherwise
    pub fn get_bin(&self, energy: f64) -> Option<usize> {
        // Energy must be >= first bin and < last bin
        if energy < self.bins[0] || energy >= *self.bins.last().unwrap() {
            return None;
        }
        
        // Binary search for the bin
        // We want to find i such that bins[i] <= energy < bins[i+1]
        let result = self.bins.binary_search_by(|&bin| {
            if bin <= energy {
                std::cmp::Ordering::Less
            } else {
                std::cmp::Ordering::Greater
            }
        });
        
        match result {
            Ok(i) => Some(i),
            Err(i) => {
                if i > 0 && i < self.bins.len() {
                    Some(i - 1)
                } else {
                    None
                }
            }
        }
    }
    
    /// Check if this filter matches a given energy (returns true if energy is in range)
    ///
    /// # Arguments
    /// * `energy` - The particle energy in eV to check
    ///
    /// # Returns
    /// `true` if the energy falls within any bin of this filter
    pub fn matches(&self, energy: f64) -> bool {
        self.get_bin(energy).is_some()
    }
    
    /// Get the number of energy bins
    pub fn num_bins(&self) -> usize {
        self.bins.len() - 1
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_energy_filter_creation() {
        let bins = vec![0.0, 1e6, 10e6, 20e6];
        let filter = EnergyFilter::new(bins.clone());
        assert_eq!(filter.bins, bins);
        assert_eq!(filter.num_bins(), 3);
    }

    #[test]
    #[should_panic(expected = "EnergyFilter requires at least 2 bin boundaries")]
    fn test_energy_filter_too_few_bins() {
        EnergyFilter::new(vec![1e6]);
    }

    #[test]
    #[should_panic(expected = "Energy bins must be in strictly ascending order")]
    fn test_energy_filter_non_ascending() {
        EnergyFilter::new(vec![1e6, 10e6, 5e6]);
    }

    #[test]
    #[should_panic(expected = "Energy bins must be in strictly ascending order")]
    fn test_energy_filter_duplicate_bins() {
        EnergyFilter::new(vec![1e6, 10e6, 10e6, 20e6]);
    }

    #[test]
    fn test_energy_filter_get_bin() {
        let filter = EnergyFilter::new(vec![0.0, 1e6, 10e6, 20e6]);
        
        assert_eq!(filter.get_bin(0.0), Some(0));
        assert_eq!(filter.get_bin(5e5), Some(0));
        assert_eq!(filter.get_bin(9.99e5), Some(0));
        
        assert_eq!(filter.get_bin(1e6), Some(1));
        assert_eq!(filter.get_bin(5e6), Some(1));
        assert_eq!(filter.get_bin(9.99e6), Some(1));
        
        assert_eq!(filter.get_bin(10e6), Some(2));
        assert_eq!(filter.get_bin(15e6), Some(2));
        assert_eq!(filter.get_bin(19.99e6), Some(2));
        
        assert_eq!(filter.get_bin(-1.0), None);
        assert_eq!(filter.get_bin(20e6), None);
        assert_eq!(filter.get_bin(25e6), None);
    }

    #[test]
    fn test_energy_filter_matches() {
        let filter = EnergyFilter::new(vec![1e6, 10e6, 20e6]);
        
        assert!(!filter.matches(0.5e6));
        assert!(filter.matches(1e6));
        assert!(filter.matches(5e6));
        assert!(filter.matches(15e6));
        assert!(!filter.matches(20e6));
        assert!(!filter.matches(25e6));
    }

    #[test]
    fn test_energy_filter_typical_neutron_bins() {
        let filter = EnergyFilter::new(vec![
            0.0,
            0.625,
            100e3,
            20e6,
        ]);
        
        assert_eq!(filter.num_bins(), 3);
        assert_eq!(filter.get_bin(0.025), Some(0));
        assert_eq!(filter.get_bin(1000.0), Some(1));
        assert_eq!(filter.get_bin(14.1e6), Some(2));
    }

    #[test]
    fn test_energy_filter_equality() {
        let filter1 = EnergyFilter::new(vec![0.0, 1e6, 10e6]);
        let filter2 = EnergyFilter::new(vec![0.0, 1e6, 10e6]);
        let filter3 = EnergyFilter::new(vec![0.0, 2e6, 10e6]);
        
        assert_eq!(filter1, filter2);
        assert_ne!(filter1, filter3);
    }

    #[test]
    fn test_energy_filter_logspace_bins() {
        let min_energy: f64 = 0.1;
        let max_energy: f64 = 20e6;
        let num_bins = 50;
        
        let log_min = min_energy.log10();
        let log_max = max_energy.log10();
        let mut bins = Vec::with_capacity(num_bins);
        for i in 0..num_bins {
            let log_value = log_min + (log_max - log_min) * (i as f64) / ((num_bins - 1) as f64);
            bins.push(10_f64.powf(log_value));
        }
        
        let filter = EnergyFilter::new(bins.clone());
        
        assert_eq!(filter.num_bins(), num_bins - 1);
        assert_eq!(filter.bins.len(), num_bins);
        assert!((filter.bins[0] - min_energy).abs() < 1e-10);
        assert!((filter.bins[num_bins - 1] - max_energy).abs() < 1e-3);
        
        assert_eq!(filter.get_bin(0.05), None);
        assert_eq!(filter.get_bin(0.1), Some(0));
        assert!(filter.get_bin(1e6).is_some());
        assert!(filter.get_bin(10e6).is_some());
        assert!(filter.get_bin(19.999e6).is_some());
        assert_eq!(filter.get_bin(25e6), None);
        
        for i in 1..bins.len() {
            assert!(filter.bins[i] > filter.bins[i-1]);
        }
        
        for i in 0..filter.num_bins() {
            let energy = (filter.bins[i] + filter.bins[i+1]) / 2.0;
            assert_eq!(filter.get_bin(energy), Some(i));
        }
    }
}
