use crate::cell::Cell;
use crate::material::Material;

/// Cell filter for tallies - filters events based on which cell they occur in
#[derive(Debug, Clone, PartialEq, Eq, Hash)]
pub struct CellFilter {
    /// The cell ID to filter on
    pub cell_id: u32,
}

impl CellFilter {
    /// Create a CellFilter from a Cell object
    ///
    /// # Arguments
    /// * `cell` - The cell object to create the filter from
    ///
    /// # Returns
    /// A new CellFilter that will match events in this cell
    ///
    /// # Panics
    /// Panics if the cell has no cell_id (None cells can't be filtered)
    pub fn new(cell: &Cell) -> Self {
        let cell_id = cell
            .cell_id
            .expect("Cannot create CellFilter for cell with no ID - assign a cell_id first");
        Self { cell_id }
    }

    /// Check if this filter matches a given cell ID
    ///
    /// # Arguments
    /// * `cell_id` - The cell ID to check against
    ///
    /// # Returns
    /// `true` if the cell ID matches this filter
    pub fn matches(&self, cell_id: u32) -> bool {
        self.cell_id == cell_id
    }

    /// Check if this filter matches a given cell object
    ///
    /// # Arguments  
    /// * `cell` - The cell object to check against
    ///
    /// # Returns
    /// `true` if the cell matches this filter
    pub fn matches_cell(&self, cell: &Cell) -> bool {
        cell.cell_id == Some(self.cell_id)
    }
}

/// Material filter for tallies - filters events based on which material they occur in
#[derive(Debug, Clone, PartialEq, Eq, Hash)]
pub struct MaterialFilter {
    /// The material ID to filter on (always Some since None materials aren't real materials)
    pub material_id: u32,
}

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
            Ok(i) => Some(i), // Exact match with bin edge
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

impl MaterialFilter {
    /// Create a MaterialFilter from a Material object
    ///
    /// # Arguments
    /// * `material` - The material object to create the filter from (must have a material_id)
    ///
    /// # Returns
    /// A new MaterialFilter that will match events in this material
    ///
    /// # Panics
    /// Panics if the material has no material_id (None materials can't be filtered)
    pub fn new(material: &Material) -> Self {
        let material_id = material.material_id.expect(
            "Cannot create MaterialFilter for material with no ID - assign a material_id first",
        );

        Self { material_id }
    }

    /// Check if this filter matches a given material ID
    ///
    /// # Arguments
    /// * `material_id` - The material ID to check against
    ///
    /// # Returns
    /// `true` if the material ID matches this filter
    pub fn matches(&self, material_id: Option<u32>) -> bool {
        material_id == Some(self.material_id)
    }

    /// Check if this filter matches a given material object
    ///
    /// # Arguments  
    /// * `material` - The material object to check against
    ///
    /// # Returns
    /// `true` if the material matches this filter
    pub fn matches_material(&self, material: &Material) -> bool {
        material.material_id == Some(self.material_id)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::cell::Cell;
    use crate::material::Material;
    use crate::region::{HalfspaceType, Region};
    use crate::surface::{BoundaryType, Surface, SurfaceKind};
    use std::sync::Arc;

    #[test]
    fn test_cell_filter_creation() {
        // Create a simple cell for testing
        let sphere = Surface {
            surface_id: Some(1),
            kind: SurfaceKind::Sphere {
                x0: 0.0,
                y0: 0.0,
                z0: 0.0,
                radius: 2.0,
            },
            boundary_type: BoundaryType::Vacuum,
        };
        let region = Region::new_from_halfspace(HalfspaceType::Below(Arc::new(sphere)));
        let cell = Cell::new(Some(42), region, Some("test_cell".to_string()), None);

        // Create filter from cell
        let filter = CellFilter::new(&cell);

        // Verify the filter was created correctly
        assert_eq!(filter.cell_id, 42);
    }

    #[test]
    fn test_cell_filter_matching() {
        // Create surfaces and regions for testing
        let sphere1 = Surface {
            surface_id: Some(1),
            kind: SurfaceKind::Sphere {
                x0: 0.0,
                y0: 0.0,
                z0: 0.0,
                radius: 2.0,
            },
            boundary_type: BoundaryType::Vacuum,
        };
        let region1 = Region::new_from_halfspace(HalfspaceType::Below(Arc::new(sphere1)));

        let sphere2 = Surface {
            surface_id: Some(2),
            kind: SurfaceKind::Sphere {
                x0: 1.0,
                y0: 1.0,
                z0: 1.0,
                radius: 3.0,
            },
            boundary_type: BoundaryType::Vacuum,
        };
        let region2 = Region::new_from_halfspace(HalfspaceType::Below(Arc::new(sphere2)));

        let cell = Cell::new(Some(42), region1, Some("test_cell".to_string()), None);

        let filter = CellFilter::new(&cell);

        // Test matching by cell ID
        assert!(filter.matches(42), "Filter should match cell ID 42");
        assert!(!filter.matches(43), "Filter should not match cell ID 43");

        // Test matching by cell object
        assert!(
            filter.matches_cell(&cell),
            "Filter should match the original cell"
        );

        // Create another cell with different ID
        let other_cell = Cell::new(Some(99), region2, Some("other_cell".to_string()), None);

        assert!(
            !filter.matches_cell(&other_cell),
            "Filter should not match different cell"
        );
    }

    #[test]
    fn test_cell_filter_equality() {
        // Create two identical cells
        let sphere = Surface {
            surface_id: Some(1),
            kind: SurfaceKind::Sphere {
                x0: 0.0,
                y0: 0.0,
                z0: 0.0,
                radius: 2.0,
            },
            boundary_type: BoundaryType::Vacuum,
        };
        let region = Region::new_from_halfspace(HalfspaceType::Below(Arc::new(sphere)));

        let cell1 = Cell::new(
            Some(42),
            region.clone(),
            Some("test_cell".to_string()),
            None,
        );

        let cell2 = Cell::new(
            Some(42),
            region.clone(),
            Some("another_name".to_string()),
            None,
        );

        let filter1 = CellFilter::new(&cell1);
        let filter2 = CellFilter::new(&cell2);

        // Filters should be equal if they have the same cell_id
        assert_eq!(
            filter1, filter2,
            "Filters with same cell_id should be equal"
        );
    }

    #[test]
    fn test_material_filter_creation() {
        let mut material = Material::with_id(123);
        material.set_name("test_material");
        let filter = MaterialFilter::new(&material);
        assert_eq!(filter.material_id, 123);
    }

    #[test]
    #[should_panic(expected = "Cannot create MaterialFilter for material with no ID")]
    fn test_material_filter_creation_no_id_panics() {
        let material = Material::new(); // No ID
        MaterialFilter::new(&material); // Should panic
    }

    #[test]
    fn test_material_filter_matching() {
        let mut material = Material::with_id(123);
        material.set_name("test_material");
        let filter = MaterialFilter::new(&material);

        assert!(
            filter.matches(Some(123)),
            "Filter should match material ID 123"
        );
        assert!(
            !filter.matches(Some(124)),
            "Filter should not match material ID 124"
        );
        assert!(
            !filter.matches(None),
            "Filter should not match None material"
        );

        assert!(
            filter.matches_material(&material),
            "Filter should match the original material"
        );

        let mut other_material = Material::with_id(456);
        other_material.set_name("other_material");
        assert!(
            !filter.matches_material(&other_material),
            "Filter should not match different material"
        );
    }

    #[test]
    fn test_material_filter_equality() {
        let mut material1 = Material::with_id(123);
        material1.set_name("material_1");
        let mut material2 = Material::with_id(123);
        material2.set_name("material_2");

        let filter1 = MaterialFilter::new(&material1);
        let filter2 = MaterialFilter::new(&material2);

        assert_eq!(
            filter1, filter2,
            "Filters with same material_id should be equal"
        );

        let mut material3 = Material::with_id(456);
        material3.set_name("material_3");
        let filter3 = MaterialFilter::new(&material3);

        assert_ne!(
            filter1, filter3,
            "Filters with different material_id should not be equal"
        );
    }

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
        
        // Test bins: [0, 1e6), [1e6, 10e6), [10e6, 20e6)
        assert_eq!(filter.get_bin(0.0), Some(0));
        assert_eq!(filter.get_bin(5e5), Some(0));
        assert_eq!(filter.get_bin(9.99e5), Some(0));
        
        assert_eq!(filter.get_bin(1e6), Some(1));
        assert_eq!(filter.get_bin(5e6), Some(1));
        assert_eq!(filter.get_bin(9.99e6), Some(1));
        
        assert_eq!(filter.get_bin(10e6), Some(2));
        assert_eq!(filter.get_bin(15e6), Some(2));
        assert_eq!(filter.get_bin(19.99e6), Some(2));
        
        // Out of range
        assert_eq!(filter.get_bin(-1.0), None);
        assert_eq!(filter.get_bin(20e6), None);
        assert_eq!(filter.get_bin(25e6), None);
    }

    #[test]
    fn test_energy_filter_matches() {
        let filter = EnergyFilter::new(vec![1e6, 10e6, 20e6]);
        
        assert!(!filter.matches(0.5e6), "Below range should not match");
        assert!(filter.matches(1e6), "Lower bound should match");
        assert!(filter.matches(5e6), "Middle of bin should match");
        assert!(filter.matches(15e6), "Second bin should match");
        assert!(!filter.matches(20e6), "Upper bound should not match");
        assert!(!filter.matches(25e6), "Above range should not match");
    }

    #[test]
    fn test_energy_filter_typical_neutron_bins() {
        // Typical neutron energy bins: thermal, epithermal, fast
        let filter = EnergyFilter::new(vec![
            0.0,      // 0 eV
            0.625,    // 0.625 eV (thermal cutoff)
            100e3,    // 100 keV (epithermal cutoff)
            20e6,     // 20 MeV (fast cutoff)
        ]);
        
        assert_eq!(filter.num_bins(), 3);
        
        // Thermal neutron (0.025 eV)
        assert_eq!(filter.get_bin(0.025), Some(0));
        
        // Epithermal neutron (1 keV = 1000 eV)
        assert_eq!(filter.get_bin(1000.0), Some(1));
        
        // Fast neutron (14.1 MeV)
        assert_eq!(filter.get_bin(14.1e6), Some(2));
    }

    #[test]
    fn test_energy_filter_equality() {
        let filter1 = EnergyFilter::new(vec![0.0, 1e6, 10e6]);
        let filter2 = EnergyFilter::new(vec![0.0, 1e6, 10e6]);
        let filter3 = EnergyFilter::new(vec![0.0, 2e6, 10e6]);
        
        assert_eq!(filter1, filter2, "Identical filters should be equal");
        assert_ne!(filter1, filter3, "Different bins should not be equal");
    }

    #[test]
    fn test_energy_filter_logspace_bins() {
        // Test with logarithmically spaced bins from 0.1 eV to 20 MeV (similar to Python example)
        let min_energy: f64 = 0.1;
        let max_energy: f64 = 20e6;
        let num_bins = 50;
        
        // Create logspace bins manually
        let log_min = min_energy.log10();
        let log_max = max_energy.log10();
        let mut bins = Vec::with_capacity(num_bins);
        for i in 0..num_bins {
            let log_value = log_min + (log_max - log_min) * (i as f64) / ((num_bins - 1) as f64);
            bins.push(10_f64.powf(log_value));
        }
        
        let filter = EnergyFilter::new(bins.clone());
        
        // Verify filter properties
        assert_eq!(filter.num_bins(), num_bins - 1);
        assert_eq!(filter.bins.len(), num_bins);
        assert!((filter.bins[0] - min_energy).abs() < 1e-10);
        assert!((filter.bins[num_bins - 1] - max_energy).abs() < 1e-3);
        
        // Test binning at various energies
        assert_eq!(filter.get_bin(0.05), None, "Below min should be None");
        assert_eq!(filter.get_bin(0.1), Some(0), "Min energy should be bin 0");
        assert!(filter.get_bin(1e6).is_some(), "1 MeV should be in range");
        assert!(filter.get_bin(10e6).is_some(), "10 MeV should be in range");
        assert!(filter.get_bin(19.999e6).is_some(), "Just below max should be in range");
        assert_eq!(filter.get_bin(25e6), None, "Above max should be None");
        
        // Verify bins are monotonically increasing
        for i in 1..bins.len() {
            assert!(filter.bins[i] > filter.bins[i-1], 
                "Bins should be strictly increasing: bins[{}]={} should be > bins[{}]={}",
                i, filter.bins[i], i-1, filter.bins[i-1]);
        }
        
        // Test that every bin is accessible
        for i in 0..filter.num_bins() {
            let energy = (filter.bins[i] + filter.bins[i+1]) / 2.0; // Midpoint of bin
            assert_eq!(filter.get_bin(energy), Some(i),
                "Energy {} should be in bin {}", energy, i);
        }
    }
}
