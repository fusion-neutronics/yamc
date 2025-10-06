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
        let cell_id = cell.cell_id
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
        let material_id = material.material_id
            .expect("Cannot create MaterialFilter for material with no ID - assign a material_id first");
        
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
            surface_id: 1,
            kind: SurfaceKind::Sphere {
                x0: 0.0,
                y0: 0.0,
                z0: 0.0,
                radius: 2.0,
            },
            boundary_type: BoundaryType::Vacuum,
        };
        let region = Region::new_from_halfspace(HalfspaceType::Below(Arc::new(sphere)));
        let cell = Cell::with_id(42, region, Some("test_cell".to_string()), None);

        // Create filter from cell
        let filter = CellFilter::new(&cell);
        
        // Verify the filter was created correctly
        assert_eq!(filter.cell_id, 42);
    }

    #[test]
    fn test_cell_filter_matching() {
        // Create surfaces and regions for testing
        let sphere1 = Surface {
            surface_id: 1,
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
            surface_id: 2,
            kind: SurfaceKind::Sphere {
                x0: 1.0,
                y0: 1.0,
                z0: 1.0,
                radius: 3.0,
            },
            boundary_type: BoundaryType::Vacuum,
        };
        let region2 = Region::new_from_halfspace(HalfspaceType::Below(Arc::new(sphere2)));
        
        let cell = Cell::with_id(42, region1, Some("test_cell".to_string()), None);

        let filter = CellFilter::new(&cell);
        
        // Test matching by cell ID
        assert!(filter.matches(42), "Filter should match cell ID 42");
        assert!(!filter.matches(43), "Filter should not match cell ID 43");
        
        // Test matching by cell object
        assert!(filter.matches_cell(&cell), "Filter should match the original cell");
        
        // Create another cell with different ID
        let other_cell = Cell::with_id(99, region2, Some("other_cell".to_string()), None);
        
        assert!(!filter.matches_cell(&other_cell), "Filter should not match different cell");
    }

    #[test]
    fn test_cell_filter_equality() {
        // Create two identical cells
        let sphere = Surface {
            surface_id: 1,
            kind: SurfaceKind::Sphere {
                x0: 0.0,
                y0: 0.0,
                z0: 0.0,
                radius: 2.0,
            },
            boundary_type: BoundaryType::Vacuum,
        };
        let region = Region::new_from_halfspace(HalfspaceType::Below(Arc::new(sphere)));
        
        let cell1 = Cell::with_id(42, region.clone(), Some("test_cell".to_string()), None);
        
        let cell2 = Cell::with_id(42, region.clone(), Some("another_name".to_string()), None);
        
        let filter1 = CellFilter::new(&cell1);
        let filter2 = CellFilter::new(&cell2);
        
        // Filters should be equal if they have the same cell_id
        assert_eq!(filter1, filter2, "Filters with same cell_id should be equal");
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
        
        assert!(filter.matches(Some(123)), "Filter should match material ID 123");
        assert!(!filter.matches(Some(124)), "Filter should not match material ID 124");
        assert!(!filter.matches(None), "Filter should not match None material");
        
        assert!(filter.matches_material(&material), "Filter should match the original material");
        
        let mut other_material = Material::with_id(456);
        other_material.set_name("other_material");
        assert!(!filter.matches_material(&other_material), "Filter should not match different material");
    }

    #[test]
    fn test_material_filter_equality() {
        let mut material1 = Material::with_id(123);
        material1.set_name("material_1");
        let mut material2 = Material::with_id(123);
        material2.set_name("material_2");
        
        let filter1 = MaterialFilter::new(&material1);
        let filter2 = MaterialFilter::new(&material2);
        
        assert_eq!(filter1, filter2, "Filters with same material_id should be equal");
        
        let mut material3 = Material::with_id(456);
        material3.set_name("material_3");
        let filter3 = MaterialFilter::new(&material3);
        
        assert_ne!(filter1, filter3, "Filters with different material_id should not be equal");
    }
}
