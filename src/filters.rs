use crate::cell::Cell;

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
    pub fn new(cell: &Cell) -> Self {
        Self {
            cell_id: cell.cell_id,
        }
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
        self.cell_id == cell.cell_id
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::cell::Cell;
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
        let cell = Cell {
            cell_id: 42,
            name: Some("test_cell".to_string()),
            region,
            material: None,
        };

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
        
        let cell = Cell {
            cell_id: 42,
            name: Some("test_cell".to_string()),
            region: region1,
            material: None,
        };

        let filter = CellFilter::new(&cell);
        
        // Test matching by cell ID
        assert!(filter.matches(42), "Filter should match cell ID 42");
        assert!(!filter.matches(43), "Filter should not match cell ID 43");
        
        // Test matching by cell object
        assert!(filter.matches_cell(&cell), "Filter should match the original cell");
        
        // Create another cell with different ID
        let other_cell = Cell {
            cell_id: 99,
            name: Some("other_cell".to_string()),
            region: region2,
            material: None,
        };
        
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
        
        let cell1 = Cell {
            cell_id: 42,
            name: Some("test_cell".to_string()),
            region: region.clone(),
            material: None,
        };
        
        let cell2 = Cell {
            cell_id: 42,
            name: Some("another_name".to_string()), // Different name, same ID
            region: region.clone(),
            material: None,
        };
        
        let filter1 = CellFilter::new(&cell1);
        let filter2 = CellFilter::new(&cell2);
        
        // Filters should be equal if they have the same cell_id
        assert_eq!(filter1, filter2, "Filters with same cell_id should be equal");
    }
}
