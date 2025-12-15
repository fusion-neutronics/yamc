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

#[cfg(test)]
mod tests {
    use super::*;
    use crate::region::{HalfspaceType, Region};
    use crate::surface::{BoundaryType, Surface, SurfaceKind};
    use std::sync::Arc;

    #[test]
    fn test_cell_filter_creation() {
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

        let filter = CellFilter::new(&cell);
        assert_eq!(filter.cell_id, 42);
    }

    #[test]
    fn test_cell_filter_matching() {
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

        assert!(filter.matches(42));
        assert!(!filter.matches(43));
        assert!(filter.matches_cell(&cell));

        let other_cell = Cell::new(Some(99), region2, Some("other_cell".to_string()), None);
        assert!(!filter.matches_cell(&other_cell));
    }

    #[test]
    fn test_cell_filter_equality() {
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

        let cell1 = Cell::new(Some(42), region.clone(), Some("test_cell".to_string()), None);
        let cell2 = Cell::new(Some(42), region.clone(), Some("another_name".to_string()), None);

        let filter1 = CellFilter::new(&cell1);
        let filter2 = CellFilter::new(&cell2);

        assert_eq!(filter1, filter2);
    }
}
