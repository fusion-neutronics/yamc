use crate::cell::Cell;
use std::collections::HashSet;

/// Geometry is a collection of cells for Monte Carlo transport

#[derive(Debug, Clone)]
pub struct Geometry {
    pub cells: Vec<Cell>,
}

impl Geometry {
    /// Create a new geometry with automatic cell ID validation and generation
    pub fn new(mut cells: Vec<Cell>) -> Result<Self, String> {
        // First, validate that all specified cell IDs are unique and positive
        let mut used_ids = HashSet::new();
        let mut cells_needing_ids = Vec::new();
        
        for (index, cell) in cells.iter().enumerate() {
            match cell.cell_id {
                None => {
                    cells_needing_ids.push(index);
                }
                Some(id) => {
                    if used_ids.contains(&id) {
                        return Err(format!("Duplicate cell_id {} found. All cell IDs must be unique.", id));
                    }
                    used_ids.insert(id);
                }
            }
        }
        
        // Generate IDs for cells that need them (cell_id == None)
        let mut next_id = 1;
        for &index in &cells_needing_ids {
            // Find next available ID
            while used_ids.contains(&next_id) {
                next_id += 1;
            }
            cells[index].set_cell_id(next_id);
            used_ids.insert(next_id);
            next_id += 1;
        }
        
        Ok(Geometry { cells })
    }

    /// Find the first cell containing the given point, or None if not found
    pub fn find_cell(&self, point: (f64, f64, f64)) -> Option<&Cell> {
        self.cells.iter().find(|cell| cell.contains(point))
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
    fn test_find_cell() {
        let s1 = Surface {
            surface_id: 1,
            kind: SurfaceKind::Sphere {
                x0: 0.0,
                y0: 0.0,
                z0: 0.0,
                radius: 2.0,
            },
            boundary_type: BoundaryType::default(),
        };
        let region = Region::new_from_halfspace(HalfspaceType::Below(Arc::new(s1)));
        let cell = Cell::with_id(1, region, Some("cell1".to_string()), None);
        let geometry = Geometry::new(vec![cell.clone()]).expect("Failed to create geometry");
        // Point inside the sphere
        assert!(geometry.find_cell((0.0, 0.0, 0.0)).is_some());
        // Point outside the sphere
        assert!(geometry.find_cell((5.0, 0.0, 0.0)).is_none());
    }

    #[test]
    fn test_cell_id_validation() {
        use crate::region::{HalfspaceType, Region};
        use crate::surface::{BoundaryType, Surface, SurfaceKind};
        
        let s1 = Surface {
            surface_id: 1,
            kind: SurfaceKind::Sphere { x0: 0.0, y0: 0.0, z0: 0.0, radius: 1.0 },
            boundary_type: BoundaryType::default(),
        };
        let region = Region::new_from_halfspace(HalfspaceType::Below(Arc::new(s1)));
        
        // Test duplicate cell IDs
        let cell1 = Cell::with_id(1, region.clone(), Some("cell1".to_string()), None);
        let cell2 = Cell::with_id(1, region.clone(), Some("cell2".to_string()), None);
        let result = Geometry::new(vec![cell1, cell2]);
        assert!(result.is_err());
        assert!(result.unwrap_err().contains("Duplicate cell_id 1"));
    }

    #[test]
    fn test_auto_id_generation() {
        use crate::region::{HalfspaceType, Region};
        use crate::surface::{BoundaryType, Surface, SurfaceKind};
        
        let s1 = Surface {
            surface_id: 1,
            kind: SurfaceKind::Sphere { x0: 0.0, y0: 0.0, z0: 0.0, radius: 1.0 },
            boundary_type: BoundaryType::default(),
        };
        let region = Region::new_from_halfspace(HalfspaceType::Below(Arc::new(s1)));
        
        // Test auto-generation of cell IDs (cell_id = None means auto-generate)
        let cell1 = Cell::new(None, region.clone(), Some("cell1".to_string()), None);
        let cell2 = Cell::new(None, region.clone(), Some("cell2".to_string()), None);
        let cell3 = Cell::with_id(5, region.clone(), Some("cell3".to_string()), None);
        
        let geometry = Geometry::new(vec![cell1, cell2, cell3]).expect("Failed to create geometry");
        
        // Check that IDs were assigned
        let ids: Vec<u32> = geometry.cells.iter().map(|c| c.cell_id.unwrap()).collect();
        assert!(ids.contains(&5)); // Explicitly set ID should remain
        // All cells should have IDs assigned
        
        // All IDs should be unique
        let mut unique_ids = ids.clone();
        unique_ids.sort();
        unique_ids.dedup();
        assert_eq!(ids.len(), unique_ids.len());
    }

    #[test]
    fn test_mixed_id_assignment() {
        use crate::region::{HalfspaceType, Region};
        use crate::surface::{BoundaryType, Surface, SurfaceKind};
        
        let s1 = Surface {
            surface_id: 1,
            kind: SurfaceKind::Sphere { x0: 0.0, y0: 0.0, z0: 0.0, radius: 1.0 },
            boundary_type: BoundaryType::default(),
        };
        let region = Region::new_from_halfspace(HalfspaceType::Below(Arc::new(s1)));
        
        // Mix of auto-generated (None) and explicit IDs
        let cell1 = Cell::new(None, region.clone(), Some("auto1".to_string()), None);
        let cell2 = Cell::with_id(3, region.clone(), Some("explicit".to_string()), None);
        let cell3 = Cell::new(None, region.clone(), Some("auto2".to_string()), None);
        
        let geometry = Geometry::new(vec![cell1, cell2, cell3]).expect("Failed to create geometry");
        
        let ids: Vec<u32> = geometry.cells.iter().map(|c| c.cell_id.unwrap()).collect();
        
        // Should contain the explicit ID
        assert!(ids.contains(&3));
        // All cells should have assigned IDs
        assert_eq!(ids.len(), 3);
        // Should have 3 unique IDs
        let mut unique_ids = ids.clone();
        unique_ids.sort();
        unique_ids.dedup();
        assert_eq!(3, unique_ids.len());
    }
}
