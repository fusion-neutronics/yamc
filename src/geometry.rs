use crate::cell::Cell;
use std::collections::HashSet;
use std::sync::Arc;

/// Geometry is a collection of cells for Monte Carlo transport

#[derive(Debug, Clone)]
pub struct Geometry {
    pub cells: Vec<Cell>,
}

impl Geometry {
    /// Create a new geometry with automatic cell ID validation and generation,
    /// and material ID validation across all materials used by cells
    pub fn new(mut cells: Vec<Cell>) -> Result<Self, String> {
        // First, validate that all specified cell IDs are unique and positive
        let mut used_cell_ids = HashSet::new();
        let mut cells_needing_ids = Vec::new();
        
        for (index, cell) in cells.iter().enumerate() {
            match cell.cell_id {
                None => {
                    cells_needing_ids.push(index);
                }
                Some(id) => {
                    if used_cell_ids.contains(&id) {
                        return Err(format!("Duplicate cell_id {} found. All cell IDs must be unique.", id));
                    }
                    used_cell_ids.insert(id);
                }
            }
        }
        
        // Generate IDs for cells that need them (cell_id == None)
        let mut next_cell_id = 1;
        for &index in &cells_needing_ids {
            // Find next available ID
            while used_cell_ids.contains(&next_cell_id) {
                next_cell_id += 1;
            }
            cells[index].set_cell_id(next_cell_id);
            used_cell_ids.insert(next_cell_id);
            next_cell_id += 1;
        }
        
        // Validate material IDs across all materials used by cells
        let mut used_material_ids = HashSet::new();
        let mut materials_needing_ids = Vec::new();
        
        for (cell_idx, cell) in cells.iter().enumerate() {
            if let Some(material_arc) = &cell.material {
                let material = material_arc.as_ref();
                match material.material_id {
                    None => {
                        materials_needing_ids.push(cell_idx);
                    }
                    Some(id) => {
                        if used_material_ids.contains(&id) {
                            return Err(format!("Duplicate material_id {} found. All material IDs must be unique across all cells.", id));
                        }
                        used_material_ids.insert(id);
                    }
                }
            }
        }
        
        // Generate IDs for materials that need them (material_id == None)
        let mut next_material_id = 1;
        for &cell_idx in &materials_needing_ids {
            if let Some(material_arc) = &cells[cell_idx].material {
                // Clone the inner Material, mutate, and re-wrap in Arc
                let mut material = (**material_arc).clone();
                while used_material_ids.contains(&next_material_id) {
                    next_material_id += 1;
                }
                material.set_material_id(next_material_id);
                used_material_ids.insert(next_material_id);
                next_material_id += 1;
                cells[cell_idx].material = Some(Arc::new(material));
            }
        }

        // Collect unique surfaces used in the geometry using Arc pointer addresses
        let mut unique_surface_ptrs = HashSet::new();
        let mut unique_surfaces = Vec::new();
        
        for cell in &cells {
            let surfaces = cell.region.surfaces_with_sense();
            for (surface, _sense) in surfaces {
                let ptr = Arc::as_ptr(&surface);
                if unique_surface_ptrs.insert(ptr) {
                    unique_surfaces.push(surface);
                }
            }
        }

        // Validate surface IDs - surfaces with IDs must be unique
        let mut used_surface_ids = HashSet::new();

        for surface in &unique_surfaces {
            if let Some(id) = surface.surface_id {
                if used_surface_ids.contains(&id) {
                    return Err(format!("Duplicate surface_id {} found. All surface IDs must be unique.", id));
                }
                used_surface_ids.insert(id);
            }
            // Surfaces with None ID are allowed - they don't need validation
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
            surface_id: Some(1),
            kind: SurfaceKind::Sphere {
                x0: 0.0,
                y0: 0.0,
                z0: 0.0,
                radius: 2.0,
            },
            boundary_type: BoundaryType::default(),
        };
        let region = Region::new_from_halfspace(HalfspaceType::Below(Arc::new(s1)));
        let cell = Cell::new(Some(1), region, Some("cell1".to_string()), None);
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
            surface_id: Some(1),
            kind: SurfaceKind::Sphere { x0: 0.0, y0: 0.0, z0: 0.0, radius: 1.0 },
            boundary_type: BoundaryType::default(),
        };
        let region = Region::new_from_halfspace(HalfspaceType::Below(Arc::new(s1)));
        
        // Test duplicate cell IDs
        let cell1 = Cell::new(Some(1), region.clone(), Some("cell1".to_string()), None);
        let cell2 = Cell::new(Some(1), region.clone(), Some("cell2".to_string()), None);
        let result = Geometry::new(vec![cell1, cell2]);
        assert!(result.is_err());
        assert!(result.unwrap_err().contains("Duplicate cell_id 1"));
    }

    #[test]
    fn test_auto_id_generation() {
        use crate::region::{HalfspaceType, Region};
        use crate::surface::{BoundaryType, Surface, SurfaceKind};
        
        let s1 = Surface {
            surface_id: Some(1),
            kind: SurfaceKind::Sphere { x0: 0.0, y0: 0.0, z0: 0.0, radius: 1.0 },
            boundary_type: BoundaryType::default(),
        };
        let region = Region::new_from_halfspace(HalfspaceType::Below(Arc::new(s1)));
        
        // Test auto-generation of cell IDs (cell_id = None means auto-generate)
        let cell1 = Cell::new(None, region.clone(), Some("cell1".to_string()), None);
        let cell2 = Cell::new(None, region.clone(), Some("cell2".to_string()), None);
        let cell3 = Cell::new(Some(5), region.clone(), Some("cell3".to_string()), None);
        
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
            surface_id: Some(1),
            kind: SurfaceKind::Sphere { x0: 0.0, y0: 0.0, z0: 0.0, radius: 1.0 },
            boundary_type: BoundaryType::default(),
        };
        let region = Region::new_from_halfspace(HalfspaceType::Below(Arc::new(s1)));
        
        // Mix of auto-generated (None) and explicit IDs
        let cell1 = Cell::new(None, region.clone(), Some("auto1".to_string()), None);
        let cell2 = Cell::new(Some(3), region.clone(), Some("explicit".to_string()), None);
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

    #[test]
    fn test_material_id_validation() {
        use crate::material::Material;
        use crate::region::{HalfspaceType, Region};
        use crate::surface::{BoundaryType, Surface, SurfaceKind};
    use std::sync::Arc;
        
        let s1 = Surface {
            surface_id: Some(1),
            kind: SurfaceKind::Sphere { x0: 0.0, y0: 0.0, z0: 0.0, radius: 1.0 },
            boundary_type: BoundaryType::default(),
        };
        let region = Region::new_from_halfspace(HalfspaceType::Below(Arc::new(s1)));
        
        // Test duplicate material IDs
    let mat1 = Arc::new(Material::with_id(10));
    let mat2 = Arc::new(Material::with_id(10)); // Same ID - should fail
        
        let cell1 = Cell::new(Some(1), region.clone(), Some("cell1".to_string()), Some(mat1));
        let cell2 = Cell::new(Some(2), region.clone(), Some("cell2".to_string()), Some(mat2));
        
        let result = Geometry::new(vec![cell1, cell2]);
        assert!(result.is_err());
        assert!(result.unwrap_err().contains("Duplicate material_id 10"));
    }

    #[test]
    fn test_material_auto_id_generation() {
        use crate::material::Material;
        use crate::region::{HalfspaceType, Region};
        use crate::surface::{BoundaryType, Surface, SurfaceKind};
    use std::sync::Arc;
        
        let s1 = Surface {
            surface_id: Some(1),
            kind: SurfaceKind::Sphere { x0: 0.0, y0: 0.0, z0: 0.0, radius: 1.0 },
            boundary_type: BoundaryType::default(),
        };
        let region = Region::new_from_halfspace(HalfspaceType::Below(Arc::new(s1)));
        
        // Materials without IDs - should get auto-assigned
    let mat1 = Arc::new(Material::new()); // material_id: None
    let mat2 = Arc::new(Material::new()); // material_id: None
        
        let cell1 = Cell::new(Some(1), region.clone(), Some("cell1".to_string()), Some(mat1.clone()));
        let cell2 = Cell::new(Some(2), region.clone(), Some("cell2".to_string()), Some(mat2.clone()));
        
        let geometry = Geometry::new(vec![cell1, cell2]).expect("Failed to create geometry");
        
        // Check that material IDs were auto-assigned
    // After geometry creation, the cell's material Arc may be updated with auto-assigned ID
    let mat1_id = geometry.cells[0].material.as_ref().unwrap().get_material_id();
    let mat2_id = geometry.cells[1].material.as_ref().unwrap().get_material_id();
    assert!(mat1_id.is_some());
    assert!(mat2_id.is_some());
    assert_ne!(mat1_id, mat2_id, "Materials should have different auto-generated IDs");
    }

    #[test]
    fn test_mixed_material_id_assignment() {
        use crate::material::Material;
        use crate::region::{HalfspaceType, Region};
        use crate::surface::{BoundaryType, Surface, SurfaceKind};
    use std::sync::Arc;
        
        let s1 = Surface {
            surface_id: Some(1),
            kind: SurfaceKind::Sphere { x0: 0.0, y0: 0.0, z0: 0.0, radius: 1.0 },
            boundary_type: BoundaryType::default(),
        };
        let region = Region::new_from_halfspace(HalfspaceType::Below(Arc::new(s1)));
        
        // Mix of explicit and auto-generated material IDs
    let mat1 = Arc::new(Material::new());        // auto-generate
    let mat2 = Arc::new(Material::with_id(5));    // explicit ID = 5
    let mat3 = Arc::new(Material::new());        // auto-generate
        
        let cell1 = Cell::new(Some(1), region.clone(), Some("cell1".to_string()), Some(mat1.clone()));
        let cell2 = Cell::new(Some(2), region.clone(), Some("cell2".to_string()), Some(mat2.clone()));
        let cell3 = Cell::new(Some(3), region.clone(), Some("cell3".to_string()), Some(mat3.clone()));
        
    let geometry = Geometry::new(vec![cell1, cell2, cell3]).expect("Failed to create geometry");
    // Check material IDs
    // After geometry creation, the cell's material Arc may be updated with auto-assigned ID
    let mat1_id = geometry.cells[0].material.as_ref().unwrap().get_material_id().unwrap();
    let mat2_id = geometry.cells[1].material.as_ref().unwrap().get_material_id().unwrap();
    let mat3_id = geometry.cells[2].material.as_ref().unwrap().get_material_id().unwrap();
    // Explicit ID should be preserved
    assert_eq!(mat2_id, 5);
    // All IDs should be unique
    let ids = vec![mat1_id, mat2_id, mat3_id];
    let mut unique_ids = ids.clone();
    unique_ids.sort();
    unique_ids.dedup();
    assert_eq!(ids.len(), unique_ids.len(), "All material IDs should be unique");
    }

    #[test]
    fn test_cells_without_materials() {
        use crate::region::{HalfspaceType, Region};
        use crate::surface::{BoundaryType, Surface, SurfaceKind};
        
        let s1 = Surface {
            surface_id: Some(1),
            kind: SurfaceKind::Sphere { x0: 0.0, y0: 0.0, z0: 0.0, radius: 1.0 },
            boundary_type: BoundaryType::default(),
        };
        let region = Region::new_from_halfspace(HalfspaceType::Below(Arc::new(s1)));
        
        // Cells without materials (void cells) should work fine
        let cell1 = Cell::new(Some(1), region.clone(), Some("void_cell1".to_string()), None);
        let cell2 = Cell::new(Some(2), region.clone(), Some("void_cell2".to_string()), None);
        
        let geometry = Geometry::new(vec![cell1, cell2]).expect("Failed to create geometry");
        assert_eq!(geometry.cells.len(), 2);
        assert_eq!(geometry.cells[0].get_cell_id(), Some(1));
        assert_eq!(geometry.cells[1].get_cell_id(), Some(2));
    }

    #[test]
    fn test_surface_id_validation() {
        use crate::region::{HalfspaceType, Region};
        use crate::surface::{BoundaryType, Surface, SurfaceKind};
        
        // Test duplicate surface IDs
        let s1 = Surface {
            surface_id: Some(10),
            kind: SurfaceKind::Sphere { x0: 0.0, y0: 0.0, z0: 0.0, radius: 1.0 },
            boundary_type: BoundaryType::default(),
        };
        let s2 = Surface {
            surface_id: Some(10), // Same ID - should fail
            kind: SurfaceKind::Sphere { x0: 2.0, y0: 0.0, z0: 0.0, radius: 1.0 },
            boundary_type: BoundaryType::default(),
        };
        
        let region1 = Region::new_from_halfspace(HalfspaceType::Below(Arc::new(s1)));
        let region2 = Region::new_from_halfspace(HalfspaceType::Below(Arc::new(s2)));
        
        let cell1 = Cell::new(Some(1), region1, Some("cell1".to_string()), None);
        let cell2 = Cell::new(Some(2), region2, Some("cell2".to_string()), None);
        
        let result = Geometry::new(vec![cell1, cell2]);
        assert!(result.is_err());
        assert!(result.unwrap_err().contains("Duplicate surface_id 10 found"));
    }

    #[test]
    fn test_surface_without_id_validation() {
        use crate::region::{HalfspaceType, Region};
        use crate::surface::{BoundaryType, Surface, SurfaceKind};
        
        // Test surface without ID - should be allowed (None IDs are valid)
        let s1 = Surface {
            surface_id: None, // No ID - should be fine
            kind: SurfaceKind::Sphere { x0: 0.0, y0: 0.0, z0: 0.0, radius: 1.0 },
            boundary_type: BoundaryType::default(),
        };
        
        let region = Region::new_from_halfspace(HalfspaceType::Below(Arc::new(s1)));
        let cell = Cell::new(Some(1), region, Some("cell1".to_string()), None);
        
        let result = Geometry::new(vec![cell]);
        assert!(result.is_ok());
    }

    #[test]
    fn test_valid_surface_ids() {
        use crate::region::{HalfspaceType, Region};
        use crate::surface::{BoundaryType, Surface, SurfaceKind};
        
        // Test valid unique surface IDs
        let s1 = Surface {
            surface_id: Some(1),
            kind: SurfaceKind::Sphere { x0: 0.0, y0: 0.0, z0: 0.0, radius: 1.0 },
            boundary_type: BoundaryType::default(),
        };
        let s2 = Surface {
            surface_id: Some(2),
            kind: SurfaceKind::Sphere { x0: 2.0, y0: 0.0, z0: 0.0, radius: 1.0 },
            boundary_type: BoundaryType::default(),
        };
        
        let region1 = Region::new_from_halfspace(HalfspaceType::Below(Arc::new(s1)));
        let region2 = Region::new_from_halfspace(HalfspaceType::Below(Arc::new(s2)));
        
        let cell1 = Cell::new(Some(1), region1, Some("cell1".to_string()), None);
        let cell2 = Cell::new(Some(2), region2, Some("cell2".to_string()), None);
        
        let geometry = Geometry::new(vec![cell1, cell2]).expect("Failed to create geometry");
        assert_eq!(geometry.cells.len(), 2);
    }
}
