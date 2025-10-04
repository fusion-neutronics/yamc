use crate::cell::Cell;

/// Geometry is a collection of cells for Monte Carlo transport

#[derive(Clone)]
pub struct Geometry {
    pub cells: Vec<Cell>,
}

impl Geometry {
    /// Find the first cell containing the given point, or None if not found
    pub fn find_cell(&self, point: (f64, f64, f64)) -> Option<&Cell> {
        self.cells.iter().find(|cell| cell.contains(point))
    }
}


#[cfg(test)]
mod tests {
    use super::*;
    use crate::cell::Cell;
    use crate::region::{Region, HalfspaceType};
    use crate::surface::{Surface, SurfaceKind, BoundaryType};
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
        let cell = Cell::new(1, region, Some("cell1".to_string()), None);
        let geometry = Geometry { cells: vec![cell.clone()] };
        // Point inside the sphere
        assert!(geometry.find_cell((0.0, 0.0, 0.0)).is_some());
        // Point outside the sphere
        assert!(geometry.find_cell((5.0, 0.0, 0.0)).is_none());
    }
}
