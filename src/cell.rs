use crate::region::Region;
use materials_for_mc::Material;
use std::sync::Arc;

/// A Cell represents a geometric region
/// This follows OpenMC's approach where cells are defined by:
/// - A region (combination of surfaces using boolean operations)
/// - A name for identification
#[derive(Clone)]
pub struct Cell {
    pub cell_id: u32,
    pub name: Option<String>,
    pub region: Region,
    pub material: Option<Material>,
}

impl Cell {
    /// Find the closest surface of this cell to a point along a direction (OpenMC/MCNP: first intersection with any region surface)
    pub fn closest_surface(&self, point: [f64; 3], direction: [f64; 3]) -> Option<std::sync::Arc<crate::surface::Surface>> {
        let mut min_dist = f64::INFINITY;
        let mut closest_surface = None;
        for (surface_arc, _sense) in self.region.surfaces_with_sense() {
            let surface: &crate::surface::Surface = surface_arc.as_ref();
            if let Some(dist) = surface.distance_to_surface(point, direction) {
                if dist > 1e-10 && dist < min_dist {
                    min_dist = dist;
                    closest_surface = Some(surface_arc.clone());
                }
            }
        }
        closest_surface
    }

    /// Compute the distance to the closest surface from a point along a direction
    pub fn distance_to_surface(&self, point: [f64; 3], direction: [f64; 3]) -> Option<f64> {
        let mut min_dist = f64::INFINITY;
        for (surface_arc, sense) in self.region.surfaces_with_sense() {
            let surface: &crate::surface::Surface = surface_arc.as_ref();
            if let Some(dist) = surface.distance_to_surface(point, direction) {
                if dist > 1e-10 && self.region.is_exit_surface((point[0], point[1], point[2]), (direction[0], direction[1], direction[2]), surface, dist, sense) {
                    if dist < min_dist {
                        min_dist = dist;
                    }
                }
            }
        }
        if min_dist < f64::INFINITY {
            Some(min_dist)
        } else {
            None
        }
    }
    /// Create a new cell with a region and optional material (fill)
    pub fn new(cell_id: u32, region: Region, name: Option<String>, material: Option<Material>) -> Self {
        Cell {
            cell_id,
            name,
            region,
            material,
        }
    }

    /// Check if a point is inside this cell's region
    pub fn contains(&self, point: (f64, f64, f64)) -> bool {
        self.region.contains(point)
    }
    pub fn material(&self) -> Option<&Material> {
        self.material.as_ref()
    }
}

#[cfg(test)]
mod tests {
// --- Surface distance tests ---
#[cfg(test)]
mod distance_tests {
    use crate::surface::{Surface, SurfaceKind};

    #[test]
    fn test_sphere_distance() {
        let sphere = Surface::new_sphere(0.0, 0.0, 0.0, 1.0, 1, None);
        // From (2,0,0) toward center
        let d = sphere.distance_to_surface([2.0, 0.0, 0.0], [-1.0, 0.0, 0.0]);
        assert!((d.unwrap() - 1.0).abs() < 1e-10);
        // From (0,0,0) outward
        let d2 = sphere.distance_to_surface([0.0, 0.0, 0.0], [1.0, 0.0, 0.0]);
        assert!((d2.unwrap() - 1.0).abs() < 1e-10);
        // No intersection
        let d3 = sphere.distance_to_surface([2.0, 0.0, 0.0], [1.0, 0.0, 0.0]);
        assert_eq!(d3, None);
        // On surface, outward
        let d4 = sphere.distance_to_surface([1.0, 0.0, 0.0], [1.0, 0.0, 0.0]);
        assert_eq!(d4, None);
    }

    #[test]
    fn test_cylinder_distance() {
        // Z-cylinder at (0,0), r=1
        let cyl = Surface::z_cylinder(0.0, 0.0, 1.0, 1, None);
        // From (2,0,0) toward center
        let d = cyl.distance_to_surface([2.0, 0.0, 0.0], [-1.0, 0.0, 0.0]);
        assert!((d.unwrap() - 1.0).abs() < 1e-10);
        // From (0,2,0) toward center
        let d2 = cyl.distance_to_surface([0.0, 2.0, 0.0], [0.0, -1.0, 0.0]);
        assert!((d2.unwrap() - 1.0).abs() < 1e-10);
        // From (0,0,0) radially outward
        let d3 = cyl.distance_to_surface([0.0, 0.0, 0.0], [1.0, 0.0, 0.0]);
        assert!((d3.unwrap() - 1.0).abs() < 1e-10);
        // No intersection
        let d4 = cyl.distance_to_surface([2.0, 0.0, 0.0], [1.0, 0.0, 0.0]);
        assert_eq!(d4, None);
        // On surface, outward
        let d5 = cyl.distance_to_surface([1.0, 0.0, 0.0], [1.0, 0.0, 0.0]);
        assert_eq!(d5, None);
    }

    #[test]
    fn test_xplane_distance() {
        let plane = Surface::x_plane(5.0, 1, None);
        // From (0,0,0) in +x direction
        let d = plane.distance_to_surface([0.0, 0.0, 0.0], [1.0, 0.0, 0.0]);
        assert_eq!(d, Some(5.0));
        // From (0,0,0) in -x direction
        let d2 = plane.distance_to_surface([0.0, 0.0, 0.0], [-1.0, 0.0, 0.0]);
        assert_eq!(d2, None);
        // From (10,0,0) in -x direction
        let d3 = plane.distance_to_surface([10.0, 0.0, 0.0], [-1.0, 0.0, 0.0]);
        assert_eq!(d3, Some(5.0));
        // Parallel direction
        let d4 = plane.distance_to_surface([0.0, 0.0, 0.0], [0.0, 1.0, 0.0]);
        assert_eq!(d4, None);
        // On plane, outward
        let d5 = plane.distance_to_surface([5.0, 0.0, 0.0], [1.0, 0.0, 0.0]);
        assert_eq!(d5, None);
    }
}
    #[test]
    fn test_cell_fill_material() {
        use materials_for_mc::Material;
        use crate::region::{Region, HalfspaceType};
        use crate::surface::{Surface, SurfaceKind, BoundaryType};
        use std::sync::Arc;

        let s1 = Surface {
            surface_id: 1,
            kind: SurfaceKind::Sphere {
                x0: 0.0,
                y0: 0.0,
                z0: 0.0,
                radius: 1.0,
            },
            boundary_type: BoundaryType::default(),
        };
        let region = Region::new_from_halfspace(HalfspaceType::Below(Arc::new(s1)));

        let mat = Material::new();
        let cell = Cell::new(1, region, Some("filled".to_string()), Some(mat.clone()));
        assert!(cell.material().is_some());
        // The default Material::new() has an empty nuclides map
        assert_eq!(cell.material().unwrap().nuclides.len(), 0);

        // Optional fill
        let cell2 = Cell::new(2, cell.region.clone(), Some("empty".to_string()), None);
        assert!(cell2.material().is_none());
    }
    #[test]
    fn test_cell_union_region() {
        // Union of two spheres
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
        let s2 = Surface {
            surface_id: 2,
            kind: SurfaceKind::Sphere {
                x0: 3.0,
                y0: 0.0,
                z0: 0.0,
                radius: 2.0,
            },
            boundary_type: BoundaryType::default(),
        };
        let region1 = Region::new_from_halfspace(HalfspaceType::Below(Arc::new(s1)));
        let region2 = Region::new_from_halfspace(HalfspaceType::Below(Arc::new(s2)));
        let region = region1.union(&region2);
    let cell = Cell::new(100, region, Some("union".to_string()), None);
        assert!(cell.contains((0.0, 0.0, 0.0))); // inside first sphere
        assert!(cell.contains((3.0, 0.0, 0.0))); // inside second sphere
        assert!(!cell.contains((6.0, 0.0, 0.0))); // outside both
    }

    #[test]
    fn test_cell_intersection_region() {
        // Intersection of two spheres
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
        let s2 = Surface {
            surface_id: 2,
            kind: SurfaceKind::Sphere {
                x0: 1.0,
                y0: 0.0,
                z0: 0.0,
                radius: 2.0,
            },
            boundary_type: BoundaryType::default(),
        };
        let region1 = Region::new_from_halfspace(HalfspaceType::Below(Arc::new(s1)));
        let region2 = Region::new_from_halfspace(HalfspaceType::Below(Arc::new(s2)));
        let region = region1.intersection(&region2);
    let cell = Cell::new(101, region, Some("intersection".to_string()), None);
        assert!(cell.contains((0.0, 0.0, 0.0))); // inside both
        assert!(cell.contains((1.0, 0.0, 0.0))); // inside both
        assert!(!cell.contains((3.0, 0.0, 0.0))); // outside both
    }

    #[test]
    fn test_cell_complement_region() {
        // Complement of a sphere
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
        let region_complement = region.complement();
    let cell = Cell::new(102, region_complement, Some("complement".to_string()), None);
        assert!(!cell.contains((0.0, 0.0, 0.0))); // inside original sphere
        assert!(cell.contains((3.0, 0.0, 0.0))); // outside original sphere
    }
    #[test]
    fn test_cell_complex_region() {
        // s1: x = 2.1, s2: x = -2.1, s3: sphere at origin, r=4.2
        let s1 = Surface {
            surface_id: 5,
            kind: SurfaceKind::Plane {
                a: 1.0,
                b: 0.0,
                c: 0.0,
                d: 2.1,
            }, // x = 2.1
            boundary_type: BoundaryType::default(),
        };
        let s2 = Surface {
            surface_id: 6,
            kind: SurfaceKind::Plane {
                a: 1.0,
                b: 0.0,
                c: 0.0,
                d: -2.1,
            }, // x = -2.1
            boundary_type: BoundaryType::default(),
        };
        let s3 = Surface {
            surface_id: 1,
            kind: SurfaceKind::Sphere {
                x0: 0.0,
                y0: 0.0,
                z0: 0.0,
                radius: 4.2,
            },
            boundary_type: BoundaryType::default(),
        };
        let region = Region::new_from_halfspace(HalfspaceType::Below(Arc::new(s1)))
            .intersection(&Region::new_from_halfspace(HalfspaceType::Above(Arc::new(
                s2,
            ))))
            .intersection(&Region::new_from_halfspace(HalfspaceType::Below(Arc::new(
                s3,
            ))));
    let cell = Cell::new(42, region, Some("complex".to_string()), None);
        // Point inside all constraints
        assert!(cell.contains((0.0, 0.0, 0.0)));
        // Point outside s1 (x > 2.1)
        assert!(!cell.contains((3.0, 0.0, 0.0)));
        // Point outside s2 (x < -2.1)
        assert!(!cell.contains((-3.0, 0.0, 0.0)));
        // Point outside sphere (r > 4.2)
        assert!(!cell.contains((0.0, 0.0, 5.0)));
    }
    use super::*;
    use crate::region::{HalfspaceType, Region};
    use crate::surface::{BoundaryType, Surface, SurfaceKind};
    use std::sync::Arc;

    #[test]
    fn test_cell_contains_simple() {
        // Sphere of radius 2 at (0,0,0)
        let sphere = Surface {
            surface_id: 1,
            kind: SurfaceKind::Sphere {
                x0: 0.0,
                y0: 0.0,
                z0: 0.0,
                radius: 2.0,
            },
            boundary_type: BoundaryType::default(),
        };
        let region = Region::new_from_halfspace(HalfspaceType::Below(Arc::new(sphere)));
    let cell = Cell::new(1, region, None, None);
        assert!(cell.contains((0.0, 0.0, 0.0)));
        assert!(!cell.contains((3.0, 0.0, 0.0)));
    }

    #[test]
    fn test_cell_union_intersection_complement() {
        // Two spheres
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
        let s2 = Surface {
            surface_id: 2,
            kind: SurfaceKind::Sphere {
                x0: 2.0,
                y0: 0.0,
                z0: 0.0,
                radius: 2.0,
            },
            boundary_type: BoundaryType::default(),
        };
        let region1 = Region::new_from_halfspace(HalfspaceType::Below(Arc::new(s1.clone())));
        let region2 = Region::new_from_halfspace(HalfspaceType::Below(Arc::new(s2.clone())));
        // Union
    let union_cell = Cell::new(2, region1.clone().union(&region2.clone()), None, None);
        assert!(union_cell.contains((0.0, 0.0, 0.0)));
        assert!(union_cell.contains((2.0, 0.0, 0.0)));
        assert!(!union_cell.contains((5.0, 0.0, 0.0)));
        // Intersection
    let intersection_cell = Cell::new(3, region1.clone().intersection(&region2.clone()), None, None);
        assert!(!intersection_cell.contains((0.0, 0.0, 0.0)));
        assert!(intersection_cell.contains((1.0, 0.0, 0.0)));
        // Complement
    let complement_cell = Cell::new(4, region1.complement(), None, None);
        assert!(!complement_cell.contains((0.0, 0.0, 0.0)));
        assert!(complement_cell.contains((5.0, 0.0, 0.0)));
    }

    #[test]
    fn test_cell_naming() {
        let sphere = Surface {
            surface_id: 1,
            kind: SurfaceKind::Sphere {
                x0: 0.0,
                y0: 0.0,
                z0: 0.0,
                radius: 2.0,
            },
            boundary_type: BoundaryType::default(),
        };
        let region = Region::new_from_halfspace(HalfspaceType::Below(Arc::new(sphere)));
    let cell = Cell::new(1, region, Some("fuel".to_string()), None);
        assert_eq!(cell.name, Some("fuel".to_string()));
    }
}
