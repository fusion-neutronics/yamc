// ...existing code...
use crate::surface::Surface;
// ...existing code...
use std::sync::Arc;

#[derive(Clone)]
pub struct Region {
    pub expr: RegionExpr,
}

#[derive(Clone)]
pub enum HalfspaceType {
    Above(Arc<Surface>),
    Below(Arc<Surface>),
}

#[derive(Clone)]
pub enum RegionExpr {
    Halfspace(HalfspaceType),
    Union(Box<RegionExpr>, Box<RegionExpr>),
    Intersection(Box<RegionExpr>, Box<RegionExpr>),
    Complement(Box<RegionExpr>),
}

// Regular Rust implementation
impl Region {

    /// Recursively collect all surfaces and their sense (true=Above, false=Below) in the region
    pub fn surfaces_with_sense(&self) -> Vec<(Arc<Surface>, bool)> {
        fn collect(expr: &RegionExpr, surfaces: &mut Vec<(Arc<Surface>, bool)>, sense: bool) {
            match expr {
                RegionExpr::Halfspace(hs) => match hs {
                    HalfspaceType::Above(surf) => surfaces.push((surf.clone(), sense)),
                    HalfspaceType::Below(surf) => surfaces.push((surf.clone(), !sense)),
                },
                RegionExpr::Union(a, b) | RegionExpr::Intersection(a, b) => {
                    collect(a, surfaces, sense);
                    collect(b, surfaces, sense);
                }
                RegionExpr::Complement(inner) => collect(inner, surfaces, !sense),
            }
        }
        let mut result = Vec::new();
        collect(&self.expr, &mut result, true);
        result
    }

    /// Check if crossing the given surface at the intersection would exit the region
    /// For unions: exit only if outside all subregions. For intersections: exit if outside any subregion.
    pub fn is_exit_surface(&self, point: (f64, f64, f64), direction: (f64, f64, f64), _surface: &Surface, dist: f64, _sense: bool) -> bool {
        let eps = 1e-8;
        let p_next = (
            point.0 + direction.0 * (dist + eps),
            point.1 + direction.1 * (dist + eps),
            point.2 + direction.2 * (dist + eps),
        );
        fn check_exit(expr: &RegionExpr, p: (f64, f64, f64)) -> bool {
            match expr {
                RegionExpr::Halfspace(hs) => match hs {
                    HalfspaceType::Above(surf) => surf.evaluate(p) <= 0.0,
                    HalfspaceType::Below(surf) => surf.evaluate(p) >= 0.0,
                },
                RegionExpr::Union(a, b) => check_exit(a, p) && check_exit(b, p), // exit union if outside all
                RegionExpr::Intersection(a, b) => check_exit(a, p) || check_exit(b, p), // exit intersection if outside any
                RegionExpr::Complement(inner) => !check_exit(inner, p),
            }
        }
        check_exit(&self.expr, p_next)
    }
    pub fn new_from_halfspace(halfspace_type: HalfspaceType) -> Self {
        Region {
            expr: RegionExpr::Halfspace(halfspace_type),
        }
    }

    pub fn intersection(&self, other: &Self) -> Self {
        Region {
            expr: RegionExpr::Intersection(
                Box::new(self.expr.clone()),
                Box::new(other.expr.clone()),
            ),
        }
    }

    pub fn union(&self, other: &Self) -> Self {
        Region {
            expr: RegionExpr::Union(Box::new(self.expr.clone()), Box::new(other.expr.clone())),
        }
    }

    pub fn complement(&self) -> Self {
        Region {
            expr: RegionExpr::Complement(Box::new(self.expr.clone())),
        }
    }

    // Updated contains method: no surface dictionary needed
    pub fn contains(&self, point: (f64, f64, f64)) -> bool {
        self.expr.evaluate_contains(point)
    }

    // Updated evaluate_contains method: no surface dictionary needed
    pub fn evaluate_contains(&self, point: (f64, f64, f64)) -> bool {
        self.expr.evaluate_contains(point)
    }

    pub fn bounding_box(&self) -> crate::bounding_box::BoundingBox {
        // Collect all axis constraints and finite bounds
        fn collect_constraints(
            expr: &RegionExpr,
            axis_lowers: &mut [f64; 3],
            axis_uppers: &mut [f64; 3],
            finite_bounds: &mut Vec<([f64; 3], [f64; 3])>,
        ) {
            match expr {
                RegionExpr::Halfspace(hs) => match hs {
                    HalfspaceType::Above(surf) => {
                        if let Some((axis, is_upper, value)) = surf.axis_constraint(true) {
                            if is_upper {
                                axis_uppers[axis] = axis_uppers[axis].min(value);
                            } else {
                                axis_lowers[axis] = axis_lowers[axis].max(value);
                            }
                        }
                        if let Some((lower, upper)) = surf.bounding_box(false) {
                            finite_bounds.push((lower, upper));
                        }
                    }
                    HalfspaceType::Below(surf) => {
                        if let Some((axis, is_upper, value)) = surf.axis_constraint(false) {
                            if is_upper {
                                axis_uppers[axis] = axis_uppers[axis].min(value);
                            } else {
                                axis_lowers[axis] = axis_lowers[axis].max(value);
                            }
                        }
                        if let Some((lower, upper)) = surf.bounding_box(true) {
                            finite_bounds.push((lower, upper));
                        }
                    }
                },
                RegionExpr::Intersection(a, b) => {
                    collect_constraints(a, axis_lowers, axis_uppers, finite_bounds);
                    collect_constraints(b, axis_lowers, axis_uppers, finite_bounds);
                }
                RegionExpr::Union(a, b) => {
                    // For union, take the union of bounds (not strict, but matches previous logic)
                    collect_constraints(a, axis_lowers, axis_uppers, finite_bounds);
                    collect_constraints(b, axis_lowers, axis_uppers, finite_bounds);
                }
                RegionExpr::Complement(inner) => {
                    // For complement, ignore constraints (could be improved)
                }
            }
        }

        let mut axis_lowers = [f64::NEG_INFINITY; 3];
        let mut axis_uppers = [f64::INFINITY; 3];
        let mut finite_bounds = Vec::new();
        collect_constraints(
            &self.expr,
            &mut axis_lowers,
            &mut axis_uppers,
            &mut finite_bounds,
        );

        // Intersect all finite bounds
        for (lower, upper) in finite_bounds {
            for i in 0..3 {
                axis_lowers[i] = axis_lowers[i].max(lower[i]);
                axis_uppers[i] = axis_uppers[i].min(upper[i]);
            }
        }

        // If any min > max, region is empty: return empty bounding box
        if axis_lowers[0] > axis_uppers[0]
            || axis_lowers[1] > axis_uppers[1]
            || axis_lowers[2] > axis_uppers[2]
        {
            return crate::bounding_box::BoundingBox::new(
                [f64::INFINITY; 3],
                [f64::NEG_INFINITY; 3],
            );
        }

        crate::bounding_box::BoundingBox::new(axis_lowers, axis_uppers)
    }
}

impl RegionExpr {
    pub fn evaluate_contains(&self, point: (f64, f64, f64)) -> bool {
        match self {
            RegionExpr::Halfspace(hs) => match hs {
                HalfspaceType::Above(surf) => surf.evaluate(point) > 0.0,
                HalfspaceType::Below(surf) => surf.evaluate(point) < 0.0,
            },
            RegionExpr::Union(a, b) => a.evaluate_contains(point) || b.evaluate_contains(point),
            RegionExpr::Intersection(a, b) => {
                a.evaluate_contains(point) && b.evaluate_contains(point)
            }
            RegionExpr::Complement(inner) => !inner.evaluate_contains(point),
        }
    }
}

#[cfg(test)]
#[test]
fn test_sphere_bb_moved_on_z_axis() {
    use crate::region::{HalfspaceType, Region};
    use crate::surface::Surface;
    use std::sync::Arc;
    // Sphere centered at (0, 0, 1) with radius 3
    let s2 = Surface::new_sphere(0.0, 0.0, 1.0, 3.0, 1, None);
    let region2 = Region::new_from_halfspace(HalfspaceType::Below(Arc::new(s2)));
    let bbox = region2.bounding_box();
    assert_eq!(bbox.lower_left, [-3.0, -3.0, -2.0]);
    assert_eq!(bbox.upper_right, [3.0, 3.0, 4.0]);
}

#[test]
fn test_sphere_with_xplanes() {
    use crate::region::{HalfspaceType, Region};
    use crate::surface::Surface;
    use std::sync::Arc;
    // XPlane at x=2.1
    let s1 = Surface::x_plane(2.1, 5, None);
    // XPlane at x=-2.1
    let s2 = Surface::x_plane(-2.1, 6, None);
    // Sphere at (0,0,0) with radius 4.2
    let s3 = Surface::new_sphere(0.0, 0.0, 0.0, 4.2, 1, None);
    // Region: x <= 2.1 & x >= -2.1 & inside sphere
    let region1 = Region::new_from_halfspace(HalfspaceType::Below(Arc::new(s1.clone())))
        .intersection(&Region::new_from_halfspace(HalfspaceType::Above(Arc::new(
            s2.clone(),
        ))))
        .intersection(&Region::new_from_halfspace(HalfspaceType::Below(Arc::new(
            s3.clone(),
        ))));
    let bbox = region1.bounding_box();
    assert_eq!(bbox.lower_left, [-2.1, -4.2, -4.2]);
    assert_eq!(bbox.upper_right, [2.1, 4.2, 4.2]);
}
mod tests {
    use super::*;
    use crate::surface::{Surface, SurfaceKind};
    use std::collections::HashMap;

    #[test]
    fn test_region_contains() {
        // Create two surfaces
        let s1 = Surface {
            surface_id: 1,
            kind: SurfaceKind::Plane {
                a: 0.0,
                b: 0.0,
                c: 1.0,
                d: -5.0,
            },
            boundary_type: crate::surface::BoundaryType::default(),
        };
        let s2 = Surface {
            surface_id: 2,
            kind: SurfaceKind::Sphere {
                x0: 0.0,
                y0: 0.0,
                z0: 0.0,
                radius: 3.0,
            },
            boundary_type: crate::surface::BoundaryType::default(),
        };

        // Map of surfaces by surface_id
        let mut surfaces = HashMap::new();
        surfaces.insert(s1.surface_id, s1.clone());
        surfaces.insert(s2.surface_id, s2.clone());

        // Build a region: inside s2 AND above s1
        let region =
            Region::new_from_halfspace(crate::region::HalfspaceType::Above(Arc::new(s1.clone())))
                .intersection(&Region::new_from_halfspace(
                    crate::region::HalfspaceType::Below(Arc::new(s2.clone())),
                ));

        // Test a point inside both
        let point = (0.0, 0.0, 0.0);
        assert!(region.contains(point));

        // Test a point outside the sphere
        let point = (0.0, 0.0, 4.0);
        assert!(!region.contains(point));
    }

    #[test]
    fn test_sphere_bounding_box() {
        // Sphere of radius 2 at (0,0,0)
        let s = Surface {
            surface_id: 1,
            kind: SurfaceKind::Sphere {
                x0: 0.0,
                y0: 0.0,
                z0: 0.0,
                radius: 2.0,
            },
            boundary_type: crate::surface::BoundaryType::default(),
        };
        let mut surfaces = HashMap::new();
        surfaces.insert(s.surface_id, s.clone());
        let region = Region::new_from_halfspace(HalfspaceType::Below(Arc::new(s.clone())));
        let bbox = region.bounding_box();
        assert_eq!(bbox.lower_left, [-2.0, -2.0, -2.0]);
        assert_eq!(bbox.upper_right, [2.0, 2.0, 2.0]);
    }

    #[test]
    fn test_box_and_sphere_bounding_box() {
        // XPlanes at x=2.1 and x=-2.1, sphere at origin with radius 4.2
        let s1 = Surface {
            surface_id: 1,
            kind: SurfaceKind::Plane {
                a: 1.0,
                b: 0.0,
                c: 0.0,
                d: 2.1,
            },
            boundary_type: crate::surface::BoundaryType::default(),
        };
        let s2 = Surface {
            surface_id: 2,
            kind: SurfaceKind::Plane {
                a: 1.0,
                b: 0.0,
                c: 0.0,
                d: -2.1,
            },
            boundary_type: crate::surface::BoundaryType::default(),
        };
        let s3 = Surface {
            surface_id: 3,
            kind: SurfaceKind::Sphere {
                x0: 0.0,
                y0: 0.0,
                z0: 0.0,
                radius: 4.2,
            },
            boundary_type: crate::surface::BoundaryType::default(),
        };
        let mut surfaces = HashMap::new();
        surfaces.insert(s1.surface_id, s1.clone());
        surfaces.insert(s2.surface_id, s2.clone());
        surfaces.insert(s3.surface_id, s3.clone());
        // Region: x >= -2.1 & x <= 2.1 & inside sphere
        let region = Region::new_from_halfspace(HalfspaceType::Above(Arc::new(s2.clone())))
            .intersection(&Region::new_from_halfspace(HalfspaceType::Below(Arc::new(
                s1.clone(),
            ))))
            .intersection(&Region::new_from_halfspace(HalfspaceType::Below(Arc::new(
                s3.clone(),
            ))));
        let bbox = region.bounding_box();
        assert_eq!(bbox.lower_left, [-2.1, -4.2, -4.2]);
        assert_eq!(bbox.upper_right, [2.1, 4.2, 4.2]);
    }

    #[test]
    fn test_zplane_bounding_box() {
        // ZPlane at z=3.5
        let s = Surface {
            surface_id: 1,
            kind: SurfaceKind::Plane {
                a: 0.0,
                b: 0.0,
                c: 1.0,
                d: 3.5,
            },
            boundary_type: crate::surface::BoundaryType::default(),
        };
        let mut surfaces = HashMap::new();
        surfaces.insert(s.surface_id, s.clone());
        // Region: z < 3.5 (Below ZPlane)
        let region = Region::new_from_halfspace(HalfspaceType::Below(Arc::new(s.clone())));
        let bbox = region.bounding_box();
        assert_eq!(bbox.lower_left[2], f64::NEG_INFINITY);
        assert_eq!(bbox.upper_right[2], 3.5);
        assert_eq!(bbox.lower_left[0], f64::NEG_INFINITY);
        assert_eq!(bbox.upper_right[0], f64::INFINITY);
        assert_eq!(bbox.lower_left[1], f64::NEG_INFINITY);
        assert_eq!(bbox.upper_right[1], f64::INFINITY);
    }

    #[test]
    fn test_xplane_bounding_box() {
        // XPlane at x=1.5
        let s = Surface {
            surface_id: 1,
            kind: SurfaceKind::Plane {
                a: 1.0,
                b: 0.0,
                c: 0.0,
                d: 1.5,
            },
            boundary_type: crate::surface::BoundaryType::default(),
        };
        let mut surfaces = HashMap::new();
        surfaces.insert(s.surface_id, s.clone());
        // Region: x < 1.5 (Below XPlane)
        let region = Region::new_from_halfspace(HalfspaceType::Below(Arc::new(s.clone())));
        let bbox = region.bounding_box();
        assert_eq!(bbox.lower_left[0], f64::NEG_INFINITY);
        assert_eq!(bbox.upper_right[0], 1.5);
        assert_eq!(bbox.lower_left[1], f64::NEG_INFINITY);
        assert_eq!(bbox.upper_right[1], f64::INFINITY);
        assert_eq!(bbox.lower_left[2], f64::NEG_INFINITY);
        assert_eq!(bbox.upper_right[2], f64::INFINITY);
    }

    #[test]
    fn test_yplane_bounding_box() {
        // YPlane at y=-2.0
        let s = Surface {
            surface_id: 1,
            kind: SurfaceKind::Plane {
                a: 0.0,
                b: 1.0,
                c: 0.0,
                d: -2.0,
            },
            boundary_type: crate::surface::BoundaryType::default(),
        };
        let mut surfaces = HashMap::new();
        surfaces.insert(s.surface_id, s.clone());
        // Region: y > -2.0 (Above YPlane)
        let region = Region::new_from_halfspace(HalfspaceType::Above(Arc::new(s.clone())));
        let bbox = region.bounding_box();
        assert_eq!(bbox.lower_left[1], -2.0);
        assert_eq!(bbox.upper_right[1], f64::INFINITY);
        assert_eq!(bbox.lower_left[0], f64::NEG_INFINITY);
        assert_eq!(bbox.upper_right[0], f64::INFINITY);
        assert_eq!(bbox.lower_left[2], f64::NEG_INFINITY);
        assert_eq!(bbox.upper_right[2], f64::INFINITY);
    }

    #[test]
    fn test_zcylinder_bounding_box() {
        // Z-cylinder at (1, 2) with radius 3
        let s = Surface {
            surface_id: 1,
            kind: SurfaceKind::Cylinder {
                axis: [0.0, 0.0, 1.0],
                origin: [1.0, 2.0, 0.0],
                radius: 3.0,
            },
            boundary_type: crate::surface::BoundaryType::default(),
        };
        // Region: inside cylinder (Below)
        let region = Region::new_from_halfspace(HalfspaceType::Below(Arc::new(s.clone())));
        let bbox = region.bounding_box();
        assert_eq!(bbox.lower_left, [-2.0, -1.0, f64::NEG_INFINITY]);
        assert_eq!(bbox.upper_right, [4.0, 5.0, f64::INFINITY]);
    }
}
