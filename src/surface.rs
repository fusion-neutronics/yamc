use crate::region::{HalfspaceType, RegionExpr};
use std::sync::Arc;

#[derive(Clone, Debug, PartialEq)]
pub enum BoundaryType {
    Transmission,
    Vacuum,
}

impl Default for BoundaryType {
    fn default() -> Self {
        BoundaryType::Vacuum
    }
}

impl BoundaryType {
    /// Parse a boundary type from a string, returning None for invalid strings
    pub fn from_str_option(s: &str) -> Option<Self> {
        match s.to_lowercase().as_str() {
            "transmission" => Some(BoundaryType::Transmission),
            "vacuum" => Some(BoundaryType::Vacuum),
            _ => None,
        }
    }
}

#[derive(Clone)]
pub struct Surface {
    pub surface_id: usize,
    pub kind: SurfaceKind,
    pub boundary_type: BoundaryType,
}

#[derive(Clone)]
pub enum SurfaceKind {
    Plane {
        a: f64,
        b: f64,
        c: f64,
        d: f64,
    },
    Sphere {
        x0: f64,
        y0: f64,
        z0: f64,
        radius: f64,
    },
    Cylinder {
        axis: [f64; 3],
        origin: [f64; 3],
        radius: f64,
    },
}

// Regular Rust implementation
impl Surface {

    /// Compute the distance from a point along a direction to the surface.
    /// Returns Some(distance) if intersection exists and distance > 0, else None.
    pub fn distance_to_surface(&self, point: [f64; 3], direction: [f64; 3]) -> Option<f64> {
        match &self.kind {
            SurfaceKind::Plane { a, b, c, d } => {
                // Plane: ax + by + cz - d = 0
                let denom = a * direction[0] + b * direction[1] + c * direction[2];
                if denom.abs() < 1e-12 {
                    // Parallel, no intersection
                    return None;
                }
                let num = d - (a * point[0] + b * point[1] + c * point[2]);
                let t = num / denom;
                if t > 0.0 {
                    Some(t)
                } else {
                    None
                }
            }
            SurfaceKind::Sphere { x0, y0, z0, radius } => {
                // Ray-sphere intersection: (p + t*v - c)·(p + t*v - c) = r^2
                let oc = [point[0] - x0, point[1] - y0, point[2] - z0];
                let a = direction[0]*direction[0] + direction[1]*direction[1] + direction[2]*direction[2];
                let b = 2.0 * (oc[0]*direction[0] + oc[1]*direction[1] + oc[2]*direction[2]);
                let c = oc[0]*oc[0] + oc[1]*oc[1] + oc[2]*oc[2] - radius*radius;
                let disc = b*b - 4.0*a*c;
                if disc < 0.0 {
                    return None;
                }
                let sqrt_disc = disc.sqrt();
                let t1 = (-b - sqrt_disc) / (2.0*a);
                let t2 = (-b + sqrt_disc) / (2.0*a);
                // Return the smallest positive t
                if t1 > 1e-12 {
                    Some(t1)
                } else if t2 > 1e-12 {
                    Some(t2)
                } else {
                    None
                }
            }
            SurfaceKind::Cylinder { axis, origin, radius } => {
                // Ray-cylinder intersection (infinite cylinder)
                // Cylinder: ( (p-c) - ((p-c)·a)a )^2 = r^2
                // Ray: p + t*v
                let p = point;
                let v = direction;
                let c = origin;
                let a_axis = axis;
                // Compute d = v - (v·a)a
                let v_dot_a = v[0]*a_axis[0] + v[1]*a_axis[1] + v[2]*a_axis[2];
                let d = [v[0] - v_dot_a*a_axis[0], v[1] - v_dot_a*a_axis[1], v[2] - v_dot_a*a_axis[2]];
                // Compute delta_p = p - c
                let delta_p = [p[0] - c[0], p[1] - c[1], p[2] - c[2]];
                let delta_p_dot_a = delta_p[0]*a_axis[0] + delta_p[1]*a_axis[1] + delta_p[2]*a_axis[2];
                let m = [delta_p[0] - delta_p_dot_a*a_axis[0], delta_p[1] - delta_p_dot_a*a_axis[1], delta_p[2] - delta_p_dot_a*a_axis[2]];
                let a_c = d[0]*d[0] + d[1]*d[1] + d[2]*d[2];
                let b_c = 2.0 * (d[0]*m[0] + d[1]*m[1] + d[2]*m[2]);
                let c_c = m[0]*m[0] + m[1]*m[1] + m[2]*m[2] - radius*radius;
                let disc = b_c*b_c - 4.0*a_c*c_c;
                if disc < 0.0 || a_c.abs() < 1e-12 {
                    return None;
                }
                let sqrt_disc = disc.sqrt();
                let t1 = (-b_c - sqrt_disc) / (2.0*a_c);
                let t2 = (-b_c + sqrt_disc) / (2.0*a_c);
                if t1 > 1e-12 {
                    Some(t1)
                } else if t2 > 1e-12 {
                    Some(t2)
                } else {
                    None
                }
            }
            _ => None,
        }
    }
    pub fn new_plane(
        a: f64,
        b: f64,
        c: f64,
        d: f64,
        surface_id: usize,
        boundary_type: Option<BoundaryType>,
    ) -> Self {
        Surface {
            surface_id,
            kind: SurfaceKind::Plane { a, b, c, d },
            boundary_type: boundary_type.unwrap_or_default(),
        }
    }

    pub fn new_sphere(
        x0: f64,
        y0: f64,
        z0: f64,
        radius: f64,
        surface_id: usize,
        boundary_type: Option<BoundaryType>,
    ) -> Self {
        Surface {
            surface_id,
            kind: SurfaceKind::Sphere { x0, y0, z0, radius },
            boundary_type: boundary_type.unwrap_or_default(),
        }
    }

    pub fn new_cylinder(
        axis: [f64; 3],
        origin: [f64; 3],
        radius: f64,
        surface_id: usize,
        boundary_type: Option<BoundaryType>,
    ) -> Self {
        Surface {
            surface_id,
            kind: SurfaceKind::Cylinder {
                axis,
                origin,
                radius,
            },
            boundary_type: boundary_type.unwrap_or_default(),
        }
    }

    pub fn x_plane(x0: f64, surface_id: usize, boundary_type: Option<BoundaryType>) -> Self {
        Self::new_plane(1.0, 0.0, 0.0, x0, surface_id, boundary_type)
    }

    pub fn y_plane(y0: f64, surface_id: usize, boundary_type: Option<BoundaryType>) -> Self {
        Self::new_plane(0.0, 1.0, 0.0, y0, surface_id, boundary_type)
    }

    pub fn z_plane(z0: f64, surface_id: usize, boundary_type: Option<BoundaryType>) -> Self {
        Self::new_plane(0.0, 0.0, 1.0, z0, surface_id, boundary_type)
    }

    /// Create a cylinder oriented along the Z axis, centered at (x0, y0), with given radius and surface_id
    pub fn z_cylinder(
        x0: f64,
        y0: f64,
        radius: f64,
        surface_id: usize,
        boundary_type: Option<BoundaryType>,
    ) -> Self {
        Self::new_cylinder(
            [0.0, 0.0, 1.0],
            [x0, y0, 0.0],
            radius,
            surface_id,
            boundary_type,
        )
    }

    /// Create a sphere with a specific boundary type
    pub fn sphere(
        x0: f64,
        y0: f64,
        z0: f64,
        radius: f64,
        surface_id: usize,
        boundary_type: Option<BoundaryType>,
    ) -> Self {
        Self::new_sphere(x0, y0, z0, radius, surface_id, boundary_type)
    }

    /// Create a cylinder with individual axis components with a specific boundary type
    pub fn cylinder(
        x0: f64,
        y0: f64,
        z0: f64,
        axis_x: f64,
        axis_y: f64,
        axis_z: f64,
        radius: f64,
        surface_id: usize,
        boundary_type: Option<BoundaryType>,
    ) -> Self {
        Self::new_cylinder(
            [axis_x, axis_y, axis_z],
            [x0, y0, z0],
            radius,
            surface_id,
            boundary_type,
        )
    }

    // Python-friendly functions that accept string boundary types
    pub fn x_plane_str(
        x0: f64,
        surface_id: usize,
        boundary_type: Option<&str>,
    ) -> Result<Self, String> {
        let boundary = match boundary_type {
            Some(s) => Some(
                BoundaryType::from_str_option(s)
                    .ok_or("boundary_type must be 'transmission' or 'vacuum'")?,
            ),
            None => None,
        };
        Ok(Self::x_plane(x0, surface_id, boundary))
    }

    pub fn y_plane_str(
        y0: f64,
        surface_id: usize,
        boundary_type: Option<&str>,
    ) -> Result<Self, String> {
        let boundary = match boundary_type {
            Some(s) => Some(
                BoundaryType::from_str_option(s)
                    .ok_or("boundary_type must be 'transmission' or 'vacuum'")?,
            ),
            None => None,
        };
        Ok(Self::y_plane(y0, surface_id, boundary))
    }

    pub fn z_plane_str(
        z0: f64,
        surface_id: usize,
        boundary_type: Option<&str>,
    ) -> Result<Self, String> {
        let boundary = match boundary_type {
            Some(s) => Some(
                BoundaryType::from_str_option(s)
                    .ok_or("boundary_type must be 'transmission' or 'vacuum'")?,
            ),
            None => None,
        };
        Ok(Self::z_plane(z0, surface_id, boundary))
    }

    pub fn sphere_str(
        x0: f64,
        y0: f64,
        z0: f64,
        radius: f64,
        surface_id: usize,
        boundary_type: Option<&str>,
    ) -> Result<Self, String> {
        let boundary = match boundary_type {
            Some(s) => Some(
                BoundaryType::from_str_option(s)
                    .ok_or("boundary_type must be 'transmission' or 'vacuum'")?,
            ),
            None => None,
        };
        Ok(Self::sphere(x0, y0, z0, radius, surface_id, boundary))
    }

    pub fn cylinder_str(
        x0: f64,
        y0: f64,
        z0: f64,
        axis_x: f64,
        axis_y: f64,
        axis_z: f64,
        radius: f64,
        surface_id: usize,
        boundary_type: Option<&str>,
    ) -> Result<Self, String> {
        let boundary = match boundary_type {
            Some(s) => Some(
                BoundaryType::from_str_option(s)
                    .ok_or("boundary_type must be 'transmission' or 'vacuum'")?,
            ),
            None => None,
        };
        Ok(Self::cylinder(
            x0, y0, z0, axis_x, axis_y, axis_z, radius, surface_id, boundary,
        ))
    }

    pub fn z_cylinder_str(
        x0: f64,
        y0: f64,
        radius: f64,
        surface_id: usize,
        boundary_type: Option<&str>,
    ) -> Result<Self, String> {
        let boundary = match boundary_type {
            Some(s) => Some(
                BoundaryType::from_str_option(s)
                    .ok_or("boundary_type must be 'transmission' or 'vacuum'")?,
            ),
            None => None,
        };
        Ok(Self::z_cylinder(x0, y0, radius, surface_id, boundary))
    }

    pub fn plane_str(
        a: f64,
        b: f64,
        c: f64,
        d: f64,
        surface_id: usize,
        boundary_type: Option<&str>,
    ) -> Result<Self, String> {
        let boundary = match boundary_type {
            Some(s) => Some(
                BoundaryType::from_str_option(s)
                    .ok_or("boundary_type must be 'transmission' or 'vacuum'")?,
            ),
            None => None,
        };
        Ok(Self::new_plane(a, b, c, d, surface_id, boundary))
    }

    /// Get the boundary type of the surface
    pub fn boundary_type(&self) -> &BoundaryType {
        &self.boundary_type
    }

    /// Set the boundary type of the surface
    pub fn set_boundary_type(&mut self, boundary_type: BoundaryType) {
        self.boundary_type = boundary_type;
    }

    pub fn evaluate(&self, point: (f64, f64, f64)) -> f64 {
        match &self.kind {
            SurfaceKind::Plane { a, b, c, d } => a * point.0 + b * point.1 + c * point.2 - d,
            SurfaceKind::Sphere { x0, y0, z0, radius } => {
                let dx = point.0 - x0;
                let dy = point.1 - y0;
                let dz = point.2 - z0;
                (dx * dx + dy * dy + dz * dz).sqrt() - radius
            }
            SurfaceKind::Cylinder {
                axis,
                origin,
                radius,
            } => {
                let v = [
                    point.0 - origin[0],
                    point.1 - origin[1],
                    point.2 - origin[2],
                ];
                let dot = v[0] * axis[0] + v[1] * axis[1] + v[2] * axis[2];
                let d = [
                    v[0] - dot * axis[0],
                    v[1] - dot * axis[1],
                    v[2] - dot * axis[2],
                ];
                (d[0] * d[0] + d[1] * d[1] + d[2] * d[2]).sqrt() - radius
            }
        }
    }

    /// Get the bounding box for this surface when used as a halfspace.
    /// For finite surfaces like spheres and cylinders, returns the bounding box bounds for the negative halfspace (inside).
    /// For axis-aligned planes, returns None as they should only contribute axis constraints.
    /// halfspace_below: true for negative halfspace (inside), false for positive halfspace (outside)
    pub fn bounding_box(&self, halfspace_below: bool) -> Option<([f64; 3], [f64; 3])> {
        match &self.kind {
            SurfaceKind::Plane { .. } => {
                // Planes contribute constraints through axis_constraint(), not bounding_box()
                None
            }
            SurfaceKind::Sphere { x0, y0, z0, radius } => {
                if halfspace_below {
                    // Inside sphere: bounded
                    Some((
                        [*x0 - *radius, *y0 - *radius, *z0 - *radius],
                        [*x0 + *radius, *y0 + *radius, *z0 + *radius],
                    ))
                } else {
                    // Outside sphere: infinite
                    None
                }
            }
            SurfaceKind::Cylinder {
                axis,
                origin,
                radius,
            } => {
                // For now, only handle Z-cylinders (axis = [0, 0, 1])
                if axis[0].abs() < 1e-10 && axis[1].abs() < 1e-10 && (axis[2] - 1.0).abs() < 1e-10 {
                    if halfspace_below {
                        // Inside cylinder: bounded in X,Y, infinite in Z
                        Some((
                            [origin[0] - *radius, origin[1] - *radius, f64::NEG_INFINITY],
                            [origin[0] + *radius, origin[1] + *radius, f64::INFINITY],
                        ))
                    } else {
                        // Outside cylinder: infinite
                        None
                    }
                } else {
                    // For non-Z-aligned cylinders, we'd need more complex bounding box calculation
                    // For now, return None (infinite bounds)
                    None
                }
            }
        }
    }

    /// Get the constraint this surface imposes on axis-aligned bounds when used as a halfspace.
    /// Returns (axis_index, is_upper_bound, value) or None if no axis constraint.
    pub fn axis_constraint(&self, halfspace_above: bool) -> Option<(usize, bool, f64)> {
        match &self.kind {
            SurfaceKind::Plane { a, b, c, d } => {
                if *a == 1.0 && *b == 0.0 && *c == 0.0 {
                    // X plane: ax + by + cz - d = x - d = 0, so x = d
                    // Above (x - d > 0, x > d) gives lower bound at d
                    // Below (x - d < 0, x < d) gives upper bound at d
                    Some((0, !halfspace_above, *d))
                } else if *a == 0.0 && *b == 1.0 && *c == 0.0 {
                    // Y plane: y - d = 0, so y = d
                    Some((1, !halfspace_above, *d))
                } else if *a == 0.0 && *b == 0.0 && *c == 1.0 {
                    // Z plane: z - d = 0, so z = d
                    Some((2, !halfspace_above, *d))
                } else {
                    None
                }
            }
            _ => None,
        }
    }
}

#[derive(Clone)]
pub struct Halfspace {
    pub expr: RegionExpr,
}

impl Halfspace {
    pub fn new_above(surface: Arc<Surface>) -> Self {
        Halfspace {
            expr: RegionExpr::Halfspace(HalfspaceType::Above(surface)),
        }
    }

    pub fn new_below(surface: Arc<Surface>) -> Self {
        Halfspace {
            expr: RegionExpr::Halfspace(HalfspaceType::Below(surface)),
        }
    }
}

#[cfg(test)]
mod tests {
    use crate::region::{HalfspaceType, Region};
    #[test]
    fn test_zcylinder_with_zplanes_bounding_box() {
        // Z-cylinder centered at (1, 2) with radius 3
        let cyl = Surface::z_cylinder(1.0, 2.0, 3.0, 1, None);
        // Z planes at z = -5 and z = 5
        let z_bottom = Surface::z_plane(-5.0, 2, None);
        let z_top = Surface::z_plane(5.0, 3, None);

        // Region: inside cylinder and between Z planes
        let region = Region::new_from_halfspace(HalfspaceType::Below(Arc::new(cyl)))
            .intersection(&Region::new_from_halfspace(HalfspaceType::Above(Arc::new(
                z_bottom,
            ))))
            .intersection(&Region::new_from_halfspace(HalfspaceType::Below(Arc::new(
                z_top,
            ))));
        let bbox = region.bounding_box();
        assert_eq!(bbox.lower_left, [-2.0, -1.0, -5.0]);
        assert_eq!(bbox.upper_right, [4.0, 5.0, 5.0]);
    }
    use super::*;

    #[test]
    fn test_plane_creation() {
        let plane = Surface::new_plane(1.0, 0.0, 0.0, 2.0, 42, None);
        match plane.kind {
            SurfaceKind::Plane { a, b, c, d } => {
                assert_eq!(a, 1.0);
                assert_eq!(b, 0.0);
                assert_eq!(c, 0.0);
                assert_eq!(d, 2.0);
            }
            _ => panic!("Not a plane"),
        }
        assert_eq!(plane.surface_id, 42);
    }

    #[test]
    fn test_sphere_creation() {
        let sphere = Surface::new_sphere(1.0, 2.0, 3.0, 5.0, 7, None);
        match sphere.kind {
            SurfaceKind::Sphere { x0, y0, z0, radius } => {
                assert_eq!(x0, 1.0);
                assert_eq!(y0, 2.0);
                assert_eq!(z0, 3.0);
                assert_eq!(radius, 5.0);
            }
            _ => panic!("Not a sphere"),
        }
        assert_eq!(sphere.surface_id, 7);
    }

    #[test]
    fn test_cylinder_creation() {
        let axis = [0.0, 1.0, 0.0];
        let origin = [1.0, 2.0, 3.0];
        let cylinder = Surface::new_cylinder(axis, origin, 2.0, 99, None);
        match cylinder.kind {
            SurfaceKind::Cylinder {
                axis: a,
                origin: o,
                radius,
            } => {
                assert_eq!(a, axis);
                assert_eq!(o, origin);
                assert_eq!(radius, 2.0);
            }
            _ => panic!("Not a cylinder"),
        }
        assert_eq!(cylinder.surface_id, 99);
    }

    #[test]
    fn test_z_cylinder_creation() {
        let zcyl = Surface::z_cylinder(1.0, 2.0, 3.0, 123, None);
        match zcyl.kind {
            SurfaceKind::Cylinder {
                axis,
                origin,
                radius,
            } => {
                assert_eq!(axis, [0.0, 0.0, 1.0]);
                assert_eq!(origin, [1.0, 2.0, 0.0]);
                assert_eq!(radius, 3.0);
            }
            _ => panic!("Not a Z cylinder"),
        }
        assert_eq!(zcyl.surface_id, 123);
    }

    #[test]
    fn test_boundary_type_default() {
        let plane = Surface::new_plane(1.0, 0.0, 0.0, 2.0, 42, None);
        assert_eq!(*plane.boundary_type(), BoundaryType::Vacuum);
    }

    #[test]
    fn test_boundary_type_vacuum() {
        let sphere = Surface::new_sphere(0.0, 0.0, 0.0, 1.0, 1, Some(BoundaryType::Vacuum));
        assert_eq!(*sphere.boundary_type(), BoundaryType::Vacuum);
    }

    #[test]
    fn test_set_boundary_type() {
        let mut cylinder = Surface::new_cylinder([0.0, 0.0, 1.0], [0.0, 0.0, 0.0], 1.0, 2, None);
        assert_eq!(*cylinder.boundary_type(), BoundaryType::Vacuum);

        cylinder.set_boundary_type(BoundaryType::Transmission);
        assert_eq!(*cylinder.boundary_type(), BoundaryType::Transmission);

        cylinder.set_boundary_type(BoundaryType::Vacuum);
        assert_eq!(*cylinder.boundary_type(), BoundaryType::Vacuum);
    }

    #[test]
    fn test_zcylinder_bounding_box() {
        let zcyl = Surface::z_cylinder(1.0, 2.0, 3.0, 123, None);

        // Test inside cylinder (halfspace_below = true)
        let bbox = zcyl.bounding_box(true);
        assert!(bbox.is_some());
        let (lower, upper) = bbox.unwrap();
        assert_eq!(lower, [-2.0, -1.0, f64::NEG_INFINITY]);
        assert_eq!(upper, [4.0, 5.0, f64::INFINITY]);

        // Test outside cylinder (halfspace_below = false)
        let bbox_outside = zcyl.bounding_box(false);
        assert!(bbox_outside.is_none()); // Outside cylinder is infinite
    }
}
