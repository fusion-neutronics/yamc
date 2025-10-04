use crate::geometry::Geometry;
use crate::materials::Materials;
use crate::source::Source;
use crate::settings::Settings;

#[derive(Debug, Clone)]
pub struct Model {
    pub geometry: Geometry,
    pub materials: Materials,
    pub source: Source,
    pub settings: Settings,
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::geometry::Geometry;
    use crate::materials::Materials;
    use crate::source::Source;
    use crate::settings::Settings;

    use crate::cell::Cell;
    use crate::region::{Region, HalfspaceType};
    use crate::surface::{Surface, SurfaceKind, BoundaryType};
    use crate::material::Material;
    use std::sync::Arc;

    #[test]
    fn test_model_construction() {
        // Sphere surface with vacuum boundary
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
        // Material with Li6 nuclide
        let mut material = Material::new();
        material.add_nuclide("Li6", 1.0).unwrap();
        let cell = Cell {
            cell_id: 1,
            name: Some("sphere_cell".to_string()),
            region,
            material: Some(material.clone()),
        };
        let geometry = Geometry { cells: vec![cell] };
        let mut materials = Materials::new();
        materials.append(material);
        let source = Source {
            position: [0.0, 0.0, 0.0],
            direction: [0.0, 0.0, 1.0],
            energy: 1e6,
        };
        let settings = Settings {
            particles: 100,
            batches: 10,
            source: source.clone(),
        };
        let model = Model {
            geometry,
            materials,
            source,
            settings,
        };
        assert_eq!(model.settings.particles, 100);
        assert_eq!(model.source.energy, 1e6);
        // Check geometry and material
        assert_eq!(model.geometry.cells.len(), 1);
        assert!(model.materials.len() > 0);
        assert!(model.materials.get(0).unwrap().nuclides.contains_key("Li6"));
    }
}