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
