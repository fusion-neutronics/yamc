pub mod bounding_box;
pub mod cell;
pub mod geometry;
pub mod particle;
pub mod source;
pub mod settings;
pub mod region;
pub mod surface;
pub mod model;
pub use bounding_box::*;
pub use cell::*;
pub use geometry::*;
pub use region::*;
pub use surface::*;
#[cfg(feature = "wasm")]
pub use wasm::material_wasm;
#[cfg(feature = "wasm")]
pub use wasm::nuclide_wasm;
#[cfg(feature = "wasm")]
pub use wasm::reaction_wasm;
#[cfg(feature = "wasm")]
pub use wasm::config_wasm;
#[cfg(feature = "wasm")]
pub use wasm::element_wasm;
#[cfg(feature = "wasm")]
pub use wasm::data_wasm;
#[cfg(feature = "wasm")]
pub use wasm::material_wasm::WasmMaterial;
#[cfg(feature = "wasm")]
pub mod wasm {
    pub mod config_wasm;
    pub mod element_wasm;
    pub mod material_wasm;
    pub mod nuclide_wasm;
    pub mod reaction_wasm;
    pub mod data_wasm;
}
mod data;
// First, import any modules and re-export the types for Rust usage
mod config;
mod element;
mod material;
mod materials;
mod nuclide;
mod reaction;
mod utilities;
mod url_cache;

pub use config::Config;
pub use element::Element;
pub use material::Material;
pub use materials::Materials;
pub use reaction::Reaction;
pub use nuclide::Nuclide;
pub use utilities::{interpolate_linear, interpolate_log_log};

// Import PyO3 items conditionally
#[cfg(feature = "pyo3")]
use pyo3::prelude::*;
#[cfg(feature = "pyo3")]
use pyo3::pymodule;
#[cfg(feature = "pyo3")]
use pyo3::wrap_pyfunction;

// Conditionally include the Python modules

// Re-export Python modules for Maturin to find

#[cfg(feature = "pyo3")]
mod python {
    pub mod config_python;
    pub mod element_python;
    pub mod material_python;
    pub mod materials_python;
    pub mod nuclide_python;
    pub mod reaction_python;
    pub mod data_python;
    pub mod bounding_box_python;
    pub mod cell_python;
    pub mod geometry_python;
    pub mod region_python;
    pub mod surface_python;
    pub mod source_python;
    pub mod settings_python;
    pub mod model_python;
    pub mod particle_python;
    pub use config_python::*;
    pub use element_python::*;
    pub use material_python::*;
    pub use materials_python::*;
    pub use nuclide_python::*;
    pub use reaction_python::*;
    pub use data_python::*;
    pub use bounding_box_python::*;
    pub use cell_python::*;
    pub use geometry_python::*;
    pub use region_python::*;
    pub use surface_python::*;
}
#[cfg(feature = "wasm")]
use wasm_bindgen::prelude::*;

#[cfg(feature = "wasm")]
#[wasm_bindgen(start)]
pub fn wasm_start() {
    console_error_panic_hook::set_once();
}

// Export WASM modules from the wasm submodule
#[cfg(feature = "wasm")]
pub use wasm::config_wasm::*;
#[cfg(feature = "wasm")]
pub use wasm::element_wasm::*;
#[cfg(feature = "wasm")]
pub use wasm::material_wasm::*;
#[cfg(feature = "wasm")]
pub use wasm::nuclide_wasm::*;
#[cfg(feature = "wasm")]
pub use wasm::reaction_wasm::*;
#[cfg(feature = "wasm")]
pub use wasm::data_wasm::*;

// If you have a main Python module entry point, update it to include PyMaterials:
#[cfg(feature = "pyo3")]
#[pymodule]
fn materials_for_mc(_py: Python<'_>, m: &PyModule) -> PyResult<()> {
    use crate::python::particle_python;
    m.add_class::<particle_python::PyParticle>()?;
    use crate::python::model_python;
    m.add_class::<model_python::PyModel>()?;
    use crate::python::material_python;
    use crate::python::materials_python;
    use crate::python::nuclide_python;
    use crate::python::reaction_python;
    use crate::python::config_python;
    use crate::python::element_python;
    use crate::python::data_python;
    use crate::python::surface_python;
    use crate::python::cell_python;
    use crate::python::region_python;
    use crate::python::geometry_python;
    use crate::python::bounding_box_python;
    use crate::python::source_python;
    use crate::python::settings_python;

    // Core materials API
    m.add_class::<material_python::PyMaterial>()?;
    m.add_class::<materials_python::PyMaterials>()?;
    m.add_class::<nuclide_python::PyNuclide>()?;
    m.add_class::<reaction_python::PyReaction>()?;
    m.add_class::<config_python::PyConfig>()?;
    m.add_class::<element_python::PyElement>()?;

    // Particle transport API
    m.add_class::<source_python::PySource>()?;
    m.add_class::<settings_python::PySettings>()?;

    // Geometry/CSG API
    m.add_class::<cell_python::PyCell>()?;
    m.add_class::<region_python::PyRegion>()?;
    m.add_class::<region_python::PyHalfspace>()?;
    m.add_class::<region_python::PyBoundingBox>()?;
    m.add_class::<surface_python::PySurface>()?;
    m.add_class::<surface_python::PyBoundaryType>()?;
    m.add_class::<geometry_python::PyGeometry>()?;

    // Functions for nuclide/data
    m.add_function(wrap_pyfunction!(nuclide_python::py_read_nuclide_from_json, m)?)?;
    m.add_function(wrap_pyfunction!(nuclide_python::clear_nuclide_cache, m)?)?;
    m.add_function(wrap_pyfunction!(data_python::natural_abundance, m)?)?;
    m.add_function(wrap_pyfunction!(data_python::element_nuclides, m)?)?;
    m.add_function(wrap_pyfunction!(data_python::element_names, m)?)?;
    m.add_function(wrap_pyfunction!(data_python::atomic_masses, m)?)?;

    // Surface constructors (OpenMC-style API)
    m.add_function(wrap_pyfunction!(surface_python::XPlane, m)?)?;
    m.add_function(wrap_pyfunction!(surface_python::YPlane, m)?)?;
    m.add_function(wrap_pyfunction!(surface_python::ZPlane, m)?)?;
    m.add_function(wrap_pyfunction!(surface_python::Sphere, m)?)?;
    m.add_function(wrap_pyfunction!(surface_python::Cylinder, m)?)?;
    m.add_function(wrap_pyfunction!(surface_python::ZCylinder, m)?)?;
    m.add_function(wrap_pyfunction!(surface_python::Plane, m)?)?;
    Ok(())
}
