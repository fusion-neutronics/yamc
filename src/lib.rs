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
    pub use config_python::*;
    pub use element_python::*;
    pub use material_python::*;
    pub use materials_python::*;
    pub use nuclide_python::*;
    pub use reaction_python::*;
    pub use data_python::*;
}
#[cfg(feature = "wasm")]
pub mod material_wasm;
#[cfg(feature = "wasm")]
pub mod nuclide_wasm;
#[cfg(feature = "wasm")]
pub mod reaction_wasm;

// WASM setup
#[cfg(feature = "wasm")]
use wasm_bindgen::prelude::*;

#[cfg(feature = "wasm")]
#[wasm_bindgen(start)]
pub fn wasm_start() {
    console_error_panic_hook::set_once();
}

// Export WASM modules
#[cfg(feature = "wasm")]
pub use config_wasm::*;
#[cfg(feature = "wasm")]
pub use element_wasm::*;
#[cfg(feature = "wasm")]
pub use material_wasm::*;
#[cfg(feature = "wasm")]
pub use nuclide_wasm::*;
#[cfg(feature = "wasm")]
pub use reaction_wasm::*;

// If you have a main Python module entry point, update it to include PyMaterials:
#[cfg(feature = "pyo3")]
#[pymodule]
fn materials_for_mc(_py: Python<'_>, m: &PyModule) -> PyResult<()> {
    use crate::python::material_python;
    use crate::python::materials_python;
    use crate::python::nuclide_python;
    use crate::python::reaction_python;
    use crate::python::config_python;
    use crate::python::element_python;
    use crate::python::data_python;

    m.add_class::<material_python::PyMaterial>()?;
    m.add_class::<materials_python::PyMaterials>()?;
    m.add_class::<nuclide_python::PyNuclide>()?;
    m.add_class::<reaction_python::PyReaction>()?; // Exposed as Reaction in Python
    m.add_class::<config_python::PyConfig>()?;
    m.add_class::<element_python::PyElement>()?;
    m.add_function(wrap_pyfunction!(
        nuclide_python::py_read_nuclide_from_json,
        m
    )?)?;
    m.add_function(wrap_pyfunction!(nuclide_python::clear_nuclide_cache, m)?)?;
    m.add_function(wrap_pyfunction!(data_python::natural_abundance, m)?)?;
    m.add_function(wrap_pyfunction!(data_python::element_nuclides, m)?)?;
    m.add_function(wrap_pyfunction!(data_python::element_names, m)?)?;
    m.add_function(wrap_pyfunction!(data_python::atomic_masses, m)?)?;
    Ok(())
}
