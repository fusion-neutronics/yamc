#[cfg(feature = "pyo3")]
use pyo3::prelude::*;
use crate::reaction::Reaction;
use crate::python::reaction_product_python::PyReactionProduct;

#[cfg(feature = "pyo3")]
#[pyclass(name = "Reaction")]
/// Lightweight reaction summary exposed to Python.
///
/// Represents a single nuclear reaction channel with reactants, products and an
/// associated energy (e.g. Q-value or representative energy in MeV).
///
/// Attributes:
///     reactants (List[str]): Names / symbols of incoming particles or nuclei.
///     products (List[str]): Names / symbols of outgoing particles or nuclei.
///     energy (float): Reaction energy value (units depend on source data, typically MeV).
#[derive(Clone, Debug)]
pub struct PyReaction {
    #[pyo3(get)]
    pub reactants: Vec<String>,
    #[pyo3(get)]
    pub products: Vec<PyReactionProduct>,
    #[pyo3(get)]
    pub energy: f64,
    #[pyo3(get)]
    pub cross_section: Vec<f64>,
    #[pyo3(get)]
    pub threshold_idx: usize,
    #[pyo3(get)]
    pub interpolation: Vec<i32>,
    #[pyo3(get)]
    pub mt_number: i32,
    #[pyo3(get)]
    pub energy_grid: Vec<f64>,
}

#[cfg(feature = "pyo3")]
#[pymethods]
impl PyReaction {
    #[new]
    pub fn new(
        reactants: Vec<String>,
        products: Vec<PyReactionProduct>,
        energy: f64,
        cross_section: Vec<f64>,
        threshold_idx: usize,
        interpolation: Vec<i32>,
        mt_number: i32,
        energy_grid: Vec<f64>,
    ) -> Self {
        PyReaction {
            reactants,
            products,
            energy,
            cross_section,
            threshold_idx,
            interpolation,
            mt_number,
            energy_grid,
        }
    }
}

#[cfg(feature = "pyo3")]
impl PyReaction {
    pub fn from_reaction(reaction: &Reaction, py: pyo3::Python) -> pyo3::PyResult<Self> {
        let reactants = vec![]; // Fill as needed
        let energy = 0.0; // Fill as needed
        let products = reaction.products.iter().map(|prod| PyReactionProduct::from_reaction_product(prod.clone(), py)).collect::<pyo3::PyResult<Vec<_>>>()?;
        Ok(PyReaction {
            reactants,
            products,
            energy,
            cross_section: reaction.cross_section.clone(),
            threshold_idx: reaction.threshold_idx,
            interpolation: reaction.interpolation.clone(),
            mt_number: reaction.mt_number,
            energy_grid: reaction.energy.clone(),
        })
    }
}
