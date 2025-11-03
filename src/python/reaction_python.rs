#[cfg(feature = "pyo3")]
use pyo3::prelude::*;
use crate::reaction::Reaction;
use crate::python::reaction_product_python::PyReactionProduct;

#[cfg(feature = "pyo3")]
#[pyclass(name = "Reaction")]
/// Lightweight reaction summary exposed to Python.
///
/// Represents a single nuclear reaction channel with its products and cross section data.
///
/// Attributes:
///     products (List[ReactionProduct]): Outgoing particles with their distributions.
///     cross_section (List[float]): Cross section values in barns.
///     threshold_idx (int): Index into parent energy grid where reaction becomes active.
///     interpolation (List[int]): Interpolation flags.
///     mt_number (int): ENDF/MT reaction identifier.
///     energy_grid (List[float]): Reaction-specific energy grid in eV.
#[derive(Clone, Debug)]
pub struct PyReaction {
    #[pyo3(get)]
    pub products: Vec<PyReactionProduct>,
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
        products: Vec<PyReactionProduct>,
        cross_section: Vec<f64>,
        threshold_idx: usize,
        interpolation: Vec<i32>,
        mt_number: i32,
        energy_grid: Vec<f64>,
    ) -> Self {
        PyReaction {
            products,
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
        let products = reaction.products.iter().map(|prod| PyReactionProduct::from_reaction_product(prod.clone(), py)).collect::<pyo3::PyResult<Vec<_>>>()?;
        Ok(PyReaction {
            products,
            cross_section: reaction.cross_section.clone(),
            threshold_idx: reaction.threshold_idx,
            interpolation: reaction.interpolation.clone(),
            mt_number: reaction.mt_number,
            energy_grid: reaction.energy.clone(),
        })
    }
}
