use pyo3::prelude::*;

#[pyclass(name = "LevelInelastic")]
pub struct PyLevelInelastic {}

#[pymethods]
impl PyLevelInelastic {
    #[new]
    fn new() -> Self {
        PyLevelInelastic {}
    }

    fn __repr__(&self) -> String {
        "LevelInelastic()".to_string()
    }
}

#[pyclass(name = "TabulatedEnergyDistribution")]
pub struct PyTabulatedEnergyDistribution {
    #[pyo3(get)]
    pub energy: Vec<f64>,
    #[pyo3(get)]
    pub energy_out: Vec<Vec<f64>>,
}

#[pymethods]
impl PyTabulatedEnergyDistribution {
    #[new]
    fn new(energy: Vec<f64>, energy_out: Vec<Vec<f64>>) -> Self {
        PyTabulatedEnergyDistribution { energy, energy_out }
    }

    fn __repr__(&self) -> String {
        format!("TabulatedEnergyDistribution(energy=[{} points], energy_out=[{} distributions])", 
                self.energy.len(), self.energy_out.len())
    }
}
