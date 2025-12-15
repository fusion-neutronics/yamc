use crate::tallies::{CellFilter, EnergyFilter, MaterialFilter};

/// Unified filter enum for tallies
#[derive(Debug, Clone, PartialEq)]
pub enum Filter {
    Cell(CellFilter),
    Material(MaterialFilter),
    Energy(EnergyFilter),
}

impl Filter {
    /// Get the type name of this filter for validation
    pub fn type_name(&self) -> &'static str {
        match self {
            Filter::Cell(_) => "CellFilter",
            Filter::Material(_) => "MaterialFilter",
            Filter::Energy(_) => "EnergyFilter",
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::tallies::{CellFilter, MaterialFilter};

    #[test]
    fn test_type_name_cell() {
        let cell_filter = CellFilter { cell_id: 42 };
        let filter = Filter::Cell(cell_filter);
        assert_eq!(filter.type_name(), "CellFilter");
    }

    #[test]
    fn test_type_name_material() {
        let material_filter = MaterialFilter { material_id: 99 };
        let filter = Filter::Material(material_filter);
        assert_eq!(filter.type_name(), "MaterialFilter");
    }
}
