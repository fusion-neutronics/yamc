use crate::material::Material;

/// Material filter for tallies - filters events based on which material they occur in
#[derive(Debug, Clone, PartialEq, Eq, Hash)]
pub struct MaterialFilter {
    /// The material ID to filter on (always Some since None materials aren't real materials)
    pub material_id: u32,
}

impl MaterialFilter {
    /// Create a MaterialFilter from a Material object
    ///
    /// # Arguments
    /// * `material` - The material object to create the filter from (must have a material_id)
    ///
    /// # Returns
    /// A new MaterialFilter that will match events in this material
    ///
    /// # Panics
    /// Panics if the material has no material_id (None materials can't be filtered)
    pub fn new(material: &Material) -> Self {
        let material_id = material.material_id.expect(
            "Cannot create MaterialFilter for material with no ID - assign a material_id first",
        );

        Self { material_id }
    }

    /// Check if this filter matches a given material ID
    ///
    /// # Arguments
    /// * `material_id` - The material ID to check against
    ///
    /// # Returns
    /// `true` if the material ID matches this filter
    pub fn matches(&self, material_id: Option<u32>) -> bool {
        material_id == Some(self.material_id)
    }

    /// Check if this filter matches a given material object
    ///
    /// # Arguments  
    /// * `material` - The material object to check against
    ///
    /// # Returns
    /// `true` if the material matches this filter
    pub fn matches_material(&self, material: &Material) -> bool {
        material.material_id == Some(self.material_id)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_material_filter_creation() {
        let mut material = Material::with_id(123);
        material.set_name("test_material");
        let filter = MaterialFilter::new(&material);
        assert_eq!(filter.material_id, 123);
    }

    #[test]
    #[should_panic(expected = "Cannot create MaterialFilter for material with no ID")]
    fn test_material_filter_creation_no_id_panics() {
        let material = Material::new();
        MaterialFilter::new(&material);
    }

    #[test]
    fn test_material_filter_matching() {
        let mut material = Material::with_id(123);
        material.set_name("test_material");
        let filter = MaterialFilter::new(&material);

        assert!(filter.matches(Some(123)));
        assert!(!filter.matches(Some(124)));
        assert!(!filter.matches(None));
        assert!(filter.matches_material(&material));

        let mut other_material = Material::with_id(456);
        other_material.set_name("other_material");
        assert!(!filter.matches_material(&other_material));
    }

    #[test]
    fn test_material_filter_equality() {
        let mut material1 = Material::with_id(123);
        material1.set_name("material_1");
        let mut material2 = Material::with_id(123);
        material2.set_name("material_2");

        let filter1 = MaterialFilter::new(&material1);
        let filter2 = MaterialFilter::new(&material2);

        assert_eq!(filter1, filter2);

        let mut material3 = Material::with_id(456);
        material3.set_name("material_3");
        let filter3 = MaterialFilter::new(&material3);

        assert_ne!(filter1, filter3);
    }
}
