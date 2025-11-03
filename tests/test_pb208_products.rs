#[cfg(test)]
mod test_pb208_products {
    use materials_for_mc::nuclide::read_nuclide_from_json;
    use std::path::Path;

    #[test]
    fn test_pb208_with_products() {
        let path = Path::new("tests/Pb208.json");
        let nuclide = read_nuclide_from_json(path, None).expect("Failed to load Pb208.json");
        
        // Check that we have reactions
        assert!(!nuclide.reactions.is_empty());
        
        // Check that some reactions have products
        let mut found_products = false;
        for temp_reactions in nuclide.reactions.values() {
            for reaction in temp_reactions.values() {
                if !reaction.products.is_empty() {
                    found_products = true;
                    println!("Reaction MT {} has {} products", 
                             reaction.mt_number, reaction.products.len());
                    
                    for (i, product) in reaction.products.iter().enumerate() {
                        println!("  Product {}: particle={:?}, emission_mode={}, distributions={}",
                                 i, product.particle, product.emission_mode, product.distribution.len());
                        
                        for (j, dist) in product.distribution.iter().enumerate() {
                            println!("    Distribution {}: {:?}", j, dist);
                        }
                    }
                    break;
                }
            }
            if found_products { break; }
        }
        
        assert!(found_products, "No products found in any reactions");
    }
}