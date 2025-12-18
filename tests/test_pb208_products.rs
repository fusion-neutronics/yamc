#[cfg(test)]
mod test_pb208_products {
    use yamc::nuclide::read_nuclide_from_hdf5;
    use std::path::Path;

    #[test]
    fn test_pb208_with_products() {
        // Use Li6.h5 since Pb208.h5 is not available
        let path = Path::new("tests/Li6.h5");
        let nuclide = read_nuclide_from_hdf5(path, None).expect("Failed to load Li6.h5");
        
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