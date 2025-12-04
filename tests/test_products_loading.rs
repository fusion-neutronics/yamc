use yamc::nuclide::read_nuclide_from_json;

#[test]
fn test_pb208_products_are_loaded() {
    // Load Pb208 data (no temperature filtering)
    let nuclide = read_nuclide_from_json("tests/Pb208.json", None).unwrap();
    
    // Check that we have reactions
    assert!(!nuclide.reactions.is_empty());
    
    // Print available temperatures
    println!("Available temperatures: {:?}", nuclide.reactions.keys().collect::<Vec<_>>());
    println!("Loaded temperatures: {:?}", nuclide.loaded_temperatures);
    println!("Available temperatures: {:?}", nuclide.available_temperatures);
    
    // Check specific temperature
    let temp = if nuclide.reactions.contains_key("294K") {
        "294K"
    } else if nuclide.reactions.contains_key("294") {
        "294"
    } else {
        // Use the first available temperature
        nuclide.reactions.keys().next().unwrap()
    };
    println!("Using temperature: {}", temp);
    
    let reactions_at_temp = &nuclide.reactions[temp];
    
    // Check specific reaction (MT=2 should have products in Pb208)
    if let Some(reaction) = reactions_at_temp.get(&2) {
        println!("MT=2 reaction has {} products", reaction.products.len());
        // Should have products for elastic scattering
        assert!(!reaction.products.is_empty());
        
        // Print product details
        for (i, product) in reaction.products.iter().enumerate() {
            println!("Product {}: {:?}", i, product.particle);
            if !product.distribution.is_empty() {
                match &product.distribution[0] {
                    yamc::reaction_product::AngleEnergyDistribution::UncorrelatedAngleEnergy { angle, energy } => {
                        println!("  Angle distribution has {} energy points", angle.energy.len());
                        match energy {
                            Some(yamc::reaction_product::EnergyDistribution::LevelInelastic { .. }) => {
                                println!("  Energy distribution: LevelInelastic");
                            }
                            Some(yamc::reaction_product::EnergyDistribution::Tabulated { energy, .. }) => {
                                println!("  Energy distribution: Tabulated with {} energy points", energy.len());
                            }
                            Some(yamc::reaction_product::EnergyDistribution::ContinuousTabular { energy, .. }) => {
                                println!("  Energy distribution: ContinuousTabular with {} energy points", energy.len());
                            }
                            None => {
                                println!("  Energy distribution: None");
                            }
                        }
                    }
                    yamc::reaction_product::AngleEnergyDistribution::KalbachMann { energy, .. } => {
                        println!("  KalbachMann distribution has {} energy points", energy.len());
                    }
                    yamc::reaction_product::AngleEnergyDistribution::CorrelatedAngleEnergy { energy, energy_out, .. } => {
                        println!("  CorrelatedAngleEnergy distribution has {} incoming energy points", energy.len());
                        println!("  With {} energy_out grids", energy_out.len());
                    }
                }
            }
        }
    } else {
        panic!("MT=2 reaction not found");
    }
    
    // Check another reaction that should have products
    println!("Available reactions: {:?}", reactions_at_temp.keys().collect::<Vec<_>>());
    
    // Count total products across all reactions
    let mut total_products = 0;
    let mut reactions_with_products = 0;
    
    for (mt, reaction) in reactions_at_temp {
        if !reaction.products.is_empty() {
            println!("MT={} has {} products", mt, reaction.products.len());
            reactions_with_products += 1;
            total_products += reaction.products.len();
        }
    }
    
    println!("Total reactions with products: {}", reactions_with_products);
    println!("Total products across all reactions: {}", total_products);
    assert!(total_products > 0, "Expected at least some products to be loaded");
}