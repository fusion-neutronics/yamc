use materials_for_mc::inelastic::*;
use materials_for_mc::particle::Particle;
use materials_for_mc::Reaction;
use rand::SeedableRng;
use rand::rngs::StdRng;

#[test]
fn test_inelastic_scatter_returns_particles() {
    let mut rng = StdRng::seed_from_u64(42);
    
    // Create a test particle with sufficient energy
    let particle = Particle::new([0.0, 0.0, 0.0], [0.0, 0.0, 1.0], 10e6); // 10 MeV
    
    // Create a test reaction without products (analytical approach)
    let reaction = Reaction {
        cross_section: vec![1.0],
        threshold_idx: 0,
        interpolation: vec![2],
        energy: vec![1e5],
        mt_number: 16, // (n,2n) reaction
        q_value: -6.0e6, // Endothermic, 6 MeV threshold
        products: vec![],
    };
    
    let awr = 9.0; // Be9 atomic weight ratio
    
    // Test the inelastic scatter function
    let outgoing_particles = inelastic_scatter(&particle, &reaction, awr, &mut rng);
    
    // For MT 16 (n,2n), we should get 2 neutrons back
    assert_eq!(outgoing_particles.len(), 2, "MT 16 (n,2n) should produce 2 neutrons");
    
    // Check that particles have reasonable energies
    for (i, p) in outgoing_particles.iter().enumerate() {
        assert!(p.energy > 0.0, "Neutron {} should have positive energy", i);
        assert!(p.energy < particle.energy, "Neutron {} energy should be less than incoming", i);
        println!("Neutron {}: energy = {:.1e} eV", i, p.energy);
    }
}

#[test]
fn test_absorption_reaction() {
    let mut rng = StdRng::seed_from_u64(42);
    
    // Create a test particle
    let particle = Particle::new([0.0, 0.0, 0.0], [0.0, 0.0, 1.0], 1e6); // 1 MeV
    
    // Create an absorption reaction (n,gamma)
    let reaction = Reaction {
        cross_section: vec![1.0],
        threshold_idx: 0,
        interpolation: vec![2],
        energy: vec![1e5],
        mt_number: 102, // (n,gamma) reaction
        q_value: 7.0e6, // Exothermic
        products: vec![],
    };
    
    let awr = 6.0; // Li6 atomic weight ratio
    
    // Test the inelastic scatter function
    let outgoing_particles = inelastic_scatter(&particle, &reaction, awr, &mut rng);
    
    // For MT 102 (n,gamma), neutron is absorbed - should get 0 neutrons back
    assert_eq!(outgoing_particles.len(), 0, "MT 102 (n,gamma) should absorb neutron - no outgoing neutrons");
}

#[test]
fn test_single_neutron_reaction() {
    let mut rng = StdRng::seed_from_u64(42);
    
    // Create a test particle with sufficient energy
    let particle = Particle::new([0.0, 0.0, 0.0], [0.0, 0.0, 1.0], 5e6); // 5 MeV
    
    // Create an inelastic scattering reaction
    let reaction = Reaction {
        cross_section: vec![1.0],
        threshold_idx: 0,
        interpolation: vec![2],
        energy: vec![1e5],
        mt_number: 51, // Discrete inelastic level
        q_value: -1.0e6, // 1 MeV excitation energy
        products: vec![],
    };
    
    let awr = 7.0; // Li7 atomic weight ratio
    
    // Test the inelastic scatter function
    let outgoing_particles = inelastic_scatter(&particle, &reaction, awr, &mut rng);
    
    // For MT 51 (inelastic), we should get 1 neutron back
    assert_eq!(outgoing_particles.len(), 1, "MT 51 (inelastic) should produce 1 neutron");
    
    // Check that particle has reasonable energy (less than incoming due to excitation)
    let outgoing_particle = &outgoing_particles[0];
    assert!(outgoing_particle.energy > 0.0, "Outgoing neutron should have positive energy");
    assert!(outgoing_particle.energy < particle.energy, "Outgoing neutron energy should be less than incoming");
    println!("Inelastic scatter: {} eV -> {} eV", particle.energy, outgoing_particle.energy);
}

#[test]
fn test_neutron_multiplicity() {
    use materials_for_mc::inelastic::get_neutron_multiplicity;
    
    // Test various MT numbers
    assert_eq!(get_neutron_multiplicity(2), 1, "Elastic should have 1 neutron");
    assert_eq!(get_neutron_multiplicity(16), 2, "(n,2n) should have 2 neutrons");
    assert_eq!(get_neutron_multiplicity(17), 3, "(n,3n) should have 3 neutrons");
    assert_eq!(get_neutron_multiplicity(37), 4, "(n,4n) should have 4 neutrons");
    assert_eq!(get_neutron_multiplicity(102), 0, "(n,gamma) should have 0 neutrons");
    assert_eq!(get_neutron_multiplicity(18), 2, "Fission should have 2 neutrons (average)");
    assert_eq!(get_neutron_multiplicity(875), 2, "(n,2n0) should have 2 neutrons");
}