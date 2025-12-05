use yamc::inelastic::*;
use yamc::particle::Particle;
use yamc::Reaction;
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
    
    // Test the inelastic scatter function
    let outgoing_particles = inelastic_scatter(&particle, &reaction, "Be9", &mut rng);
    
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
#[should_panic(expected = "Unsupported inelastic reaction MT number: 102")]
fn test_unsupported_mt_number_panics() {
    use yamc::inelastic::get_inelastic_neutron_multiplicity;
    
    // MT 102 (n,gamma) is not an inelastic reaction, should panic
    get_inelastic_neutron_multiplicity(102);
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
    
    // Test the inelastic scatter function  
    let outgoing_particles = inelastic_scatter(&particle, &reaction, "Li6", &mut rng);
    
    // For MT 51 (inelastic), we should get 1 neutron back
    assert_eq!(outgoing_particles.len(), 1, "MT 51 (inelastic) should produce 1 neutron");
    
    // Check that particle has reasonable energy (less than incoming due to excitation)
    let outgoing_particle = &outgoing_particles[0];
    assert!(outgoing_particle.energy > 0.0, "Outgoing neutron should have positive energy");
    assert!(outgoing_particle.energy < particle.energy, "Outgoing neutron energy should be less than incoming");
    println!("Inelastic scatter: {} eV -> {} eV", particle.energy, outgoing_particle.energy);
}

#[test]
fn test_inelastic_neutron_multiplicity() {
    use yamc::inelastic::get_inelastic_neutron_multiplicity;
    
    // Test inelastic reaction MT numbers only
    assert_eq!(get_inelastic_neutron_multiplicity(4), 1, "(n,n') should have 1 neutron");
    assert_eq!(get_inelastic_neutron_multiplicity(16), 2, "(n,2n) should have 2 neutrons");
    assert_eq!(get_inelastic_neutron_multiplicity(17), 3, "(n,3n) should have 3 neutrons");
    assert_eq!(get_inelastic_neutron_multiplicity(37), 4, "(n,4n) should have 4 neutrons");
    assert_eq!(get_inelastic_neutron_multiplicity(51), 1, "Inelastic level 1 should have 1 neutron");
    assert_eq!(get_inelastic_neutron_multiplicity(60), 1, "Inelastic level 10 should have 1 neutron");
    assert_eq!(get_inelastic_neutron_multiplicity(90), 1, "Inelastic level 40 should have 1 neutron");
    assert_eq!(get_inelastic_neutron_multiplicity(91), 1, "Continuum inelastic should have 1 neutron");
    assert_eq!(get_inelastic_neutron_multiplicity(22), 1, "(n,n'alpha) should have 1 neutron");
    assert_eq!(get_inelastic_neutron_multiplicity(28), 1, "(n,n'p) should have 1 neutron");
}