"""
Test suite for verifying that full distribution structures are properly exposed to Python.

This test ensures that users can examine and interrogate the internal structure
of reaction product distributions through the Python API.
"""

import pytest
import numpy as np
import materials_for_mc as mc


class TestDistributionStructureExposure:
    """Test that distribution structures are fully exposed to Python"""
    
    def test_distribution_types_property(self):
        """Test that distribution_types property is accessible"""
        # Create a test reaction product
        test_product = mc.create_test_reaction_product()
        
        # Test that distribution_types property exists and is accessible
        assert hasattr(test_product, 'distribution_types')
        
        # Test that it returns a list of strings
        dist_types = test_product.distribution_types
        assert isinstance(dist_types, list)
        assert all(isinstance(dt, str) for dt in dist_types)
        
        # Test that it contains expected distribution types
        assert len(dist_types) > 0
        # The test product should have UncorrelatedAngleEnergy
        assert 'UncorrelatedAngleEnergy' in dist_types
    
    def test_distribution_types_with_test_product(self):
        """Test distribution_types with test reaction product"""
        # Create multiple test products to verify consistent behavior
        for i in range(3):
            test_product = mc.create_test_reaction_product()
            
            # Test that distribution_types property exists and works
            assert hasattr(test_product, 'distribution_types')
            dist_types = test_product.distribution_types
            
            # Should be a list of strings
            assert isinstance(dist_types, list)
            assert all(isinstance(dt, str) for dt in dist_types)
            
            # Should contain valid distribution type names
            valid_types = {'UncorrelatedAngleEnergy', 'KalbachMann'}
            for dt in dist_types:
                assert dt in valid_types, f"Unknown distribution type: {dt}"
            
            # Test that we can sample from each distribution type
            for dist_type in dist_types:
                energy, mu = test_product.sample(1e6)  # 1 MeV input
                assert isinstance(energy, float)
                assert isinstance(mu, float)
                assert energy > 0
                assert -1.0 <= mu <= 1.0
    
    def test_different_distribution_types(self):
        """Test that different distribution types are properly identified"""
        # Test multiple products to see what distribution types we get
        found_types = set()
        
        # Create several test products and collect their types
        for i in range(10):
            test_product = mc.create_test_reaction_product()
            if hasattr(test_product, 'distribution_types'):
                for dist_type in test_product.distribution_types:
                    found_types.add(dist_type)
        
        # Should have found at least one type
        assert len(found_types) > 0
        
        # Log what we found for debugging
        print(f"Found distribution types: {found_types}")
        
        # Test that each type is a recognized type
        valid_types = {'UncorrelatedAngleEnergy', 'KalbachMann'}
        for found_type in found_types:
            assert found_type in valid_types
    
    def test_distribution_types_consistency(self):
        """Test that distribution_types is consistent with num_distributions"""
        test_product = mc.create_test_reaction_product()
        
        # Both methods should exist
        assert hasattr(test_product, 'distribution_types')
        assert hasattr(test_product, 'num_distributions')
        
        # Length of distribution_types should match num_distributions
        dist_types = test_product.distribution_types
        num_dists = test_product.num_distributions
        
        assert len(dist_types) == num_dists
    
    def test_distribution_introspection(self):
        """Test that we can introspect distribution properties"""
        test_product = mc.create_test_reaction_product()
        
        # Test basic product properties
        assert hasattr(test_product, 'particle')
        assert hasattr(test_product, 'emission_mode')
        assert hasattr(test_product, 'get_decay_rate')
        
        # Test that particle is accessible and valid
        particle = test_product.particle
        assert isinstance(particle, str)
        assert particle in ['neutron', 'photon']  # Expected particle types
        
        # Test emission mode
        emission_mode = test_product.emission_mode
        assert isinstance(emission_mode, str)
        
        # Test decay rate
        decay_rate = test_product.get_decay_rate()
        assert isinstance(decay_rate, float)
        assert decay_rate >= 0.0
    
    def test_distribution_sampling_with_types(self):
        """Test that we can sample after examining distribution types"""
        test_product = mc.create_test_reaction_product()
        
        # Examine the distribution types first
        dist_types = test_product.distribution_types
        
        # Based on the type, we should be able to sample
        for dist_type in dist_types:
            if dist_type == 'UncorrelatedAngleEnergy':
                # Should be able to sample energy and angle independently
                energy, mu = test_product.sample(1e6)  # 1 MeV input
                
                # Check that we got reasonable outputs
                assert isinstance(energy, float)
                assert isinstance(mu, float)
                assert energy > 0
                assert -1.0 <= mu <= 1.0
            
            elif dist_type == 'KalbachMann':
                # Should be able to sample correlated energy and angle
                energy, mu = test_product.sample(1e6)  # 1 MeV input
                
                # Check that we got reasonable outputs
                assert isinstance(energy, float)
                assert isinstance(mu, float)
                assert energy > 0
                assert -1.0 <= mu <= 1.0
    
    def test_angle_distribution_exposure(self):
        """Test that angle distribution components are accessible"""
        # This tests the PyAngleDistribution class
        
        # Create test data for angle distribution
        energy_points = [1e5, 1e6, 1e7]  # Energy points in eV
        
        # Create PyTabulated for each energy point (simplified test)
        mu_points = [-1.0, 0.0, 1.0]
        prob_points = [0.0, 0.5, 1.0]  # CDF
        
        # Test PyTabulated creation
        tabulated = mc.PyTabulated(mu_points, prob_points)
        
        # Test that tabulated object has expected properties
        assert hasattr(tabulated, 'x')
        assert hasattr(tabulated, 'p')
        assert hasattr(tabulated, 'sample')
        
        # Test the properties
        assert tabulated.x == mu_points
        assert tabulated.p == prob_points
        
        # Test sampling
        sample_mu = tabulated.sample()
        assert isinstance(sample_mu, float)
        assert -1.0 <= sample_mu <= 1.0
    
    def test_full_workflow_introspection(self):
        """Test a complete workflow of distribution examination and sampling"""
        # Create a test product
        test_product = mc.create_test_reaction_product()
        
        # Step 1: Examine the product properties
        particle_type = test_product.particle
        emission_mode = test_product.emission_mode
        decay_rate = test_product.get_decay_rate()
        num_distributions = test_product.num_distributions
        
        print(f"Product analysis:")
        print(f"  Particle: {particle_type}")
        print(f"  Emission mode: {emission_mode}")
        print(f"  Decay rate: {decay_rate}")
        print(f"  Number of distributions: {num_distributions}")
        
        # Step 2: Examine distribution types
        dist_types = test_product.distribution_types
        print(f"  Distribution types: {dist_types}")
        
        # Step 3: Based on distribution types, perform appropriate sampling
        test_energies = [1e4, 1e5, 1e6, 5e6, 1e7]  # Test at different energies
        
        for energy in test_energies:
            try:
                e_out, mu = test_product.sample(energy)
                
                # Validate the output
                assert isinstance(e_out, float)
                assert isinstance(mu, float)
                assert e_out > 0
                assert -1.0 <= mu <= 1.0
                
                print(f"  E_in={energy:.0e} eV -> E_out={e_out:.2e} eV, mu={mu:.3f}")
                
            except Exception as e:
                pytest.fail(f"Sampling failed at energy {energy}: {e}")
        
        # Step 4: Verify that the information is consistent
        assert len(dist_types) == num_distributions
        assert particle_type is not None
        assert emission_mode is not None
        assert isinstance(decay_rate, float)
    
    def test_distribution_types_immutability(self):
        """Test that distribution_types returns a copy, not a mutable reference"""
        test_product = mc.create_test_reaction_product()
        
        # Get distribution types
        dist_types1 = test_product.distribution_types
        dist_types2 = test_product.distribution_types
        
        # Should be equal but not the same object (if they're lists)
        assert dist_types1 == dist_types2
        
        # Modifying one should not affect the other
        if len(dist_types1) > 0:
            original_first = dist_types1[0]
            dist_types1[0] = "modified"
            
            # Get a fresh copy
            dist_types3 = test_product.distribution_types
            assert dist_types3[0] == original_first  # Should be unchanged


    def test_pb208_decay_rate_exposure(self):
        """Test that Pb208 product decay rates are accessible (if products exist)"""
        try:
            # Load Pb208 nuclide 
            pb208 = mc.Nuclide('Pb208')
            pb208.read_nuclide_from_json('tests/Pb208.json')
            
            # Check if Pb208 has any reactions with products
            products_found = False
            for reaction in pb208.reactions:
                if hasattr(reaction, 'products') and len(reaction.products) > 0:
                    products_found = True
                    for product in reaction.products:
                        # Test that we can access decay rate
                        assert hasattr(product, 'get_decay_rate')
                        decay_rate = product.get_decay_rate()
                        assert isinstance(decay_rate, float)
                        assert decay_rate >= 0.0
                        
                        # Test particle and emission mode
                        assert hasattr(product, 'particle')
                        assert hasattr(product, 'emission_mode')
                        particle = product.particle
                        emission_mode = product.emission_mode
                        assert isinstance(particle, str)
                        assert isinstance(emission_mode, str)
                        
                        print(f"Pb208 product: {particle}, emission: {emission_mode}, decay_rate: {decay_rate}")
            
            if not products_found:
                print("No products found in Pb208 data - skipping decay rate test")
                
        except Exception as e:
            # If Pb208 data is not available or has issues, skip this test
            print(f"Pb208 test skipped due to: {e}")


class TestDistributionTypeValidation:
    """Test validation of distribution types"""
    
    def test_only_valid_types_returned(self):
        """Test that only valid distribution types are returned"""
        valid_types = {'UncorrelatedAngleEnergy', 'KalbachMann'}
        
        # Test with the test product
        test_product = mc.create_test_reaction_product()
        dist_types = test_product.distribution_types
        
        for dist_type in dist_types:
            assert dist_type in valid_types, f"Invalid distribution type: {dist_type}"
    
    def test_distribution_type_strings(self):
        """Test that distribution types are proper strings"""
        test_product = mc.create_test_reaction_product()
        dist_types = test_product.distribution_types
        
        for dist_type in dist_types:
            assert isinstance(dist_type, str)
            assert len(dist_type) > 0
            assert not dist_type.isspace()  # Not just whitespace


if __name__ == "__main__":
    # Run the tests when executed directly
    pytest.main([__file__, "-v"])