"""
Test suite for reaction product sampling functionality.

Tests both the Rust implementation via PyO3 bindings and
the Python interface for nuclear reaction product sampling.
"""

import pytest
import numpy as np
import materials_for_mc as mc


class TestTabulatedSampling:
    """Test tabulated distribution sampling"""
    
    def test_tabulated_creation(self):
        """Test creating a tabulated distribution"""
        x = [-1.0, 0.0, 1.0]
        p = [0.2, 0.6, 1.0]  # CDF values
        
        tab = mc.PyTabulated(x, p)
        
        assert tab.x == x
        assert tab.p == p
    
    def test_tabulated_sampling(self):
        """Test sampling from tabulated distribution"""
        # Create a simple distribution - uniform over [-1, 1]
        x = [-1.0, 0.0, 1.0]
        p = [0.0, 0.5, 1.0]  # CDF for uniform distribution
        
        tab = mc.PyTabulated(x, p)
        
        # Sample multiple times
        samples = [tab.sample() for _ in range(1000)]
        
        # Check all samples are within bounds
        assert all(-1.0 <= s <= 1.0 for s in samples)
        
        # Check rough statistics (probabilistic test)
        mean = np.mean(samples)
        assert abs(mean) < 0.2  # Should be roughly centered for uniform distribution
    
    def test_cdf_conversion(self):
        """Test PDF to CDF conversion"""
        # PDF values
        x = [-1.0, -0.5, 0.0, 0.5, 1.0]
        p = [0.1, 0.2, 0.4, 0.2, 0.1]  # PDF
        
        pdf = mc.PyTabulated(x, p)
        cdf = pdf.to_cdf()
        
        # CDF should be normalized to 1.0
        assert abs(cdf.p[-1] - 1.0) < 1e-10
        
        # CDF should be monotonically increasing
        for i in range(1, len(cdf.p)):
            assert cdf.p[i] >= cdf.p[i-1]


class TestIsotropicSampling:
    """Test isotropic scattering utilities"""
    
    def test_sample_isotropic(self):
        """Test isotropic mu sampling"""
        samples = [mc.sample_isotropic() for _ in range(1000)]
        
        # All samples should be valid cosines
        assert all(-1.0 <= s <= 1.0 for s in samples)
        
        # Should be roughly uniform (mean ~ 0)
        mean = np.mean(samples)
        assert abs(mean) < 0.1


class TestSampleTabulated:
    """Test convenience function for tabulated sampling"""
    
    def test_sample_tabulated_function(self):
        """Test the sample_tabulated convenience function"""
        x = [-1.0, 0.0, 1.0]
        p = [0.25, 0.5, 1.0]  # CDF
        
        samples = [mc.sample_tabulated(x, p) for _ in range(1000)]
        
        # Check bounds
        assert all(-1.0 <= s <= 1.0 for s in samples)
        
        # Check distribution roughly
        mean = np.mean(samples)
        assert abs(mean) < 0.2


class TestReactionProduct:
    """Test ReactionProduct sampling"""
    
    def test_create_test_reaction_product(self):
        """Test creating a test reaction product"""
        product = mc.create_test_reaction_product()
        
        assert product.particle == "neutron"
        assert product.emission_mode == "prompt"
        assert product.get_decay_rate() == 0.0
        assert product.num_distributions > 0
    
    def test_reaction_product_sampling(self):
        """Test sampling from reaction product"""
        product = mc.create_test_reaction_product()
        
        incoming_energy = 14e6  # 14 MeV
        
        # Sample multiple times
        results = [product.sample(incoming_energy) for _ in range(100)]
        
        for e_out, mu in results:
            # Energy should be reasonable (elastic scattering)
            assert e_out == incoming_energy  # No energy change
            
            # Mu should be valid cosine
            assert -1.0 <= mu <= 1.0
    
    def test_reaction_product_properties(self):
        """Test reaction product property methods"""
        product = mc.create_test_reaction_product()
        
        # Test particle type checking
        assert product.is_particle_type("neutron")
        assert not product.is_particle_type("photon")
        
        # Test emission mode
        assert product.is_prompt()
        assert not product.is_delayed()
    
    def test_multiple_sampling(self):
        """Test sampling multiple particles from one product"""
        product = mc.create_test_reaction_product()
        
        results = product.sample_multiple(14e6)
        
        # Should get at least one result
        assert len(results) >= 1
        
        for e_out, mu in results:
            assert e_out > 0
            assert -1.0 <= mu <= 1.0


class TestSamplingStatistics:
    """Test statistical properties of sampling"""
    
    def test_angular_distribution_statistics(self):
        """Test that angular distributions have correct statistical properties"""
        product = mc.create_test_reaction_product()
        
        # Sample many mu values
        n_samples = 10000
        mu_samples = [product.sample(14e6)[1] for _ in range(n_samples)]
        
        # Check mean and standard deviation are reasonable for the test distribution
        mean_mu = np.mean(mu_samples)
        std_mu = np.std(mu_samples)
        
        # The test distribution should be roughly centered and not too narrow
        assert abs(mean_mu) < 0.1  # Roughly isotropic
        assert std_mu > 0.3       # Not too peaked
    
    def test_energy_conservation(self):
        """Test energy conservation for elastic scattering"""
        product = mc.create_test_reaction_product()
        
        energies = [1e5, 1e6, 14e6, 1e8]  # Various incoming energies
        
        for incoming_energy in energies:
            e_out, mu = product.sample(incoming_energy)
            
            # For the test product (elastic), energy should be conserved
            assert abs(e_out - incoming_energy) < 1e-10
    
    def test_sampling_reproducibility_with_seed(self):
        """Test that sampling is reproducible when using the same conditions"""
        # Note: Python sampling uses thread_rng internally, so we can't easily seed it
        # This test just checks that sampling produces different results on repeated calls
        product = mc.create_test_reaction_product()
        
        results1 = [product.sample(14e6) for _ in range(10)]
        results2 = [product.sample(14e6) for _ in range(10)]
        
        # Results should be different (very unlikely to be identical)
        assert results1 != results2


class TestErrorHandling:
    """Test error handling and edge cases"""
    
    def test_empty_tabulated_distribution(self):
        """Test handling of empty distributions"""
        # Empty arrays should not crash
        empty_tab = mc.PyTabulated([], [])
        
        # Should return 0.0 for empty distribution
        result = empty_tab.sample()
        assert result == 0.0
    
    def test_single_point_distribution(self):
        """Test single-point distribution"""
        single_tab = mc.PyTabulated([0.5], [1.0])
        
        # Should always return the single value
        for _ in range(10):
            assert single_tab.sample() == 0.5


class TestPhysicsValidation:
    """Test that sampling results are physically reasonable"""
    
    def test_mu_bounds(self):
        """Test that all mu values are valid cosines"""
        product = mc.create_test_reaction_product()
        
        # Test at various energies
        energies = np.logspace(4, 8, 20)  # 10 keV to 100 MeV
        
        for energy in energies:
            for _ in range(10):  # Multiple samples per energy
                e_out, mu = product.sample(energy)
                
                # Mu must be a valid cosine
                assert -1.0 <= mu <= 1.0, f"Invalid mu={mu} at E={energy}"
    
    def test_energy_positivity(self):
        """Test that outgoing energies are positive"""
        product = mc.create_test_reaction_product()
        
        energies = [1e3, 1e6, 14e6, 1e8]
        
        for energy in energies:
            for _ in range(10):
                e_out, mu = product.sample(energy)
                assert e_out > 0, f"Non-positive energy {e_out} at E={energy}"


class TestIntegration:
    """Integration tests with real nuclear data (if available)"""
    
    def test_with_nuclide_data(self):
        """Test sampling integration with nuclide data"""
        # This would test integration with actual nuclide reactions
        # For now, just test that our interfaces are compatible
        
        product = mc.create_test_reaction_product()
        
        # Simulate typical neutron transport energies
        transport_energies = [
            2.53e-8,  # Thermal
            1e-6,     # Epithermal
            1e-3,     # Intermediate
            1e0,      # Fast
            14e6,     # Fusion
        ]
        
        for energy in transport_energies:
            e_out, mu = product.sample(energy)
            
            # Basic physics checks
            assert e_out > 0
            assert -1.0 <= mu <= 1.0


if __name__ == "__main__":
    pytest.main([__file__, "-v"])