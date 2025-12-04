#!/usr/bin/env python3
"""Tests for tally functionality."""

import pytest
import yamc as mc


class TestTally:
    """Test tally creation and basic properties."""
    
    def test_tally_creation(self):
        """Test creating a new tally."""
        tally = mc.Tally()
        assert tally.scores == []  # Default score
        assert tally.name is None
        assert tally.id is None
        assert tally.units == ""
        assert tally.mean == []
        
    def test_tally_score_assignment(self):
        """Test setting tally score as single integer."""
        tally = mc.Tally()
        
        # Test absorption
        tally.scores = [101]
        assert tally.scores == [101]
        
        # Test elastic scattering  
        tally.scores = [2]
        assert tally.scores == [2]
        
        # Test fission
        tally.scores = [18]
        assert tally.scores == [18]
        
    def test_tally_score_type_validation(self):
        """Test that tally score only accepts integers or string scores."""
        tally = mc.Tally()
        
        # These should fail
        with pytest.raises(TypeError):
            tally.scores = '101'  # Should be list, not string
            
        with pytest.raises(TypeError):
            tally.scores = 101  # Should be list, not int

        with pytest.raises(TypeError):
            tally.scores = 101.0  # Should be list, not float
    
    def test_tally_flux_score_string(self):
        """Test setting flux score using string."""
        tally = mc.Tally()
        
        # Test flux score with string
        tally.scores = ['flux']
        assert tally.scores == ['flux']

        # Test mixed scores
        tally.scores = [101, 'flux']
        assert tally.scores == [101, 'flux']
        
        # Test invalid string
        with pytest.raises(ValueError):
            tally.scores = ['invalid_score']
    
    def test_tally_flux_score_constant(self):
        """Test that FLUX_SCORE constant is accessible."""
        # Can use constant directly
        tally = mc.Tally()
        tally.scores = ['flux']
        assert tally.scores == ['flux']
            
    def test_tally_name_and_id(self):
        """Test setting tally name and ID."""
        tally = mc.Tally()
        
        # Test name
        tally.name = "Absorption Tally"
        assert tally.name == "Absorption Tally"
        
        # Test ID
        tally.id = 42
        assert tally.id == 42
        
        # Test None values
        tally.name = None
        tally.id = None
        assert tally.name is None
        assert tally.id is None
        
    def test_tally_repr(self):
        """Test tally string representation."""
        tally = mc.Tally()
        tally.scores = [101]
        tally.name = "Test Tally"
        tally.id = 1
        
        repr_str = repr(tally)
        assert "Tally(" in repr_str
        assert "scores=[101]" in repr_str
        assert 'Some("Test Tally")' in repr_str  # Rust Option format
        assert "Some(1)" in repr_str


class TestTallySimulation:
    """Test tallies in actual Monte Carlo simulation."""
    
    @pytest.fixture
    def simple_model(self):
        """Create a simple model for testing tallies."""
        # Create sphere surface with vacuum boundary
        sphere = mc.Sphere(
            surface_id=1,
            x0=0.0,
            y0=0.0,
            z0=0.0,
            r=2.0,
            boundary_type='vacuum',
        )
        region = -sphere
        
        # Create material with Li6
        material = mc.Material()
        material.add_nuclide("Li6", 1.0)
        material.set_density("g/cm3", 5.5)
        material.read_nuclides_from_json({"Li6": "tests/Li6.json"})
        
        # Create cell
        cell = mc.Cell(
            cell_id=1,
            name="sphere_cell",
            region=region,
            fill=material,
        )
        geometry = mc.Geometry(cells=[cell])
        
        # Create source and settings
        source = mc.IndependentSource(space=[0.0, 0.0, 0.0], angle=mc.stats.Monodirectional(reference_uvw=[0.0, 0.0, 1.0]), energy=1e6)
        settings = mc.Settings(particles=100, batches=5, source=source)
        
        return geometry, settings
    
    def test_simulation_with_absorption_tally(self, simple_model):
        """Test simulation with absorption tally."""
        geometry, settings = simple_model
        
        # Create absorption tally
        absorption_tally = mc.Tally()
        absorption_tally.scores = [101]  # MT 101 = absorption
        absorption_tally.name = "Absorption Tally"
        tallies = [absorption_tally]
        
        # Create and run model
        model = mc.Model(geometry=geometry, settings=settings, tallies=tallies)
        model.run()
        
        # Check absorption tally
        assert absorption_tally.name == "Absorption Tally"
        assert absorption_tally.n_batches == 5
        assert absorption_tally.particles_per_batch == 100
        
        # Statistics should be calculated
        if absorption_tally.total_count()[0] > 0:
            assert absorption_tally.std_dev[0] >= 0.0
            assert absorption_tally.rel_error[0] >= 0.0
    
    def test_simulation_with_multiple_tallies(self, simple_model):
        """Test simulation with multiple tallies."""
        geometry, settings = simple_model
        
        # Create multiple tallies
        absorption_tally = mc.Tally()
        absorption_tally.scores = [101]  # Absorption
        absorption_tally.name = "Absorption Events"
        
        elastic_tally = mc.Tally()
        elastic_tally.scores = [2]  # Elastic scattering
        elastic_tally.name = "Elastic Scattering Events"
        
        tallies = [absorption_tally, elastic_tally]
        
        # Create and run model
        model = mc.Model(geometry=geometry, settings=settings, tallies=tallies)
        model.run()

        
        # Check all tallies
        assert absorption_tally.name == "Absorption Events"
        assert elastic_tally.name == "Elastic Scattering Events"
        
        # All should have proper batch information
        for tally_result in model.tallies:
            assert tally_result.n_batches == 5
            assert tally_result.particles_per_batch == 100
            assert len(tally_result.batch_data[0]) == 5
            
    def test_simulation_without_user_tallies(self, simple_model):
        """Test simulation with only leakage tally (no user tallies)."""
        geometry, settings = simple_model
        
        # Create model with no user tallies
        model = mc.Model(geometry=geometry, settings=settings, tallies=[])
        model.run()

    def test_tally_statistics_consistency(self, simple_model):
        """Test that tally statistics are consistent and reasonable."""
        geometry, settings = simple_model
        
        # Create absorption tally
        absorption_tally = mc.Tally()
        absorption_tally.scores = [101]
        absorption_tally.name = "Statistics Test"
        tallies = [absorption_tally]
        
        # Run simulation
        model = mc.Model(geometry=geometry, settings=settings, tallies=tallies)
        model.run()
        
        # Test statistics consistency
        assert absorption_tally.n_batches == len(absorption_tally.batch_data[0])
        assert absorption_tally.total_count()[0] == sum(absorption_tally.batch_data[0])

        # If we have results, test statistical relationships
        if absorption_tally.total_count()[0] > 0:
            # Mean should be total divided by (batches * particles_per_batch)
            expected_mean = absorption_tally.total_count()[0] / (
                absorption_tally.n_batches * absorption_tally.particles_per_batch
            )
            assert abs(absorption_tally.mean[0] - expected_mean) < 1e-10

            # Relative error should be std_dev / mean (if mean > 0)
            if absorption_tally.mean[0] > 0:
                expected_rel_error = absorption_tally.std_dev[0] / absorption_tally.mean[0]
                assert abs(absorption_tally.rel_error[0] - expected_rel_error) < 1e-10


class TestFluxTally:
    """Test flux tallying functionality."""
    
    @pytest.fixture
    def flux_test_model(self):
        """Create a simple model for testing flux tallies."""
        # Create sphere surface with vacuum boundary
        sphere = mc.Sphere(
            surface_id=1,
            x0=0.0,
            y0=0.0,
            z0=0.0,
            r=5.0,
            boundary_type='vacuum',
        )
        region = -sphere
        
        # Create material with Li6
        material = mc.Material()
        material.add_nuclide("Li6", 1.0)
        material.set_density("g/cm3", 0.46)
        material.read_nuclides_from_json({"Li6": "tests/Li6.json"})
        
        # Create cell
        cell = mc.Cell(
            cell_id=1,
            name="sphere_cell",
            region=region,
            fill=material,
        )
        geometry = mc.Geometry(cells=[cell])
        
        # Create source and settings - point source at center
        source = mc.IndependentSource(
            space=[0.0, 0.0, 0.0],
            angle=mc.stats.Isotropic(),
            energy=14.1e6  # 14.1 MeV
        )
        settings = mc.Settings(particles=1000, batches=10, source=source)
        
        return geometry, settings
    
    def test_flux_tally_using_string(self, flux_test_model):
        """Test flux tally using 'flux' string."""
        geometry, settings = flux_test_model
        
        # Create flux tally using string
        flux_tally = mc.Tally()
        flux_tally.scores = ['flux']
        flux_tally.name = "Flux Tally"
        tallies = [flux_tally]
        
        # Create and run model
        model = mc.Model(geometry=geometry, settings=settings, tallies=tallies)
        model.run()
        
        # Check flux tally results
        assert flux_tally.name == "Flux Tally"
        assert flux_tally.n_batches == 10
        assert flux_tally.particles_per_batch == 1000
        assert len(flux_tally.mean) == 1
        
        # Flux should be positive
        assert flux_tally.mean[0] > 0.0
        
        # Should have statistics
        assert flux_tally.std_dev[0] >= 0.0
        assert flux_tally.rel_error[0] >= 0.0
    
    def test_flux_tally_using_constant(self, flux_test_model):
        """Test flux tally using FLUX_SCORE constant."""
        geometry, settings = flux_test_model
        
        # Create flux tally using constant
        flux_tally = mc.Tally()
        flux_tally.scores = ['flux']
        flux_tally.name = "Flux with Constant"
        tallies = [flux_tally]
        
        # Create and run model
        model = mc.Model(geometry=geometry, settings=settings, tallies=tallies)
        model.run()
        
        # Check results
        assert flux_tally.mean[0] > 0.0
    
    def test_flux_and_reaction_tally_mixed(self, flux_test_model):
        """Test tally with both flux and reaction scores."""
        geometry, settings = flux_test_model
        
        # Create tally with flux and absorption
        tally = mc.Tally()
        tally.scores = ['flux', 101]  # flux and absorption
        tally.name = "Mixed Tally"
        tallies = [tally]
        
        # Create and run model
        model = mc.Model(geometry=geometry, settings=settings, tallies=tallies)
        model.run()
        
        # Should have two scores
        assert len(tally.mean) == 2
        assert tally.scores == ['flux', 101]
        
        # Flux (index 0) should be positive
        assert tally.mean[0] > 0.0
        
        # Absorption might be zero or positive
        assert tally.mean[1] >= 0.0


class TestTallyIntegration:
    """Integration tests for tally system."""
    
    def test_tally_display_output(self):
        """Test that tally display output is reasonable."""
        tally = mc.Tally()
        tally.scores = [101]
        tally.name = "Test Display"

        # Test string output doesn't crash
        str_output = str(tally)
        assert "Test Display" in str_output
        assert "Mean:" in str_output
        
    def test_model_constructor_with_tallies(self):
        """Test that model accepts tallies in constructor."""
        # Create minimal geometry
        sphere = mc.Sphere(surface_id=1, x0=0.0, y0=0.0, z0=0.0, r=1.0, boundary_type='vacuum')
        region = -sphere
        material = mc.Material()
        material.add_nuclide("Li6", 1.0)
        material.set_density("g/cm3", 1.0)
        material.read_nuclides_from_json({"Li6": "tests/Li6.json"})
        cell = mc.Cell(cell_id=1, name="test", region=region, fill=material)
        geometry = mc.Geometry(cells=[cell])
        
        # Create source and settings
        source = mc.IndependentSource(space=[0.0, 0.0, 0.0], angle=mc.stats.Monodirectional([0.0,0.0,1.0]), energy=1e6)
        settings = mc.Settings(particles=10, batches=1, source=source)
        
        # Create tally
        tally = mc.Tally()
        tally.scores = [101]
        tallies = [tally]
        
        # Model should accept tallies
        model = mc.Model(geometry=geometry, settings=settings, tallies=tallies)
        assert model is not None
        
        # Should be able to access tallies
        assert len(model.tallies) == 1
        assert model.tallies[0].scores == [101]