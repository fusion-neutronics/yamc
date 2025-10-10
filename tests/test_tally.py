#!/usr/bin/env python3
"""Tests for tally functionality."""

import pytest
import materials_for_mc as mc


class TestTally:
    """Test tally creation and basic properties."""
    
    def test_tally_creation(self):
        """Test creating a new tally."""
        tally = mc.Tally()
        assert tally.score == 0  # Default score
        assert tally.name is None
        assert tally.id is None
        assert tally.units == ""
        assert tally.mean == 0.0
        
    def test_tally_score_assignment(self):
        """Test setting tally score as single integer."""
        tally = mc.Tally()
        
        # Test absorption
        tally.score = 101
        assert tally.score == 101
        assert isinstance(tally.score, int)
        
        # Test elastic scattering  
        tally.score = 2
        assert tally.score == 2
        
        # Test fission
        tally.score = 18
        assert tally.score == 18
        
    def test_tally_score_type_validation(self):
        """Test that tally score only accepts integers."""
        tally = mc.Tally()
        
        # These should fail
        with pytest.raises(TypeError):
            tally.score = [101, 2]  # No lists
            
        with pytest.raises(TypeError):
            tally.score = "101"  # No strings
            
        with pytest.raises(TypeError):
            tally.score = 101.0  # No floats
            
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
        tally.score = 101
        tally.name = "Test Tally"
        tally.id = 1
        
        repr_str = repr(tally)
        assert "Tally(" in repr_str
        assert "score=101" in repr_str
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
            boundary_type='Vacuum',
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
        absorption_tally.score = 101  # MT 101 = absorption
        absorption_tally.name = "Absorption Tally"
        tallies = [absorption_tally]
        
        # Create and run model
        model = mc.Model(geometry=geometry, settings=settings, tallies=tallies)
        results = model.run()
        
        # Should return leakage tally + user tallies
        assert len(results) == 2
        leakage_tally, absorption_result = results
        
        # Check leakage tally
        assert leakage_tally.name == "Leakage"
        assert leakage_tally.units == "particles"
        assert leakage_tally.n_batches == 5
        assert leakage_tally.particles_per_batch == 100
        
        # Check absorption tally
        assert absorption_result.name == "Absorption Tally"
        assert absorption_result.units == "events"
        assert absorption_result.n_batches == 5
        assert absorption_result.particles_per_batch == 100
        
        # Results should be reasonable (non-negative, finite)
        assert leakage_tally.mean >= 0.0
        assert leakage_tally.total_count() >= 0
        assert absorption_result.mean >= 0.0
        assert absorption_result.total_count() >= 0
        
        # Statistics should be calculated
        if absorption_result.total_count() > 0:
            assert absorption_result.std_dev >= 0.0
            assert absorption_result.rel_error >= 0.0
    
    def test_simulation_with_multiple_tallies(self, simple_model):
        """Test simulation with multiple tallies."""
        geometry, settings = simple_model
        
        # Create multiple tallies
        absorption_tally = mc.Tally()
        absorption_tally.score = 101  # Absorption
        absorption_tally.name = "Absorption Events"
        
        elastic_tally = mc.Tally()
        elastic_tally.score = 2  # Elastic scattering
        elastic_tally.name = "Elastic Scattering Events"
        
        tallies = [absorption_tally, elastic_tally]
        
        # Create and run model
        model = mc.Model(geometry=geometry, settings=settings, tallies=tallies)
        results = model.run()
        
        # Should return leakage + 2 user tallies
        assert len(results) == 3
        leakage_tally, absorption_result, elastic_result = results
        
        # Check all tallies
        assert leakage_tally.name == "Leakage"
        assert absorption_result.name == "Absorption Events"
        assert elastic_result.name == "Elastic Scattering Events"
        
        # All should have proper batch information
        for tally_result in results:
            assert tally_result.n_batches == 5
            assert tally_result.particles_per_batch == 100
            assert len(tally_result.batch_data) == 5
            
    def test_simulation_without_user_tallies(self, simple_model):
        """Test simulation with only leakage tally (no user tallies)."""
        geometry, settings = simple_model
        
        # Create model with no user tallies
        model = mc.Model(geometry=geometry, settings=settings, tallies=[])
        results = model.run()
        
        # Should only return leakage tally
        assert len(results) == 1
        leakage_tally = results[0]
        
        assert leakage_tally.name == "Leakage"
        assert leakage_tally.units == "particles"
        assert leakage_tally.n_batches == 5
        assert leakage_tally.particles_per_batch == 100
        
    def test_tally_statistics_consistency(self, simple_model):
        """Test that tally statistics are consistent and reasonable."""
        geometry, settings = simple_model
        
        # Create absorption tally
        absorption_tally = mc.Tally()
        absorption_tally.score = 101
        absorption_tally.name = "Statistics Test"
        tallies = [absorption_tally]
        
        # Run simulation
        model = mc.Model(geometry=geometry, settings=settings, tallies=tallies)
        results = model.run()
        
        _, absorption_result = results
        
        # Test statistics consistency
        assert absorption_result.n_batches == len(absorption_result.batch_data)
        assert absorption_result.total_count() == sum(absorption_result.batch_data)
        
        # If we have results, test statistical relationships
        if absorption_result.total_count() > 0:
            # Mean should be total divided by (batches * particles_per_batch)
            expected_mean = absorption_result.total_count() / (
                absorption_result.n_batches * absorption_result.particles_per_batch
            )
            assert abs(absorption_result.mean - expected_mean) < 1e-10
            
            # Relative error should be std_dev / mean (if mean > 0)
            if absorption_result.mean > 0:
                expected_rel_error = absorption_result.std_dev / absorption_result.mean
                assert abs(absorption_result.rel_error - expected_rel_error) < 1e-10


class TestTallyIntegration:
    """Integration tests for tally system."""
    
    def test_tally_display_output(self):
        """Test that tally display output is reasonable."""
        tally = mc.Tally()
        tally.score = 101
        tally.name = "Test Display"
        
        # Test string output doesn't crash
        str_output = str(tally)
        assert "Test Display" in str_output
        assert "Mean:" in str_output
        
    def test_model_constructor_with_tallies(self):
        """Test that model accepts tallies in constructor."""
        # Create minimal geometry
        sphere = mc.Sphere(surface_id=1, x0=0.0, y0=0.0, z0=0.0, r=1.0, boundary_type='Vacuum')
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
        tally.score = 101
        tallies = [tally]
        
        # Model should accept tallies
        model = mc.Model(geometry=geometry, settings=settings, tallies=tallies)
        assert model is not None
        
        # Should be able to access tallies
        assert len(model.tallies) == 1
        assert model.tallies[0].score == 101