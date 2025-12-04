"""
Comprehensive pytest tests for yaml Python API

Tests the OpenMC-compatible IndependentSource and stats module functionality.
"""

import pytest
import yaml as mc


class TestIndependentSource:
    """Test cases for IndependentSource class."""

    def test_parameterless_constructor(self):
        """Test that IndependentSource() works without arguments with correct defaults."""
        source = mc.IndependentSource()
        
        # Check default values
        assert source.space == [0.0, 0.0, 0.0]
        assert source.energy == 14.06e6
        
        # Sample should work
        particle = source.sample()
        assert particle.position == [0.0, 0.0, 0.0]
        assert particle.energy == 14.06e6
        assert particle.alive is True
        
        # Direction should be normalized (isotropic by default)
        direction = particle.direction
        magnitude = sum(d * d for d in direction) ** 0.5
        assert abs(magnitude - 1.0) < 1e-10

    def test_space_property(self):
        """Test space property getter and setter."""
        source = mc.IndependentSource()
        
        # Test getter
        assert source.space == [0.0, 0.0, 0.0]
        
        # Test setter
        source.space = [1.0, 2.0, 3.0]
        assert source.space == [1.0, 2.0, 3.0]
        
        # Test particle sampling reflects new space
        particle = source.sample()
        assert particle.position == [1.0, 2.0, 3.0]

    def test_energy_property(self):
        """Test energy property getter and setter.""" 
        source = mc.IndependentSource()
        
        # Test getter (default)
        assert source.energy == 14.06e6
        
        # Test setter
        source.energy = 2.5e6
        assert source.energy == 2.5e6
        
        # Test particle sampling reflects new energy
        particle = source.sample()
        assert particle.energy == 2.5e6

    def test_set_angle_isotropic(self):
        """Test setting angle to isotropic using property."""
        source = mc.IndependentSource()
        
        # Should be isotropic by default, but set it explicitly using property
        source.angle = mc.stats.Isotropic()
        
        # Sample multiple particles and check they have different directions
        directions = []
        for _ in range(10):
            particle = source.sample()
            direction = particle.direction
            
            # Each direction should be normalized
            magnitude = sum(d * d for d in direction) ** 0.5
            assert abs(magnitude - 1.0) < 1e-10
            
            directions.append(tuple(direction))
        
        # With isotropic sampling, very unlikely to get duplicates
        unique_directions = set(directions)
        assert len(unique_directions) > 1

    def test_set_angle_monodirectional(self):
        """Test setting angle to monodirectional using property."""
        source = mc.IndependentSource()
        
        # Set monodirectional using property
        reference_direction = [0.0, 0.0, 1.0]
        source.angle = mc.stats.Monodirectional(reference_direction)
        
        # All sampled particles should have the same direction
        for _ in range(10):
            particle = source.sample()
            assert particle.direction == reference_direction

    def test_angle_switching(self):
        """Test switching between isotropic and monodirectional angles."""
        source = mc.IndependentSource()
        
        # Start with monodirectional
        source.angle = mc.stats.Monodirectional([1.0, 0.0, 0.0])
        particle1 = source.sample()
        assert particle1.direction == [1.0, 0.0, 0.0]
        
        # Switch to isotropic
        source.angle = mc.stats.Isotropic()
        particle2 = source.sample()
        particle3 = source.sample()
        
        # Directions should be normalized but likely different
        for p in [particle2, particle3]:
            magnitude = sum(d * d for d in p.direction) ** 0.5
            assert abs(magnitude - 1.0) < 1e-10
        
        # Very unlikely to be the same
        assert particle2.direction != particle3.direction

    def test_repr(self):
        """Test string representation."""
        source = mc.IndependentSource()
        repr_str = repr(source)
        
        # Should contain space and energy info
        assert "IndependentSource" in repr_str
        assert "space=[0.0, 0.0, 0.0]" in repr_str
        assert "energy=14060000" in repr_str


class TestStatsModule:
    """Test cases for the stats submodule."""

    def test_stats_module_exists(self):
        """Test that stats submodule exists and contains expected classes."""
        assert hasattr(mc, 'stats')
        assert hasattr(mc.stats, 'Isotropic')
        assert hasattr(mc.stats, 'Monodirectional')

    def test_old_direct_access_fails(self):
        """Test that old direct access to distributions fails."""
        with pytest.raises(AttributeError):
            mc.Isotropic()
        
        with pytest.raises(AttributeError):
            mc.Monodirectional([1.0, 0.0, 0.0])


class TestIsotropic:
    """Test cases for Isotropic distribution."""

    def test_construction(self):
        """Test Isotropic construction."""
        iso = mc.stats.Isotropic()
        assert iso is not None

    def test_sample(self):
        """Test Isotropic sampling."""
        iso = mc.stats.Isotropic()
        
        # Sample multiple directions
        directions = []
        for _ in range(100):
            direction = iso.sample()
            
            # Should be 3D
            assert len(direction) == 3
            
            # Should be normalized
            magnitude = sum(d * d for d in direction) ** 0.5
            assert abs(magnitude - 1.0) < 1e-10
            
            directions.append(tuple(direction))
        
        # Should have many unique directions (isotropic)
        unique_directions = set(directions)
        assert len(unique_directions) > 50  # Conservative check for randomness

    def test_repr(self):
        """Test Isotropic string representation."""
        iso = mc.stats.Isotropic()
        repr_str = repr(iso)
        assert "Isotropic" in repr_str


class TestMonodirectional:
    """Test cases for Monodirectional distribution."""

    def test_construction(self):
        """Test Monodirectional construction."""
        mono = mc.stats.Monodirectional([1.0, 0.0, 0.0])
        assert mono is not None

    def test_sample_consistency(self):
        """Test that Monodirectional always returns the same direction."""
        reference = [0.0, 1.0, 0.0]
        mono = mc.stats.Monodirectional(reference)
        
        # All samples should be identical
        for _ in range(10):
            direction = mono.sample()
            assert direction == reference

    def test_different_directions(self):
        """Test Monodirectional with different reference directions."""
        test_directions = [
            [1.0, 0.0, 0.0],
            [0.0, 1.0, 0.0], 
            [0.0, 0.0, 1.0],
            [-1.0, 0.0, 0.0],
            [0.5773502691896257, 0.5773502691896257, 0.5773502691896257]  # Normalized (1,1,1)
        ]
        
        for ref_dir in test_directions:
            mono = mc.stats.Monodirectional(ref_dir)
            sample = mono.sample()
            assert sample == ref_dir

    def test_repr(self):
        """Test Monodirectional string representation."""
        mono = mc.stats.Monodirectional([1.0, 0.0, 0.0])
        repr_str = repr(mono)
        assert "Monodirectional" in repr_str
        assert "[1.0, 0.0, 0.0]" in repr_str


class TestIntegration:
    """Integration tests for the complete API."""

    def test_openmc_workflow(self):
        """Test a complete OpenMC-style workflow."""
        # Create source with defaults
        source = mc.IndependentSource()
        
        # Modify properties
        source.space = [0.0, 0.0, -10.0]
        source.energy = 14.1e6
        
        # Set different angular distributions
        source.angle = mc.stats.Monodirectional([0.0, 0.0, 1.0])
        particle1 = source.sample()
        assert particle1.position == [0.0, 0.0, -10.0]
        assert particle1.energy == 14.1e6
        assert particle1.direction == [0.0, 0.0, 1.0]
        
        # Switch to isotropic
        source.angle = mc.stats.Isotropic()
        particle2 = source.sample()
        assert particle2.position == [0.0, 0.0, -10.0]
        assert particle2.energy == 14.1e6
        # Direction should be different (very likely)
        assert particle2.direction != [0.0, 0.0, 1.0]

    def test_direct_distribution_usage(self):
        """Test using distributions directly outside of source."""
        # Create distributions
        iso = mc.stats.Isotropic()
        mono = mc.stats.Monodirectional([1.0, 0.0, 0.0])
        
        # Sample from each
        iso_samples = [iso.sample() for _ in range(5)]
        mono_samples = [mono.sample() for _ in range(5)]
        
        # Isotropic should vary
        unique_iso = set(tuple(s) for s in iso_samples)
        assert len(unique_iso) > 1
        
        # Monodirectional should be constant
        unique_mono = set(tuple(s) for s in mono_samples)
        assert len(unique_mono) == 1
        assert mono_samples[0] == [1.0, 0.0, 0.0]

    def test_flexible_constructor(self):
        """Test the new flexible constructor with keyword arguments."""
        # Test all combinations of keyword arguments
        
        # Only space
        source1 = mc.IndependentSource(space=[1.0, 2.0, 3.0])
        assert source1.space == [1.0, 2.0, 3.0]
        assert source1.energy == 14.06e6  # default
        
        # Only energy
        source2 = mc.IndependentSource(energy=2e6)
        assert source2.space == [0.0, 0.0, 0.0]  # default
        assert source2.energy == 2e6
        
        # Only angle
        angle_dist = mc.stats.Monodirectional([0.0, 0.0, 1.0])
        source3 = mc.IndependentSource(angle=angle_dist)
        particle = source3.sample()
        assert source3.space == [0.0, 0.0, 0.0]  # default
        assert source3.energy == 14.06e6  # default
        assert particle.direction == [0.0, 0.0, 1.0]
        
        # Space and energy
        source4 = mc.IndependentSource(space=[5.0, 0.0, -2.0], energy=3e6)
        assert source4.space == [5.0, 0.0, -2.0]
        assert source4.energy == 3e6
        
        # Energy and angle
        source5 = mc.IndependentSource(
            energy=1e6, 
            angle=mc.stats.Monodirectional([1.0, 0.0, 0.0])
        )
        assert source5.energy == 1e6
        particle5 = source5.sample()
        assert particle5.direction == [1.0, 0.0, 0.0]
        
        # All three arguments (original failing case)
        source6 = mc.IndependentSource(
            space=[0.0, 0.0, 0.0], 
            angle=mc.stats.Monodirectional([0.0, 0.0, 1.0]), 
            energy=1e6
        )
        assert source6.space == [0.0, 0.0, 0.0]
        assert source6.energy == 1e6
        particle6 = source6.sample()
        assert particle6.direction == [0.0, 0.0, 1.0]

    def test_backwards_compatibility_failures(self):
        """Test that old positional constructor patterns fail appropriately."""
        # Old constructor with positional arguments should fail
        with pytest.raises(TypeError):
            mc.IndependentSource([0.0, 0.0, 0.0], [0.0, 0.0, 1.0], 1e6)


if __name__ == "__main__":
    pytest.main([__file__])