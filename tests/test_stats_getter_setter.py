#!/usr/bin/env python3
"""
Test getter/setter functionality for stats distributions.
"""
import yamc as mc


def test_monodirectional_getter_setter():
    """Test that Monodirectional has both getter and setter for reference_uvw."""
    # Create with initial direction
    mono = mc.stats.Monodirectional([1.0, 0.0, 0.0])
    assert mono.reference_uvw == [1.0, 0.0, 0.0]
    
    # Test setter
    mono.reference_uvw = [0.0, 1.0, 0.0]
    assert mono.reference_uvw == [0.0, 1.0, 0.0]
    
    # Test another direction
    mono.reference_uvw = [0.0, 0.0, 1.0]
    assert mono.reference_uvw == [0.0, 0.0, 1.0]
    
    # Test that samples return the set direction
    sample1 = mono.sample()
    sample2 = mono.sample()
    assert sample1 == sample2 == [0.0, 0.0, 1.0]


def test_isotropic_has_no_settable_properties():
    """Test that Isotropic doesn't have settable properties (as expected)."""
    iso = mc.stats.Isotropic()
    
    # Should not have reference_uvw property
    assert not hasattr(iso, 'reference_uvw')
    
    # Should still sample properly
    sample1 = iso.sample()
    sample2 = iso.sample()
    assert len(sample1) == 3
    assert len(sample2) == 3
    # Samples should be normalized
    import math
    mag1 = math.sqrt(sum(x*x for x in sample1))
    mag2 = math.sqrt(sum(x*x for x in sample2))
    assert abs(mag1 - 1.0) < 1e-10
    assert abs(mag2 - 1.0) < 1e-10


def test_setter_consistency_with_samples():
    """Test that setting reference_uvw affects sampling consistently."""
    mono = mc.stats.Monodirectional([1.0, 0.0, 0.0])
    
    # Initial samples should match reference_uvw
    for _ in range(5):
        sample = mono.sample()
        assert sample == mono.reference_uvw
    
    # Change direction and test again
    mono.reference_uvw = [-1.0, 0.0, 0.0]
    for _ in range(5):
        sample = mono.sample()
        assert sample == [-1.0, 0.0, 0.0]
        assert sample == mono.reference_uvw


if __name__ == "__main__":
    test_monodirectional_getter_setter()
    test_isotropic_has_no_settable_properties()
    test_setter_consistency_with_samples()
    print("All getter/setter tests passed!")