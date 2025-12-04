#!/usr/bin/env python3
"""
Test getter/setter functionality specifically for IndependentSource.angle property.
"""
import yamc as mc


def test_angle_getter_setter_consistency():
    """Test that angle getter/setter work consistently."""
    source = mc.IndependentSource()
    
    # Default should be Isotropic
    angle = source.angle
    assert isinstance(angle, mc.stats.Isotropic)
    assert str(angle) == "Isotropic()"
    
    # Set to Monodirectional and verify getter returns correct object
    mono = mc.stats.Monodirectional([1.0, 0.0, 0.0])
    source.angle = mono
    
    retrieved_angle = source.angle
    assert isinstance(retrieved_angle, mc.stats.Monodirectional)
    assert retrieved_angle.reference_uvw == [1.0, 0.0, 0.0]
    assert str(retrieved_angle) == "Monodirectional(reference_uvw=[1.0, 0.0, 0.0])"
    
    # Set to Isotropic and verify
    iso = mc.stats.Isotropic()
    source.angle = iso
    
    retrieved_angle = source.angle
    assert isinstance(retrieved_angle, mc.stats.Isotropic)
    assert str(retrieved_angle) == "Isotropic()"


def test_angle_property_consistency():
    """Test that angle property works consistently across multiple sources."""
    source1 = mc.IndependentSource()
    source2 = mc.IndependentSource()
    
    # Set same angle on both using property
    mono = mc.stats.Monodirectional([0.0, 1.0, 0.0])
    source1.angle = mono
    source2.angle = mono
    
    # Both should produce same samples
    p1 = source1.sample()
    p2 = source2.sample()
    assert p1.direction == p2.direction
    
    # And getters should return equivalent objects
    angle1 = source1.angle
    angle2 = source2.angle
    assert isinstance(angle1, mc.stats.Monodirectional)
    assert isinstance(angle2, mc.stats.Monodirectional)
    assert angle1.reference_uvw == angle2.reference_uvw


def test_angle_getter_returns_independent_objects():
    """Test that angle getter returns independent objects (not shared references)."""
    source = mc.IndependentSource()
    source.angle = mc.stats.Monodirectional([1.0, 0.0, 0.0])
    
    # Get angle object twice
    angle1 = source.angle
    angle2 = source.angle
    
    # Should be equivalent but independent objects
    assert angle1.reference_uvw == angle2.reference_uvw
    assert angle1 is not angle2  # Different objects
    
    # Modifying one shouldn't affect the other or the source
    angle1.reference_uvw = [0.0, 1.0, 0.0]
    
    # Source and angle2 should be unchanged
    assert source.angle.reference_uvw == [1.0, 0.0, 0.0]
    assert angle2.reference_uvw == [1.0, 0.0, 0.0]


def test_all_properties_have_getters_setters():
    """Test that all properties consistently have both getters and setters."""
    source = mc.IndependentSource()
    
    # Test space property
    assert hasattr(source, 'space')  # getter
    original_space = source.space
    source.space = [1.0, 2.0, 3.0]  # setter
    assert source.space == [1.0, 2.0, 3.0]
    
    # Test energy property  
    assert hasattr(source, 'energy')  # getter
    original_energy = source.energy
    source.energy = 2e6  # setter
    assert source.energy == 2e6
    
    # Test angle property
    assert hasattr(source, 'angle')  # getter
    original_angle = source.angle
    source.angle = mc.stats.Monodirectional([0.0, 0.0, 1.0])  # setter
    retrieved_angle = source.angle
    assert isinstance(retrieved_angle, mc.stats.Monodirectional)
    assert retrieved_angle.reference_uvw == [0.0, 0.0, 1.0]


if __name__ == "__main__":
    test_angle_getter_setter_consistency()
    test_angle_property_vs_method_consistency()
    test_angle_getter_returns_independent_objects()
    test_all_properties_have_getters_setters()
    print("All IndependentSource angle property tests passed!")