import constructive_solid_geometry_for_mc as csg4mc

def test_sphere_bb_moved_on_z_axis():
    s2 = csg4mc.Sphere(x0=0, y0=0, z0=1, r=3, surface_id=1)
    region2 = -s2
    bb = region2.bounding_box()
    assert bb.lower_left == [-3.0, -3.0, -2.0]
    assert bb.upper_right == [3.0, 3.0, 4.0]

def test_sphere_with_xplanes():
    s1 = csg4mc.XPlane(x0=2.1, surface_id=5)
    s2 = csg4mc.XPlane(x0=-2.1, surface_id=6)
    s3 = csg4mc.Sphere(x0=0, y0=0, z0=0, r=4.2, surface_id=1)

    region1 = -s1 & +s2 & -s3
    assert region1.contains((0, 0, 0))
    bb = region1.bounding_box()
    assert bb.lower_left == [-2.1, -4.2, -4.2]
    assert bb.upper_right == [2.1, 4.2, 4.2]

def test_zcylinder_with_zplanes():
    # Create a Z-cylinder centered at (1, 2) with radius 3
    cyl = csg4mc.ZCylinder(x0=1.0, y0=2.0, r=3.0, surface_id=1)
    # Create Z planes to bound the cylinder in Z direction
    z_bottom = csg4mc.ZPlane(z0=-5.0, surface_id=2)
    z_top = csg4mc.ZPlane(z0=5.0, surface_id=3)
    
    # Region inside cylinder and between the Z planes
    region = -cyl & +z_bottom & -z_top
    
    # Test that points are contained as expected
    assert region.contains((1.0, 2.0, 0.0))  # Center of cylinder
    assert region.contains((3.0, 2.0, 0.0))  # On cylinder surface in +X
    assert region.contains((1.0, 4.0, 0.0))  # On cylinder surface in +Y
    assert not region.contains((5.0, 2.0, 0.0))  # Outside cylinder
    assert not region.contains((1.0, 2.0, 6.0))  # Above Z plane
    
    # Get bounding box - should be bounded by planes in X, Y and Z planes in Z
    bb = region.bounding_box()
    # X bounds: cylinder center (1) ± radius (3) = [-2, 4]
    # Y bounds: cylinder center (2) ± radius (3) = [-1, 5] 
    # Z bounds: between the Z planes = [-5, 5]
    assert bb.lower_left == [-2.0, -1.0, -5.0]
    assert bb.upper_right == [4.0, 5.0, 5.0]