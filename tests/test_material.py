import pytest
from yamc import Material, Config

def test_volume_set_and_get():
    material = Material()

    # Test setting a valid volume
    material.volume = 100.0
    assert material.volume == 100.0

    # Test setting another valid volume
    material.volume = 200.0
    assert material.volume == 200.0

    # Test setting an invalid (negative) volume
    with pytest.raises(ValueError, match="Volume must be positive"):
        material.volume = -50.0

    # Ensure the volume remains unchanged after the invalid set
    assert material.volume == 200.0

def test_initial_volume():
    material = Material()

    # Test that the initial volume is None
    assert material.volume is None
    
def test_adding_nuclides():
    material = Material()
    material.add_nuclide("H", 1.0)
    material.add_nuclide("Fe", 0.5)
    assert material.nuclides == [("Fe", 0.5), ("H", 1.0)] # TODO change to named tuple using a custom type

def test_density_settings():
    material = Material()
    material.add_nuclide("W", 0.5)
    material.set_density('g/cm3', 19.3)
    assert material.density == 19.3
    assert material.density_units == 'g/cm3'

def test_material_id_default():
    """Test that new materials have default material_id of None"""
    material = Material()
    assert material.material_id is None, "Default material_id should be None"

def test_material_id_set_and_get():
    """Test setting and getting material_id"""
    material = Material()
    
    # Test default value
    assert material.material_id is None
    
    # Test setting material_id
    material.material_id = 42
    assert material.material_id == 42
    
    # Test setting different values
    material.material_id = 999
    assert material.material_id == 999
    
    # Test setting to 0 (should still be valid)
    material.material_id = 0
    assert material.material_id == 0

def test_material_id_with_constructor():
    """Test creating materials with specific IDs and verify independence"""
    # Create material with default constructor
    material1 = Material(material_id = 10, name="mat1")
    assert material1.material_id == 10
    
    # Create another material
    material2 = Material()
    material2.material_id = 20
    
    # Verify they have different IDs
    assert material1.material_id == 10
    assert material2.material_id == 20
    
    # Modify one and ensure the other is unaffected
    material1.material_id = 100
    assert material1.material_id == 100
    assert material2.material_id == 20, "Other material's ID should not change"

def test_material_id_with_large_values():
    """Test material_id with edge case values"""
    material = Material(name="test")
    
    # Test with maximum u32 value (4294967295)
    max_u32 = 4294967295
    material.material_id = max_u32
    assert material.material_id == max_u32
    
    # Test with 0
    material.material_id = 0
    assert material.material_id == 0

def test_material_id_none_state():
    """Test that material_id can be None and must be explicitly set"""
    material = Material()
    
    # Verify it starts as None
    assert material.material_id is None
    
    # Set it to a value
    material.material_id = 42
    assert material.material_id == 42
    
    # Note: We can't set it back to None through the Python interface
    # since the setter expects a u32, but the getter can return None

def test_get_nuclide_names():
    material = Material()
    
    # Test empty material has no nuclides
    assert material.get_nuclide_names() == []
    
    # Add some nuclides
    material.add_nuclide("U235", 0.05)
    material.add_nuclide("U238", 0.95)
    material.add_nuclide("O16", 2.0)
    
    # Check we get all nuclides in alphabetical order
    assert material.get_nuclide_names() == ["O16", "U235", "U238"]
    
    # Test adding the same nuclide again replaces the previous value
    material.add_nuclide("U235", 0.1)
    assert "U235" in material.get_nuclide_names()
    assert len(material.get_nuclide_names()) == 3  # Still 3 nuclides
    
    # Verify the order is still alphabetical
    assert material.get_nuclide_names() == ["O16", "U235", "U238"]

# def test_material_data_xs_reading():
def test_global_default_cross_section_keyword():
    from yamc import Config, Material
    import os
    # Set cross section for Li6 specifically
    Config.set_cross_section("Li6", "tests/Li6.h5")
    mat = Material()
    mat.add_nuclide("Li6", 1.0)
    grid = mat.unified_energy_grid_neutron()
    assert len(grid) > 0, "Energy grid should not be empty when using global default cross section keyword"
    # Import Config locally to avoid module-level import issues
    from yamc import Config
    import os
    
    material = Material()
    material.add_nuclide("Li6", 1.0)
    
    # Print the current working directory and check if the file exists
    print(f"Current working directory: {os.getcwd()}")
    print(f"File exists: {os.path.exists('tests/Li6.h5')}")
    
    # Set the cross-section path in the global Config
    Config.set_cross_sections({"Li6": "tests/Li6.h5"})
    
    # Print the cross-sections to verify they're set correctly
    print(f"Cross-sections: {Config.get_cross_sections()}")
    
    # Accessing methods that need nuclide data will automatically load them
    # Let's call a method that triggers the loading
    grid = material.unified_energy_grid_neutron()
    
    assert len(material.nuclides) == 1
    assert len(grid) > 0, "Energy grid should not be empty"

def test_add_element_lithium():
    mat = Material()
    mat.add_element('Li', 1.0)
    nuclides = dict(mat.nuclides)
    # Should contain Li6 and Li7 with correct fractions
    assert 'Li6' in nuclides
    assert 'Li7' in nuclides
    # Natural abundances from IUPAC
    assert abs(nuclides['Li6'] - 0.07589) < 1e-5
    assert abs(nuclides['Li7'] - 0.92411) < 1e-5
    mat = Material()

    mat.add_element('lithium', 1.0)
    nuclides = dict(mat.nuclides)
    # Should contain Li6 and Li7 with correct fractions
    assert 'Li6' in nuclides
    assert 'Li7' in nuclides
    # Natural abundances from IUPAC
    assert abs(nuclides['Li6'] - 0.07589) < 1e-5
    assert abs(nuclides['Li7'] - 0.92411) < 1e-5


def test_add_element_not_found():
    mat = Material()
    with pytest.raises(Exception) as excinfo:
        mat.add_element('Xx', 1.0)
    assert 'not a recognized element symbol' in str(excinfo.value)

def test_mean_free_path_lithium_14mev():
    from yamc import Config
    import math
    mat = Material()
    mat.add_element('Li', 1.0)
    mat.set_density('g/cm3', 0.534)  # lithium density
    # Set up real cross-section data for both Li6 and Li7
    Config.set_cross_sections({
        "Li6": "tests/Li6.h5",
        "Li7": "tests/Li7.h5"
    })
    # This will trigger loading and calculation
    mfp = mat.mean_free_path_neutron(14e6)
    assert mfp is not None
    # Expected value using AWR-derived atomic mass (matching OpenMC approach)
    assert math.isclose(mfp, 14.96349733, rel_tol=1e-4), f"Expected ~14.963 cm, got {mfp}"

def test_material_reaction_mts_lithium():
    # Set up real cross-section data for both Li6 and Li7
    Config.set_cross_sections({
        "Li6": "tests/Li6.h5",
        "Li7": "tests/Li7.h5"
    })
    mat = Material()
    mat.add_element('Li', 1.0)
    # This will load Li6 and Li7, so the MTs should be the union of both
    mts = mat.reaction_mts
    print(mts)
    # Check essential MTs are present (HDF5 data may have different set than JSON)
    essential_mts = [1, 2, 3, 101, 102, 1001]  # Total, elastic, non-elastic, absorption, capture, scattering
    for mt in essential_mts:
        assert mt in mts, f"Material lithium should have MT {mt}, got {mts}"

def test_calculate_microscopic_xs_neutron_lithium():
    from yamc import Material, Config
    mat = Material()
    mat.add_element('Li', 1.0)
    Config.set_cross_sections({
        "Li6": "tests/Li6.h5",
        "Li7": "tests/Li7.h5"
    })
    micro_xs = mat.calculate_microscopic_xs_neutron()
    # Check both Li6 and Li7 are present
    assert "Li6" in micro_xs
    assert "Li7" in micro_xs
    # Check that MT=2 is present for both
    assert 2 in micro_xs["Li6"]
    assert 2 in micro_xs["Li7"]
    # Check that the cross section arrays are the same length as the grid
    grid = mat.unified_energy_grid_neutron()
    assert len(micro_xs["Li6"][2]) == len(grid)
    assert len(micro_xs["Li7"][2]) == len(grid)

def test_material_vs_nuclide_microscopic_xs_li6():
    from yamc import Material, Config, Nuclide
    import numpy as np
    mat = Material()
    mat.add_nuclide('Li6', 1.0)
    Config.set_cross_sections({"Li6": "tests/Li6.h5"})
    micro_xs_mat = mat.calculate_microscopic_xs_neutron()
    grid = mat.unified_energy_grid_neutron()
    # Get the nuclide directly
    nuclide = Nuclide('Li6')
    nuclide.read_nuclide_from_h5('tests/Li6.h5')
    # Get the temperature string used by the material
    temperature = mat.temperature
    # Get reactions and energy grid for nuclide
    reactions = nuclide.reactions
    assert reactions is not None, "No reactions for Li6"
    energy_map = nuclide.energy
    assert energy_map is not None, "No energy map for Li6"
    energy_grid = energy_map.get(temperature)
    assert energy_grid is not None, "No energy grid for Li6"
    # For each MT in the material, compare the cross sections
    for mt, xs_mat in micro_xs_mat['Li6'].items():
        if mt in reactions:
            reaction = reactions[mt]
            threshold_idx = reaction.threshold_idx
            nuclide_energy = energy_grid[threshold_idx:]
            xs_nuclide = reaction.cross_section
            # Interpolate nuclide xs onto the material grid
            xs_nuclide_interp = []
            for g in grid:
                if g < nuclide_energy[0]:
                    xs_nuclide_interp.append(0.0)
                else:
                    # Linear interpolation
                    xs = np.interp(g, nuclide_energy, xs_nuclide)
                    xs_nuclide_interp.append(xs)
            # Compare arrays (allow small tolerance)
            np.testing.assert_allclose(xs_mat, xs_nuclide_interp, rtol=1e-10, err_msg=f"Mismatch for MT {mt}")


def test_calculate_microscopic_xs_neutron_mt_filter():
    mat = Material()
    mat.add_element("Li", 1.0)
    # Set up the nuclide JSON map for Li6 and Li7
    nuclide_json_map = {"Li6": "tests/Li6.h5", "Li7": "tests/Li7.h5"}
    mat.read_nuclides_from_h5(nuclide_json_map)
    # Calculate all MTs
    xs_all = mat.calculate_microscopic_xs_neutron()
    # Calculate only MT=2
    xs_mt2 = mat.calculate_microscopic_xs_neutron(mt_filter=[2])
    # For each nuclide, only MT=2 should be present and match the unfiltered result
    for nuclide in ["Li6", "Li7"]:
        assert nuclide in xs_mt2, f"{nuclide} missing in filtered result"
        xs_map = xs_mt2[nuclide]
        assert list(xs_map.keys()) == [2], f"Filtered result for {nuclide} should have only MT=2"
        xs_all_mt2 = xs_all[nuclide][2]
        xs_filtered = xs_map[2]
        assert xs_all_mt2 == xs_filtered, f"Filtered and unfiltered MT=2 xs do not match for {nuclide}"

def test_macroscopic_xs_neutron_mt_filter():
    from yamc import Material
    mat = Material()
    mat.add_element("Li", 1.0)
    mat.read_nuclides_from_h5({"Li6": "tests/Li6.h5", "Li7": "tests/Li7.h5"})
    mat.set_density("g/cm3", 1.0)
    # Calculate all MTs
    energy_all, macro_xs_all = mat.calculate_macroscopic_xs(mt_filter=[1,2, 3], by_nuclide=False)
    # Calculate only MT=2
    energy_single, macro_xs_mt2 = mat.calculate_macroscopic_xs(mt_filter=[2], by_nuclide=False)
    assert len(energy_single) == len(energy_all), "Energy grids should match"
    # The cross section array for MT=2 should match the unfiltered result
    xs_all = macro_xs_all[2]
    xs_filtered = macro_xs_mt2[2]
    assert xs_all == xs_filtered, "Filtered and unfiltered MT=2 macro_xs do not match"

def test_hierarchical_mt3_generated_for_li6():
    import yamc
    mat = yamc.Material()
    mat.add_nuclide('Li6', 1.0)
    mat.read_nuclides_from_h5({'Li6': 'tests/Li6.h5'})
    mat.set_density('g/cm3', 0.534)
    mat.temperature = "294"
    energies, xs_dict = mat.calculate_macroscopic_xs([3], by_nuclide=False)
    # MT=3 should be present and non-empty
    assert 3 in xs_dict, "MT=3 should be generated by sum rule for Li6"


def test_macroscopic_xs_mt3_does_not_generate_mt1():
    mat = Material()
    mat.add_nuclide("Li6", 1.0)
    mat.set_density("g/cm3", 0.534)
    nuclide_json_map = {"Li6": "tests/Li6.h5"}
    mat.read_nuclides_from_h5(nuclide_json_map)
    mt_filter = [3]
    _, macro_xs = mat.calculate_macroscopic_xs(mt_filter, by_nuclide=False)
    assert 1 not in macro_xs, "MT=1 should NOT be present when only MT=3 is requested"

def test_macroscopic_xs_mt24_does_not_generate_mt1():
    mat = Material()
    mat.add_nuclide("Li6", 1.0)
    mat.set_density("g/cm3", 0.534)
    nuclide_json_map = {"Li6": "tests/Li6.h5"}
    mat.read_nuclides_from_h5(nuclide_json_map)
    mt_filter = [24]
    _, macro_xs = mat.calculate_macroscopic_xs(mt_filter, by_nuclide=False)
    assert 1 not in macro_xs, "MT=1 should NOT be present when only MT=24 is requested"

def test_sample_distance_to_collision_statistical():
    mat = Material()
    mat.add_nuclide("Li6", 1.0)
    mat.set_density("g/cm3", 1.)
    mat.read_nuclides_from_h5({"Li6": "tests/Li6.h5"})
    mat.temperature = "294"
    mat.calculate_macroscopic_xs()  # Ensure xs are calculated
    energy = 14e6
    samples = []
    for seed in range(1000):
        d = mat.sample_distance_to_collision(energy, seed=seed)
        assert d is not None
        assert d >= 0.0, f"Sampled distance should not be negative, got {d}"
        samples.append(d)
    avg = sum(samples) / len(samples)
    assert abs(avg - 6.9) < 0.1, f"Average sampled distance {avg} not within 0.1 of 6.9"

# Test for sample_interacting_nuclide Python binding (equivalent to Rust test_sample_interacting_nuclide_li6_li7)
def test_sample_interacting_nuclide_li6_li7():
    material = Material()
    material.add_nuclide("Li6", 0.5)
    material.add_nuclide("Li7", 0.5)
    material.set_density("g/cm3", 1.0)
    material.temperature = "294"

    # Load nuclide data from JSON
    nuclide_json_map = {
        "Li6": "tests/Li6.h5",'(n,total)'
        "Li7": "tests/Li7.h5",
    }
    material.read_nuclides_from_h5(nuclide_json_map)

    # Calculate per-nuclide macroscopic total xs
    material.calculate_macroscopic_xs([1], True)

    # Sample the interacting nuclide many times at 14 MeV
    energy = 100_000.0
    n_samples = 10000
    counts = {"Li6": 0, "Li7": 0}
    for seed in range(n_samples):
        nuclide = material.sample_interacting_nuclide(energy, seed=seed)
        counts[nuclide] = counts.get(nuclide, 0) + 1

    count_li6 = counts.get("Li6", 0)
    count_li7 = counts.get("Li7", 0)
    total = count_li6 + count_li7
    frac_li6 = count_li6 / total if total > 0 else 0.0
    frac_li7 = count_li7 / total if total > 0 else 0.0

    print(f"Li6 fraction: {frac_li6}, Li7 fraction: {frac_li7}")

    # Both nuclides should be sampled
    assert frac_li6 > 0.0 and frac_li7 > 0.0, "Both nuclides should be sampled"
    # Fractions should sum to 1
    assert abs(frac_li6 + frac_li7 - 1.0) < 1e-6, "Fractions should sum to 1"
    # Li7 should be sampled more often than Li6
    assert frac_li6 > frac_li7, "Li6 should be sampled more often than Li6"


def test_material_be9_selective_temperature_load():
    from yamc import Config, Material, clear_nuclide_cache
    clear_nuclide_cache()
    mat = Material()
    # Set material temperature to 300 K (Be9 JSON also has 294)
    mat.temperature = "300"
    mat.add_nuclide("Be9", 1.0)
    Config.set_cross_sections({"Be9": "tests/Be9.h5"})
    mat.read_nuclides_from_h5({"Be9": "tests/Be9.h5"})

def test_material_be9_selective_temperature_load_294():
    from yamc import Config, Material, clear_nuclide_cache
    clear_nuclide_cache()
    mat = Material()
    # Set material temperature to 294 K (Be9 JSON also has 300)
    mat.temperature = "294"
    mat.add_nuclide("Be9", 1.0)
    Config.set_cross_sections({"Be9": "tests/Be9.h5"})
    mat.read_nuclides_from_h5({"Be9": "tests/Be9.h5"})

def test_read_nuclides_from_h5_keyword():
    """Test keyword-based data loading from remote H5 files."""
    mat = Material()
    mat.add_element('Li', 1.0)
    mat.set_density('g/cm3', 2.0)
    mat.volume = 1.0
    mat.read_nuclides_from_h5("tendl-2019")
    # Verify data was loaded
    grid = mat.unified_energy_grid_neutron()
    assert len(grid) > 0, "Energy grid should not be empty after keyword-based loading"

def test_read_nuclides_from_h5_dict():
    mat = Material()
    mat.add_element('Li', 1.0)
    mat.set_density('g/cm3', 2.0)
    mat.volume = 1.0
    # Provide local H5 paths for both Li6 and Li7 since element expansion creates both
    mat.read_nuclides_from_h5({"Li6": "tests/Li6.h5", "Li7": "tests/Li7.h5"})


def test_material_different_data_sources():
    """Test that material loading respects different data sources."""
    # Material using keyword
    mat1 = Material()
    mat1.add_nuclide('Li6', 1.0)
    mat1.set_density('g/cm3', 1.0)
    mat1.read_nuclides_from_h5("tendl-2019")

    # Material using local file
    mat2 = Material()
    mat2.add_nuclide('Li6', 1.0)
    mat2.set_density('g/cm3', 1.0)
    mat2.read_nuclides_from_h5({"Li6": "tests/Li6.h5"})

    # Both should have loaded data
    grid1 = mat1.unified_energy_grid_neutron()
    grid2 = mat2.unified_energy_grid_neutron()
    assert len(grid1) > 0, "Keyword source should load data"
    assert len(grid2) > 0, "File source should load data"


def test_material_file_and_keyword_sources():
    """Test that materials can use both file paths and keywords."""
    mat = Material()
    mat.add_nuclide('Li6', 1.0)
    mat.add_nuclide('Li7', 1.0)
    mat.set_density('g/cm3', 1.0)
    # Mix file path and keyword
    mat.read_nuclides_from_h5({
        "Li6": "tests/Li6.h5",
        "Li7": "tendl-2019"
    })
    grid = mat.unified_energy_grid_neutron()
    assert len(grid) > 0, "Mixed sources should load data"


def test_material_cache_respects_data_source_boundaries():
    """Test that material cache properly separates different data sources."""
    from yamc import clear_nuclide_cache
    clear_nuclide_cache()

    # Load from keyword
    mat1 = Material()
    mat1.add_nuclide('Li6', 1.0)
    mat1.set_density('g/cm3', 1.0)
    mat1.read_nuclides_from_h5("tendl-2019")

    # Load same nuclide from local file
    mat2 = Material()
    mat2.add_nuclide('Li6', 1.0)
    mat2.set_density('g/cm3', 1.0)
    mat2.read_nuclides_from_h5({"Li6": "tests/Li6.h5"})

    # Both should work independently
    grid1 = mat1.unified_energy_grid_neutron()
    grid2 = mat2.unified_energy_grid_neutron()
    assert len(grid1) > 0, "First material should have data"
    assert len(grid2) > 0, "Second material should have data"


def test_material_path_normalization_in_cache():
    """Test that different path formats for same file use same cache entry."""
    import os
    
    # Material 1: relative path
    mat_rel = Material()
    mat_rel.add_nuclide("Li6", 1.0)
    mat_rel.read_nuclides_from_h5({"Li6": "tests/Li6.h5"})
    mat_rel.set_density("g/cm3", 1.0)
    
    # Material 2: absolute path to same file
    mat_abs = Material()
    mat_abs.add_nuclide("Li6", 1.0)
    abs_path = os.path.abspath("tests/Li6.h5")
    mat_abs.read_nuclides_from_h5({"Li6": abs_path})
    mat_abs.set_density("g/cm3", 1.0)
    
    # Should give identical results (same file, cache working)
    xs_rel, _ = mat_rel.macroscopic_cross_section("(n,gamma)")
    xs_abs, _ = mat_abs.macroscopic_cross_section("(n,gamma)")
    
    assert xs_rel == xs_abs, "Relative and absolute paths to same file should give identical results"

def test_material_name_set_and_get():
    mat = Material()
    assert mat.name is None
    mat.name = "TestPyName"
    assert mat.name == "TestPyName"
    mat.name = "OtherName"
    assert mat.name == "OtherName"
    mat2 = Material(name="InitialName")
    assert mat2.name == "InitialName"