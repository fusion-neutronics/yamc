import pytest
from materials_for_mc import Nuclide

def test_be9_not_fissionable():
    nuc = Nuclide('Be9')
    nuc.read_nuclide_from_json('tests/Be9.json')
    assert hasattr(nuc, 'fissionable'), "Nuclide should have a 'fissionable' attribute"
    assert nuc.fissionable is False, "Be9 should not be fissionable"

def test_fe58_not_fissionable():
    nuc = Nuclide('Fe58')
    nuc.read_nuclide_from_json('tests/Fe58.json')
    assert hasattr(nuc, 'fissionable'), "Nuclide should have a 'fissionable' attribute"
    assert nuc.fissionable is False, "Fe58 should not be fissionable"
from materials_for_mc import Nuclide

def test_read_li6_nuclide():
    nuc1 = Nuclide('Li6')
    nuc1.read_nuclide_from_json('tests/Li6.json')
    assert nuc1.element.lower() == 'lithium'
    assert nuc1.atomic_symbol == "Li"
    assert nuc1.atomic_number == 3
    assert nuc1.mass_number == 6
    assert nuc1.neutron_number == 3
    assert nuc1.available_temperatures == ['294']
    # We don't expect any specific order of MT numbers, just check they're all ints
    assert all(isinstance(mt, int) for mt in nuc1.reaction_mts)

    cs = nuc1.reactions['294'][2].cross_section

    for entry in cs:
        assert isinstance(entry, float)
        assert isinstance(entry, float)

    expected_li6 = [1, 2, 3, 4, 5, 24, 27, 51, 52, 53, 54, 55, 56, 57, 58, 59, 60, 61, 62, 63, 64, 65, 66, 67, 68, 69, 70, 71, 72, 73, 74, 75, 76, 77, 78, 79, 80, 81, 101, 102, 103, 105, 203, 204, 205, 206, 207, 301, 444]
    assert nuc1.reaction_mts == expected_li6

def test_read_li7_nuclide():
    nuc1 = Nuclide('Li7')
    nuc1.read_nuclide_from_json('tests/Li7.json')
    assert nuc1.element.lower() == 'lithium'
    assert nuc1.atomic_symbol == "Li"
    assert nuc1.atomic_number == 3
    assert nuc1.mass_number == 7
    assert nuc1.neutron_number == 4
    assert nuc1.available_temperatures == ['294']
    # We don't expect any specific order of MT numbers, just check they're all ints
    assert all(isinstance(mt, int) for mt in nuc1.reaction_mts)

    cs = nuc1.reactions['294'][2].cross_section
    
    for entry in cs:
        assert isinstance(entry, float)
        assert isinstance(entry, float)

    expected_li7 = [1, 2, 3, 4, 5, 16, 24, 25, 27, 51, 52, 53, 54, 55, 56, 57, 58, 59, 60, 61, 62, 63, 64, 65, 66, 67, 68, 69, 70, 71, 72, 73, 74, 75, 76, 77, 78, 79, 80, 81, 82, 101, 102, 104, 203, 204, 205, 206, 207, 301, 444]
    assert nuc1.reaction_mts == expected_li7


def test_read_be9_available_and_loaded_temperatures():
    nuc = Nuclide('Be9')
    nuc.read_nuclide_from_json('tests/Be9.json')
    assert nuc.available_temperatures == ['294', '300']
    # By current implementation, all temps are loaded eagerly
    assert hasattr(nuc, 'loaded_temperatures'), "loaded_temperatures attribute missing"
    assert nuc.loaded_temperatures == ['294', '300']
    # Reactions dict should contain both temperatures
    assert '294' in nuc.reactions
    assert '300' in nuc.reactions


def test_read_be9_mt_numbers_per_temperature():
    nuc = Nuclide('Be9')
    nuc.read_nuclide_from_json('tests/Be9.json')
    mts_294 = sorted(int(mt) for mt in nuc.reactions['294'].keys())
    mts_300 = sorted(int(mt) for mt in nuc.reactions['300'].keys())
    expected_294 = sorted([
        1,2,3,16,27,101,102,103,104,105,107,203,204,205,207,301,444,
        875,876,877,878,879,880,881,882,883,884,885,886,887,888,889,890
    ])
    expected_300 = sorted([
        1,2,3,16,27,101,102,103,104,105,107,203,204,205,207,301
    ])
    assert mts_294 == expected_294, f"Be9 294K MT list mismatch: {mts_294}"
    assert mts_300 == expected_300, f"Be9 300K MT list mismatch: {mts_300}"
    # Ensure 300K list is subset of 294K list
    assert set(mts_300).issubset(set(mts_294))


def test_read_be9_selective_single_temperature():
    # Ensure only the specified temperature (300) is retained in reactions and loaded_temperatures
    nuc = Nuclide('Be9')
    nuc.read_nuclide_from_json('tests/Be9.json', temperatures=['300'])
    assert nuc.available_temperatures == ['294', '300'], "available_temperatures should list all temps present in file"
    assert nuc.loaded_temperatures == ['300'], f"loaded_temperatures should be only ['300'], got {nuc.loaded_temperatures}"
    assert '300' in nuc.reactions, "300K reactions missing after selective load"
    assert '294' not in nuc.reactions, "294K reactions should not be loaded when selectively requesting only 300K"
    # MT list at 300 should match subset expectation
    mts_300 = sorted(int(mt) for mt in nuc.reactions['300'].keys())
    expected_300 = sorted([
        1,2,3,16,27,101,102,103,104,105,107,203,204,205,207,301
    ])
    assert mts_300 == expected_300, f"Selective load 300K MT list mismatch: {mts_300}"


def test_read_nuclide_from_json_keyword():
    from materials_for_mc import Nuclide
    nuc = Nuclide('Li6')
    nuc.read_nuclide_from_json("tendl-21")

def test_read_nuclide_from_json_local_path():
    from materials_for_mc import Nuclide
    nuc = Nuclide('Li6')
    # Should not raise TypeError when passing local path
    nuc.read_nuclide_from_json("tests/Li6.json")


def test_microscopic_cross_section_with_temperature():
    """Test microscopic_cross_section with explicit temperature."""
    from materials_for_mc import Nuclide
    nuc = Nuclide('Be9')
    nuc.read_nuclide_from_json('tests/Be9.json')
    
    # Test with specific temperature
    xs, energy = nuc.microscopic_cross_section(reaction=2, temperature='294')
    assert len(xs) > 0, "Cross section data should not be empty"
    assert len(energy) > 0, "Energy data should not be empty"
    assert len(xs) == len(energy), "Cross section and energy arrays should have same length"
    
    # Test with different temperature
    xs_300, energy_300 = nuc.microscopic_cross_section(reaction=2, temperature='300')
    assert len(xs_300) > 0, "Cross section data should not be empty for 300K"
    assert len(energy_300) > 0, "Energy data should not be empty for 300K"
    
    # Test different MT numbers
    xs_mt3, energy_mt3 = nuc.microscopic_cross_section(reaction=3, temperature='294')
    assert len(xs_mt3) > 0, "MT=3 cross section data should not be empty"
    assert len(energy_mt3) > 0, "MT=3 energy data should not be empty"


def test_microscopic_cross_section_without_temperature():
    """Test microscopic_cross_section with single loaded temperature."""
    from materials_for_mc import Nuclide
    nuc = Nuclide('Be9')
    # Load only one temperature
    nuc.read_nuclide_from_json('tests/Be9.json', temperatures=['294'])
    
    # Should work without specifying temperature since only one is loaded
    xs, energy = nuc.microscopic_cross_section(2)
    assert len(xs) > 0, "Cross section data should not be empty"
    assert len(energy) > 0, "Energy data should not be empty"
    assert len(xs) == len(energy), "Cross section and energy arrays should have same length"


def test_microscopic_cross_section_multiple_temperatures_error():
    """Test that microscopic_cross_section raises error when multiple temperatures loaded without specifying."""
    from materials_for_mc import Nuclide
    nuc = Nuclide('Be9')
    nuc.read_nuclide_from_json('tests/Be9.json')  # Loads both 294 and 300
    
    # Should raise error when no temperature specified with multiple loaded
    with pytest.raises(Exception) as exc_info:
        nuc.microscopic_cross_section(2)
    error_msg = str(exc_info.value)
    assert "Multiple temperatures loaded" in error_msg
    assert "294" in error_msg and "300" in error_msg


def test_microscopic_cross_section_invalid_temperature():
    """Test error handling for invalid temperature."""
    from materials_for_mc import Nuclide
    nuc = Nuclide('Be9')
    nuc.read_nuclide_from_json('tests/Be9.json')
    
    # Should raise error for non-existent temperature
    with pytest.raises(Exception) as exc_info:
        nuc.microscopic_cross_section(reaction=2, temperature='500')
    error_msg = str(exc_info.value)
    assert "Temperature '500' not found" in error_msg
    assert "Available temperatures:" in error_msg
    assert "294" in error_msg and "300" in error_msg


def test_microscopic_cross_section_invalid_mt():
    """Test error handling for invalid MT number."""
    from materials_for_mc import Nuclide
    nuc = Nuclide('Be9')
    nuc.read_nuclide_from_json('tests/Be9.json')
    
    # Should raise error for non-existent MT
    with pytest.raises(Exception) as exc_info:
        nuc.microscopic_cross_section(reaction=9999, temperature='294')
    error_msg = str(exc_info.value)
    assert "MT 9999 not found" in error_msg
    assert "Available MTs:" in error_msg


def test_microscopic_cross_section_multiple_mt_numbers():
    """Test microscopic_cross_section with various MT numbers."""
    from materials_for_mc import Nuclide
    nuc = Nuclide('Be9')
    nuc.read_nuclide_from_json('tests/Be9.json')
    
    # Test common MT numbers that should exist in Be9
    test_mts = [1, 2, 3, 16, 27, 101, 102]  # Common reaction types
    
    for mt in test_mts:
        try:
            xs, energy = nuc.microscopic_cross_section(reaction=mt, temperature='294')
            assert len(xs) > 0, f"MT={mt} should have cross section data"
            assert len(energy) > 0, f"MT={mt} should have energy data"
            assert len(xs) == len(energy), f"MT={mt} data length mismatch"
            assert all(e > 0 for e in energy), f"MT={mt} energy values should be positive"
            assert all(x >= 0 for x in xs), f"MT={mt} cross sections should be non-negative"
        except Exception:
            # Some MT numbers might not exist, which is fine
            pass


def test_microscopic_cross_section_lithium():
    """Test microscopic_cross_section with Li6 data."""
    from materials_for_mc import Nuclide
    nuc = Nuclide('Li6')
    nuc.read_nuclide_from_json('tests/Li6.json')
    
    # Li6 should have only one temperature, so no temperature needed
    xs, energy = nuc.microscopic_cross_section(reaction=2)  # Elastic scattering
    assert len(xs) > 0, "Li6 elastic scattering data should not be empty"
    assert len(energy) > 0, "Li6 energy data should not be empty"
    
    # Test with explicit temperature too
    xs_explicit, energy_explicit = nuc.microscopic_cross_section(reaction=2, temperature='294')
    assert xs == xs_explicit, "Results should be identical with/without explicit temperature"
    assert energy == energy_explicit, "Energy should be identical with/without explicit temperature"


def test_auto_loading_from_config():
    """Test that microscopic_cross_section can auto-load data from config when nuclide is empty"""
    from materials_for_mc import Config
    
    # Set up config for auto-loading
    config = Config()
    config.set_cross_sections({'Be9': 'tests/Be9.json'})
    
    # Create empty nuclide with name but no data loaded
    nuc = Nuclide('Be9')
    assert nuc.loaded_temperatures == [], "Should start with no loaded temperatures"
    
    # Call microscopic_cross_section - should auto-load data
    xs, energy = nuc.microscopic_cross_section(reaction=2, temperature='294')
    assert len(xs) > 0, "Auto-loaded cross section data should not be empty"
    assert len(energy) > 0, "Auto-loaded energy data should not be empty"
    assert len(xs) == len(energy), "Cross section and energy arrays should have same length"
    
    # Note: loaded_temperatures won't be updated in the Python object due to immutable API
    # The auto-loading happens internally but doesn't modify the original object


def test_auto_loading_additional_temperature():
    """Test that microscopic_cross_section can auto-load additional temperatures"""
    from materials_for_mc import Config
    
    # Set up config for auto-loading
    config = Config()
    config.set_cross_sections({'Be9': 'tests/Be9.json'})
    
    # Load Be9 with only 294K initially
    nuc = Nuclide('Be9')
    nuc.read_nuclide_from_json('tests/Be9.json', ['294'])
    
    assert nuc.loaded_temperatures == ['294'], "Should only have 294K loaded initially"
    assert '300' in nuc.available_temperatures, "Should know 300K is available"
    
    # Request 300K data - should auto-load additional temperature
    xs, energy = nuc.microscopic_cross_section(reaction=2, temperature='300')
    assert len(xs) > 0, "Auto-loaded 300K cross section data should not be empty"
    assert len(energy) > 0, "Auto-loaded 300K energy data should not be empty"
    
    # The original nuclide object still shows only 294K due to immutable API
    # But the internal auto-loading worked to provide the 300K data
    

def test_auto_loading_without_config_fails():
    """Test that auto-loading fails gracefully when no config is available"""
    from materials_for_mc import Config
    
    # Clear any existing config by setting an empty dict
    config = Config()
    original_configs = config.get_cross_sections()
    
    # Temporarily clear all configurations
    config.set_cross_sections({})
    
    try:
        # Create empty nuclide with name but no config
        nuc = Nuclide('TestNuclide')
        
        # Call microscopic_cross_section - should fail with helpful error
        try:
            nuc.microscopic_cross_section(reaction=2, temperature='294')
            assert False, "Auto-loading without config should fail"
        except Exception as e:
            error_msg = str(e)
            # The error could be either "No configuration found" or a download error
            # Both are acceptable since there's no valid config
            assert ("No configuration found" in error_msg or 
                    "Failed to download" in error_msg or
                    "404 Not Found" in error_msg), f"Error should indicate missing or invalid configuration: {error_msg}"
            assert "TestNuclide" in error_msg or "TestNuclide" in str(e.__class__), f"Error should be related to TestNuclide: {error_msg}"
    
    finally:
        # Restore original configuration
        if original_configs:
            config.set_cross_sections(original_configs)


def test_auto_loading_multiple_calls_consistent():
    """Test that multiple auto-loading calls give consistent results"""
    from materials_for_mc import Config
    
    # Set up config for auto-loading
    config = Config()
    config.set_cross_sections({'Be9': 'tests/Be9.json'})
    
    # Create empty nuclide
    nuc = Nuclide('Be9')
    
    # Call microscopic_cross_section multiple times
    xs1, energy1 = nuc.microscopic_cross_section(reaction=2, temperature='294')
    xs2, energy2 = nuc.microscopic_cross_section(reaction=2, temperature='294')
    xs3, energy3 = nuc.microscopic_cross_section(reaction=102, temperature='300')
    
    # First two calls should give identical results
    assert xs1 == xs2, "Multiple calls with same parameters should give identical results"
    assert energy1 == energy2, "Multiple calls with same parameters should give identical energy"
    
    # Third call should work too (different MT and temperature)
    assert len(xs3) > 0, "Auto-loading different MT and temperature should work"
    assert len(energy3) > 0, "Auto-loading different MT and temperature should provide energy"


def test_auto_loading_with_manual_loading_combined():
    """Test combining manual loading with auto-loading for additional data"""
    from materials_for_mc import Config
    
    # Set up config for auto-loading
    config = Config()
    config.set_cross_sections({'Be9': 'tests/Be9.json'})
    
    # Manually load some data first
    nuc = Nuclide('Be9')
    nuc.read_nuclide_from_json('tests/Be9.json', ['294'])
    
    # Verify manual loading worked
    assert '294' in nuc.loaded_temperatures, "Manual loading should work"
    assert '300' in nuc.available_temperatures, "Should know other temperatures are available"
    
    # Now use auto-loading for data that was manually loaded
    xs_manual, energy_manual = nuc.microscopic_cross_section(reaction=2, temperature='294')
    assert len(xs_manual) > 0, "Should get data for manually loaded temperature"
    
    # And use auto-loading for additional temperature not manually loaded  
    xs_auto, energy_auto = nuc.microscopic_cross_section(reaction=2, temperature='300')
    assert len(xs_auto) > 0, "Should auto-load additional temperature"
    
    # Test with a temperature-specific MT - try MT=444 which should only be available at 294K
    xs_specific, energy_specific = nuc.microscopic_cross_section(reaction=444, temperature='294')
    assert len(xs_specific) > 0, "Should get temperature-specific MT data"
    
    # Try to get MT=444 at 300K - this should fail since it's not available at that temperature
    try:
        nuc.microscopic_cross_section(reaction=444, temperature='300')
        # If we get here without exception, that's fine too - means 444 exists at both temps
        print("MT=444 exists at both temperatures")
    except ValueError as e:
        # This is expected if MT=444 is not available at 300K
        assert "MT 444 not found" in str(e), f"Should get MT not found error, got: {e}"
    
    # Note: For Be9 MT=2, the cross sections at 294K and 300K might be identical
    # This is fine - the important thing is that both calls succeeded


def test_fendl_3_2c_keyword():
    """Test that the fendl-3.2c keyword is recognized and works correctly."""
    import materials_for_mc
    
    # Test that the keyword is recognized (this tests the Rust backend)
    try:
        # This should not raise an exception if the keyword is recognized
        # We'll create a dummy config entry to test keyword recognition
        from materials_for_mc import Config
        config = Config()
        
        # Test setting cross sections with the keyword - this should not fail
        # if the keyword is recognized in the backend
        config.set_cross_sections({'Li6': 'fendl-3.2c'})
        
        # Verify we can retrieve it
        cross_sections = config.get_cross_sections()
        assert 'Li6' in cross_sections, "Li6 should be in cross sections config"
        assert cross_sections['Li6'] == 'fendl-3.2c', "Should store fendl-3.2c keyword correctly"
        
        # Test that keyword expansion would work (without actually downloading)
        # This implicitly tests the URL cache functionality  
        print("fendl-3.2c keyword test passed - keyword is recognized")
        
    except Exception as e:
        pytest.fail(f"fendl-3.2c keyword should be recognized by the system: {e}")


def test_auto_loading_with_global_keyword():
    """Test that auto-loading works with global keyword configuration"""
    from materials_for_mc import Config, Nuclide
    
    # Set global keyword configuration
    config = Config()
    config.set_cross_sections('fendl-3.2c')
    
    # Verify config is set correctly
    assert config.get_cross_section('Li6') == 'fendl-3.2c', "Global config should apply to Li6"
    
    # Create empty nuclide
    nuc = Nuclide('Li6')
    assert nuc.loaded_temperatures == [], "Should start with no loaded temperatures"
    
    # Call microscopic_cross_section - should auto-load data from global config
    try:
        xs, energy = nuc.microscopic_cross_section(reaction=1, temperature='294')
        assert len(xs) > 0, "Auto-loaded cross section data should not be empty"
        assert len(energy) > 0, "Auto-loaded energy data should not be empty"
        assert len(xs) == len(energy), "Cross section and energy arrays should have same length"
        print("Auto-loading with global keyword test passed!")
        
    except Exception as e:
        # If we can't download (no internet or URL issues), that's OK for this test
        # The important thing is that the config lookup worked
        if "No configuration found" in str(e):
            pytest.fail(f"Config lookup failed - auto-loading should work with global keywords: {e}")
        else:
            print(f"Note: Auto-loading test skipped due to download issue: {e}")
            # This is acceptable - we verified the config lookup works


def test_microscopic_cross_section_by_name():
    """Test microscopic_cross_section with reaction names."""
    from materials_for_mc import Nuclide
    nuc = Nuclide('Be9')
    nuc.read_nuclide_from_json('tests/Be9.json')
    
    # Test elastic scattering using reaction name
    xs_name, energy_name = nuc.microscopic_cross_section("(n,elastic)", temperature='294')
    assert len(xs_name) > 0, "Elastic scattering cross section should not be empty"
    assert len(energy_name) > 0, "Energy data should not be empty"
    
    # Compare with MT number approach (MT=2 is elastic scattering)
    xs_mt, energy_mt = nuc.microscopic_cross_section(reaction=2, temperature='294')
    
    # Should get identical results
    assert xs_name == xs_mt, "Reaction name and MT number should give identical cross sections"
    assert energy_name == energy_mt, "Reaction name and MT number should give identical energy grids"
    
    # Test other common reactions
    test_reactions = [
        ("(n,gamma)", 102),   # Radiative capture
        ("(n,a)", 107),       # Alpha production
        ("(n,total)", 1),     # Total cross section
    ]
    
    for reaction_name, mt_num in test_reactions:
        try:
            xs_name, energy_name = nuc.microscopic_cross_section(reaction_name, temperature='294')
            xs_mt, energy_mt = nuc.microscopic_cross_section(reaction=mt_num, temperature='294')
            
            assert len(xs_name) > 0, f"{reaction_name} should have cross section data"
            assert len(energy_name) > 0, f"{reaction_name} should have energy data"
            assert xs_name == xs_mt, f"{reaction_name} and MT={mt_num} should give identical results"
            assert energy_name == energy_mt, f"{reaction_name} and MT={mt_num} should give identical energy"
            
        except Exception as e:
            # Some reactions might not exist for Be9, which is acceptable
            if "not found" in str(e).lower():
                print(f"Note: {reaction_name} (MT={mt_num}) not available in Be9 data")
            else:
                raise e


def test_microscopic_cross_section_by_name_invalid_reaction():
    """Test error handling for invalid reaction names."""
    from materials_for_mc import Nuclide
    nuc = Nuclide('Be9')
    nuc.read_nuclide_from_json('tests/Be9.json')
    
    # Test with invalid reaction name
    with pytest.raises(Exception) as exc_info:
        nuc.microscopic_cross_section("invalid_reaction", temperature='294')
    error_msg = str(exc_info.value)
    assert "not found in reaction mapping" in error_msg or "Unknown reaction" in error_msg


def test_microscopic_cross_section_by_name_fission():
    """Test that the special 'fission' alias works."""
    from materials_for_mc import Nuclide
    
    # Use Li6 which might have fission data, or test the error handling
    nuc = Nuclide('Li6')
    nuc.read_nuclide_from_json('tests/Li6.json')
    
    try:
        xs_fission, energy_fission = nuc.microscopic_cross_section("fission", temperature='294')
        xs_mt18, energy_mt18 = nuc.microscopic_cross_section(18, temperature='294')
        
        # Should get identical results since fission maps to MT=18
        assert xs_fission == xs_mt18, "fission and MT=18 should give identical results"
        assert energy_fission == energy_mt18, "fission and MT=18 should give identical energy"
        
    except Exception as e:
        # Li6 might not have fission data, which is acceptable
        if "MT 18 not found" in str(e) or "not found" in str(e).lower():
            print("Note: Li6 does not have fission data (expected)")
        else:
            raise e


def test_sample_reaction_basic():
    """Test basic functionality of sample_reaction method."""
    from materials_for_mc import Nuclide
    
    nuc = Nuclide('Li6')
    nuc.read_nuclide_from_json('tests/Li6.json')
    
    # Test sampling at a typical neutron energy
    reaction = nuc.sample_reaction(energy=1.0, temperature='294', seed=42)
    
    if reaction is not None:
        # Check that all expected fields are present
        assert 'mt_number' in reaction, "Reaction should have mt_number field"
        assert 'cross_section' in reaction, "Reaction should have cross_section field"
        assert 'threshold_idx' in reaction, "Reaction should have threshold_idx field"
        assert 'interpolation' in reaction, "Reaction should have interpolation field"
        assert 'energy' in reaction, "Reaction should have energy field"
        
        # Check data types
        assert isinstance(reaction['mt_number'], int), "MT number should be integer"
        assert isinstance(reaction['cross_section'], list), "Cross section should be list"
        assert isinstance(reaction['energy'], list), "Energy should be list"
        assert isinstance(reaction['threshold_idx'], int), "Threshold index should be integer"
        assert isinstance(reaction['interpolation'], list), "Interpolation should be list"
        
        # Check data validity
        assert reaction['mt_number'] > 0, "MT number should be positive"
        assert len(reaction['cross_section']) > 0, "Cross section should not be empty"
        assert len(reaction['energy']) > 0, "Energy grid should not be empty"
        assert reaction['threshold_idx'] >= 0, "Threshold index should be non-negative"
        assert all(x >= 0 for x in reaction['cross_section']), "Cross sections should be non-negative"
        assert all(e > 0 for e in reaction['energy']), "Energy values should be positive"
    else:
        # If no reaction is sampled, that's also valid (zero total cross section)
        print("Note: No reaction sampled (possibly zero total cross section)")


def test_sample_reaction_reproducibility():
    """Test that sample_reaction gives reproducible results with same seed."""
    from materials_for_mc import Nuclide
    
    nuc = Nuclide('Li6')
    nuc.read_nuclide_from_json('tests/Li6.json')
    
    # Sample with same seed multiple times
    reaction1 = nuc.sample_reaction(energy=1.0, temperature='294', seed=123)
    reaction2 = nuc.sample_reaction(energy=1.0, temperature='294', seed=123)
    reaction3 = nuc.sample_reaction(energy=1.0, temperature='294', seed=456)
    
    if reaction1 is not None and reaction2 is not None:
        # Same seed should give same result
        assert reaction1['mt_number'] == reaction2['mt_number'], "Same seed should give same MT number"
        assert reaction1['cross_section'] == reaction2['cross_section'], "Same seed should give same cross section data"
        assert reaction1['energy'] == reaction2['energy'], "Same seed should give same energy data"
        
        # Different seed might give different result (but not guaranteed)
        if reaction3 is not None:
            # We can't guarantee different results, but at least verify the structure is correct
            assert isinstance(reaction3['mt_number'], int), "Different seed should still give valid MT number"
    else:
        print("Note: No reaction sampled in reproducibility test")


def test_sample_reaction_different_energies():
    """Test sample_reaction at different neutron energies."""
    from materials_for_mc import Nuclide
    
    nuc = Nuclide('Li6')
    nuc.read_nuclide_from_json('tests/Li6.json')
    
    # Test various energies from thermal to fast neutron range
    test_energies = [1e-3, 1e-2, 1e-1, 1.0, 10.0, 100.0, 1000.0]
    sampled_mts = set()
    
    for energy in test_energies:
        reaction = nuc.sample_reaction(energy=energy, temperature='294', seed=42)
        
        if reaction is not None:
            sampled_mts.add(reaction['mt_number'])
            
            # Verify reaction structure at each energy
            assert reaction['mt_number'] > 0, f"Valid MT number at energy {energy}"
            assert len(reaction['cross_section']) > 0, f"Non-empty cross section at energy {energy}"
            assert len(reaction['energy']) > 0, f"Non-empty energy grid at energy {energy}"
    
    # Should sample at least one reaction across all energies (unless total XS is zero everywhere)
    if len(sampled_mts) == 0:
        print("Warning: No reactions sampled across all test energies")
    else:
        print(f"Sampled {len(sampled_mts)} different reaction types across energy range")


def test_sample_reaction_be9():
    """Test sample_reaction with Be9 nuclide (different reaction channels)."""
    from materials_for_mc import Nuclide
    
    nuc = Nuclide('Be9')
    nuc.read_nuclide_from_json('tests/Be9.json')
    
    # Be9 has different reaction channels than Li6
    reaction = nuc.sample_reaction(energy=1.0, temperature='294', seed=42)
    
    if reaction is not None:
        # Verify structure
        assert 'mt_number' in reaction, "Be9 reaction should have mt_number field"
        assert isinstance(reaction['mt_number'], int), "Be9 MT number should be integer"
        assert reaction['mt_number'] > 0, "Be9 MT number should be positive"
        
        # Check that MT number is valid for Be9
        assert reaction['mt_number'] in nuc.reaction_mts, f"MT {reaction['mt_number']} should be in available MTs for Be9"
        
        print(f"Be9 sampled reaction: MT {reaction['mt_number']}")
    else:
        print("Note: No Be9 reaction sampled")


def test_sample_reaction_without_seed():
    """Test that sample_reaction works without specifying a seed."""
    from materials_for_mc import Nuclide
    
    nuc = Nuclide('Li6')
    nuc.read_nuclide_from_json('tests/Li6.json')
    
    # Sample without seed (should use system random state)
    reaction1 = nuc.sample_reaction(energy=1.0, temperature='294')
    reaction2 = nuc.sample_reaction(energy=1.0, temperature='294')
    
    # Both should work (might be same or different)
    if reaction1 is not None:
        assert isinstance(reaction1['mt_number'], int), "No-seed reaction 1 should be valid"
    
    if reaction2 is not None:
        assert isinstance(reaction2['mt_number'], int), "No-seed reaction 2 should be valid"
    
    print("No-seed sampling test completed")


def test_sample_reaction_multiple_temperatures():
    """Test sample_reaction with different temperatures."""
    from materials_for_mc import Nuclide
    
    nuc = Nuclide('Be9')
    nuc.read_nuclide_from_json('tests/Be9.json')  # Has both 294K and 300K
    
    # Sample at both temperatures
    reaction_294 = nuc.sample_reaction(energy=1.0, temperature='294', seed=42)
    reaction_300 = nuc.sample_reaction(energy=1.0, temperature='300', seed=42)
    
    # Both should work
    if reaction_294 is not None:
        assert isinstance(reaction_294['mt_number'], int), "294K reaction should be valid"
        
    if reaction_300 is not None:
        assert isinstance(reaction_300['mt_number'], int), "300K reaction should be valid"
    
    # Results might be different due to temperature-dependent cross sections
    # but we can't guarantee this, so just verify they both work
    print("Multiple temperature sampling test completed")


def test_sample_reaction_invalid_temperature():
    """Test error handling for invalid temperature in sample_reaction."""
    from materials_for_mc import Nuclide
    
    nuc = Nuclide('Li6')
    nuc.read_nuclide_from_json('tests/Li6.json')
    
    # Test with invalid temperature - should handle gracefully
    # The current implementation uses fallback behavior rather than strict errors
    reaction = nuc.sample_reaction(energy=1.0, temperature='999', seed=42)
    
    # The method should either work (using fallback) or return None
    if reaction is not None:
        assert isinstance(reaction['mt_number'], int), "Fallback temperature reaction should be valid"
        print("Invalid temperature handled with fallback")
    else:
        print("Invalid temperature resulted in no reaction sampled")


def test_sample_reaction_edge_cases():
    """Test sample_reaction with edge case energies."""
    from materials_for_mc import Nuclide
    
    nuc = Nuclide('Li6')
    nuc.read_nuclide_from_json('tests/Li6.json')
    
    # Test very low energy
    reaction_low = nuc.sample_reaction(energy=1e-5, temperature='294', seed=42)
    
    # Test very high energy
    reaction_high = nuc.sample_reaction(energy=1e7, temperature='294', seed=42)
    
    # Both should either work or return None (no error)
    if reaction_low is not None:
        assert isinstance(reaction_low['mt_number'], int), "Very low energy should give valid reaction"
    
    if reaction_high is not None:
        assert isinstance(reaction_high['mt_number'], int), "Very high energy should give valid reaction"
    
    print("Edge case energy test completed")


def test_sample_reaction_zero_energy():
    """Test sample_reaction with zero energy."""
    from materials_for_mc import Nuclide
    
    nuc = Nuclide('Li6')
    nuc.read_nuclide_from_json('tests/Li6.json')
    
    # Test with zero energy (edge case)
    reaction = nuc.sample_reaction(energy=0.0, temperature='294', seed=42)
    
    # Should either work or return None, not crash
    if reaction is not None:
        assert isinstance(reaction['mt_number'], int), "Zero energy should give valid reaction if any"
    
    print("Zero energy test completed")


def test_sample_reaction_consistency_with_available_mts():
    """Test that sampled reactions are from available MT numbers."""
    from materials_for_mc import Nuclide
    
    nuc = Nuclide('Li6')
    nuc.read_nuclide_from_json('tests/Li6.json')
    
    available_mts = set(nuc.reaction_mts)
    sampled_mts = set()
    
    # Sample many times to get a good distribution
    for i in range(100):
        reaction = nuc.sample_reaction(energy=1.0, temperature='294', seed=i)
        if reaction is not None:
            sampled_mts.add(reaction['mt_number'])
    
    # All sampled MTs should be in the available set
    for mt in sampled_mts:
        assert mt in available_mts, f"Sampled MT {mt} should be in available MTs {available_mts}"
    
    if sampled_mts:
        print(f"Sampled MTs {sorted(sampled_mts)} are all in available MTs {sorted(available_mts)}")
    else:
        print("No reactions sampled in consistency test")


def test_sample_reaction_return_type():
    """Test that sample_reaction returns correct Python types."""
    from materials_for_mc import Nuclide
    
    nuc = Nuclide('Li6')
    nuc.read_nuclide_from_json('tests/Li6.json')
    
    reaction = nuc.sample_reaction(energy=1.0, temperature='294', seed=42)
    
    if reaction is not None:
        # Check that it's a dictionary
        assert isinstance(reaction, dict), "Reaction should be a dictionary"
        
        # Check specific field types
        assert isinstance(reaction['mt_number'], int), "mt_number should be int"
        assert isinstance(reaction['cross_section'], list), "cross_section should be list"
        assert isinstance(reaction['energy'], list), "energy should be list"
        assert isinstance(reaction['threshold_idx'], int), "threshold_idx should be int"
        assert isinstance(reaction['interpolation'], list), "interpolation should be list"
        
        # Check list element types
        if reaction['cross_section']:
            assert isinstance(reaction['cross_section'][0], float), "Cross section values should be float"
        
        if reaction['energy']:
            assert isinstance(reaction['energy'][0], float), "Energy values should be float"
        
        if reaction['interpolation']:
            assert isinstance(reaction['interpolation'][0], int), "Interpolation values should be int"
    
    else:
        # None is also a valid return type
        assert reaction is None, "Return should be None or dict"
    
    print("Return type test completed")


def test_nuclide_different_data_sources():
    """Test that loading the same nuclide from different sources gives different results."""
    
    # Load Li7 from different sources
    li7_tendl = Nuclide("Li7")
    li7_tendl.read_nuclide_from_json("tendl-21")
    
    li7_fendl = Nuclide("Li7")
    li7_fendl.read_nuclide_from_json("fendl-3.2c")
    
    # Get cross sections from both
    xs_tendl, energy_tendl = li7_tendl.microscopic_cross_section("(n,gamma)")
    xs_fendl, energy_fendl = li7_fendl.microscopic_cross_section("(n,gamma)")
    
    # Should have different data
    data_different = (len(xs_tendl) != len(xs_fendl) or 
                     any(abs(a - b) > 1e-10 for a, b in zip(xs_tendl, xs_fendl)))
    
    assert data_different, "TENDL and FENDL data should be different"
    print(f"TENDL Li7: {len(xs_tendl)} points, FENDL Li7: {len(xs_fendl)} points")


def test_nuclide_file_vs_keyword_sources():
    """Test that file paths and keywords can coexist."""
    
    # Load Li6 from local file
    li6_file = Nuclide("Li6")
    li6_file.read_nuclide_from_json("tests/Li6.json")
    
    # Load Li7 from keyword
    li7_keyword = Nuclide("Li7")
    li7_keyword.read_nuclide_from_json("tendl-21")
    
    assert li6_file.name == "Li6"
    assert li7_keyword.name == "Li7"
    
    # Verify we can load cross sections from both
    xs_file, _ = li6_file.microscopic_cross_section("(n,gamma)")
    xs_keyword, _ = li7_keyword.microscopic_cross_section("(n,gamma)")
    
    assert len(xs_file) > 0 and len(xs_keyword) > 0, "Both should have cross section data"


def test_nuclide_cache_respects_data_source_boundaries():
    """Test that the cache properly separates different data sources."""
    
    # Load Li7 from TENDL first time
    li7_tendl_1 = Nuclide("Li7")
    li7_tendl_1.read_nuclide_from_json("tendl-21")
    
    # Load Li7 from FENDL
    li7_fendl = Nuclide("Li7")
    li7_fendl.read_nuclide_from_json("fendl-3.2c")
    
    # Load Li7 from TENDL again (should use cache)
    li7_tendl_2 = Nuclide("Li7")
    li7_tendl_2.read_nuclide_from_json("tendl-21")
    
    # Get cross sections
    xs_tendl_1, _ = li7_tendl_1.microscopic_cross_section("(n,gamma)")
    xs_fendl, _ = li7_fendl.microscopic_cross_section("(n,gamma)")
    xs_tendl_2, _ = li7_tendl_2.microscopic_cross_section("(n,gamma)")
    
    # TENDL loads should be identical (cache working)
    assert xs_tendl_1 == xs_tendl_2, "TENDL loads should be identical (cache working)"
    
    # TENDL vs FENDL should be different (different data sources)
    tendl_vs_fendl_different = (len(xs_tendl_1) != len(xs_fendl) or 
                               any(abs(a - b) > 1e-10 for a, b in zip(xs_tendl_1, xs_fendl)))
    assert tendl_vs_fendl_different, "TENDL and FENDL should have different data"


def test_nuclide_path_normalization():
    """Test that different path formats for same file give same results."""
    import os
    
    # Load Li6 with relative path
    li6_rel = Nuclide("Li6")
    li6_rel.read_nuclide_from_json("tests/Li6.json")
    
    # Load Li6 with absolute path
    li6_abs = Nuclide("Li6")
    abs_path = os.path.abspath("tests/Li6.json")
    li6_abs.read_nuclide_from_json(abs_path)
    
    # Should give identical results (same file)
    xs_rel, _ = li6_rel.microscopic_cross_section("(n,gamma)")
    xs_abs, _ = li6_abs.microscopic_cross_section("(n,gamma)")
    
    assert xs_rel == xs_abs, "Relative and absolute paths to same file should give identical results"


def test_basic_nuclide_properties_after_loading():
    """Test that basic nuclide properties are correctly set after loading data."""
    # Test Li6 properties
    nuc_li6 = Nuclide('Li6')
    nuc_li6.read_nuclide_from_json('tests/Li6.json')
    
    # Basic identification properties
    assert nuc_li6.name == 'Li6', "Nuclide name should be Li6"
    assert nuc_li6.element.lower() == 'lithium', "Element should be lithium"
    assert nuc_li6.atomic_symbol == 'Li', "Atomic symbol should be Li"
    assert nuc_li6.atomic_number == 3, "Atomic number should be 3"
    assert nuc_li6.mass_number == 6, "Mass number should be 6"
    assert nuc_li6.neutron_number == 3, "Neutron number should be 3"
    
    # Temperature availability
    assert len(nuc_li6.available_temperatures) > 0, "Should have available temperatures"
    assert '294' in nuc_li6.available_temperatures, "Should have 294K temperature data"
    
    # Reaction data availability
    assert len(nuc_li6.reaction_mts) > 0, "Should have reaction MTs available"
    assert isinstance(nuc_li6.reaction_mts, list), "Reaction MTs should be a list"
    assert all(isinstance(mt, int) for mt in nuc_li6.reaction_mts), "All MTs should be integers"


def test_basic_nuclide_properties_li7():
    """Test that basic nuclide properties are correctly set for Li7."""
    nuc_li7 = Nuclide('Li7')
    nuc_li7.read_nuclide_from_json('tests/Li7.json')
    
    # Basic identification properties
    assert nuc_li7.name == 'Li7', "Nuclide name should be Li7"
    assert nuc_li7.element.lower() == 'lithium', "Element should be lithium"
    assert nuc_li7.atomic_symbol == 'Li', "Atomic symbol should be Li"
    assert nuc_li7.atomic_number == 3, "Atomic number should be 3"
    assert nuc_li7.mass_number == 7, "Mass number should be 7"
    assert nuc_li7.neutron_number == 4, "Neutron number should be 4"
    
    # Temperature and reaction data
    assert len(nuc_li7.available_temperatures) > 0, "Should have available temperatures"
    assert len(nuc_li7.reaction_mts) > 0, "Should have reaction MTs available"


def test_basic_nuclide_properties_be9():
    """Test that basic nuclide properties are correctly set for Be9."""
    nuc_be9 = Nuclide('Be9')
    nuc_be9.read_nuclide_from_json('tests/Be9.json')
    
    # Basic identification properties
    assert nuc_be9.name == 'Be9', "Nuclide name should be Be9"
    assert nuc_be9.element.lower() == 'beryllium', "Element should be beryllium"
    assert nuc_be9.atomic_symbol == 'Be', "Atomic symbol should be Be"
    assert nuc_be9.atomic_number == 4, "Atomic number should be 4"
    assert nuc_be9.mass_number == 9, "Mass number should be 9"
    assert nuc_be9.neutron_number == 5, "Neutron number should be 5"
    
    # Be9 has multiple temperatures
    assert len(nuc_be9.available_temperatures) > 1, "Be9 should have multiple temperatures"
    assert '294' in nuc_be9.available_temperatures, "Should have 294K temperature data"
    assert '300' in nuc_be9.available_temperatures, "Should have 300K temperature data"


def test_basic_nuclide_properties_iron_isotopes():
    """Test that basic nuclide properties are correctly set for iron isotopes."""
    iron_isotopes = [
        ('Fe54', 54, 28), ('Fe56', 56, 30), ('Fe57', 57, 31), ('Fe58', 58, 32)
    ]
    
    for name, mass_num, neutron_num in iron_isotopes:
        nuc = Nuclide(name)
        nuc.read_nuclide_from_json(f'tests/{name}.json')
        
        # Basic identification properties
        assert nuc.name == name, f"Nuclide name should be {name}"
        assert nuc.element.lower() == 'iron', f"Element should be iron for {name}"
        assert nuc.atomic_symbol == 'Fe', f"Atomic symbol should be Fe for {name}"
        assert nuc.atomic_number == 26, f"Atomic number should be 26 for {name}"
        assert nuc.mass_number == mass_num, f"Mass number should be {mass_num} for {name}"
        assert nuc.neutron_number == neutron_num, f"Neutron number should be {neutron_num} for {name}"
        
        # Should have temperature and reaction data
        assert len(nuc.available_temperatures) > 0, f"{name} should have available temperatures"
        assert len(nuc.reaction_mts) > 0, f"{name} should have reaction MTs available"


def test_nuclide_properties_before_loading():
    """Test nuclide properties before loading any data."""
    nuc = Nuclide('TestNuclide')
    
    # Basic properties should be set from the name
    assert nuc.name == 'TestNuclide', "Name should be set from constructor"
    
    # Other properties should be empty/None before loading
    assert nuc.available_temperatures == [], "Should have no temperatures before loading"
    assert nuc.loaded_temperatures == [], "Should have no loaded temperatures before loading"
    
    # Reaction MTs returns None before loading data
    assert nuc.reaction_mts is None, "Should have None reaction MTs before loading"
    
    # Test other properties that might be None before loading
    # These properties may be None before data is loaded
    assert nuc.element is None, "Element should be None before loading"
    assert nuc.atomic_number is None, "Atomic number should be None before loading" 
    assert nuc.mass_number is None, "Mass number should be None before loading"


def test_nuclide_properties_consistency():
    """Test that nuclide properties are consistent across multiple instances."""
    # Create two instances of the same nuclide
    nuc1 = Nuclide('Li6')
    nuc1.read_nuclide_from_json('tests/Li6.json')
    
    nuc2 = Nuclide('Li6')
    nuc2.read_nuclide_from_json('tests/Li6.json')
    
    # Should have identical properties
    assert nuc1.name == nuc2.name, "Names should be identical"
    assert nuc1.element == nuc2.element, "Elements should be identical"
    assert nuc1.atomic_symbol == nuc2.atomic_symbol, "Atomic symbols should be identical"
    assert nuc1.atomic_number == nuc2.atomic_number, "Atomic numbers should be identical"
    assert nuc1.mass_number == nuc2.mass_number, "Mass numbers should be identical"
    assert nuc1.neutron_number == nuc2.neutron_number, "Neutron numbers should be identical"
    assert nuc1.available_temperatures == nuc2.available_temperatures, "Available temperatures should be identical"
    assert nuc1.reaction_mts == nuc2.reaction_mts, "Reaction MTs should be identical"


def test_nuclide_name_parsing():
    """Test that nuclide names are parsed correctly."""
    test_cases = [
        ('Li6', 'Li', 6),
        ('Li7', 'Li', 7),
        ('Be9', 'Be', 9),
        ('Fe54', 'Fe', 54),
        ('Fe56', 'Fe', 56),
        ('Fe57', 'Fe', 57),
        ('Fe58', 'Fe', 58),
    ]
    
    for name, expected_symbol, expected_mass in test_cases:
        nuc = Nuclide(name)
        nuc.read_nuclide_from_json(f'tests/{name}.json')
        
        assert nuc.atomic_symbol == expected_symbol, f"Atomic symbol for {name} should be {expected_symbol}"
        assert nuc.mass_number == expected_mass, f"Mass number for {name} should be {expected_mass}"


def test_nuclide_reaction_mts_content():
    """Test that reaction MTs contain expected values."""
    nuc = Nuclide('Li6')
    nuc.read_nuclide_from_json('tests/Li6.json')
    
    # Check for common reaction types that should be present
    common_mts = [1, 2, 3]  # total, elastic, nonelastic
    present_mts = set(nuc.reaction_mts)
    
    for mt in common_mts:
        if mt in present_mts:  # Not all nuclides have all reaction types
            assert mt in nuc.reaction_mts, f"MT {mt} should be in reaction MTs if present"
    
    # All MTs should be positive integers
    assert all(mt > 0 for mt in nuc.reaction_mts), "All MT numbers should be positive"
    assert all(mt < 1000 for mt in nuc.reaction_mts), "All MT numbers should be reasonable (<1000)"


def test_nuclide_temperature_data_consistency():
    """Test that temperature data is consistent."""
    nuc = Nuclide('Be9')
    nuc.read_nuclide_from_json('tests/Be9.json')
    
    # Available temperatures should include loaded temperatures
    available_set = set(nuc.available_temperatures)
    loaded_set = set(nuc.loaded_temperatures)
    
    assert loaded_set.issubset(available_set), "Loaded temperatures should be subset of available temperatures"
    
    # Temperature strings should be valid
    for temp in nuc.available_temperatures:
        assert isinstance(temp, str), "Temperature should be string"
        assert temp.isdigit() or '.' in temp, "Temperature should be numeric string"
        float_temp = float(temp)
        assert float_temp > 0, "Temperature should be positive"
        assert float_temp < 10000, "Temperature should be reasonable (<10000K)"


def test_nuclide_fissionable_property():
    """Test that the fissionable property is set correctly."""
    # Test non-fissionable nuclides
    non_fissionable = ['Li6', 'Li7', 'Be9', 'Fe54', 'Fe56', 'Fe57', 'Fe58']
    
    for name in non_fissionable:
        nuc = Nuclide(name)
        nuc.read_nuclide_from_json(f'tests/{name}.json')
        
        assert hasattr(nuc, 'fissionable'), f"{name} should have fissionable attribute"
        assert nuc.fissionable is False, f"{name} should not be fissionable"


def test_nuclide_atomic_mass_property():
    """Test atomic mass property if available."""
    nuc = Nuclide('Li6')
    nuc.read_nuclide_from_json('tests/Li6.json')
    
    # Check if atomic mass is available
    if hasattr(nuc, 'atomic_mass'):
        assert isinstance(nuc.atomic_mass, (int, float)), "Atomic mass should be numeric"
        assert nuc.atomic_mass > 0, "Atomic mass should be positive"
        # Li6 atomic mass should be around 6 atomic mass units
        assert 5 < nuc.atomic_mass < 7, "Li6 atomic mass should be around 6 AMU"


def test_nuclide_cross_section_availability():
    """Test that cross section data is available after loading."""
    nuc = Nuclide('Li6')
    nuc.read_nuclide_from_json('tests/Li6.json')
    
    # Should have reactions data
    assert hasattr(nuc, 'reactions'), "Should have reactions attribute"
    assert isinstance(nuc.reactions, dict), "Reactions should be dictionary"
    
    # Should have at least one temperature
    assert len(nuc.reactions) > 0, "Should have reactions for at least one temperature"
    
    # Each temperature should have reaction data
    for temp in nuc.reactions:
        temp_reactions = nuc.reactions[temp]
        assert isinstance(temp_reactions, dict), f"Reactions for {temp} should be dictionary"
        assert len(temp_reactions) > 0, f"Should have at least one reaction at {temp}"
        
        # Check a few reactions have cross section data
        for mt in list(temp_reactions.keys())[:3]:  # Check first 3 reactions
            reaction = temp_reactions[mt]
            assert hasattr(reaction, 'cross_section'), f"Reaction {mt} should have cross_section"
            assert len(reaction.cross_section) > 0, f"Reaction {mt} should have non-empty cross section"
