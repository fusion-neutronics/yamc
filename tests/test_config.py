import pytest
import yamc

@pytest.fixture(autouse=True)
def clear_config():
    """Clear config state before each test"""
    yamc.Config.clear()
    yield
    # No cleanup needed after test

def test_set_cross_sections_with_dict():
    """Test setting cross sections with a dictionary"""
    cross_sections = {
        "Li6": "tendl-2019",
        "Li7": "../../tests/Li7.h5"
    }
    yamc.Config.set_cross_sections(cross_sections)
    
    # Verify the cross sections were set
    assert yamc.Config.get_cross_section("Li6") == "tendl-2019"
    assert yamc.Config.get_cross_section("Li7") == "../../tests/Li7.h5"

def test_set_cross_sections_with_string():
    """Test setting cross sections with a global keyword string"""
    yamc.Config.set_cross_sections("tendl-2019")
    
    # Any nuclide should now return the global keyword
    assert yamc.Config.get_cross_section("Fe56") == "tendl-2019"
    assert yamc.Config.get_cross_section("Li6") == "tendl-2019"

def test_set_cross_section_single_nuclide_path():
    """Test setting a single nuclide with a file path"""
    yamc.Config.set_cross_section("Fe56", "../../tests/Fe56.h5")
    assert yamc.Config.get_cross_section("Fe56") == "../../tests/Fe56.h5"

def test_set_cross_section_single_nuclide_keyword():
    """Test setting a single nuclide with a keyword"""
    yamc.Config.set_cross_section("Fe56", "tendl-2019")
    assert yamc.Config.get_cross_section("Fe56") == "tendl-2019"

def test_set_cross_section_global_keyword():
    """Test setting a global keyword using set_cross_section"""
    yamc.Config.set_cross_section("tendl-2019")
    assert yamc.Config.get_cross_section("Li6") == "tendl-2019"
    assert yamc.Config.get_cross_section("Fe56") == "tendl-2019"

def test_set_cross_sections_invalid_type():
    """Test that set_cross_sections raises TypeError for invalid input"""
    with pytest.raises(TypeError):
        yamc.Config.set_cross_sections(123)  # Invalid type

def test_mixed_global_and_specific_config():
    """Test mixing global default with specific nuclide overrides"""
    # Set global default to TENDL
    yamc.Config.set_cross_sections("tendl-2019")
    
    # Override specific nuclides to FENDL
    yamc.Config.set_cross_section("Fe56", "fendl-3.1d")
    yamc.Config.set_cross_section("Li6", "tests/Li6.h5")
    
    # Check that specific overrides work
    assert yamc.Config.get_cross_section("Fe56") == "fendl-3.1d"
    assert yamc.Config.get_cross_section("Li6") == "tests/Li6.h5"
    
    # Check that other nuclides fall back to global default
    assert yamc.Config.get_cross_section("Be9") == "tendl-2019"
    assert yamc.Config.get_cross_section("U235") == "tendl-2019"

def test_explicit_path_override_in_nuclide_loading():
    """Test that explicit paths in read_nuclide_from_h5 work correctly"""
    # Set up global config with local HDF5 files
    yamc.Config.set_cross_sections({
        "Li6": "tests/Li6.h5",
        "Be9": "tests/Be9.h5"
    })

    # Create nuclide and load from explicit path
    li6_explicit = yamc.Nuclide("Li6")
    li6_explicit.read_nuclide_from_h5("tests/Li6.h5")

    # Should load successfully
    assert li6_explicit.name == "Li6"

    # Verify it loaded
    assert len(li6_explicit.available_temperatures) > 0

# def test_set_cross_section_invalid_keyword():
#     """Test that set_cross_section raises error for invalid keyword"""
#     import pytest
#     # This should raise a panic that gets converted to a Python exception
#     with pytest.raises(Exception) as exc_info:
#         yamc.Config.set_cross_section("invalid-keyword")
#     assert "Invalid cross section keyword" in str(exc_info.value)