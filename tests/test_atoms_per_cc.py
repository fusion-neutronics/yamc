import pytest
from yamc import Material, Config


def test_get_atoms_per_barn_cm():
    # Configure cross-sections to load nuclide data
    Config.set_cross_sections({"Li6": "tests/Li6.h5", "Li7": "tests/Li7.h5"})

    # Test with Li isotopes - proper atomic mass calculation
    material = Material()
    material.add_nuclide("Li6", 0.5)
    material.add_nuclide("Li7", 0.5)
    material.set_density("g/cm3", 1.0)
    # Load nuclide data (required for atomic mass from HDF5)
    material.read_nuclides_from_h5({"Li6": "tests/Li6.h5", "Li7": "tests/Li7.h5"})

    atoms = material.get_atoms_per_barn_cm()
    assert len(atoms) == 2, "Should have 2 nuclides in the dict"

    # Verify both nuclides have positive atom densities
    assert atoms["Li6"] > 0, "Li6 atoms/barn-cm should be positive"
    assert atoms["Li7"] > 0, "Li7 atoms/barn-cm should be positive"
    # The values should be approximately equal since fractions are equal
    assert atoms["Li6"] == pytest.approx(atoms["Li7"], rel=0.1), "Li6 and Li7 should have similar atom densities"

    # Test with different density units (kg/m³)
    material = Material()
    material.add_nuclide("Li6", 1.0)
    material.set_density("kg/m3", 1000.0)  # 1000 kg/m³ = 1 g/cm³
    material.read_nuclides_from_h5({"Li6": "tests/Li6.h5"})

    atoms_kg_m3 = material.get_atoms_per_barn_cm()
    assert len(atoms_kg_m3) == 1, "Should have 1 nuclide in the dict"

    # Compare with same material using g/cm³
    material_g_cm3 = Material()
    material_g_cm3.add_nuclide("Li6", 1.0)
    material_g_cm3.set_density("g/cm3", 1.0)
    material_g_cm3.read_nuclides_from_h5({"Li6": "tests/Li6.h5"})

    atoms_g_cm3 = material_g_cm3.get_atoms_per_barn_cm()

    # Both should give the same result since the densities are equivalent
    assert atoms_kg_m3["Li6"] == pytest.approx(atoms_g_cm3["Li6"], rel=0.01), "Different density units should give consistent results"
