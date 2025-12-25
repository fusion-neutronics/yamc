import yamc
yamc.Config.set_cross_sections({
    "Be9": "tests/Be9.h5",
    "Fe54": "tests/Fe54.h5",
    "Fe56": "tests/Fe56.h5",
    "Fe57": "tests/Fe57.h5",
    "Fe58": "tests/Fe58.h5",
    "Li6": "tests/Li6.h5",
    "Li7": "tests/Li7.h5",
})

mat1 = yamc.Material(material_id=1, name='Test Material 1')
mat1.add_nuclide('Li6', 0.5)
mat1.add_nuclide('Li7', 0.5)
mat1.add_nuclide('Fe56', 1.0)
mat1.add_nuclide('Be9', 1.0)
mat1.set_density('g/cm3', 2.0)
print(mat1)

# Demonstrate partial override: provide only one mapping (others pulled from Config)
# mat1.read_nuclides_from_h5({"Li7": "tests/Li7.h5"})

xs, energy = mat1.macroscopic_cross_section(reaction='(n,gamma)')
xs, energy = mat1.macroscopic_cross_section(reaction=1)
print(xs)