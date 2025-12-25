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
mat1.set_density('g/cm3', 2.0)
print(mat1)

mat2 = yamc.Material(material_id=2)
mat2.add_nuclide('Li7', 0.5)
mat2.set_density('g/cm3', 2.0)
print(mat2)

mat3 = yamc.Material(name='Test Material 3')
mat3.add_nuclide('Fe56', 1.0)
mat3.set_density('g/cm3', 2.0)
print(mat3)

mat4 = yamc.Material()
mat4.add_nuclide('Be9', 1.0)
mat4.set_density('g/cm3', 2.0)
print(mat4)
