import yamc

yamc.Config.set_cross_sections("tendl-2019")

mat1 = yamc.Material()
mat1.add_element('Li', 0.5)
mat1.set_density('g/cm3', 2.0)
mat1.volume = 4.2
mat1.read_nuclides_from_json("tendl-2019")
print(mat1)
