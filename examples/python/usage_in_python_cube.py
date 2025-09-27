import constructive_solid_geometry_for_mc as csg4mc


s1 = csg4mc.XPlane(x0=2.1, surface_id=5)
s2 = csg4mc.XPlane(x0=-2.1, surface_id=6)
s3 = csg4mc.Sphere(x0=0, y0=0, z0=0, r=4.2, surface_id=1)
# s3 = csg4mc.Cylinder(x0=0, y0=0, z0=0, axis_x=0, axis_y=0, axis_z=1, r=1, surface_id=2)

surfaces_dict = {s1.id: s1, s2.id: s2, s3.id: s3}

region1 = -s1 & +s2 & -s3
inside = region1.contains((0, 0, 0))
print("Point inside region1?", inside)
print(region1.bounding_box())

# region2 = +s1 & -s2 & -s3
# inside = region2.contains((0, 0, 0), surfaces_dict)
# print("Point inside region1?", inside)
# print(region2.bounding_box(surfaces_dict))


import numpy as np

results = []
for y in np.linspace(-10, 10, 11):
    for x in np.linspace(-10, 10, 11):
        contains = region1.contains((x, y, 0.))
        # print(f"Point ({x}, {y}, 0) inside region2? {contains}")
        results.append(int(contains))

results_np = np.array(results).reshape((11, 11))
print(results_np)

# import openmc
# import matplotlib.pyplot as plt

# s1 = openmc.XPlane(x0=2.1)
# s2 = openmc.XPlane(x0=-2.1)
# s3 = openmc.Sphere(x0=0., y0=0., z0=0., r=4.2)

# region1 = -s1 & +s2 & -s3
# print(region1.bounding_box)

# region1.plot(basis='xy')
# plt.show()