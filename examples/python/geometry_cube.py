import numpy as np
import yaml as mc


s1 = mc.XPlane(x0=2.1, surface_id=5)
s2 = mc.XPlane(x0=-2.1, surface_id=6)
s3 = mc.Sphere(x0=0, y0=0, z0=0, r=4.2, surface_id=1)
# s3 = mc.Cylinder(x0=0, y0=0, z0=0, axis_x=0, axis_y=0, axis_z=1, r=1, surface_id=2)

surfaces_dict = {s1.surface_id: s1, s2.surface_id: s2, s3.surface_id: s3}

region1 = -s1 & +s2 & -s3
inside = region1.contains((0, 0, 0))
print("Point inside region1?", inside)
print(region1.bounding_box())

results = []
for y in np.linspace(-10, 10, 11):
    for x in np.linspace(-10, 10, 11):
        contains = region1.contains((x, y, 0.))
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




s1 = mc.Plane(0, 0, 1, 5)
s2 = mc.Sphere(x0=0,y0=0,z0=1, r=3, surface_id=1)
s3 = mc.Cylinder(x0=0, y0=0, z0=0, axis_x=0, axis_y=0, axis_z=1, r=1, surface_id=2)

region1 = -s1 & +s2 | ~(-s3)

inside = region1.contains((0, 0, 0))

print("Point inside region1?", inside)

s4 = mc.XPlane(1.0, surface_id=10)
region2 = -s2

inside = region2.contains((0, 0, 0))

print("Point inside region2?", inside)

bb = region2.bounding_box()
print("Bounding box of region2:", bb.lower_left, bb.upper_right)

print(f'Bounding box center {bb.center}')

print(f"bb width {bb.width}")


results = []
for x in np.linspace(
    bb.lower_left[0], bb.upper_right[0], 10
):
    for y in np.linspace(
        bb.lower_left[1], bb.upper_right[1], 10
    ):
        contains = region2.contains((x, y, 0))
        # print(f"Point ({x}, {y}, 0) inside region2? {contains}")
        results.append(int(contains))

results_np = np.array(results).reshape((10, 10))
print(results_np)


# plt.imshow(results_np)
# plt.show()