#!/usr/bin/env python3
"""
Simple example plotting the Li6 inelastic cross sections.
"""

import yamc
import matplotlib.pyplot as plt

li6_nuc = yamc.Nuclide("Li6")
li6_nuc.read_nuclide_from_json("tests/Li6.h5")

xs_4, energy_4 = li6_nuc.microscopic_cross_section(4, "294")
xs_57, energy_57 = li6_nuc.microscopic_cross_section(57, "294")

plt.figure(figsize=(10, 6))
plt.plot(energy_4, xs_4, label="4", linewidth=1.5)

# found with nuclide.reaction_mts 
mts = [
    51, 52, 53, 54, 55, 56, 57, 58, 59, 60, 61, 62, 63, 64, 65,
    66, 67, 68, 69, 70, 71, 72, 73, 74, 75, 76, 77, 78, 79, 80, 81,
]
for mt in mts:
    xs, energy = li6_nuc.microscopic_cross_section(mt, "294")
    plt.plot(energy, xs, label=str(mt), linewidth=1.0, alpha=0.5)

plt.xlabel("Energy (eV)")
plt.ylabel("Cross Section (barns)")
plt.title("Li6 inelastic Cross Section Comparison")
plt.legend()
plt.grid(True, alpha=0.3)
plt.tight_layout()

# Save and show
plt.savefig("li6_comparison.png", dpi=157, bbox_inches="tight")
print("Plot saved as li6_comparison.png")
