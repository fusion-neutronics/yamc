#!/usr/bin/env python3
"""
Test script to demonstrate the new q_value functionality.
This script shows how to access reaction Q-values from the nuclear data.
"""

import yamc

def test_q_values():
    """Test q_value functionality with Li6 and Be9 nuclides."""
    
    print("Testing Q-value functionality in yamc")
    print("=" * 50)
    
    # Test Li6 nuclide
    print("\n1. Loading Li6 nuclide...")
    li6 = yamc.Nuclide()
    li6.read_nuclide_from_json('tests/Li6.h5')
    
    print(f"   Name: {li6.name}")
    print(f"   Available temperatures: {li6.available_temperatures}")
    
    # Get reactions and show q_values
    reactions = li6.reactions['294']  # Get reactions at 294K
    print(f"   Number of reactions: {len(reactions)}")
    
    print("\n   Li6 Q-values for common reactions:")
    common_reactions = {
        1: "Total",
        2: "Elastic scattering",
        102: "(n,gamma)",
        103: "(n,p)", 
        105: "(n,t)",
        24: "(n,2n+alpha)"
    }
    
    for mt, description in common_reactions.items():
        if mt in reactions:
            q_val = reactions[mt].q_value
            print(f"   MT={mt:3d} {description:20s}: {q_val:10.1f} eV")
    
    # Test Be9 nuclide
    print("\n2. Loading Be9 nuclide...")
    be9 = yamc.Nuclide()
    be9.read_nuclide_from_json('tests/Be9.h5')
    
    print(f"   Name: {be9.name}")
    print(f"   Available temperatures: {be9.available_temperatures}")
    
    # Get reactions and show q_values
    be9_reactions = be9.reactions['294']  # Get reactions at 294K
    print(f"   Number of reactions: {len(be9_reactions)}")
    
    print("\n   Be9 Q-values for common reactions:")
    common_reactions_be9 = {
        1: "Total",
        2: "Elastic scattering", 
        3: "Non-elastic",
        16: "(n,2n)",
        102: "(n,gamma)",
        103: "(n,p)",
        104: "(n,d)",
        105: "(n,t)"
    }
    
    for mt, description in common_reactions_be9.items():
        if mt in be9_reactions:
            q_val = be9_reactions[mt].q_value
            print(f"   MT={mt:3d} {description:20s}: {q_val:10.1f} eV")
    
    print("\n3. Demonstrating access to other reaction properties...")
    # Show that other properties are still available
    mt2_reaction = reactions[2]  # Elastic scattering
    print(f"   MT=2 properties:")
    print(f"   - MT number: {mt2_reaction.mt_number}")
    print(f"   - Q-value: {mt2_reaction.q_value} eV")
    print(f"   - Cross section points: {len(mt2_reaction.cross_section)}")
    print(f"   - Energy grid points: {len(mt2_reaction.energy_grid)}")
    print(f"   - Number of products: {len(mt2_reaction.products)}")
    
    print("\nQ-value functionality test completed successfully!")
    # Test passed if we reach this point

if __name__ == "__main__":
    test_q_values()