#!/usr/bin/env python
"""
Analyze poloidal (换向) loop current in tokamak wall.

This script:
1. Identifies which NODESET corresponds to poloidal loop
2. Extracts hole potential from restart files  
3. Computes loop current I = phi_h / mu_0
4. Plots current vs time
"""
import h5py
import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path

mu0 = 4 * np.pi * 1e-7

print("="*60)
print("Poloidal Loop Current Analysis")
print("="*60)

# Step 1: Identify poloidal hole
print("\nStep 1: Identifying poloidal vs toroidal holes...")

with h5py.File('tokamak_mesh_EHL2_66holes.h5', 'r') as f:
    coords = f['mesh/R'][:]
    
    for i in [65, 66]:
        key = f'mesh/NODESET{i:04d}'
        indices = f[key][:] - 1
        hole_coords = coords[indices]
        
        R = np.sqrt(hole_coords[:,0]**2 + hole_coords[:,1]**2)
        Z = hole_coords[:,2]
        phi = np.arctan2(hole_coords[:,1], hole_coords[:,0])
        phi = np.unwrap(phi)
        
        R_range = R.max() - R.min()
        Z_range = Z.max() - Z.min()
        phi_range = phi.max() - phi.min()
        
        print(f"\nNODESET{i:04d} ({len(indices)} vertices):")
        print(f"  R range: [{R.min():.3f}, {R.max():.3f}] m (Δ={R_range:.3f})")
        print(f"  Z range: [{Z.min():.3f}, {Z.max():.3f}] m (Δ={Z_range:.3f})")
        print(f"  Phi range: [{phi.min():.3f}, {phi.max():.3f}] rad (Δ={phi_range:.3f})")
        
        if phi_range < 1.0:
            print(f"  → This is the POLOIDAL loop (绕小环)")
            poloidal_hole_num = i
        elif phi_range > 5.0:
            print(f"  → This is the TOROIDAL loop (绕大环)")

print(f"\n✓ Identified: POLOIDAL hole = NODESET{poloidal_hole_num:04d}")

# Step 2: Find hole DOF index
print("\nStep 2: Finding hole DOF index...")

rst_files = sorted(Path('.').glob('pThinCurr_*.rst'))
if len(rst_files) == 0:
    print("ERROR: No restart files found.")
    exit(1)

with h5py.File(rst_files[0], 'r') as f:
    potential = f['potential'][:]
    nelems = len(potential)

nholes = 66
poloidal_dof = nelems - nholes + (poloidal_hole_num - 1)
print(f"  Total DOFs: {nelems}")
print(f"  Poloidal hole DOF index: {poloidal_dof}")

# Step 3: Extract current time series
print("\nStep 3: Extracting poloidal loop current...")

times, currents = [], []
for rst_file in sorted(rst_files):
    with h5py.File(rst_file, 'r') as f:
        t_str = rst_file.stem.split('_')[-1]
        t = int(t_str) * 1e-3
        potential = f['potential'][:]
        phi_hole = potential[poloidal_dof]
        I_hole = phi_hole / mu0
        times.append(t)
        currents.append(I_hole)

times = np.array(times)
currents = np.array(currents)
print(f"✓ Extracted {len(times)} time points")

# Step 4: Plot
print("\nStep 4: Plotting results...")
fig, ax = plt.subplots(figsize=(10, 5))
ax.plot(times, currents, 'b-', linewidth=2)
ax.set_xlabel('Time [s]')
ax.set_ylabel('Current [A]')
ax.set_title(f'Poloidal Loop Current (NODESET{poloidal_hole_num:04d})')
ax.grid(True, alpha=0.3)
ax.axhline(0, color='k', linestyle='--', linewidth=0.5)
plt.tight_layout()
plt.savefig('poloidal_loop_current.png', dpi=150)
print("  Saved: poloidal_loop_current.png")

print("\n" + "="*60)
print(f"Summary: Peak current = {np.abs(currents).max():.2f} A")
print("="*60)
