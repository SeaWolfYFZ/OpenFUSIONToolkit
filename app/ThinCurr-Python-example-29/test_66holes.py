#!/usr/bin/env python
"""
Test script to generate mesh with 66 holes (all 64 ports + 2 homology bases)
and test if this resolves the oscillation issue in the skipped port.
"""

import os
import sys
import h5py
import numpy as np

# Modify the script to include all 64 ports as holes
print("=" * 60)
print("Step 1: Generating mesh with 66 holes")
print("=" * 60)

# Run the modified ThinCurr_compute_holes.py
os.system(
    "python ThinCurr_compute_holes.py "
    + "--in_file=tokamak_mesh_EHL2.h5 "
    + "--out_file=tokamak_mesh_EHL2_66holes.h5 "
    + "--plot_final=False "
    + "--debug"
)

# Check the result
print("\n" + "=" * 60)
print("Step 2: Checking generated NODESETs")
print("=" * 60)

with h5py.File("tokamak_mesh_EHL2_66holes.h5", "r") as f:
    num_nodesets = f["mesh/NUM_NODESETS"][0]
    print(f"Number of NODESETs: {num_nodesets}")

    # List all NODESET sizes
    print("\nNODESET details:")
    for i in range(1, num_nodesets + 1):
        key = f"mesh/NODESET{i:04d}"
        if key in f:
            size = len(f[key])
            print(f"  NODESET{i:04d}: {size} vertices")

    # Check for SIDESETs
    if "mesh/NUM_SIDESETS" in f:
        num_sidesets = f["mesh/NUM_SIDESETS"][0]
        print(f"\nNumber of SIDESETs: {num_sidesets}")
    else:
        print("\nNumber of SIDESETs: 0 (no closures)")

print("\n" + "=" * 60)
print("Step 3: Running ThinCurr simulation with 66-hole mesh")
print("=" * 60)

# Now test with ThinCurr
from OpenFUSIONToolkit import OFT_env
from OpenFUSIONToolkit.ThinCurr import ThinCurr

myOFT = OFT_env(nthreads=28)
tw = ThinCurr(myOFT)

print(f"\nSetting up model with 66-hole mesh...")
tw.setup_model(mesh_file="tokamak_mesh_EHL2_66holes.h5", xml_filename="oft_in.xml")

print(f"\nModel summary:")
print(f"  nholes = {tw.nholes}")
print(f"  np_active = {tw.np_active}")
print(f"  nelems = {tw.nelems}")

# Compute L matrix
print("\nComputing L matrix (HODLR)...")
tw.compute_Lmat(use_hodlr=True, cache_file="DATA_HOLDR_L_66holes.save")

# Compute R matrix
print("\nComputing R matrix...")
tw.compute_Rmat()


# Load current waveform
def load_time_current_txt(path):
    times_ms = []
    currents_a = []
    with open(path, "r") as f:
        for line in f:
            stripped = line.strip()
            if not stripped or stripped.startswith("#") or "Time" in stripped:
                continue
            parts = stripped.split()
            if len(parts) >= 2:
                times_ms.append(float(parts[0]))
                currents_a.append(float(parts[1]))
    return times_ms, currents_a


times_ms, currents_a = load_time_current_txt("Central_Solenoid_Current.txt")
coil_currs = np.array([[1e-3 * i, currents_a[i]] for i in range(len(currents_a))])

# Run short simulation to test
print("\nRunning time-domain simulation (100 steps)...")
dt = 1.0e-3
nsteps_test = 100

tw.run_td(dt, nsteps_test, status_freq=10, coil_currs=coil_currs)

print("\n" + "=" * 60)
print("Test complete!")
print("=" * 60)
print("\nNext steps:")
print("1. Check jumper hist file for oscillation patterns")
print("2. If successful, run full simulation")
print("3. Compare with 65-hole results")
