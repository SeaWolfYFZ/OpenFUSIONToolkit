# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Build Commands

```bash
# Configure and build (using config_cmake.sh)
bash config_cmake.sh && cd build_release && make -j && make install

# Alternative: manual cmake configuration
mkdir build && cd build
cmake -DCMAKE_BUILD_TYPE=Release -DOFT_BUILD_PYTHON=ON -DOFT_BUILD_TESTS=ON ../src
make -j && make install
```

### Running Tests

```bash
# Run all tests (from build directory)
cd build_release && make test

# Run all tests including slow ones
make test_full

# Run specific pytest tests
cd build_release/tests && pytest -m "not slow" physics/test_ThinCurr.py
```

## Project Architecture

OpenFUSIONToolkit (OFT) is a Fortran-based finite element framework for plasma and fusion research. It uses **Doxygen** for code documentation.

### Core Source Directories (`src/`)

| Directory | Purpose |
|-----------|---------|
| `base/` | Core utilities, I/O, memory management |
| `grid/` | Mesh handling (triangular, tetrahedral, hexahedral) |
| `lin_alg/` | Linear algebra (sparse/dense matrices, solvers) |
| `fem/` | Finite element infrastructure |
| `physics/` | Physics modules (ThinCurr, TokaMaker, MUG, Marklin) |
| `python/` | Python bindings via f2py |
| `docs/` | Doxygen documentation sources |

### Physics Modules

- **ThinCurr**: Thin-wall eddy current modeling (primary focus of this repo)
  - Key files: `thin_wall.F90`, `thin_wall_solvers.F90`, `thin_wall_hodlr.F90`
- **TokaMaker**: Axisymmetric MHD equilibria
  - Key files: `grad_shaf*.F90`
- **MUG**: Extended MHD simulations
- **Marklin**: 3D force-free equilibria

### Python Package Structure

```
src/python/OpenFUSIONToolkit/
├── ThinCurr/    # ThinCurr Python interface
├── TokaMaker/   # TokaMaker Python interface
├── Marklin/     # Marklin Python interface
└── _core.py     # Core utilities
```

## Key Technical Concepts

### ThinCurr: Mesh Topology

ThinCurr models surface currents using a scalar potential: `J_s = ∇χ × n̂`. This requires special handling for multiply-connected geometries:

- **"Holes"**: Special DOFs for topologically non-trivial closed loops (e.g., poloidal/toroidal directions in a torus). Stored at indices `np_active+1` to `np_active+nholes` in solution vectors.
- **"Closures"**: Gauge-fixing vertices removed from the active set for closed surfaces.
- **Data structures**: `hole_mesh` type in `thin_wall.F90`; `kfh/lfh` arrays for face-hole interactions.

### ThinCurr: HODLR Matrix Compression

For large models, ThinCurr uses Hierarchical Off-Diagonal Low-Rank approximation:
- Dense L matrix assembly is O(N²); HODLR reduces to O(N log N)
- Binary tree spatial partitioning (`oft_tw_block`, `oft_tw_level` types)
- ACA+ (Adaptive Cross Approximation) for off-diagonal blocks
- Implementation: `thin_wall_hodlr.F90`

## Code Style Notes

- **Fortran comments**: Use `!` for regular comments, `!!` for Doxygen documentation
- **Doxygen format**: `!>` for file/module-level docs, `!<` for member documentation
- **LaTeX in docs**: Use `\f[...\f]` for display equations, `\f$...\f$` for inline

## CMake Build Options

Key options in `config_cmake.sh`:
- `OFT_BUILD_PYTHON`: Build Python bindings
- `OFT_BUILD_TESTS`: Enable test targets
- `OFT_BUILD_DOCS`: Build Doxygen documentation
- `OFT_USE_OpenMP`: Enable OpenMP parallelization
- `OFT_USE_MPI`: Enable MPI support

## External Dependencies

The project uses several external libraries (paths in `config_cmake.sh`):
- HDF5: File I/O
- METIS: Mesh partitioning
- ARPACK: Eigenvalue solvers
- OpenBLAS: BLAS/LAPACK
- FoX: XML parsing