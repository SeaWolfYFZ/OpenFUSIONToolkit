# ThinCurr C++ Example

C++ example for OpenFUSIONToolkit ThinCurr module with improved project structure.

## Project Structure

```
ThinCurr-C-example-11/
├── CMakeLists.txt              # CMake build configuration
├── config_example.in           # Example configuration file (Fortran namelist)
├── include/                    # Header files
│   ├── oft_capi.h              # OFT C API declarations
│   └── error_utils.hpp         # Error handling utilities
├── src/                        # Source files
│   └── main.cpp                # Main program
├── oft_in.xml                  # XML configuration file
├── tokamak_mesh_holes_16.h5    # Mesh data file
└── build/                      # Build directory (generated)
    └── bin/
        └── tokamak-test        # Compiled executable
```

## Build Instructions

### CMake (Recommended)
```bash
# Or manually
mkdir build && cd build
cmake ..
make
```

## Running the Example

### Prepare Configuration File
Copy necessary files (including mesh data, XML configuration, Fortran namelist configuration)
```bash
cp ./tokamak_mesh_holes_16.h5 ./oft_in.xml ./oftcppin ./tokamak-visualization.ipynb  ./build/bin
# Edit oftcppin as needed
```

### Run the Simulation
```bash
# Ensure required files are present:
# - tokamak_mesh_holes_16.h5  (mesh data)
# - oft_in.xml                (XML configuration)
# - oftcppin                  (Fortran namelist configuration, hardcoded name)

# Run the simulation (no command line arguments needed)
cd build/bin && ./tokamak-test
```

## Key Improvements

1. **Header/Source Separation**: Clear separation of interface and implementation
2. **CMake Build System**: Standardized, portable build configuration
3. **Modular Design**: Independent components for better maintainability
4. **User-Provided Configuration**: External Fortran namelist files (hardcoded as `oftcppin`) instead of automatic temporary files
5. **Simplified Usage**: No command line arguments needed - just create the `oftcppin` file and run the executable
6. **Error Handling**: Consistent error checking and reporting

## Files Description

### Configuration Files
- `config_example.in`: Example Fortran namelist configuration file with runtime, mesh, and HODLR options
- `oft_in.xml`: XML configuration file defining resistivity, coils, and other physical parameters
- `tokamak_mesh_holes_16.h5`: HDF5 mesh data file

### Header Files
- `include/oft_capi.h`: All `extern "C"` function declarations from OFT Fortran modules
- `include/error_utils.hpp`: Utility functions for error checking and logging

### Source Files
- `src/main.cpp`: Main program with simulation workflow (expects configuration file named `oftcppin`)

### Build Files
- `CMakeLists.txt`: CMake configuration with proper library linking

## Dependencies

- OpenFUSIONToolkit library: `liboftpy.so`
- C++11 compatible compiler
- CMake 3.10 or higher (for CMake build)

## Notes

- This example requires a Fortran namelist configuration file named `oftcppin` (hardcoded filename)
- The `config_example.in` file shows the required format with common options (copy it to `oftcppin`)
- All OFT function signatures remain unchanged for compatibility
- Compared to ThinCurr-C-example-10, this version uses user-provided configuration files (hardcoded as `oftcppin`) instead of automatic temporary file creation
- Ensure all required data files (HDF5 mesh, XML config) are present in the working directory
- No command line arguments are needed - just run the executable after creating the `oftcppin` file