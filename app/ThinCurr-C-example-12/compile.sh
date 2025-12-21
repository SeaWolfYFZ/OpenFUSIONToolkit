#!/bin/bash
# ThinCurr C++ Example Build Script
# Two build methods available:
# 1. Direct compilation (original method)
# 2. CMake build (recommended)

set -e

OFT_INSTALL_DIR="/home/yfz/project/OpenFUSIONToolkit-251219/install_release"

echo "=== ThinCurr C++ Example Build Script ==="
echo ""

if [ "$1" = "direct" ] || [ "$1" = "-d" ]; then
    echo "Building with direct compilation..."
    echo "g++ -std=c++11 -o tokamak-test-10 tokamak-test-10.cpp -L$OFT_INSTALL_DIR/bin -Wl,-rpath,$OFT_INSTALL_DIR/bin -loftpy"
    g++ -std=c++11 -o tokamak-test-10 tokamak-test-10.cpp -L$OFT_INSTALL_DIR/bin -Wl,-rpath,$OFT_INSTALL_DIR/bin -loftpy
    echo "Direct build complete: tokamak-test-10"
else
    echo "Building with CMake (recommended)..."
    echo ""

    if [ ! -d "build" ]; then
        echo "Creating build directory..."
        mkdir -p build
    fi

    cd build
    echo "Running CMake..."
    cmake ..

    echo ""
    echo "Compiling..."
    make -j$(nproc)

    echo ""
    echo "CMake build complete!"
    echo "Executable: build/bin/tokamak-test-10"
    echo ""
    echo "To run: ./build/bin/tokamak-test-10"
fi
