#!/bin/bash
# AthenaK Restart Reader Build Script

set -e  # Exit on any error

echo "=== AthenaK Restart Reader Build Script ==="
echo

# Default values
BUILD_TYPE="Release"
ENABLE_MPI="ON"
CLEAN_BUILD=false
INSTALL_PREFIX=""

# Parse command line arguments
while [[ $# -gt 0 ]]; do
    case $1 in
        --debug)
            BUILD_TYPE="Debug"
            shift
            ;;
        --no-mpi)
            ENABLE_MPI="OFF"
            shift
            ;;
        --clean)
            CLEAN_BUILD=true
            shift
            ;;
        --prefix)
            INSTALL_PREFIX="$2"
            shift 2
            ;;
        --help|-h)
            echo "Usage: $0 [OPTIONS]"
            echo "Options:"
            echo "  --debug      Build in Debug mode (default: Release)"
            echo "  --no-mpi     Disable MPI support (default: enabled)"
            echo "  --clean      Clean build directory first"
            echo "  --prefix DIR Set installation prefix"
            echo "  --help       Show this help"
            exit 0
            ;;
        *)
            echo "Unknown option: $1"
            echo "Use --help for usage information"
            exit 1
            ;;
    esac
done

# Display build configuration
echo "Build Configuration:"
echo "  Build Type: $BUILD_TYPE"
echo "  MPI Support: $ENABLE_MPI"
echo "  Clean Build: $CLEAN_BUILD"
if [[ -n "$INSTALL_PREFIX" ]]; then
    echo "  Install Prefix: $INSTALL_PREFIX"
fi
echo

# Clean build directory if requested
if [[ "$CLEAN_BUILD" == true ]]; then
    echo "Cleaning build directory..."
    rm -rf build
fi

# Create build directory
mkdir -p build
cd build

# Configure CMake
echo "Configuring with CMake..."
CMAKE_ARGS=(
    "-DCMAKE_BUILD_TYPE=$BUILD_TYPE"
    "-DENABLE_MPI=$ENABLE_MPI"
)

if [[ -n "$INSTALL_PREFIX" ]]; then
    CMAKE_ARGS+=("-DCMAKE_INSTALL_PREFIX=$INSTALL_PREFIX")
fi

cmake .. "${CMAKE_ARGS[@]}"

# Build
echo "Building..."
make -j$(nproc)

# Success message
echo
echo "=== Build Complete ==="
echo "Artifacts built:"
echo "  restart_reader         - Main CLI"
echo "  librestart_reader_lib.a - Static library (internal use)"
echo
echo "To install: make install"
echo "To test: ./restart_reader <restart_file.rst>"
if [[ "$ENABLE_MPI" == "ON" ]]; then
    echo "To test MPI: mpirun -np 4 ./restart_reader <restart_file.rst>"
fi
