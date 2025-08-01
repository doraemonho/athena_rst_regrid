# AthenaK Restart Reader

A utility for reading and analyzing AthenaK restart files with MPI support.

## Features

- ✅ Full compatibility with AthenaK restart file format
- ✅ MPI parallel processing (matches AthenaK's MPI implementation)
- ✅ Automatic load balancing across MPI ranks
- ✅ Support for all physics modules (MHD, Hydro, Turbulence, etc.)
- ✅ Detailed mesh and physics data display

## Quick Start

### Build (Default: Release with MPI)
```bash
./build.sh
```

### Build Options
```bash
./build.sh --debug      # Debug build
./build.sh --no-mpi     # Serial build (no MPI)
./build.sh --clean      # Clean build
./build.sh --help       # Show all options
```

### Usage

**Serial:**
```bash
./build/restart_reader path/to/restart.rst
```

**MPI Parallel:**
```bash
mpirun -np 4 ./build/restart_reader path/to/restart.rst
```

## Build Requirements

- CMake 3.10+
- C++17 compiler
- MPI implementation (OpenMPI, MPICH, etc.) - optional

## Output

The tool displays:
- Restart file metadata (time, cycle, mesh info)
- Physics module configuration  
- MPI distribution information
- Sample physics data from each rank
- Data validation checks

## Files

| File | Description |
|------|-------------|
| `restart_reader.cpp/hpp` | Main restart reader implementation |
| `io_wrapper.cpp/hpp` | MPI-aware file I/O wrapper |
| `mpi_distribution.cpp/hpp` | Load balancing and MPI distribution |
| `physics_reader.cpp/hpp` | Physics data reading and parsing |
| `parameter_parser.cpp/hpp` | Parameter parsing utilities |
| `main.cpp` | Command-line interface |
| `debug_mpi.cpp` | MPI debugging utility |

## MPI Implementation

The restart reader uses **identical MPI patterns** to AthenaK:
- Collective `MPI_File_open()` operations
- Rank 0 reads data, broadcasts to all ranks
- Same data distribution and load balancing algorithms
- Compatible with all AthenaK restart files