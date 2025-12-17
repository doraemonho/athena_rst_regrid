# AthenaK Restart Regridding Tool (Standalone)

`restart_reader` is a small C++17 CLI for working with AthenaK restart files (`.rst`).
It supports:

- Read/inspect a restart file (prints mesh + sample physics values).
- Upscale a restart by 2× in each active dimension and write a new `.rst`.
- Downsample by 2× in each active dimension and write an Athena-style `.bin`.

MPI is optional. When enabled, MeshBlocks are distributed across ranks using the same
distribution logic as AthenaK.

## Build

### Script (recommended)

```bash
./build.sh            # Release + MPI (default)
./build.sh --debug
./build.sh --no-mpi
./build.sh --clean
```

### Manual (CMake)

```bash
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release -DENABLE_MPI=ON
cmake --build build -j
```

### Requirements

- CMake 3.10+
- C++17 compiler
- MPI implementation (optional; enabled by default)

## Usage

### Read-only (inspect)

```bash
./build/restart_reader path/to/restart.rst
```

### Upscale to restart (`.rst`)

```bash
./build/restart_reader input.rst --upscale output.rst
```

### Downsample to binary (`.bin`)

```bash
./build/restart_reader input.rst --downsample-bin output.bin
```

### MPI examples

```bash
mpirun -np 4 ./build/restart_reader input.rst
mpirun -np 4 ./build/restart_reader input.rst --upscale output.rst
mpirun -np 4 ./build/restart_reader input.rst --downsample-bin output.bin
```

## Downsampling notes (what you get)

- Input requirements:
  - MHD data must be present (`<mhd>` block); EOS must be `ideal` or `isothermal`.
  - Global mesh and per-MeshBlock cell counts must be divisible by 2 in each active
    dimension.
  - Assumes a fully refined mesh where each parent has all `2^D` children and fine gids
    are ordered by child index.
- Output variables (always 8, written as `float`): `dens`, `vel1`, `vel2`, `vel3`,
  `press`, `Bcc1`, `Bcc2`, `Bcc3`.

Algorithm details: `src/DOWNSAMPLING.md`.

## Source layout

- Core reader: `src/restart_reader.*`, `src/physics_reader.*`, `src/io_wrapper.*`,
  `src/mpi_distribution.*`, `src/parameter_parser.*`
- Upscaling to `.rst`: `src/upscaler.*`, `src/prolongation.hpp`, `src/restart_writer.*`
- Downsampling to `.bin`: `src/downsampler.*`, `src/binary_writer.*`

## Relationship to the in-tree version

This directory is the standalone snapshot of the in-tree tool under
`../athenak_restart_multiphase/restart_reader/`.
