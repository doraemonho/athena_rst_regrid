# AthenaK Restart Regridding Tool (Standalone)

`restart_reader` is a small C++17 CLI for working with AthenaK restart files (`.rst`).
It supports:

- Read/inspect a restart file (prints mesh + sample physics values).
- Upscale a restart by 2× in each active dimension and write a new `.rst`.
- Downsample by 2× or 4× in each active dimension and write an Athena-style `.bin`.

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
./build/restart_reader input.rst --downsample-bin output.bin --downsample-factor 4
./build/restart_reader input.rst --downsample-turb-bin turb.bin
./build/restart_reader input.rst --downsample-bin mhd.bin --downsample-turb-bin turb.bin \
  --downsample-factor 4
```

### MPI examples

```bash
mpirun -np 4 ./build/restart_reader input.rst
mpirun -np 4 ./build/restart_reader input.rst --upscale output.rst
mpirun -np 4 ./build/restart_reader input.rst --downsample-bin output.bin
mpirun -np 4 ./build/restart_reader input.rst --downsample-bin output.bin \
  --downsample-factor 4
mpirun -np 4 ./build/restart_reader input.rst --downsample-turb-bin turb.bin
```

## Downsampling notes (what you get)

- Input requirements:
  - MHD data must be present (`<mhd>` block); EOS must be `ideal` or `isothermal`.
  - Global mesh and per-MeshBlock cell counts must be divisible by the downsample
    factor (`2` or `4`) in each active dimension.
  - Assumes a fully refined mesh where each coarse MeshBlock combines `f^D` fine
    MeshBlocks (where `f` is the downsample factor and `D` is 1D/2D/3D), and fine gids
    are ordered by Morton/Z-order child index (see `src/DOWNSAMPLING.md`).
- Output variables (AthenaK `mhd_w_bcc` ordering, written as `float`):
  - `ideal` EOS: `dens`, `velx`, `vely`, `velz`, `eint`, `bcc1`, `bcc2`, `bcc3`.
  - `isothermal` EOS: `dens`, `velx`, `vely`, `velz`, `bcc1`, `bcc2`, `bcc3`.
- Turbulence output (separate `.bin`, AthenaK `turb_force` ordering, written as `float`):
  - `force1`, `force2`, `force3`.

Algorithm details: `src/DOWNSAMPLING.md`.

## Source layout

- Core reader: `src/restart_reader.*`, `src/physics_reader.*`, `src/io_wrapper.*`,
  `src/mpi_distribution.*`, `src/parameter_parser.*`
- Upscaling to `.rst`: `src/upscaler.*`, `src/prolongation.hpp`, `src/restart_writer.*`
- Downsampling to `.bin`: `src/downsampler.*`, `src/binary_writer.*`

## Relationship to the in-tree version

This directory is the standalone snapshot of the in-tree tool under
`../athenak_restart_multiphase/restart_reader/`.
