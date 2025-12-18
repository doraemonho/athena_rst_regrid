# MPI usage notes

The main README for this project is `../README.md`. This file highlights MPI-specific
details that are easy to forget.

## Build with MPI

MPI is enabled by default.

```bash
cd ..
./build.sh
```

To disable MPI (serial build):

```bash
cd ..
./build.sh --no-mpi
```

Manual CMake (equivalent to the script):

```bash
cd ..
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release -DENABLE_MPI=ON
cmake --build build -j
```

## Run

```bash
cd ..
mpirun -np 4 ./build/restart_reader /path/to/restart.rst
mpirun -np 4 ./build/restart_reader /path/to/restart.rst --upscale out.rst
mpirun -np 4 ./build/restart_reader /path/to/restart.rst --downsample-bin out.bin
mpirun -np 4 ./build/restart_reader /path/to/restart.rst --downsample-bin out.bin \
  --downsample-factor 4
```

## Implementation details (high level)

- File I/O uses MPI-IO (`MPI_File_open`, `MPI_File_read_at[_all]`,
  `MPI_File_write_at[_all]`) via `src/io_wrapper.*`.
- `RestartReader` reads some metadata on rank 0 (parameter string, header, MeshBlock
  lists) and broadcasts it to all ranks; physics payloads are read using rank-local and
  collective `*_at[_all]` calls to match AthenaK’s layout.
- Upscaling/downsampling may redistribute data across ranks (each constructs its own
  `MPIDistribution` for the output layout).

## Troubleshooting

- Don’t run with more ranks than MeshBlocks in the restart file.
- Some WSL2 MPI setups can fail to spawn ranks correctly; if `mpirun` reports every
  process as rank 0/1, try running on native Linux/HPC instead.
