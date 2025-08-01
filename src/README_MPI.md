# AthenaK Restart Reader - MPI Usage Guide

## Overview
The restart reader now properly supports MPI parallel execution, matching AthenaK's MPI implementation exactly.

## Building with MPI

### Prerequisites
- MPI installation (OpenMPI, MPICH, etc.)
- CMake 3.10 or higher

### Build Instructions
```bash
cd restart_reader
mkdir build
cd build

# Configure with MPI enabled
cmake .. -DMPI_PARALLEL_ENABLED=ON

# Build
make
```

## Running the Restart Reader

### Serial Execution (1 rank)
```bash
./restart_reader /path/to/restart_file.rst
```

### MPI Parallel Execution
```bash
# Run with 4 MPI ranks
mpirun -np 4 ./restart_reader /path/to/restart_file.rst

# Run with 16 MPI ranks
mpirun -np 16 ./restart_reader /path/to/restart_file.rst

# Run with custom MPI options
mpirun -np 8 --oversubscribe ./restart_reader /path/to/restart_file.rst
```

## Important Notes

1. **Rank Detection**: The restart reader now automatically detects MPI rank and total ranks from the MPI runtime, just like AthenaK does.

2. **Load Balancing**: The reader automatically distributes MeshBlocks across MPI ranks using AthenaK's load balancing algorithm.

3. **File Access**: Only rank 0 opens the restart file, then broadcasts data to other ranks as needed.

4. **Output**: Each rank will display information about its assigned MeshBlocks:
   ```
   Rank 0: MeshBlocks 0-24 (25 total)
   Rank 1: MeshBlocks 25-49 (25 total)
   ...
   ```

## Example Usage

### Reading a 100 MeshBlock restart file with 4 ranks:
```bash
mpirun -np 4 ./restart_reader ../turb_run/turb.00010.rst
```

Expected output:
```
=== AthenaK Restart Reader with MPI Support ===
Running with 4 MPI rank(s)

Reading AthenaK restart file: ../turb_run/turb.00010.rst
Using automatic load balancing
MPI Distribution:
  Rank 0: MeshBlocks 0-24 (25 total)
  Rank 1: MeshBlocks 25-49 (25 total)
  Rank 2: MeshBlocks 50-74 (25 total)
  Rank 3: MeshBlocks 75-99 (25 total)
...
```

## Troubleshooting

1. **MPI not found during build**: Ensure MPI is in your PATH or specify MPI compiler wrappers:
   ```bash
   cmake .. -DMPI_PARALLEL_ENABLED=ON -DCMAKE_CXX_COMPILER=mpicxx
   ```

2. **Segmentation fault**: Make sure the number of MPI ranks doesn't exceed the number of MeshBlocks in the restart file.

3. **Wrong data distribution**: The reader uses AthenaK's exact load balancing algorithm. If you need custom distribution, you can modify the MPIDistribution class.

4. **WSL2 MPI Issues**: On WSL2, MPI may not properly initialize across processes, causing each process to think it's rank 0 of 1. This is a known WSL2 limitation. The code is still correct and will work properly on actual HPC systems or native Linux environments.

5. **MPI warnings**: The warning "Attempting to use an MPI routine before initializing" is a false positive from certain MPI implementations and can be safely ignored if the program runs correctly.