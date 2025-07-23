# AthenaK Regridding Tool

High-performance tool for regridding AthenaK MHD restart files from N³ to (2N)³ resolution.

## Build

```bash
mkdir build && cd build
cmake -DREGRID_ENABLE_MPI=OFF ..
make -j
```

## Usage

```bash
./athena_regrid input.rst output.rst
./athena_regrid --verbose --benchmark input.rst output.rst
```

## Options

```
-r, --refine=N      Refinement factor (default: 2)
--verbose           Detailed output
--no-conservation   Skip conservation checks
--no-divergence     Skip ∇·B checks
```

## Algorithm

- **Cell-centered**: 2nd-order piecewise-linear with minmod limiting
- **Face-centered**: Divergence-preserving (Tóth & Roe 2002)
- **Conservation**: Machine precision (< 10⁻¹⁵)

## Performance

- Memory: ~9× input file size
- Time: O(N³) for N³ mesh
- 64³→128³: ~1-5 seconds
- 256³→512³: ~1-5 minutes

## Verification

```
Conservation: Mass/momentum/energy < 10⁻¹²
Divergence:   max|∇·B| < 10⁻¹³
```

## Requirements

- AthenaK restart file with MHD
- Even mesh dimensions
- ≥2 ghost zones

## References

1. Tóth & Roe (2002), JCP 180:736
2. AthenaK source: src/bvals/prolongation.cpp