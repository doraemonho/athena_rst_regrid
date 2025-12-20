# Downsampling (2× / 4× coarsening) in `github_version_regridding/src`

This document explains, step-by-step, how the standalone tool downsamples an
AthenaK restart (`.rst`) by a factor of **2** or **4** in each active dimension
and writes an Athena-style binary output (`.bin`).

The implementation lives in:

- `github_version_regridding/src/downsampler.hpp`
- `github_version_regridding/src/downsampler.cpp`
- `github_version_regridding/src/binary_writer.hpp`
- `github_version_regridding/src/binary_writer.cpp`

## 0) Entry point (CLI → downsampler)

The CLI is in `github_version_regridding/src/main.cpp`.

Downsampling is enabled by:

- `--downsample-bin <mhd.bin>`: MHD `mhd_w_bcc` output
- `--downsample-turb-bin <turb.bin>`: turbulence forcing `turb_force` output

Both accept `--downsample-factor {2|4}` (default: `2`).

```cpp
// main.cpp (simplified)
} else if (downsample_mhd_mode || downsample_turb_mode) {
  Downsampler downsampler(reader, downsample_factor);
  if (downsample_mhd_mode) downsampler.DownsampleToBinary(output_mhd_filename);
  if (downsample_turb_mode) {
    downsampler.DownsampleTurbulenceToBinary(output_turb_filename);
  }
}
```

So the call stack is roughly:

1. `RestartReader::ReadRestartFile(...)` reads the `.rst` and fills in mesh/physics data.
2. `Downsampler::DownsampleToBinary(mhd.bin)` performs MHD coarsening + writes `.bin`.
3. `Downsampler::DownsampleTurbulenceToBinary(turb.bin)` performs turbulence coarsening
   + writes `.bin`.

## 1) Preconditions (what the downsampler assumes)

The downsampler is intentionally narrow/specific.

Common requirements (MHD + turbulence outputs):

1. **Global mesh and per-MeshBlock cell counts must be divisible by the downsample
   factor** (in every active dimension).
2. **The restart must represent a globally refined mesh** with *complete sibling
   groups*:
   - `nmb_total` must be divisible by `f^D` where `f` is the downsample factor and `D`
     is 1D/2D/3D.
   - fine MeshBlocks are assumed to be ordered in contiguous sibling groups in
     the expected child order (explained below).
3. **`root_level >= log2(f)`** since downsampling reduces `root_level` by that many
   refinement levels (1 for `f=2`, 2 for `f=4`).

Output-specific requirements:

- For MHD output (`DownsampleToBinary()`):
  - MHD must be present (`config.has_mhd == true`).
  - EOS must be `ideal` or `isothermal`.
  - `ideal` requires energy (`config.nmhd >= 5`).
- For turbulence output (`DownsampleTurbulenceToBinary()`):
  - turbulence forcing must be present (`config.has_turbulence == true`).

If any of these are violated, the corresponding downsample routine returns `false`.

## 2) Step-by-step logic inside `DownsampleToBinary()`

### Step 2.1: Determine dimensionality and sibling-group size

The constructor decides whether the dataset is 1D/2D/3D using the *global* mesh
indices from the restart header:

```cpp
// downsampler.cpp (simplified)
multi_d_ = (mesh_indcs.nx2 > 1);
three_d_ = (mesh_indcs.nx3 > 1);
nfine_per_coarse_ = downsample_factor;
if (multi_d_) nfine_per_coarse_ *= downsample_factor;
if (three_d_) nfine_per_coarse_ *= downsample_factor;
```

Interpretation:

- 1D: `f` fine MeshBlocks combine into 1 coarse MeshBlock.
- 2D: `f^2` fine MeshBlocks combine into 1 coarse MeshBlock.
- 3D: `f^3` fine MeshBlocks combine into 1 coarse MeshBlock.

Examples:

- `f=2`: 2 / 4 / 8
- `f=4`: 4 / 16 / 64

### Step 2.2: Parse EOS settings (and floors) from the restart parameter string

`Downsampler::ParseEOS()` reads `eos`, `gamma`, `iso_sound_speed`, and the
AthenaK EOS floor parameters from the restart’s text parameter block (via
`ParameterParser`):

- `dfloor`, `pfloor`, `tfloor`, `sfloor`, `tmax`

Defaults match AthenaK’s `EquationOfState` constructor:

- floors default to `FLT_MIN`
- `tmax` defaults to `FLT_MAX`

```cpp
// downsampler.cpp
eos_ = "ideal";
auto eos_opt = ParameterParser::ExtractStringInBlock(params, "mhd", "eos");
if (eos_opt.has_value()) {
  eos_ = *eos_opt;
}

auto gamma_opt = ParameterParser::ExtractRealInBlock(params, "mhd", "gamma");
if (!gamma_opt.has_value()) {
  gamma_opt = ParameterParser::ExtractRealInBlock(params, "hydro", "gamma");
}
if (gamma_opt.has_value()) {
  gamma_ = *gamma_opt;
}

auto cs_opt = ParameterParser::ExtractRealInBlock(params, "mhd", "iso_sound_speed");
if (cs_opt.has_value()) {
  iso_sound_speed_ = *cs_opt;
}
```

### Step 2.3: Compute the coarse mesh configuration (geometry + counts)

`ComputeCoarseMeshConfig()` constructs the *coarsened* mesh description:

- `dx*` is multiplied by `f` in each active dimension.
- global `nx*` is divided by `f` in each active dimension.
- `coarse_root_level = fine_root_level - log2(f)`.
- `coarse_nmb_total = fine_nmb_total / nfine_per_coarse`.

```cpp
// downsampler.cpp (simplified)
coarse_mesh_size_.dx1 = fine.dx1 * f;
coarse_mesh_size_.dx2 = fine.dx2 * (multi_d_ ? f : 1.0);
coarse_mesh_size_.dx3 = fine.dx3 * (three_d_ ? f : 1.0);

coarse_mesh_indcs_.nx1 = fine_indcs.nx1 / f;
coarse_mesh_indcs_.nx2 = multi_d_ ? fine_indcs.nx2 / f : 1;
coarse_mesh_indcs_.nx3 = three_d_ ? fine_indcs.nx3 / f : 1;

coarse_nmb_total_ = reader_.GetNMBTotal() / nfine_per_coarse_;
coarse_root_level_ = reader_.GetRootLevel() - log2(f);
```

Important detail: **the per-MeshBlock cell counts are *not* changed**. Instead,
each coarse MeshBlock is formed by merging `f^D` fine MeshBlocks and filling the
coarse MeshBlock with `f^D` *sub-chunks* (each sub-chunk is `1/f`-resolution in
each active dimension).

### Step 2.4: Build the coarse logical locations (and validate fine ordering)

`BuildCoarseLogicalLocations()` consumes `reader_.GetLlocEachMB()` (logical
location per fine gid) and constructs a list of parent (coarse) logical
locations. It enforces that:

- every sibling group has a consistent parent,
- all siblings are at the same refinement level,
- the sibling order matches a Morton/Z-order child index for the `f^D` children.

Child indexing convention:

```cpp
// f = 2^L, where L = log2(f)
// D = number of active dimensions (1/2/3)
//
// child_index bit layout (LSB-first, interleaved by dimension):
//   bit (b*D + d) == (ox[d] >> b) & 1
//
// Example (3D, f=4 => L=2, D=3):
//   bit 0 = ox1 bit0, bit 1 = ox2 bit0, bit 2 = ox3 bit0
//   bit 3 = ox1 bit1, bit 4 = ox2 bit1, bit 5 = ox3 bit1
//
// Offsets are recovered by de-interleaving the bits back into ox1/ox2/ox3.
```

So for `f=2`, this reduces to the original 3-bit child index encoding. For `f=4`,
each offset is 0/1/2/3 and the order matches the restart file’s Morton ordering.

Parent logical location is computed by integer-dividing the fine logical
coordinates by `f` and decrementing the level by `log2(f)`:

```cpp
// downsampler.cpp
parent.level = fine.level - log2(f);
parent.lx1   = fine.lx1 / f;
parent.lx2   = fine.lx2 / f;
parent.lx3   = fine.lx3 / f;
```

The validation check that enforces ordering is effectively:

```cpp
// downsampler.cpp (conceptual)
child.lx1 == f * parent.lx1 + ox1
child.lx2 == f * parent.lx2 + ox2
child.lx3 == f * parent.lx3 + ox3
```

### Step 2.5: Create a *new* MPI distribution for coarse MeshBlocks

Even if the fine restart was read with some distribution, the downsampler
constructs a separate distribution for the smaller coarse set:

```cpp
MPIDistribution coarse_dist(nranks, coarse_nmb_total_);
coarse_dist.SetupDistribution();
```

Each rank knows the coarse gid range it owns:

- `coarse_gid_start = coarse_dist.GetStartingMeshBlockID(my_rank)`
- `coarse_nmb_thisrank = coarse_dist.GetNumMeshBlocks(my_rank)`

### Step 2.6: For each local fine MeshBlock, compute one downsampled “chunk”

The main loop (in `DownsampleToBinary()`) iterates over fine MeshBlocks owned by
the current rank:

```cpp
for (int lf = 0; lf < fine_nmb_thisrank; ++lf) {
  int fine_gid   = fine_gid_start + lf;
  int coarse_gid = fine_gid / nfine_per_coarse_;
  int child_idx  = fine_gid % nfine_per_coarse_;

  std::vector<float> chunk;
  ComputeDownsampledChunk(lf, &chunk);
  // ... send or locally insert the chunk into the owning coarse MeshBlock ...
}
```

#### What a “chunk” is

`chunk` is **one child’s contribution** to its parent coarse MeshBlock:

- It has **1/f** the cell count in each active dimension:
  - `nx1_sub = nx1/f`, `nx2_sub = nx2/f` (if 2D/3D), `nx3_sub = nx3/f` (if 3D)
- It stores output variables in *var-major* order, matching AthenaK’s
  `mhd_w_bcc` binary output:
  - `ideal` EOS (8 vars): `dens, velx, vely, velz, eint, bcc1, bcc2, bcc3`
  - `isothermal` EOS (7 vars): `dens, velx, vely, velz, bcc1, bcc2, bcc3`

For turbulence output (`DownsampleTurbulenceToBinary()`), the chunk uses the same
layout but stores `nforce` variables (typically 3) matching AthenaK `turb_force`:
`force1, force2, force3`.

Layout:

```text
chunk = [var0 subgrid][var1 subgrid]...[var(N-1) subgrid]
each subgrid is flattened (i,j,k) via CellIndex(i,j,k,nx1_sub,nx2_sub)
```

### Step 2.7: Assemble each coarse MeshBlock by inserting its child chunks

Each fine MeshBlock computes `coarse_gid` and identifies the rank that owns that
coarse gid (using `FindRankForGid(...)` + the coarse distribution’s gid ranges).

If the coarse MeshBlock is owned locally, the child chunk is copied directly
into the correct sub-region of the coarse MeshBlock’s output buffer using child
offsets:

```cpp
// downsampler.cpp (simplified)
int ox1 = ChildOffsetX1(child_index);
int ox2 = multi_d_ ? ChildOffsetX2(child_index) : 0;
int ox3 = three_d_ ? ChildOffsetX3(child_index) : 0;

int i_off = ox1 * nx1_sub;
int j_off = ox2 * nx2_sub;
int k_off = ox3 * nx3_sub;

out_var[CellIndex(i_off + i, j_off + j, k_off + k, nx1, nx2)] =
  chunk_var[CellIndex(i, j, k, nx1_sub, nx2_sub)];
```

If the owning coarse rank is remote, the rank packs the chunk into a send buffer
as:

1. a per-message `int32_t n_chunks`
2. repeated records of:
   - `ChunkHeader{ coarse_gid, child_index }`
   - `float data[out_vars * sub_cells]`

Then it sends the message via `MPI_Isend`.

The owning rank receives remote chunks in a loop using `MPI_Probe`/`MPI_Recv`,
unpacks them, and performs the same insertion.

### Step 2.8: Validate completeness

After all local + remote chunks are inserted, each rank verifies that every
owned coarse MeshBlock received exactly `nfine_per_coarse_` children.

### Step 2.9: Update parameters and write Athena `.bin`

The parameter string embedded in the output is updated so that `<mesh>` has the
coarsened global sizes (`nx1/nx2/nx3`) using:

- `ParameterParser::ReplaceMeshNx(...)` in
  `github_version_regridding/src/parameter_parser.cpp`

Finally, `BinaryWriter::WriteBinaryFile(...)` writes:

- a text preheader describing variable names and sizes,
- the updated parameter string,
- one fixed-size binary record per coarse MeshBlock containing:
  - indices + logical location,
  - physical bounds,
  - raw `float` array for the `out_vars` variables on that MeshBlock.

## 3) Step-by-step logic inside `ComputeDownsampledChunk()`

`ComputeDownsampledChunk(local_fine_mb, &chunk)` performs the *math* of the
downsampling for one fine MeshBlock.

### Step 3.1: Compute sub-grid sizes

Given a fine MeshBlock of size `(nx1,nx2,nx3)`, the sub-chunk size is
`(nx1/f, nx2/f, nx3/f)` in active dimensions.

### Step 3.2: Restrict face-centered magnetic fields (area averages)

The restart stores magnetic field on faces (`B1f/B2f/B3f`).

The downsampler first constructs coarse face fields (`b1c/b2c/b3c`) on the
sub-grid by *averaging the fine faces that cover the same coarse face area*.

Example: restricting `B1f` to `b1c` (x1-faces). In 3D a coarse x1-face is covered
by an `f×f` block of fine x1-faces in the (x2,x3) plane:

```cpp
// downsampler.cpp (3D case)
Real sum = 0.0;
for (int kk = 0; kk < f; ++kk) {
  for (int jj = 0; jj < f; ++jj) {
    sum += b1f[FaceIndexX1(fi, fj + jj, fk + kk, nout1p1, nout2)];
  }
}
sum *= 1.0 / (f*f);   // average over f^2 fine faces
b1c[FaceIndexX1(i, j, k, nx1_sub + 1, nx2_sub)] = sum;
```

2D and 1D collapse the averaging in the inactive dimensions.

### Step 3.3: Restrict cell-centered conserved quantities (volume averages)

The restart’s MHD cell-centered array is read as var-major blocks:

```text
mhd[var][cell]  where cell = CellIndex(i,j,k,nout1,nout2)
```

For each coarse cell `(i,j,k)` in the sub-grid, the downsampler averages the
corresponding `f^D` fine cells:

- 1D: average `f` cells
- 2D: average `f^2` cells
- 3D: average `f^3` cells

Example (3D, `f=2`):

```cpp
rho = 0.125 * (rho000 + rho100 + rho010 + rho110 + rho001 + rho101 + rho011 + rho111);
m1  = 0.125 * (m1000 + ...);
m2  = 0.125 * (...);
m3  = 0.125 * (...);
eng = 0.125 * (...); // only if energy exists
```

For `f=4`, the same idea applies but with `1/64` and a 4×4×4 set of fine cells.

### Step 3.4: Convert to primitive variables + compute internal energy

The output is **not** written as conserved variables; the chunk stores:

- density `rho`
- velocity components `v = m/rho` (written as `velx/vely/velz`)
- internal energy density `eint` (**ideal EOS only**)
- cell-centered magnetic field `bcc1/bcc2/bcc3`

Cell-centered magnetic fields are computed by averaging the two neighboring
coarse faces around the cell center:

```cpp
bx = 0.5 * (b1c[b1_left] + b1c[b1_right]);
by = 0.5 * (b2c[b2_bot]  + b2c[b2_top]);
bz = 0.5 * (b3c[b3_bot]  + b3c[b3_top]);
```

Internal energy calculation (ideal EOS only):

- **Ideal EOS**:
  - apply density floor: `rho = max(rho, dfloor)` (momentum/energy unchanged)
  - kinetic energy `ke = 0.5*(m·m)/rho`
  - magnetic energy `me = 0.5*(B·B)`
  - internal energy `eint = eng - ke - me`
  - apply EOS floors (matching AthenaK `SingleC2P_IdealMHD`):
    - energy floor via `pfloor`
    - temperature floor via `tfloor`
    - temperature ceiling via `tmax`
    - entropy floor via `sfloor`
  - write the floored `eint` to the output

Pressure is **not** written; if needed, compute it from the output as:

- `ideal` EOS: `press = (gamma - 1) * eint`
- `isothermal` EOS: `press = cs^2 * rho`

## 4) What changes in the output (data + metadata)

Compared to the input restart:

- **Resolution**: `nx*` is divided by `f` in each active dimension.
- **Cell sizes**: `dx*` is multiplied by `f` in each active dimension.
- **MeshBlocks**: total MeshBlock count is divided by `f^D`, and each coarse
  MeshBlock is formed by merging `f^D` fine MeshBlocks.
- **Variables written** (AthenaK `mhd_w_bcc` ordering):
  - `ideal` EOS (8 vars): `dens, velx, vely, velz, eint, bcc1, bcc2, bcc3`
  - `isothermal` EOS (7 vars): `dens, velx, vely, velz, bcc1, bcc2, bcc3`
- **Turbulence variables written** (separate `.bin`, AthenaK `turb_force` ordering):
  - `force1, force2, force3`
- **Parameter string**: `<mesh>` `nx1/nx2/nx3` is rewritten to match the coarse mesh.
- **Root level**: decremented by `log2(f)`.

## 5) What changed in *this codebase* to support downsampling

Relative to the in-tree reader copy in
`athenak_restart_multiphase/restart_reader/src/`, the standalone repo adds:

- New implementation files:
  - `github_version_regridding/src/downsampler.cpp` + `.hpp`
  - `github_version_regridding/src/binary_writer.cpp` + `.hpp`
- CLI wiring in `github_version_regridding/src/main.cpp`:
  - adds `--downsample-bin` / `--downsample-turb-bin` / `--downsample-factor`
  - calls `Downsampler::DownsampleToBinary()` and/or
    `Downsampler::DownsampleTurbulenceToBinary()`
- Extended parameter utilities in `github_version_regridding/src/parameter_parser.*`:
  - `ExtractStringInBlock(...)`, `ExtractRealInBlock(...)` (EOS parsing)
  - `ReplaceMeshNx(...)` (rewrite `<mesh>` sizes for the coarse output)
