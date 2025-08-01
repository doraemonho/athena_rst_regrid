#ifndef COMMON_HPP_
#define COMMON_HPP_

#include <cstdint>

// Type definitions to match AthenaK
using Real = double;
using IOWrapperSizeT = std::uint64_t;

// MPI type definitions for Real
#if MPI_PARALLEL_ENABLED
#include <mpi.h>
#if SINGLE_PRECISION_ENABLED
#define MPI_ATHENA_REAL MPI_FLOAT
#else
#define MPI_ATHENA_REAL MPI_DOUBLE
#endif
#endif

// Physical size in a Mesh or MeshBlock
struct RegionSize {
    Real x1min, x2min, x3min;
    Real x1max, x2max, x3max;
    Real dx1, dx2, dx3;       // (uniform) grid spacing
};

// Cell indices and number of active and ghost cells in a Mesh or MeshBlock
struct RegionIndcs {
    int ng;                       // number of ghost cells
    int nx1, nx2, nx3;            // number of active cells (not including ghost zones)
    int is, ie, js, je, ks, ke;   // indices of ACTIVE cells
    int cnx1, cnx2, cnx3;         // number of active coarse cells (not including gzs)
    int cis, cie, cjs, cje, cks, cke;  // indices of ACTIVE coarse cells
};

// Logical location and level of MeshBlock
struct LogicalLocation {
    std::int32_t lx1, lx2, lx3, level;
};

// Random number generator state structure (for turbulence)
#define NTAB 32
struct RNG_State {
    int64_t idum;
    int64_t idum2;
    int64_t iy;
    int64_t iv[NTAB];
    int iset;
    double gset;
};

// Physics configuration
struct PhysicsConfig {
    bool has_hydro = false;
    bool has_mhd = false;
    bool has_radiation = false;
    bool has_turbulence = false;
    bool has_z4c = false;
    bool has_adm = false;
    
    int nhydro = 0;
    int nmhd = 0;
    int nrad = 0;
    int nforce = 3;
    int nz4c = 0;
    int nadm = 0;
    int nco = 0;  // number of compact objects
};

// Validation tolerance
constexpr Real TOLERANCE = 1e-6;

#endif // COMMON_HPP_