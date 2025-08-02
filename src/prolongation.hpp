#ifndef PROLONGATION_HPP_
#define PROLONGATION_HPP_
//========================================================================================
// Prolongation operators adapted from AthenaK source code
// Reference: /mnt/g/athenak_restart_multiphase/src/mesh/prolongation.hpp
// 
// Standalone implementation without Kokkos dependencies
//========================================================================================

#include <cmath>
#include <algorithm>
#include "common.hpp"

// Helper macro for sign function
#define SIGN(x) ((x) > 0.0 ? 1.0 : ((x) < 0.0 ? -1.0 : 0.0))

//----------------------------------------------------------------------------------------
//! \fn ProlongCCStandalone()
//! \brief 2nd-order (piecewise-linear) prolongation operator for cell-centered variables
//! This is a standalone version of AthenaK's ProlongCC without Kokkos
//!
//! Interpolates from coarse cell at (k,j,i) to 8 fine cells starting at (fk,fj,fi)
//! For 2D problems, only 4 fine cells are filled (fk+0)
//! For 1D problems, only 2 fine cells are filled (fk+0, fj+0)

inline void ProlongCCStandalone(
    const Real* coarse_data,      // Input: coarse grid data
    Real* fine_data,              // Output: fine grid data  
    const int ck, const int cj, const int ci,    // Coarse cell indices
    const int fk, const int fj, const int fi,    // Fine cell starting indices
    const int cnx1, const int cnx2, const int cnx3,  // Coarse grid dimensions
    const int fnx1, const int fnx2, const int fnx3,  // Fine grid dimensions
    const bool multi_d, const bool three_d)          // Dimensionality flags
{
    // Helper lambda to access 3D array as 1D
    auto idx3d = [](int k, int j, int i, int nx1, int nx2) {
        return k * nx2 * nx1 + j * nx1 + i;
    };
    
    // Get value at coarse cell
    Real val = coarse_data[idx3d(ck, cj, ci, cnx1, cnx2)];
    
    // Calculate x1-gradient using min-mod limiter
    Real dl = 0.0, dr = 0.0, dvar1 = 0.0;
    if (ci > 0 && ci < cnx1-1) {
        dl = val - coarse_data[idx3d(ck, cj, ci-1, cnx1, cnx2)];
        dr = coarse_data[idx3d(ck, cj, ci+1, cnx1, cnx2)] - val;
        dvar1 = 0.125 * (SIGN(dl) + SIGN(dr)) * std::min(std::abs(dl), std::abs(dr));
    }
    
    // Calculate x2-gradient using min-mod limiter
    Real dvar2 = 0.0;
    if (multi_d && cj > 0 && cj < cnx2-1) {
        dl = val - coarse_data[idx3d(ck, cj-1, ci, cnx1, cnx2)];
        dr = coarse_data[idx3d(ck, cj+1, ci, cnx1, cnx2)] - val;
        dvar2 = 0.125 * (SIGN(dl) + SIGN(dr)) * std::min(std::abs(dl), std::abs(dr));
    }
    
    // Calculate x3-gradient using min-mod limiter
    Real dvar3 = 0.0;
    if (three_d && ck > 0 && ck < cnx3-1) {
        dl = val - coarse_data[idx3d(ck-1, cj, ci, cnx1, cnx2)];
        dr = coarse_data[idx3d(ck+1, cj, ci, cnx1, cnx2)] - val;
        dvar3 = 0.125 * (SIGN(dl) + SIGN(dr)) * std::min(std::abs(dl), std::abs(dr));
    }
    
    // Interpolate to the finer grid
    fine_data[idx3d(fk, fj, fi, fnx1, fnx2)] = val - dvar1 - dvar2 - dvar3;
    fine_data[idx3d(fk, fj, fi+1, fnx1, fnx2)] = val + dvar1 - dvar2 - dvar3;
    
    if (multi_d) {
        fine_data[idx3d(fk, fj+1, fi, fnx1, fnx2)] = val - dvar1 + dvar2 - dvar3;
        fine_data[idx3d(fk, fj+1, fi+1, fnx1, fnx2)] = val + dvar1 + dvar2 - dvar3;
    }
    
    if (three_d) {
        fine_data[idx3d(fk+1, fj, fi, fnx1, fnx2)] = val - dvar1 - dvar2 + dvar3;
        fine_data[idx3d(fk+1, fj, fi+1, fnx1, fnx2)] = val + dvar1 - dvar2 + dvar3;
        fine_data[idx3d(fk+1, fj+1, fi, fnx1, fnx2)] = val - dvar1 + dvar2 + dvar3;
        fine_data[idx3d(fk+1, fj+1, fi+1, fnx1, fnx2)] = val + dvar1 + dvar2 + dvar3;
    }
}

//----------------------------------------------------------------------------------------
//! \fn ProlongFCSharedX1FaceStandalone()
//! \brief 2nd-order prolongation for face-centered variables on shared X1-faces

inline void ProlongFCSharedX1FaceStandalone(
    const Real* coarse_x1f,       // Input: coarse x1-face data
    Real* fine_x1f,               // Output: fine x1-face data
    const int ck, const int cj, const int ci,    // Coarse face indices
    const int fk, const int fj, const int fi,    // Fine face starting indices
    const int cnx1, const int cnx2, const int cnx3,  // Coarse grid dimensions
    const int fnx1, const int fnx2, const int fnx3,  // Fine grid dimensions
    const bool multi_d, const bool three_d)
{
    // Helper lambda to access 3D face array as 1D
    // Note: x1-faces have dimension (nx3, nx2, nx1+1)
    auto idx3d_x1f = [](int k, int j, int i, int nx1, int nx2) {
        return k * nx2 * (nx1+1) + j * (nx1+1) + i;
    };
    
    Real val = coarse_x1f[idx3d_x1f(ck, cj, ci, cnx1, cnx2)];
    
    // Prolongate by interpolating in x2/x3 directions
    Real dvar2 = 0.0;
    if (multi_d && cj > 0 && cj < cnx2-1) {
        Real dl = val - coarse_x1f[idx3d_x1f(ck, cj-1, ci, cnx1, cnx2)];
        Real dr = coarse_x1f[idx3d_x1f(ck, cj+1, ci, cnx1, cnx2)] - val;
        dvar2 = 0.125 * (SIGN(dl) + SIGN(dr)) * std::min(std::abs(dl), std::abs(dr));
    }
    
    Real dvar3 = 0.0;
    if (three_d && ck > 0 && ck < cnx3-1) {
        Real dl = val - coarse_x1f[idx3d_x1f(ck-1, cj, ci, cnx1, cnx2)];
        Real dr = coarse_x1f[idx3d_x1f(ck+1, cj, ci, cnx1, cnx2)] - val;
        dvar3 = 0.125 * (SIGN(dl) + SIGN(dr)) * std::min(std::abs(dl), std::abs(dr));
    }
    
    fine_x1f[idx3d_x1f(fk, fj, fi, fnx1, fnx2)] = val - dvar2 - dvar3;
    if (multi_d) {
        fine_x1f[idx3d_x1f(fk, fj+1, fi, fnx1, fnx2)] = val + dvar2 - dvar3;
    }
    if (three_d) {
        fine_x1f[idx3d_x1f(fk+1, fj, fi, fnx1, fnx2)] = val - dvar2 + dvar3;
        fine_x1f[idx3d_x1f(fk+1, fj+1, fi, fnx1, fnx2)] = val + dvar2 + dvar3;
    }
}

//----------------------------------------------------------------------------------------
//! \fn ProlongFCSharedX2FaceStandalone()
//! \brief 2nd-order prolongation for face-centered variables on shared X2-faces

inline void ProlongFCSharedX2FaceStandalone(
    const Real* coarse_x2f,       // Input: coarse x2-face data
    Real* fine_x2f,               // Output: fine x2-face data
    const int ck, const int cj, const int ci,    // Coarse face indices
    const int fk, const int fj, const int fi,    // Fine face starting indices
    const int cnx1, const int cnx2, const int cnx3,  // Coarse grid dimensions
    const int fnx1, const int fnx2, const int fnx3,  // Fine grid dimensions
    const bool three_d)
{
    // Helper lambda to access 3D face array as 1D
    // Note: x2-faces have dimension (nx3, nx2+1, nx1)
    auto idx3d_x2f = [](int k, int j, int i, int nx1, int nx2) {
        return k * (nx2+1) * nx1 + j * nx1 + i;
    };
    
    Real val = coarse_x2f[idx3d_x2f(ck, cj, ci, cnx1, cnx2)];
    
    // Prolongate by interpolating in x1/x3 directions
    Real dl = 0.0, dr = 0.0, dvar1 = 0.0;
    if (ci > 0 && ci < cnx1-1) {
        dl = val - coarse_x2f[idx3d_x2f(ck, cj, ci-1, cnx1, cnx2)];
        dr = coarse_x2f[idx3d_x2f(ck, cj, ci+1, cnx1, cnx2)] - val;
        dvar1 = 0.125 * (SIGN(dl) + SIGN(dr)) * std::min(std::abs(dl), std::abs(dr));
    }
    
    Real dvar3 = 0.0;
    if (three_d && ck > 0 && ck < cnx3-1) {
        dl = val - coarse_x2f[idx3d_x2f(ck-1, cj, ci, cnx1, cnx2)];
        dr = coarse_x2f[idx3d_x2f(ck+1, cj, ci, cnx1, cnx2)] - val;
        dvar3 = 0.125 * (SIGN(dl) + SIGN(dr)) * std::min(std::abs(dl), std::abs(dr));
    }
    
    fine_x2f[idx3d_x2f(fk, fj, fi, fnx1, fnx2)] = val - dvar1 - dvar3;
    fine_x2f[idx3d_x2f(fk, fj, fi+1, fnx1, fnx2)] = val + dvar1 - dvar3;
    if (three_d) {
        fine_x2f[idx3d_x2f(fk+1, fj, fi, fnx1, fnx2)] = val - dvar1 + dvar3;
        fine_x2f[idx3d_x2f(fk+1, fj, fi+1, fnx1, fnx2)] = val + dvar1 + dvar3;
    }
}

//----------------------------------------------------------------------------------------
//! \fn ProlongFCSharedX3FaceStandalone()
//! \brief 2nd-order prolongation for face-centered variables on shared X3-faces

inline void ProlongFCSharedX3FaceStandalone(
    const Real* coarse_x3f,       // Input: coarse x3-face data
    Real* fine_x3f,               // Output: fine x3-face data
    const int ck, const int cj, const int ci,    // Coarse face indices
    const int fk, const int fj, const int fi,    // Fine face starting indices
    const int cnx1, const int cnx2, const int cnx3,  // Coarse grid dimensions
    const int fnx1, const int fnx2, const int fnx3,  // Fine grid dimensions
    const bool multi_d)
{
    // Helper lambda to access 3D face array as 1D
    // Note: x3-faces have dimension (nx3+1, nx2, nx1)
    auto idx3d_x3f = [](int k, int j, int i, int nx1, int nx2) {
        return k * nx2 * nx1 + j * nx1 + i;
    };
    
    Real val = coarse_x3f[idx3d_x3f(ck, cj, ci, cnx1, cnx2)];
    
    // Prolongate by interpolating in x1/x2 directions
    Real dl = 0.0, dr = 0.0, dvar1 = 0.0;
    if (ci > 0 && ci < cnx1-1) {
        dl = val - coarse_x3f[idx3d_x3f(ck, cj, ci-1, cnx1, cnx2)];
        dr = coarse_x3f[idx3d_x3f(ck, cj, ci+1, cnx1, cnx2)] - val;
        dvar1 = 0.125 * (SIGN(dl) + SIGN(dr)) * std::min(std::abs(dl), std::abs(dr));
    }
    
    Real dvar2 = 0.0;
    if (multi_d && cj > 0 && cj < cnx2-1) {
        dl = val - coarse_x3f[idx3d_x3f(ck, cj-1, ci, cnx1, cnx2)];
        dr = coarse_x3f[idx3d_x3f(ck, cj+1, ci, cnx1, cnx2)] - val;
        dvar2 = 0.125 * (SIGN(dl) + SIGN(dr)) * std::min(std::abs(dl), std::abs(dr));
    }
    
    fine_x3f[idx3d_x3f(fk, fj, fi, fnx1, fnx2)] = val - dvar1 - dvar2;
    fine_x3f[idx3d_x3f(fk, fj, fi+1, fnx1, fnx2)] = val + dvar1 - dvar2;
    if (multi_d) {
        fine_x3f[idx3d_x3f(fk, fj+1, fi, fnx1, fnx2)] = val - dvar1 + dvar2;
        fine_x3f[idx3d_x3f(fk, fj+1, fi+1, fnx1, fnx2)] = val + dvar1 + dvar2;
    }
}

//----------------------------------------------------------------------------------------
//! \fn ProlongFCInternalStandalone()
//! \brief 2nd-order prolongation for face-centered variables on internal edges
//! of new fine cells within one coarse cell using divergence-preserving interpolation
//! Based on Toth & Roe, JCP 180, 736 (2002)

inline void ProlongFCInternalStandalone(
    Real* fine_x1f, Real* fine_x2f, Real* fine_x3f,    // Fine face data arrays
    const int fk, const int fj, const int fi,          // Fine cell starting indices
    const int fnx1, const int fnx2, const int fnx3,    // Fine grid dimensions
    const bool three_d)
{
    // Helper lambdas to access face arrays
    auto idx3d_x1f = [](int k, int j, int i, int nx1, int nx2) {
        return k * nx2 * (nx1+1) + j * (nx1+1) + i;
    };
    auto idx3d_x2f = [](int k, int j, int i, int nx1, int nx2) {
        return k * (nx2+1) * nx1 + j * nx1 + i;
    };
    auto idx3d_x3f = [](int k, int j, int i, int nx1, int nx2) {
        return k * nx2 * nx1 + j * nx1 + i;
    };
    
    // Prolongate internal fields in 3D
    if (three_d) {
        Real Uxx = 0.0, Vyy = 0.0, Wzz = 0.0;
        Real Uxyz = 0.0, Vxyz = 0.0, Wxyz = 0.0;
        
        for (int jj = 0; jj < 2; jj++) {
            int jsgn = 2*jj - 1;
            int fjj = fj + jj, fjp = fj + 2*jj;
            for (int ii = 0; ii < 2; ii++) {
                int isgn = 2*ii - 1;
                int fii = fi + ii, fip = fi + 2*ii;
                
                Uxx += isgn * (jsgn * (fine_x2f[idx3d_x2f(fk, fjp, fii, fnx1, fnx2)] + 
                                      fine_x2f[idx3d_x2f(fk+1, fjp, fii, fnx1, fnx2)]) +
                              (fine_x3f[idx3d_x3f(fk+2, fjj, fii, fnx1, fnx2)] - 
                               fine_x3f[idx3d_x3f(fk, fjj, fii, fnx1, fnx2)]));
                
                Vyy += jsgn * ((fine_x3f[idx3d_x3f(fk+2, fjj, fii, fnx1, fnx2)] - 
                               fine_x3f[idx3d_x3f(fk, fjj, fii, fnx1, fnx2)]) +
                              isgn * (fine_x1f[idx3d_x1f(fk, fjj, fip, fnx1, fnx2)] + 
                                     fine_x1f[idx3d_x1f(fk+1, fjj, fip, fnx1, fnx2)]));
                
                Wzz += isgn * (fine_x1f[idx3d_x1f(fk+1, fjj, fip, fnx1, fnx2)] - 
                              fine_x1f[idx3d_x1f(fk, fjj, fip, fnx1, fnx2)]) +
                       jsgn * (fine_x2f[idx3d_x2f(fk+1, fjp, fii, fnx1, fnx2)] - 
                              fine_x2f[idx3d_x2f(fk, fjp, fii, fnx1, fnx2)]);
                
                Uxyz += isgn * jsgn * (fine_x1f[idx3d_x1f(fk+1, fjj, fip, fnx1, fnx2)] - 
                                      fine_x1f[idx3d_x1f(fk, fjj, fip, fnx1, fnx2)]);
                Vxyz += isgn * jsgn * (fine_x2f[idx3d_x2f(fk+1, fjp, fii, fnx1, fnx2)] - 
                                      fine_x2f[idx3d_x2f(fk, fjp, fii, fnx1, fnx2)]);
                Wxyz += isgn * jsgn * (fine_x3f[idx3d_x3f(fk+2, fjj, fii, fnx1, fnx2)] - 
                                      fine_x3f[idx3d_x3f(fk, fjj, fii, fnx1, fnx2)]);
            }
        }
        
        Uxx *= 0.125;  Vyy *= 0.125;  Wzz *= 0.125;
        Uxyz *= 0.0625; Vxyz *= 0.0625; Wxyz *= 0.0625;
        
        // Update internal x1-faces
        fine_x1f[idx3d_x1f(fk, fj, fi+1, fnx1, fnx2)] = 
            0.5 * (fine_x1f[idx3d_x1f(fk, fj, fi, fnx1, fnx2)] + 
                   fine_x1f[idx3d_x1f(fk, fj, fi+2, fnx1, fnx2)]) + Uxx - Vxyz - Wxyz;
        fine_x1f[idx3d_x1f(fk, fj+1, fi+1, fnx1, fnx2)] = 
            0.5 * (fine_x1f[idx3d_x1f(fk, fj+1, fi, fnx1, fnx2)] + 
                   fine_x1f[idx3d_x1f(fk, fj+1, fi+2, fnx1, fnx2)]) + Uxx - Vxyz + Wxyz;
        fine_x1f[idx3d_x1f(fk+1, fj, fi+1, fnx1, fnx2)] = 
            0.5 * (fine_x1f[idx3d_x1f(fk+1, fj, fi, fnx1, fnx2)] + 
                   fine_x1f[idx3d_x1f(fk+1, fj, fi+2, fnx1, fnx2)]) + Uxx + Vxyz - Wxyz;
        fine_x1f[idx3d_x1f(fk+1, fj+1, fi+1, fnx1, fnx2)] = 
            0.5 * (fine_x1f[idx3d_x1f(fk+1, fj+1, fi, fnx1, fnx2)] + 
                   fine_x1f[idx3d_x1f(fk+1, fj+1, fi+2, fnx1, fnx2)]) + Uxx + Vxyz + Wxyz;
        
        // Update internal x2-faces
        fine_x2f[idx3d_x2f(fk, fj+1, fi, fnx1, fnx2)] = 
            0.5 * (fine_x2f[idx3d_x2f(fk, fj, fi, fnx1, fnx2)] + 
                   fine_x2f[idx3d_x2f(fk, fj+2, fi, fnx1, fnx2)]) + Vyy - Uxyz - Wxyz;
        fine_x2f[idx3d_x2f(fk, fj+1, fi+1, fnx1, fnx2)] = 
            0.5 * (fine_x2f[idx3d_x2f(fk, fj, fi+1, fnx1, fnx2)] + 
                   fine_x2f[idx3d_x2f(fk, fj+2, fi+1, fnx1, fnx2)]) + Vyy - Uxyz + Wxyz;
        fine_x2f[idx3d_x2f(fk+1, fj+1, fi, fnx1, fnx2)] = 
            0.5 * (fine_x2f[idx3d_x2f(fk+1, fj, fi, fnx1, fnx2)] + 
                   fine_x2f[idx3d_x2f(fk+1, fj+2, fi, fnx1, fnx2)]) + Vyy + Uxyz - Wxyz;
        fine_x2f[idx3d_x2f(fk+1, fj+1, fi+1, fnx1, fnx2)] = 
            0.5 * (fine_x2f[idx3d_x2f(fk+1, fj, fi+1, fnx1, fnx2)] + 
                   fine_x2f[idx3d_x2f(fk+1, fj+2, fi+1, fnx1, fnx2)]) + Vyy + Uxyz + Wxyz;
        
        // Update internal x3-faces
        fine_x3f[idx3d_x3f(fk+1, fj, fi, fnx1, fnx2)] = 
            0.5 * (fine_x3f[idx3d_x3f(fk+2, fj, fi, fnx1, fnx2)] + 
                   fine_x3f[idx3d_x3f(fk, fj, fi, fnx1, fnx2)]) + Wzz - Uxyz - Vxyz;
        fine_x3f[idx3d_x3f(fk+1, fj, fi+1, fnx1, fnx2)] = 
            0.5 * (fine_x3f[idx3d_x3f(fk+2, fj, fi+1, fnx1, fnx2)] + 
                   fine_x3f[idx3d_x3f(fk, fj, fi+1, fnx1, fnx2)]) + Wzz - Uxyz + Vxyz;
        fine_x3f[idx3d_x3f(fk+1, fj+1, fi, fnx1, fnx2)] = 
            0.5 * (fine_x3f[idx3d_x3f(fk+2, fj+1, fi, fnx1, fnx2)] + 
                   fine_x3f[idx3d_x3f(fk, fj+1, fi, fnx1, fnx2)]) + Wzz + Uxyz - Vxyz;
        fine_x3f[idx3d_x3f(fk+1, fj+1, fi+1, fnx1, fnx2)] = 
            0.5 * (fine_x3f[idx3d_x3f(fk+2, fj+1, fi+1, fnx1, fnx2)] + 
                   fine_x3f[idx3d_x3f(fk, fj+1, fi+1, fnx1, fnx2)]) + Wzz + Uxyz + Vxyz;
        
    // Prolongate internal fields in 2D
    } else {
        Real tmp1 = 0.25 * (fine_x2f[idx3d_x2f(fk, fj+2, fi+1, fnx1, fnx2)] - 
                           fine_x2f[idx3d_x2f(fk, fj, fi+1, fnx1, fnx2)] -
                           fine_x2f[idx3d_x2f(fk, fj+2, fi, fnx1, fnx2)] + 
                           fine_x2f[idx3d_x2f(fk, fj, fi, fnx1, fnx2)]);
        Real tmp2 = 0.25 * (fine_x1f[idx3d_x1f(fk, fj, fi, fnx1, fnx2)] - 
                           fine_x1f[idx3d_x1f(fk, fj, fi+2, fnx1, fnx2)] -
                           fine_x1f[idx3d_x1f(fk, fj+1, fi, fnx1, fnx2)] + 
                           fine_x1f[idx3d_x1f(fk, fj+1, fi+2, fnx1, fnx2)]);
        
        fine_x1f[idx3d_x1f(fk, fj, fi+1, fnx1, fnx2)] = 
            0.5 * (fine_x1f[idx3d_x1f(fk, fj, fi, fnx1, fnx2)] + 
                   fine_x1f[idx3d_x1f(fk, fj, fi+2, fnx1, fnx2)]) + tmp1;
        fine_x1f[idx3d_x1f(fk, fj+1, fi+1, fnx1, fnx2)] = 
            0.5 * (fine_x1f[idx3d_x1f(fk, fj+1, fi, fnx1, fnx2)] + 
                   fine_x1f[idx3d_x1f(fk, fj+1, fi+2, fnx1, fnx2)]) + tmp1;
        fine_x2f[idx3d_x2f(fk, fj+1, fi, fnx1, fnx2)] = 
            0.5 * (fine_x2f[idx3d_x2f(fk, fj, fi, fnx1, fnx2)] + 
                   fine_x2f[idx3d_x2f(fk, fj+2, fi, fnx1, fnx2)]) + tmp2;
        fine_x2f[idx3d_x2f(fk, fj+1, fi+1, fnx1, fnx2)] = 
            0.5 * (fine_x2f[idx3d_x2f(fk, fj, fi+1, fnx1, fnx2)] + 
                   fine_x2f[idx3d_x2f(fk, fj+2, fi+1, fnx1, fnx2)]) + tmp2;
    }
}

#endif // PROLONGATION_HPP_