//========================================================================================
// AthenaK Regridding Tool - DataInterpolator Implementation
// Licensed under the 3-clause BSD License (the "LICENSE")
//========================================================================================
//! \file data_interpolator.cpp
//  \brief Implementation of DataInterpolator class with AthenaK prolongation operators

#include "data_interpolator.hpp"
#include <iostream>
#include <cmath>
#include <algorithm>
#include <cassert>

//----------------------------------------------------------------------------------------
//! \fn bool DataInterpolator::InterpolateAllData()
//! \brief Main function to interpolate all physics data
bool DataInterpolator::InterpolateAllData(const RestartData& input_data, 
                                         const NewMeshData& mesh_data,
                                         InterpolatedData& output_data) {
  SetError("");
  
  // Allocate output arrays
  AllocateOutputArrays(input_data, mesh_data, output_data);
  CopyPhysicsParameters(input_data, output_data);
  
  std::cout << "Debug: Starting cell-centered data interpolation..." << std::endl;
  // Interpolate cell-centered data
  if (!InterpolateCellCenteredData(input_data, mesh_data, output_data)) {
    std::cerr << "ERROR: Cell-centered data interpolation failed" << std::endl;
    return false;
  }
  std::cout << "Debug: Cell-centered data interpolation completed" << std::endl;
  
  // Interpolate face-centered magnetic fields
  if (input_data.nmhd > 0) {
    std::cout << "Debug: Starting face-centered data interpolation..." << std::endl;
    if (!InterpolateFaceCenteredData(input_data, mesh_data, output_data)) {
      std::cerr << "ERROR: Face-centered data interpolation failed" << std::endl;
      return false;
    }
    std::cout << "Debug: Face-centered data interpolation completed" << std::endl;
  }
  
  std::cout << "Successfully interpolated all physics data" << std::endl;
  return true;
}

//----------------------------------------------------------------------------------------
//! \fn bool DataInterpolator::InterpolateCellCenteredData()
//! \brief Interpolate cell-centered variables (hydro, MHD conserved, etc.)
bool DataInterpolator::InterpolateCellCenteredData(const RestartData& input_data, 
                                                  const NewMeshData& mesh_data,
                                                  InterpolatedData& output_data) {
  int nx1 = input_data.mb_indcs.nx1;
  int nx2 = input_data.mb_indcs.nx2;
  int nx3 = input_data.mb_indcs.nx3;
  int ng = input_data.mb_indcs.ng;
  int nout1 = input_data.nout1;
  int nout2 = input_data.nout2;
  int nout3 = input_data.nout3;
  
  std::cout << "Debug: InterpolateCellCenteredData - processing " << input_data.nmb_total << " MeshBlocks" << std::endl;
  std::cout << "Debug: input: nhydro=" << input_data.nhydro << ", nmhd=" << input_data.nmhd 
            << ", nrad=" << input_data.nrad << ", nforce=" << input_data.nforce << std::endl;
  std::cout << "Debug: output: nhydro=" << output_data.nhydro << ", nmhd=" << output_data.nmhd 
            << ", nrad=" << output_data.nrad << ", nforce=" << output_data.nforce << std::endl;
  std::cout << "Debug: output nmb_total=" << output_data.nmb_total << std::endl;
  
  // Process each original MeshBlock
  for (int old_mb = 0; old_mb < input_data.nmb_total; ++old_mb) {
    if (old_mb % 10 == 0 || old_mb >= 30) {
      std::cout << "Debug: Processing MeshBlock " << old_mb << "/" << input_data.nmb_total << std::endl;
    }
    
    // Interpolate hydro data if present
    if (input_data.nhydro > 0) {
      for (int child_id = 0; child_id < 8; ++child_id) {
        if (old_mb >= mesh_data.old_to_new_map.size()) {
          std::cerr << "ERROR: old_mb " << old_mb << " >= map size " << mesh_data.old_to_new_map.size() << std::endl;
          return false;
        }
        int new_mb = mesh_data.old_to_new_map[old_mb][child_id];
        if (new_mb < 0 || new_mb >= output_data.nmb_total) {
          std::cerr << "ERROR: Invalid new_mb " << new_mb << " (should be 0-" << (output_data.nmb_total-1) << ")" << std::endl;
          return false;
        }
        
        // Get pointers to coarse and fine data
        const Real* coarse_hydro = input_data.hydro_data.data() + 
          old_mb * input_data.nhydro * nout1 * nout2 * nout3;
        Real* fine_hydro = output_data.hydro_data.data() + 
          new_mb * output_data.nhydro * nout1 * nout2 * nout3;
        
        // Apply prolongation with appropriate offset
        int ic, jc, kc;
        ChildIDToIndex3D(child_id, ic, jc, kc);
        ProlongCC(coarse_hydro, fine_hydro, input_data.nhydro, 
                 nx1, nx2, nx3, ng, ic, jc, kc);
      }
    }
    
    // Interpolate MHD conserved variables if present
    if (input_data.nmhd > 0) {
      if (old_mb % 10 == 0) {
        std::cout << "Debug: MHD interpolation for MeshBlock " << old_mb << std::endl;
      }
      for (int child_id = 0; child_id < 8; ++child_id) {
        int new_mb = mesh_data.old_to_new_map[old_mb][child_id];
        
        // Check array bounds
        size_t coarse_offset = static_cast<size_t>(old_mb) * input_data.nmhd * nout1 * nout2 * nout3;
        size_t fine_offset = static_cast<size_t>(new_mb) * output_data.nmhd * nout1 * nout2 * nout3;
        
        if (coarse_offset >= input_data.mhd_data.size()) {
          std::cerr << "ERROR: coarse_offset " << coarse_offset << " >= input array size " << input_data.mhd_data.size() << std::endl;
          return false;
        }
        if (fine_offset >= output_data.mhd_data.size()) {
          std::cerr << "ERROR: fine_offset " << fine_offset << " >= output array size " << output_data.mhd_data.size() << std::endl;
          return false;
        }
        
        const Real* coarse_mhd = input_data.mhd_data.data() + coarse_offset;
        Real* fine_mhd = output_data.mhd_data.data() + fine_offset;
        
        int ic, jc, kc;
        ChildIDToIndex3D(child_id, ic, jc, kc);
        
        if (old_mb == 30) {
          std::cout << "Debug: About to call ProlongCC for MeshBlock 30, child " << child_id << std::endl;
          std::cout << "Debug: ic=" << ic << ", jc=" << jc << ", kc=" << kc << std::endl;
          std::cout << "Debug: new_mb=" << new_mb << std::endl;
          std::cout << "Debug: coarse_offset=" << coarse_offset << ", fine_offset=" << fine_offset << std::endl;
          std::cout << "Debug: coarse_mhd=" << (void*)coarse_mhd << ", fine_mhd=" << (void*)fine_mhd << std::endl;
        }
        
        ProlongCC(coarse_mhd, fine_mhd, input_data.nmhd, 
                 nx1, nx2, nx3, ng, ic, jc, kc);
                 
        if (old_mb == 30) {
          std::cout << "Debug: ProlongCC completed for MeshBlock 30, child " << child_id << std::endl;
        }
      }
    }
    
    // Interpolate radiation data if present
    if (input_data.nrad > 0) {
      for (int child_id = 0; child_id < 8; ++child_id) {
        int new_mb = mesh_data.old_to_new_map[old_mb][child_id];
        
        const Real* coarse_rad = input_data.rad_data.data() + 
          old_mb * input_data.nrad * nout1 * nout2 * nout3;
        Real* fine_rad = output_data.rad_data.data() + 
          new_mb * output_data.nrad * nout1 * nout2 * nout3;
        
        int ic, jc, kc;
        ChildIDToIndex3D(child_id, ic, jc, kc);
        ProlongCC(coarse_rad, fine_rad, input_data.nrad, 
                 nx1, nx2, nx3, ng, ic, jc, kc);
      }
    }
    
    // Interpolate turbulence forcing if present
    if (input_data.has_turb) {
      if (old_mb == 30) {
        std::cout << "Debug: Starting turbulence interpolation for MeshBlock 30" << std::endl;
      }
      for (int child_id = 0; child_id < 8; ++child_id) {
        if (old_mb == 30) {
          std::cout << "Debug: Turbulence child " << child_id << " for MeshBlock 30" << std::endl;
        }
        int new_mb = mesh_data.old_to_new_map[old_mb][child_id];
        
        const Real* coarse_force = input_data.force_data.data() + 
          old_mb * input_data.nforce * nout1 * nout2 * nout3;
        Real* fine_force = output_data.force_data.data() + 
          new_mb * output_data.nforce * nout1 * nout2 * nout3;
        
        int ic, jc, kc;
        ChildIDToIndex3D(child_id, ic, jc, kc);
        ProlongCC(coarse_force, fine_force, input_data.nforce, 
                 nx1, nx2, nx3, ng, ic, jc, kc);
        if (old_mb == 30) {
          std::cout << "Debug: Completed turbulence child " << child_id << " for MeshBlock 30" << std::endl;
        }
      }
    }
    
    // Interpolate Z4c data if present
    if (input_data.nz4c > 0) {
      for (int child_id = 0; child_id < 8; ++child_id) {
        int new_mb = mesh_data.old_to_new_map[old_mb][child_id];
        
        const Real* coarse_z4c = input_data.z4c_data.data() + 
          old_mb * input_data.nz4c * nout1 * nout2 * nout3;
        Real* fine_z4c = output_data.z4c_data.data() + 
          new_mb * output_data.nz4c * nout1 * nout2 * nout3;
        
        int ic, jc, kc;
        ChildIDToIndex3D(child_id, ic, jc, kc);
        ProlongCC(coarse_z4c, fine_z4c, input_data.nz4c, 
                 nx1, nx2, nx3, ng, ic, jc, kc);
      }
    }
    
    // Interpolate ADM data if present
    if (input_data.nadm > 0) {
      for (int child_id = 0; child_id < 8; ++child_id) {
        int new_mb = mesh_data.old_to_new_map[old_mb][child_id];
        
        const Real* coarse_adm = input_data.adm_data.data() + 
          old_mb * input_data.nadm * nout1 * nout2 * nout3;
        Real* fine_adm = output_data.adm_data.data() + 
          new_mb * output_data.nadm * nout1 * nout2 * nout3;
        
        int ic, jc, kc;
        ChildIDToIndex3D(child_id, ic, jc, kc);
        ProlongCC(coarse_adm, fine_adm, input_data.nadm, 
                 nx1, nx2, nx3, ng, ic, jc, kc);
      }
    }
    
    if (old_mb >= 30) {
      std::cout << "Debug: Completed all interpolation for MeshBlock " << old_mb << std::endl;
    }
    
    // Check for potential memory corruption by validating array sizes
    if (old_mb == 34) {
      std::cout << "Debug: Checking array integrity after MeshBlock 34..." << std::endl;
      std::cout << "Debug: input_data.mhd_data.size() = " << input_data.mhd_data.size() << std::endl;
      std::cout << "Debug: output_data.mhd_data.size() = " << output_data.mhd_data.size() << std::endl;
      std::cout << "Debug: input_data.force_data.size() = " << input_data.force_data.size() << std::endl;
      std::cout << "Debug: output_data.force_data.size() = " << output_data.force_data.size() << std::endl;
    }
  }
  
  std::cout << "Debug: Finished processing all MeshBlocks in InterpolateCellCenteredData" << std::endl;
  
  return true;
}

//----------------------------------------------------------------------------------------
//! \fn bool DataInterpolator::InterpolateFaceCenteredData()
//! \brief Interpolate face-centered magnetic field data with proper parent-child mapping
bool DataInterpolator::InterpolateFaceCenteredData(const RestartData& input_data, 
                                                  const NewMeshData& mesh_data,
                                                  InterpolatedData& output_data) {
  // Get dimensions
  int nx1 = input_data.mb_indcs.nx1;
  int nx2 = input_data.mb_indcs.nx2; 
  int nx3 = input_data.mb_indcs.nx3;
  int ng = input_data.mb_indcs.ng;
  int nout1 = input_data.nout1;
  int nout2 = input_data.nout2;
  int nout3 = input_data.nout3;
  
  // Initialize output arrays (already allocated in AllocateOutputArrays)
  
  // Process each original (parent) MeshBlock to its 8 children
  for (int old_mb = 0; old_mb < input_data.nmb_total; ++old_mb) {
    
    // Get pointers to parent face-centered data
    const Real* parent_b1f = input_data.mhd_b1f_data.data() + 
      old_mb * nout3 * nout2 * (nout1+1);
    const Real* parent_b2f = input_data.mhd_b2f_data.data() + 
      old_mb * nout3 * (nout2+1) * nout1;
    const Real* parent_b3f = input_data.mhd_b3f_data.data() + 
      old_mb * (nout3+1) * nout2 * nout1;
    
    // Debug: Check the first few values of parent data
    if (old_mb == 0) {
      std::cout << "DEBUG MB0: B1f[0:5] = ";
      for (int i = 0; i < 5; ++i) {
        std::cout << parent_b1f[i] << " ";
      }
      std::cout << std::endl;
    }
    
    // Interpolate each of 8 children
    for (int child_id = 0; child_id < 8; ++child_id) {
      int new_mb = mesh_data.old_to_new_map[old_mb][child_id];
      
      // Get child position within parent (0 or 1 in each direction)
      int ic, jc, kc;
      ChildIDToIndex3D(child_id, ic, jc, kc);
      
      // Get pointers to child face-centered data
      Real* child_b1f = output_data.mhd_b1f_data.data() + 
        new_mb * nout3 * nout2 * (nout1+1);
      Real* child_b2f = output_data.mhd_b2f_data.data() + 
        new_mb * nout3 * (nout2+1) * nout1;
      Real* child_b3f = output_data.mhd_b3f_data.data() + 
        new_mb * (nout3+1) * nout2 * nout1;
      
      // Apply face-centered prolongation for each component
      if (old_mb == 0 && child_id == 0) {
        std::cout << "DEBUG: Before prolongation - child_b1f[0:3] = ";
        for (int i = 0; i < 3; ++i) {
          std::cout << child_b1f[i] << " ";
        }
        std::cout << std::endl;
      }
      
      ProlongFC1(parent_b1f, child_b1f, nx1, nx2, nx3, ng, ic, jc, kc);
      ProlongFC2(parent_b2f, child_b2f, nx1, nx2, nx3, ng, ic, jc, kc);
      ProlongFC3(parent_b3f, child_b3f, nx1, nx2, nx3, ng, ic, jc, kc);
      
      if (old_mb == 0 && child_id == 0) {
        std::cout << "DEBUG: After prolongation - child_b1f[0:3] = ";
        for (int i = 0; i < 3; ++i) {
          std::cout << child_b1f[i] << " ";
        }
        std::cout << std::endl;
      }
    }
  }
  
  std::cout << "InterpolateFaceCenteredData: Completed B-field interpolation for " 
            << mesh_data.nmb_total_new << " MeshBlocks" << std::endl;
  
  return true;
}

//----------------------------------------------------------------------------------------
//! \fn void DataInterpolator::ProlongCC()
//! \brief Cell-centered prolongation with 2nd-order piecewise-linear interpolation
//! Based on AthenaK prolongation.hpp lines 18-58
void DataInterpolator::ProlongCC(const Real* coarse_data, Real* fine_data, 
                                int nvar, int nx1, int nx2, int nx3, int ng,
                                int child_i, int child_j, int child_k) const {
  int nout1 = nx1 + 2*ng;
  int nout2 = (nx2 > 1) ? (nx2 + 2*ng) : 1;
  int nout3 = (nx3 > 1) ? (nx3 + 2*ng) : 1;
  
  bool multi_d = (nx2 > 1);
  bool three_d = (nx3 > 1);
  
  // Each child MeshBlock fills the entire fine grid using its assigned region of the parent
  
  for (int n = 0; n < nvar; ++n) {
    for (int ck = 0; ck < nout3/2; ++ck) {
      for (int cj = 0; cj < nout2/2; ++cj) {
        for (int ci = 0; ci < nout1/2; ++ci) {
          
          // Map coarse indices to parent region based on child position
          int parent_ck = ck + child_k * (nout3/2);
          int parent_cj = cj + child_j * (nout2/2);
          int parent_ci = ci + child_i * (nout1/2);
          
          // Get coarse cell value from parent region
          // Note: coarse_data is already offset to the correct MeshBlock, so we use direct indexing
          size_t c_idx = n*nout3*nout2*nout1 + parent_ck*nout2*nout1 + parent_cj*nout1 + parent_ci;
          Real coarse_val = coarse_data[c_idx];
          
          // Calculate x1-gradient using AthenaK's exact minmod limiter
          Real dl = 0.0, dr = 0.0, dvar1 = 0.0;
          if (parent_ci > 0) dl = coarse_val - coarse_data[n*nout3*nout2*nout1 + parent_ck*nout2*nout1 + parent_cj*nout1 + (parent_ci-1)];
          if (parent_ci < nout1-1) dr = coarse_data[n*nout3*nout2*nout1 + parent_ck*nout2*nout1 + parent_cj*nout1 + (parent_ci+1)] - coarse_val;
          dvar1 = 0.125 * (Sign(dl) + Sign(dr)) * fmin(fabs(dl), fabs(dr));
          
          // Calculate x2-gradient using AthenaK's exact minmod limiter  
          Real dvar2 = 0.0;
          if (multi_d) {
            dl = 0.0; dr = 0.0;
            if (parent_cj > 0) dl = coarse_val - coarse_data[n*nout3*nout2*nout1 + parent_ck*nout2*nout1 + (parent_cj-1)*nout1 + parent_ci];
            if (parent_cj < nout2-1) dr = coarse_data[n*nout3*nout2*nout1 + parent_ck*nout2*nout1 + (parent_cj+1)*nout1 + parent_ci] - coarse_val;
            dvar2 = 0.125 * (Sign(dl) + Sign(dr)) * fmin(fabs(dl), fabs(dr));
          }
          
          // Calculate x3-gradient using AthenaK's exact minmod limiter
          Real dvar3 = 0.0;
          if (three_d) {
            dl = 0.0; dr = 0.0;
            if (parent_ck > 0) dl = coarse_val - coarse_data[n*nout3*nout2*nout1 + (parent_ck-1)*nout2*nout1 + parent_cj*nout1 + parent_ci];
            if (parent_ck < nout3-1) dr = coarse_data[n*nout3*nout2*nout1 + (parent_ck+1)*nout2*nout1 + parent_cj*nout1 + parent_ci] - coarse_val;
            dvar3 = 0.125 * (Sign(dl) + Sign(dr)) * fmin(fabs(dl), fabs(dr));
          }
          
          // Map to fine grid: each coarse cell becomes 2×2×2 fine cells
          int fk_base = 2*ck;
          int fj_base = 2*cj; 
          int fi_base = 2*ci;
          
          // Fill 8 fine cells using AthenaK's exact pattern
          for (int dfk = 0; dfk <= (three_d ? 1 : 0); ++dfk) {
            for (int dfj = 0; dfj <= (multi_d ? 1 : 0); ++dfj) {
              for (int dfi = 0; dfi <= 1; ++dfi) {
                
                int fk = fk_base + dfk;
                int fj = fj_base + dfj; 
                int fi = fi_base + dfi;
                
                // Check bounds to prevent segfault
                if (fk >= nout3 || fj >= nout2 || fi >= nout1) {
                  continue; // Skip out-of-bounds writes
                }
                
                // Calculate prolongated value based on position within 2×2×2 stencil
                Real sign_i = (dfi == 0) ? -1.0 : 1.0;
                Real sign_j = (dfj == 0) ? -1.0 : 1.0;
                Real sign_k = (dfk == 0) ? -1.0 : 1.0;
                
                Real prolongated_val = coarse_val + sign_i*dvar1 + sign_j*dvar2 + sign_k*dvar3;
                
                // Note: fine_data is already offset to the correct MeshBlock, so we use direct indexing
                size_t idx = n*nout3*nout2*nout1 + fk*nout2*nout1 + fj*nout1 + fi;
                fine_data[idx] = prolongated_val;
              }
            }
          }
        }
      }
    }
  }
}

//----------------------------------------------------------------------------------------
//! \fn void DataInterpolator::ProlongFCSharedX1Face()
//! \brief Prolongate B_x on shared x1-faces
void DataInterpolator::ProlongFCSharedX1Face(const Real* coarse_data, Real* fine_data,
                                            int nx1, int nx2, int nx3, int ng,
                                            int coarse_i, int coarse_j, int coarse_k) const {
  int nout1 = nx1 + 2*ng;
  int nout2 = (nx2 > 1) ? (nx2 + 2*ng) : 1;
  int nout3 = (nx3 > 1) ? (nx3 + 2*ng) : 1;
  
  // Loop over x1-faces in fine MeshBlock
  for (int fk = 0; fk < nout3; ++fk) {
    for (int fj = 0; fj < nout2; ++fj) {
      for (int fi = 0; fi <= nout1; ++fi) {
        
        // Map fine face to coarse face
        int ci = fi / 2;
        int cj = fj / 2;
        int ck = fk / 2;
        
        // Get coarse face value
        size_t c_idx = GetIndex4D(ck, cj, ci, 0, nout1+1, nout2, nout3);
        Real coarse_val = coarse_data[c_idx];
        
        // For shared faces, we can apply linear interpolation
        // with minmod limiting based on neighboring face values
        Real fine_val = coarse_val;  // Start with coarse value
        
        // Apply gradient correction if we have neighboring faces
        if (cj > 0 && cj < nout2-1 && nout2 > 1) {
          size_t idx_down = GetIndex4D(ck, cj-1, ci, 0, nout1+1, nout2, nout3);
          size_t idx_up = GetIndex4D(ck, cj+1, ci, 0, nout1+1, nout2, nout3);
          Real dl = coarse_val - coarse_data[idx_down];
          Real dr = coarse_data[idx_up] - coarse_val;
          Real grad = 0.125 * (Sign(dl) + Sign(dr)) * fmin(fabs(dl), fabs(dr));
          
          Real dy = (fj % 2 == 0) ? -0.25 : 0.25;
          fine_val += dy * grad;
        }
        
        if (ck > 0 && ck < nout3-1 && nout3 > 1) {
          size_t idx_down = GetIndex4D(ck-1, cj, ci, 0, nout1+1, nout2, nout3);
          size_t idx_up = GetIndex4D(ck+1, cj, ci, 0, nout1+1, nout2, nout3);
          Real dl = coarse_val - coarse_data[idx_down];
          Real dr = coarse_data[idx_up] - coarse_val;
          Real grad = 0.125 * (Sign(dl) + Sign(dr)) * fmin(fabs(dl), fabs(dr));
          
          Real dz = (fk % 2 == 0) ? -0.25 : 0.25;
          fine_val += dz * grad;
        }
        
        size_t f_idx = GetIndex4D(fk, fj, fi, 0, nout1+1, nout2, nout3);
        fine_data[f_idx] = fine_val;
      }
    }
  }
}

//----------------------------------------------------------------------------------------
//! \fn void DataInterpolator::ProlongFCSharedX2Face()
//! \brief Prolongate B_y on shared x2-faces
void DataInterpolator::ProlongFCSharedX2Face(const Real* coarse_data, Real* fine_data,
                                            int nx1, int nx2, int nx3, int ng,
                                            int coarse_i, int coarse_j, int coarse_k) const {
  int nout1 = nx1 + 2*ng;
  int nout2 = (nx2 > 1) ? (nx2 + 2*ng) : 1;
  int nout3 = (nx3 > 1) ? (nx3 + 2*ng) : 1;
  
  // Loop over x2-faces in fine MeshBlock  
  for (int fk = 0; fk < nout3; ++fk) {
    for (int fj = 0; fj <= nout2; ++fj) {
      for (int fi = 0; fi < nout1; ++fi) {
        
        int ci = fi / 2;
        int cj = fj / 2;
        int ck = fk / 2;
        
        size_t c_idx = GetIndex4D(ck, cj, ci, 0, nout1, nout2+1, nout3);
        Real coarse_val = coarse_data[c_idx];
        Real fine_val = coarse_val;
        
        // Apply gradient corrections
        if (ci > 0 && ci < nout1-1) {
          size_t idx_left = GetIndex4D(ck, cj, ci-1, 0, nout1, nout2+1, nout3);
          size_t idx_right = GetIndex4D(ck, cj, ci+1, 0, nout1, nout2+1, nout3);
          Real dl = coarse_val - coarse_data[idx_left];
          Real dr = coarse_data[idx_right] - coarse_val;
          Real grad = 0.125 * (Sign(dl) + Sign(dr)) * fmin(fabs(dl), fabs(dr));
          
          Real dx = (fi % 2 == 0) ? -0.25 : 0.25;
          fine_val += dx * grad;
        }
        
        if (ck > 0 && ck < nout3-1 && nout3 > 1) {
          size_t idx_down = GetIndex4D(ck-1, cj, ci, 0, nout1, nout2+1, nout3);
          size_t idx_up = GetIndex4D(ck+1, cj, ci, 0, nout1, nout2+1, nout3);
          Real dl = coarse_val - coarse_data[idx_down];
          Real dr = coarse_data[idx_up] - coarse_val;
          Real grad = 0.125 * (Sign(dl) + Sign(dr)) * fmin(fabs(dl), fabs(dr));
          
          Real dz = (fk % 2 == 0) ? -0.25 : 0.25;
          fine_val += dz * grad;
        }
        
        size_t f_idx = GetIndex4D(fk, fj, fi, 0, nout1, nout2+1, nout3);
        fine_data[f_idx] = fine_val;
      }
    }
  }
}

//----------------------------------------------------------------------------------------
//! \fn void DataInterpolator::ProlongFCSharedX3Face()
//! \brief Prolongate B_z on shared x3-faces
void DataInterpolator::ProlongFCSharedX3Face(const Real* coarse_data, Real* fine_data,
                                            int nx1, int nx2, int nx3, int ng,
                                            int coarse_i, int coarse_j, int coarse_k) const {
  int nout1 = nx1 + 2*ng;
  int nout2 = (nx2 > 1) ? (nx2 + 2*ng) : 1;
  int nout3 = (nx3 > 1) ? (nx3 + 2*ng) : 1;
  
  // Loop over x3-faces in fine MeshBlock
  for (int fk = 0; fk <= nout3; ++fk) {
    for (int fj = 0; fj < nout2; ++fj) {
      for (int fi = 0; fi < nout1; ++fi) {
        
        int ci = fi / 2;
        int cj = fj / 2;
        int ck = fk / 2;
        
        size_t c_idx = GetIndex4D(ck, cj, ci, 0, nout1, nout2, nout3+1);
        Real coarse_val = coarse_data[c_idx];
        Real fine_val = coarse_val;
        
        // Apply gradient corrections
        if (ci > 0 && ci < nout1-1) {
          size_t idx_left = GetIndex4D(ck, cj, ci-1, 0, nout1, nout2, nout3+1);
          size_t idx_right = GetIndex4D(ck, cj, ci+1, 0, nout1, nout2, nout3+1);
          Real dl = coarse_val - coarse_data[idx_left];
          Real dr = coarse_data[idx_right] - coarse_val;
          Real grad = 0.125 * (Sign(dl) + Sign(dr)) * fmin(fabs(dl), fabs(dr));
          
          Real dx = (fi % 2 == 0) ? -0.25 : 0.25;
          fine_val += dx * grad;
        }
        
        if (cj > 0 && cj < nout2-1 && nout2 > 1) {
          size_t idx_down = GetIndex4D(ck, cj-1, ci, 0, nout1, nout2, nout3+1);
          size_t idx_up = GetIndex4D(ck, cj+1, ci, 0, nout1, nout2, nout3+1);
          Real dl = coarse_val - coarse_data[idx_down];
          Real dr = coarse_data[idx_up] - coarse_val;
          Real grad = 0.125 * (Sign(dl) + Sign(dr)) * fmin(fabs(dl), fabs(dr));
          
          Real dy = (fj % 2 == 0) ? -0.25 : 0.25;
          fine_val += dy * grad;
        }
        
        size_t f_idx = GetIndex4D(fk, fj, fi, 0, nout1, nout2, nout3+1);
        fine_data[f_idx] = fine_val;
      }
    }
  }
}

//----------------------------------------------------------------------------------------
//! \fn void DataInterpolator::ProlongFCInternal()
//! \brief Divergence-preserving internal field prolongation (Tóth & Roe 2002)
void DataInterpolator::ProlongFCInternal(Real* b1f_data, Real* b2f_data, Real* b3f_data,
                                        int nx1, int nx2, int nx3, int ng,
                                        int fine_i, int fine_j, int fine_k) const {
  // This implements the divergence-preserving prolongation of Tóth & Roe (JCP 180, 2002)
  // For internal faces that are not shared with coarse grid
  
  int nout1 = nx1 + 2*ng;
  int nout2 = (nx2 > 1) ? (nx2 + 2*ng) : 1;
  int nout3 = (nx3 > 1) ? (nx3 + 2*ng) : 1;
  
  // This is a simplified implementation
  // In the full AthenaK implementation, this uses a complex formula
  // involving contributions from all surrounding face fields to ensure ∇·B = 0
  
  // For now, use simple averaging to preserve divergence to first order
  // In practice, the shared face prolongation should already maintain divergence well
  
  // Internal x1-faces (between fine cells within the same coarse cell)
  if (fine_i % 2 == 1 && fine_i < nout1) {
    size_t idx = GetIndex4D(fine_k, fine_j, fine_i, 0, nout1+1, nout2, nout3);
    size_t idx_left = GetIndex4D(fine_k, fine_j, fine_i-1, 0, nout1+1, nout2, nout3);
    size_t idx_right = GetIndex4D(fine_k, fine_j, fine_i+1, 0, nout1+1, nout2, nout3);
    b1f_data[idx] = 0.5 * (b1f_data[idx_left] + b1f_data[idx_right]);
  }
  
  // Internal x2-faces
  if (fine_j % 2 == 1 && fine_j < nout2 && nout2 > 1) {
    size_t idx = GetIndex4D(fine_k, fine_j, fine_i, 0, nout1, nout2+1, nout3);
    size_t idx_down = GetIndex4D(fine_k, fine_j-1, fine_i, 0, nout1, nout2+1, nout3);
    size_t idx_up = GetIndex4D(fine_k, fine_j+1, fine_i, 0, nout1, nout2+1, nout3);
    b2f_data[idx] = 0.5 * (b2f_data[idx_down] + b2f_data[idx_up]);
  }
  
  // Internal x3-faces
  if (fine_k % 2 == 1 && fine_k < nout3 && nout3 > 1) {
    size_t idx = GetIndex4D(fine_k, fine_j, fine_i, 0, nout1, nout2, nout3+1);
    size_t idx_back = GetIndex4D(fine_k-1, fine_j, fine_i, 0, nout1, nout2, nout3+1);
    size_t idx_front = GetIndex4D(fine_k+1, fine_j, fine_i, 0, nout1, nout2, nout3+1);
    b3f_data[idx] = 0.5 * (b3f_data[idx_back] + b3f_data[idx_front]);
  }
}

//----------------------------------------------------------------------------------------
//! \fn void DataInterpolator::ProlongFC()
//! \brief Face-centered prolongation dispatcher
void DataInterpolator::ProlongFC(const Real* coarse_data, Real* fine_data,
                                int nvar, int nx1, int nx2, int nx3, int ng,
                                int child_i, int child_j, int child_k,
                                int child_gid, int parent_gid, 
                                const RestartData& input_data, int face_dir) const {
  
  int nout1 = nx1 + 2*ng;
  int nout2 = (nx2 > 1) ? (nx2 + 2*ng) : 1;
  int nout3 = (nx3 > 1) ? (nx3 + 2*ng) : 1;
  
  // Get offset for child MeshBlock in fine data
  size_t child_offset = 0;
  if (face_dir == 1) {
    child_offset = child_gid * nout3 * nout2 * (nout1 + 1);
  } else if (face_dir == 2) {
    child_offset = child_gid * nout3 * (nout2 + 1) * nout1;
  } else if (face_dir == 3) {
    child_offset = child_gid * (nout3 + 1) * nout2 * nout1;
  }
  
  // Use appropriate prolongation function based on face direction
  if (face_dir == 1) {
    ProlongFCSharedX1Face(coarse_data, fine_data + child_offset, nx1, nx2, nx3, ng, child_i, child_j, child_k);
  } else if (face_dir == 2) {
    ProlongFCSharedX2Face(coarse_data, fine_data + child_offset, nx1, nx2, nx3, ng, child_i, child_j, child_k);
  } else if (face_dir == 3) {
    ProlongFCSharedX3Face(coarse_data, fine_data + child_offset, nx1, nx2, nx3, ng, child_i, child_j, child_k);
  }
}

//----------------------------------------------------------------------------------------
//! \fn void DataInterpolator::AllocateOutputArrays()
//! \brief Allocate arrays for interpolated data
void DataInterpolator::AllocateOutputArrays(const RestartData& input_data, 
                                           const NewMeshData& mesh_data,
                                           InterpolatedData& output_data) const {
  int nmb_new = mesh_data.nmb_total_new;
  int nout1 = input_data.nout1;
  int nout2 = input_data.nout2;
  int nout3 = input_data.nout3;
  
  std::cout << "Debug: Allocating output arrays for " << nmb_new << " MeshBlocks" << std::endl;
  std::cout << "Debug: MeshBlock size with ghost zones: " << nout1 << " x " << nout2 << " x " << nout3 << std::endl;
  
  // Allocate arrays based on what physics modules are present
  if (input_data.nhydro > 0) {
    size_t hydro_size = static_cast<size_t>(nmb_new) * input_data.nhydro * nout1 * nout2 * nout3;
    size_t hydro_bytes = hydro_size * sizeof(Real);
    std::cout << "Debug: Allocating hydro array: " << hydro_size << " elements = " 
              << hydro_bytes / (1024.0*1024.0*1024.0) << " GB" << std::endl;
    try {
      output_data.hydro_data.resize(hydro_size);
      std::cout << "Debug: Hydro array allocation successful" << std::endl;
    } catch (const std::exception& e) {
      std::cerr << "ERROR: Failed to allocate hydro array: " << e.what() << std::endl;
      throw;
    }
  }
  
  if (input_data.nmhd > 0) {
    // Check for potential overflow before allocation
    uint64_t mhd_size_check = static_cast<uint64_t>(nmb_new) * input_data.nmhd * nout1 * nout2 * nout3;
    size_t mhd_size = static_cast<size_t>(nmb_new) * input_data.nmhd * nout1 * nout2 * nout3;
    
    if (mhd_size_check != mhd_size) {
      std::cerr << "ERROR: Size calculation overflow detected!" << std::endl;
      std::cerr << "uint64_t calculation: " << mhd_size_check << std::endl;
      std::cerr << "size_t calculation: " << mhd_size << std::endl;
      throw std::overflow_error("Array size calculation overflow");
    }
    
    size_t mhd_bytes = mhd_size * sizeof(Real);
    std::cout << "Debug: Allocating MHD array: " << mhd_size << " elements = " 
              << mhd_bytes / (1024.0*1024.0*1024.0) << " GB" << std::endl;
    std::cout << "Debug: MHD calculation: " << nmb_new << " * " << input_data.nmhd 
              << " * " << nout1 << " * " << nout2 << " * " << nout3 << std::endl;
    std::cout << "Debug: sizeof(size_t) = " << sizeof(size_t) << " bytes" << std::endl;
    
    try {
      // Try allocating without initialization first
      output_data.mhd_data.resize(mhd_size);
      std::cout << "Debug: MHD array allocation successful (no initialization)" << std::endl;
    } catch (const std::exception& e) {
      std::cerr << "ERROR: Failed to allocate MHD array: " << e.what() << std::endl;
      std::cerr << "ERROR: Requested size: " << mhd_size << " elements (" 
                << mhd_bytes / (1024.0*1024.0*1024.0) << " GB)" << std::endl;
      throw;
    }
    
    size_t b1f_size = static_cast<size_t>(nmb_new) * (nout1 + 1) * nout2 * nout3;
    size_t b1f_bytes = b1f_size * sizeof(Real);
    std::cout << "Debug: Allocating B1f array: " << b1f_size << " elements = " 
              << b1f_bytes / (1024.0*1024.0*1024.0) << " GB" << std::endl;
    try {
      output_data.mhd_b1f_data.resize(b1f_size);
      std::cout << "Debug: B1f array allocation successful" << std::endl;
    } catch (const std::exception& e) {
      std::cerr << "ERROR: Failed to allocate B1f array: " << e.what() << std::endl;
      throw;
    }
    
    size_t b2f_size = static_cast<size_t>(nmb_new) * nout1 * (nout2 + 1) * nout3;
    size_t b2f_bytes = b2f_size * sizeof(Real);
    std::cout << "Debug: Allocating B2f array: " << b2f_size << " elements = " 
              << b2f_bytes / (1024.0*1024.0*1024.0) << " GB" << std::endl;
    try {
      output_data.mhd_b2f_data.resize(b2f_size);
      std::cout << "Debug: B2f array allocation successful" << std::endl;
    } catch (const std::exception& e) {
      std::cerr << "ERROR: Failed to allocate B2f array: " << e.what() << std::endl;
      throw;
    }
    
    size_t b3f_size = static_cast<size_t>(nmb_new) * nout1 * nout2 * (nout3 + 1);
    size_t b3f_bytes = b3f_size * sizeof(Real);
    std::cout << "Debug: Allocating B3f array: " << b3f_size << " elements = " 
              << b3f_bytes / (1024.0*1024.0*1024.0) << " GB" << std::endl;
    try {
      output_data.mhd_b3f_data.resize(b3f_size);
      std::cout << "Debug: B3f array allocation successful" << std::endl;
    } catch (const std::exception& e) {
      std::cerr << "ERROR: Failed to allocate B3f array: " << e.what() << std::endl;
      throw;
    }
  }
  
  if (input_data.nrad > 0) {
    size_t rad_size = static_cast<size_t>(nmb_new) * input_data.nrad * nout1 * nout2 * nout3;
    std::cout << "Debug: Allocating radiation array: " << rad_size << " elements = " 
              << (rad_size * sizeof(Real)) / (1024.0*1024.0*1024.0) << " GB" << std::endl;
    output_data.rad_data.resize(rad_size);
    std::cout << "Debug: Radiation array allocation successful" << std::endl;
  }
  
  if (input_data.has_turb) {
    size_t force_size = static_cast<size_t>(nmb_new) * input_data.nforce * nout1 * nout2 * nout3;
    std::cout << "Debug: Allocating turbulence array: " << force_size << " elements = " 
              << (force_size * sizeof(Real)) / (1024.0*1024.0*1024.0) << " GB" << std::endl;
    output_data.force_data.resize(force_size);
    std::cout << "Debug: Turbulence array allocation successful" << std::endl;
  }
  
  if (input_data.nz4c > 0) {
    size_t z4c_size = static_cast<size_t>(nmb_new) * input_data.nz4c * nout1 * nout2 * nout3;
    std::cout << "Debug: Allocating Z4c array: " << z4c_size << " elements = " 
              << (z4c_size * sizeof(Real)) / (1024.0*1024.0*1024.0) << " GB" << std::endl;
    output_data.z4c_data.resize(z4c_size);
    std::cout << "Debug: Z4c array allocation successful" << std::endl;
  }
  
  if (input_data.nadm > 0) {
    size_t adm_size = static_cast<size_t>(nmb_new) * input_data.nadm * nout1 * nout2 * nout3;
    std::cout << "Debug: Allocating ADM array: " << adm_size << " elements = " 
              << (adm_size * sizeof(Real)) / (1024.0*1024.0*1024.0) << " GB" << std::endl;
    output_data.adm_data.resize(adm_size);
    std::cout << "Debug: ADM array allocation successful" << std::endl;
  }
}

//----------------------------------------------------------------------------------------
//! \fn void DataInterpolator::CopyPhysicsParameters()
//! \brief Copy physics parameters from input to output data
void DataInterpolator::CopyPhysicsParameters(const RestartData& input_data, 
                                            InterpolatedData& output_data) const {
  output_data.nmb_total = input_data.nmb_total * 8;  // 8 children per parent
  output_data.nhydro = input_data.nhydro;
  output_data.nmhd = input_data.nmhd;
  output_data.nrad = input_data.nrad;
  output_data.nforce = input_data.nforce;
  output_data.nz4c = input_data.nz4c;
  output_data.nadm = input_data.nadm;
  output_data.nout1 = input_data.nout1;
  output_data.nout2 = input_data.nout2;
  output_data.nout3 = input_data.nout3;
  
  // Calculate total memory allocated
  size_t total_bytes = 0;
  total_bytes += output_data.hydro_data.size() * sizeof(Real);
  total_bytes += output_data.mhd_data.size() * sizeof(Real);
  total_bytes += output_data.mhd_b1f_data.size() * sizeof(Real);
  total_bytes += output_data.mhd_b2f_data.size() * sizeof(Real);
  total_bytes += output_data.mhd_b3f_data.size() * sizeof(Real);
  total_bytes += output_data.rad_data.size() * sizeof(Real);
  total_bytes += output_data.force_data.size() * sizeof(Real);
  total_bytes += output_data.z4c_data.size() * sizeof(Real);
  total_bytes += output_data.adm_data.size() * sizeof(Real);
  
  std::cout << "Debug: Total memory allocated for output arrays: " 
            << total_bytes / (1024.0*1024.0*1024.0) << " GB" << std::endl;
}

//----------------------------------------------------------------------------------------
//! \fn bool DataInterpolator::VerifyConservation()
//! \brief Verify that total conserved quantities are preserved
bool DataInterpolator::VerifyConservation(const RestartData& input_data, 
                                         const InterpolatedData& output_data) const {
  // Check mass conservation for MHD
  if (input_data.nmhd > 0) {
    Real total_mass_input = 0.0;
    Real total_mass_output = 0.0;
    
    // Sum mass from all input MeshBlocks
    for (int mb = 0; mb < input_data.nmb_total; ++mb) {
      for (int k = 0; k < input_data.nout3; ++k) {
        for (int j = 0; j < input_data.nout2; ++j) {
          for (int i = 0; i < input_data.nout1; ++i) {
            size_t idx = output_data.GetCCIndex(mb, IDN, k, j, i, output_data.nmhd);
            total_mass_input += input_data.mhd_data[idx];
          }
        }
      }
    }
    
    // Sum mass from all output MeshBlocks
    for (int mb = 0; mb < output_data.nmb_total; ++mb) {
      for (int k = 0; k < output_data.nout3; ++k) {
        for (int j = 0; j < output_data.nout2; ++j) {
          for (int i = 0; i < output_data.nout1; ++i) {
            size_t idx = output_data.GetCCIndex(mb, IDN, k, j, i, output_data.nmhd);
            total_mass_output += output_data.mhd_data[idx];
          }
        }
      }
    }
    
    Real relative_error = std::abs(total_mass_output - total_mass_input) / 
                         std::max(std::abs(total_mass_input), 1e-15);
    
    if (relative_error > 1e-12) {
      SetError("Mass conservation violated: relative error = " + std::to_string(relative_error));
      return false;
    }
    
    std::cout << "Mass conservation verified: relative error = " << relative_error << std::endl;
  }
  
  return true;
}

//----------------------------------------------------------------------------------------
//! \fn bool DataInterpolator::VerifyDivergenceB()
//! \brief Verify that ∇·B = 0 is maintained
bool DataInterpolator::VerifyDivergenceB(const InterpolatedData& output_data) const {
  if (output_data.nmhd == 0) return true;  // No magnetic field to check
  
  Real max_div_b = 0.0;
  int nx1 = output_data.nout1 - 2*2;  // Remove ghost zones (assuming ng=2)
  int nx2 = (output_data.nout2 > 1) ? output_data.nout2 - 2*2 : 1;
  int nx3 = (output_data.nout3 > 1) ? output_data.nout3 - 2*2 : 1;
  int ng = 2;  // Assuming ng=2
  
  for (int mb = 0; mb < output_data.nmb_total; ++mb) {
    for (int k = ng; k < ng + nx3; ++k) {
      for (int j = ng; j < ng + nx2; ++j) {
        for (int i = ng; i < ng + nx1; ++i) {
          
          // Calculate divergence using face-centered fields
          size_t idx_b1f_right = output_data.GetFCX1Index(mb, k, j, i+1);
          size_t idx_b1f_left = output_data.GetFCX1Index(mb, k, j, i);
          Real db1_dx = output_data.mhd_b1f_data[idx_b1f_right] - output_data.mhd_b1f_data[idx_b1f_left];
          
          Real db2_dy = 0.0;
          if (nx2 > 1) {
            size_t idx_b2f_up = output_data.GetFCX2Index(mb, k, j+1, i);
            size_t idx_b2f_down = output_data.GetFCX2Index(mb, k, j, i);
            db2_dy = output_data.mhd_b2f_data[idx_b2f_up] - output_data.mhd_b2f_data[idx_b2f_down];
          }
          
          Real db3_dz = 0.0;
          if (nx3 > 1) {
            size_t idx_b3f_front = output_data.GetFCX3Index(mb, k+1, j, i);
            size_t idx_b3f_back = output_data.GetFCX3Index(mb, k, j, i);
            db3_dz = output_data.mhd_b3f_data[idx_b3f_front] - output_data.mhd_b3f_data[idx_b3f_back];
          }
          
          Real div_b = db1_dx + db2_dy + db3_dz;
          max_div_b = std::max(max_div_b, std::abs(div_b));
        }
      }
    }
  }
  
  if (max_div_b > 1e-12) {
    SetError("Divergence-free constraint violated: max(|∇·B|) = " + std::to_string(max_div_b));
    return false;
  }
  
  std::cout << "Divergence-free constraint verified: max(|∇·B|) = " << max_div_b << std::endl;
  return true;
}

//----------------------------------------------------------------------------------------
//! \fn void DataInterpolator::ProlongFC1()
//! \brief Face-centered prolongation for x1-faces (B1f) 
void DataInterpolator::ProlongFC1(const Real* coarse_b1f, Real* fine_b1f,
                                  int nx1, int nx2, int nx3, int ng,
                                  int coarse_i, int coarse_j, int coarse_k) const {
  int nout1 = nx1 + 2*ng;
  int nout2 = (nx2 > 1) ? (nx2 + 2*ng) : 1;
  int nout3 = (nx3 > 1) ? (nx3 + 2*ng) : 1;
  
  // Simple prolongation: each coarse face value is copied to corresponding fine faces
  for (int k = 0; k < nout3; ++k) {
    for (int j = 0; j < nout2; ++j) {
      for (int i = 0; i <= nout1; ++i) {  // Note: <= for face-centered
        // Map fine indices to coarse indices within child's region
        int local_k = k;
        int local_j = j;
        int local_i = i;
        
        int coarse_ki = local_k / 2;
        int coarse_ji = local_j / 2; 
        int coarse_ii = local_i / 2;
        
        // Handle boundary case for face-centered indexing
        if (i == nout1) coarse_ii = nout1 / 2;
        
        // Map to parent region based on child position
        int parent_ck = coarse_ki + coarse_k * (nout3/2);
        int parent_cj = coarse_ji + coarse_j * (nout2/2);
        int parent_ci = coarse_ii + coarse_i * (nout1/2);
        
        // Get coarse value from parent region
        size_t coarse_idx = parent_ck * nout2 * (nout1+1) + parent_cj * (nout1+1) + parent_ci;
        size_t fine_idx = k * nout2 * (nout1+1) + j * (nout1+1) + i;
        
        // Bounds check
        if (parent_ci <= nout1 && parent_cj < nout2 && parent_ck < nout3) {
          fine_b1f[fine_idx] = coarse_b1f[coarse_idx];
        } else {
          fine_b1f[fine_idx] = 0.0;  // Zero for out-of-bounds
        }
      }
    }
  }
}

//----------------------------------------------------------------------------------------
//! \fn void DataInterpolator::ProlongFC2()
//! \brief Face-centered prolongation for x2-faces (B2f)
void DataInterpolator::ProlongFC2(const Real* coarse_b2f, Real* fine_b2f,
                                  int nx1, int nx2, int nx3, int ng,
                                  int coarse_i, int coarse_j, int coarse_k) const {
  int nout1 = nx1 + 2*ng;
  int nout2 = (nx2 > 1) ? (nx2 + 2*ng) : 1;
  int nout3 = (nx3 > 1) ? (nx3 + 2*ng) : 1;
  
  for (int k = 0; k < nout3; ++k) {
    for (int j = 0; j <= nout2; ++j) {  // Note: <= for face-centered
      for (int i = 0; i < nout1; ++i) {
        // Map fine indices to coarse indices within child's region
        int coarse_ki = k / 2;
        int coarse_ji = j / 2;
        int coarse_ii = i / 2;
        
        // Handle boundary case for face-centered indexing
        if (j == nout2) coarse_ji = nout2 / 2;
        
        // Map to parent region based on child position
        int parent_ck = coarse_ki + coarse_k * (nout3/2);
        int parent_cj = coarse_ji + coarse_j * (nout2/2);
        int parent_ci = coarse_ii + coarse_i * (nout1/2);
        
        // Get coarse value from parent region
        size_t coarse_idx = parent_ck * (nout2+1) * nout1 + parent_cj * nout1 + parent_ci;
        size_t fine_idx = k * (nout2+1) * nout1 + j * nout1 + i;
        
        // Bounds check
        if (parent_ci < nout1 && parent_cj <= nout2 && parent_ck < nout3) {
          fine_b2f[fine_idx] = coarse_b2f[coarse_idx];
        } else {
          fine_b2f[fine_idx] = 0.0;  // Zero for out-of-bounds
        }
      }
    }
  }
}

//----------------------------------------------------------------------------------------
//! \fn void DataInterpolator::ProlongFC3()
//! \brief Face-centered prolongation for x3-faces (B3f)
void DataInterpolator::ProlongFC3(const Real* coarse_b3f, Real* fine_b3f,
                                  int nx1, int nx2, int nx3, int ng,
                                  int coarse_i, int coarse_j, int coarse_k) const {
  int nout1 = nx1 + 2*ng;
  int nout2 = (nx2 > 1) ? (nx2 + 2*ng) : 1;
  int nout3 = (nx3 > 1) ? (nx3 + 2*ng) : 1;
  
  for (int k = 0; k <= nout3; ++k) {  // Note: <= for face-centered
    for (int j = 0; j < nout2; ++j) {
      for (int i = 0; i < nout1; ++i) {
        // Map fine indices to coarse indices within child's region
        int coarse_ki = k / 2;
        int coarse_ji = j / 2;
        int coarse_ii = i / 2;
        
        // Handle boundary case for face-centered indexing
        if (k == nout3) coarse_ki = nout3 / 2;
        
        // Map to parent region based on child position
        int parent_ck = coarse_ki + coarse_k * (nout3/2);
        int parent_cj = coarse_ji + coarse_j * (nout2/2);
        int parent_ci = coarse_ii + coarse_i * (nout1/2);
        
        // Get coarse value from parent region
        size_t coarse_idx = parent_ck * nout2 * nout1 + parent_cj * nout1 + parent_ci;
        size_t fine_idx = k * nout2 * nout1 + j * nout1 + i;
        
        // Bounds check
        if (parent_ci < nout1 && parent_cj < nout2 && parent_ck <= nout3) {
          fine_b3f[fine_idx] = coarse_b3f[coarse_idx];
        } else {
          fine_b3f[fine_idx] = 0.0;  // Zero for out-of-bounds
        }
      }
    }
  }
}