#include "upscaler.hpp"
#include "restart_writer.hpp"
#include "prolongation.hpp"
#include <iostream>
#include <cstring>
#include <limits>

#if MPI_PARALLEL_ENABLED
#include <mpi.h>
#endif

Upscaler::Upscaler(RestartReader& reader) : reader_(reader) {
    CalculateFineMeshConfiguration();
    GenerateFineLogicalLocations();
}

void Upscaler::CalculateFineMeshConfiguration() {
    // Get coarse mesh configuration
    const RegionSize& coarse_mesh_size = reader_.GetMeshSize();
    const RegionIndcs& coarse_mesh_indcs = reader_.GetMeshIndcs();
    const RegionIndcs& coarse_mb_indcs = reader_.GetMBIndcs();
    
    // Double the resolution for fine mesh
    fine_mesh_size_ = coarse_mesh_size;
    fine_mesh_size_.dx1 = coarse_mesh_size.dx1 / 2.0;
    fine_mesh_size_.dx2 = coarse_mesh_size.dx2 / 2.0;
    fine_mesh_size_.dx3 = coarse_mesh_size.dx3 / 2.0;
    
    // Update mesh indices (double the number of cells)
    fine_mesh_indcs_ = coarse_mesh_indcs;
    fine_mesh_indcs_.nx1 = 2 * coarse_mesh_indcs.nx1;
    fine_mesh_indcs_.nx2 = 2 * coarse_mesh_indcs.nx2;
    fine_mesh_indcs_.nx3 = 2 * coarse_mesh_indcs.nx3;
    
    // Active cell indices: is/js/ks start at ng, ie/je/ke are calculated from nx
    // Based on AthenaK mesh.cpp lines 289-310
    fine_mesh_indcs_.is = coarse_mesh_indcs.ng;  // = ng
    fine_mesh_indcs_.ie = fine_mesh_indcs_.is + fine_mesh_indcs_.nx1 - 1;
    fine_mesh_indcs_.js = coarse_mesh_indcs.ng;  // = ng
    fine_mesh_indcs_.je = fine_mesh_indcs_.js + fine_mesh_indcs_.nx2 - 1;
    fine_mesh_indcs_.ks = coarse_mesh_indcs.ng;  // = ng
    fine_mesh_indcs_.ke = fine_mesh_indcs_.ks + fine_mesh_indcs_.nx3 - 1;
    
    // MeshBlock indices remain the same (same number of cells per meshblock)
    fine_mb_indcs_ = coarse_mb_indcs;
    
    // Each coarse meshblock becomes 2/4/8 fine meshblocks based on dimensionality
    int nfine_per_coarse = 2;  // 1D
    if (coarse_mesh_indcs.nx2 > 1) nfine_per_coarse = 4;  // 2D
    if (coarse_mesh_indcs.nx3 > 1) nfine_per_coarse = 8;  // 3D
    fine_nmb_total_ = reader_.GetNMBTotal() * nfine_per_coarse;
    
    if (reader_.GetMyRank() == 0) {
        std::cout << "\nFine mesh configuration:" << std::endl;
        std::cout << "  Total meshblocks: " << reader_.GetNMBTotal() << " -> " << fine_nmb_total_ << std::endl;
        std::cout << "  Mesh resolution: " << coarse_mesh_indcs.nx1 << "³ -> " 
                  << fine_mesh_indcs_.nx1 << "³" << std::endl;
        std::cout << "  Meshblock size remains: " << coarse_mb_indcs.nx1 << "³" << std::endl;
    }
}

void Upscaler::GenerateFineLogicalLocations() {
    const std::vector<LogicalLocation>& coarse_llocs = reader_.GetLlocEachMB();
    const std::vector<float>& coarse_costs = reader_.GetCostEachMB();
    
    fine_lloc_eachmb_.resize(fine_nmb_total_);
    fine_cost_eachmb_.resize(fine_nmb_total_);
    
    // Determine number of fine blocks per coarse block (for 3D it's 8)
    const RegionIndcs& mesh_indcs = reader_.GetMeshIndcs();
    int nfine_per_coarse = 8;  // For 3D case
    
    // For each coarse meshblock, create 8 fine meshblocks (3D case)
    for (int cmb = 0; cmb < reader_.GetNMBTotal(); cmb++) {
        const LogicalLocation& cloc = coarse_llocs[cmb];
        float ccost = coarse_costs[cmb] / static_cast<float>(nfine_per_coarse);  // Distribute cost evenly
        
        // Create 8 fine blocks (2x2x2) for each coarse block
        for (int k = 0; k < 2; k++) {
            for (int j = 0; j < 2; j++) {
                for (int i = 0; i < 2; i++) {
                    int fmb = cmb * 8 + k * 4 + j * 2 + i;
                    
                    // Fine logical location has level+1 and doubled indices
                    fine_lloc_eachmb_[fmb].level = cloc.level + 1;
                    fine_lloc_eachmb_[fmb].lx1 = 2 * cloc.lx1 + i;
                    fine_lloc_eachmb_[fmb].lx2 = 2 * cloc.lx2 + j;
                    fine_lloc_eachmb_[fmb].lx3 = 2 * cloc.lx3 + k;
                    
                    fine_cost_eachmb_[fmb] = ccost;
                }
            }
        }
    }
}

bool Upscaler::UpscaleMeshBlock(int coarse_mb_id, int fine_mb_start_id,
                               const Real* coarse_data, Real* fine_data,
                               int nvars, bool is_face_centered, int face_dir) {
    const RegionIndcs& mb_indcs = reader_.GetMBIndcs();
    int cnx1 = mb_indcs.nx1;
    int cnx2 = mb_indcs.nx2;
    int cnx3 = mb_indcs.nx3;
    int cng = mb_indcs.ng;
    
    // Total cells including ghost zones
    int cnx1_tot = cnx1 + 2*cng;
    int cnx2_tot = (cnx2 > 1) ? cnx2 + 2*cng : 1;
    int cnx3_tot = (cnx3 > 1) ? cnx3 + 2*cng : 1;
    
    bool multi_d = (cnx2 > 1);
    bool three_d = (cnx3 > 1);
    
    // Process each variable
    for (int v = 0; v < nvars; v++) {
        // Get pointer to coarse data for this variable
        const Real* coarse_var = coarse_data + v * cnx3_tot * cnx2_tot * cnx1_tot;
        
        // For each of the 8 fine meshblocks
        for (int fmb = 0; fmb < 8; fmb++) {
            int fk_offset = (fmb / 4) * cnx3;
            int fj_offset = ((fmb / 2) % 2) * cnx2;
            int fi_offset = (fmb % 2) * cnx1;
            
            // Get pointer to fine data for this meshblock and variable
            // Use size_t for all calculations to avoid integer overflow
            size_t cells_per_mb_size = static_cast<size_t>(cnx3_tot) * cnx2_tot * cnx1_tot;
            size_t mb_offset = static_cast<size_t>(fine_mb_start_id + fmb) * nvars * cells_per_mb_size;
            size_t var_offset = static_cast<size_t>(v) * cells_per_mb_size;
            size_t fine_offset = mb_offset + var_offset;
            Real* fine_var = fine_data + fine_offset;
            
            // DEBUG: Print pointer calculations
            if (v == 0 && fmb == 0) {
                std::cout << "\nDEBUG: Meshblock pointer calculations:" << std::endl;
                std::cout << "  fine_mb_start_id=" << fine_mb_start_id << " fmb=" << fmb << std::endl;
                std::cout << "  nvars=" << nvars << " cnx1_tot=" << cnx1_tot << " cnx2_tot=" << cnx2_tot << " cnx3_tot=" << cnx3_tot << std::endl;
                std::cout << "  cells_per_mb_size=" << cells_per_mb_size << std::endl;
                std::cout << "  mb_offset=" << mb_offset << " var_offset=" << var_offset << std::endl;
                std::cout << "  fine_offset=" << fine_offset << std::endl;
                std::cout << "  fine_data base=" << (void*)fine_data << " fine_var=" << (void*)fine_var << std::endl;
                std::cout << "  Total fine data size=" << (static_cast<size_t>(fine_nmb_total_) * nvars * cells_per_mb_size) << std::endl;
            }
            
            // Prolongate cells from coarse to fine meshblock
            // The approach is: for each fine meshblock, determine which coarse cells
            // contribute to it and where they map in the fine grid
            
            // Calculate which portion of the coarse meshblock maps to this fine meshblock
            // Each fine meshblock gets half of the coarse cells in each direction
            int ci_start = cng + (fi_offset > 0 ? cnx1/2 : 0);
            int ci_end = cng + (fi_offset > 0 ? cnx1 : cnx1/2);
            int cj_start = cng + (multi_d && fj_offset > 0 ? cnx2/2 : 0);
            int cj_end = cng + (multi_d ? (fj_offset > 0 ? cnx2 : cnx2/2) : 1);
            int ck_start = cng + (three_d && fk_offset > 0 ? cnx3/2 : 0);
            int ck_end = cng + (three_d ? (fk_offset > 0 ? cnx3 : cnx3/2) : 1);
            
            // DEBUG: Print loop bounds for first meshblock
            if (v == 0 && fmb == 0) {
                std::cout << "\nDEBUG: Loop bounds for fmb=" << fmb << ":" << std::endl;
                std::cout << "  ci: " << ci_start << " to " << ci_end << std::endl;
                std::cout << "  cj: " << cj_start << " to " << cj_end << std::endl;
                std::cout << "  ck: " << ck_start << " to " << ck_end << std::endl;
            }
            
            // Loop over coarse cells that map to this fine meshblock
            for (int k = ck_start; k < ck_end; k++) {
                for (int j = cj_start; j < cj_end; j++) {
                    for (int i = ci_start; i < ci_end; i++) {
                        // Calculate where this coarse cell maps in the fine meshblock
                        // The mapping is: fine_index = 2 * (coarse_index - coarse_start) + fine_ng
                        int fi = cng + 2 * (i - ci_start);
                        int fj = cng + (multi_d ? 2 * (j - cj_start) : 0);
                        int fk = cng + (three_d ? 2 * (k - ck_start) : 0);
                        
                        // Apply prolongation operator
                        ProlongCCStandalone(coarse_var, fine_var,
                                          k, j, i, fk, fj, fi,
                                          cnx1_tot, cnx2_tot, cnx3_tot,
                                          cnx1_tot, cnx2_tot, cnx3_tot,
                                          multi_d, three_d);
                    }
                }
            }
            
            // Fill ghost zones by copying from neighboring cells (simple approach)
            // NOTE: This is a simplified placeholder. Proper ghost zone filling would require:
            // 1. Information about neighboring meshblocks
            // 2. Periodic/reflecting/outflow boundary conditions at domain edges
            // 3. Proper interpolation for AMR boundaries
            // For restart files, ghost zones may not be critical if the simulation
            // will immediately update them based on boundary conditions
            for (int k = 0; k < cnx3_tot; k++) {
                for (int j = 0; j < cnx2_tot; j++) {
                    for (int i = 0; i < cnx1_tot; i++) {
                        if (i < cng || i >= cnx1 + cng ||
                            j < cng || j >= cnx2 + cng ||
                            k < cng || k >= cnx3 + cng) {
                            // Simple copy from nearest active cell
                            int ci = std::max(cng, std::min(i, cnx1 + cng - 1));
                            int cj = std::max(cng, std::min(j, cnx2 + cng - 1));
                            int ck = std::max(cng, std::min(k, cnx3 + cng - 1));
                            
                            int src_idx = ck * cnx2_tot * cnx1_tot + cj * cnx1_tot + ci;
                            int dst_idx = k * cnx2_tot * cnx1_tot + j * cnx1_tot + i;
                            fine_var[dst_idx] = fine_var[src_idx];
                        }
                    }
                }
            }
        }
    }
    
    return true;
}

bool Upscaler::UpscaleFaceCenteredData(int coarse_mb_id, int fine_mb_start_id,
                                      const Real* coarse_x1f, const Real* coarse_x2f, const Real* coarse_x3f,
                                      Real* fine_x1f, Real* fine_x2f, Real* fine_x3f) {
    const RegionIndcs& mb_indcs = reader_.GetMBIndcs();
    int cnx1 = mb_indcs.nx1;
    int cnx2 = mb_indcs.nx2;
    int cnx3 = mb_indcs.nx3;
    int cng = mb_indcs.ng;
    
    // Total cells including ghost zones
    int cnx1_tot = cnx1 + 2*cng;
    int cnx2_tot = (cnx2 > 1) ? cnx2 + 2*cng : 1;
    int cnx3_tot = (cnx3 > 1) ? cnx3 + 2*cng : 1;
    
    bool multi_d = (cnx2 > 1);
    bool three_d = (cnx3 > 1);
    
    // Face array dimensions - CRITICAL: These must match AthenaK's face array layout
    int x1f_size = cnx3_tot * cnx2_tot * (cnx1_tot + 1);
    int x2f_size = cnx3_tot * (cnx2_tot + (multi_d ? 1 : 0)) * cnx1_tot;
    int x3f_size = (cnx3_tot + (three_d ? 1 : 0)) * cnx2_tot * cnx1_tot;
    
    // Inner cell bounds (without ghost zones)
    int cis = cng, cie = cng + cnx1 - 1;
    int cjs = cng, cje = cng + cnx2 - 1;
    int cks = cng, cke = cng + cnx3 - 1;
    
    if (reader_.GetMyRank() == 0 && coarse_mb_id == 0) {
        std::cout << "DEBUG UPSCALE: cnx1=" << cnx1 << " cnx2=" << cnx2 << " cnx3=" << cnx3 << " cng=" << cng << std::endl;
        std::cout << "DEBUG UPSCALE: cnx1_tot=" << cnx1_tot << " cnx2_tot=" << cnx2_tot << " cnx3_tot=" << cnx3_tot << std::endl;
        std::cout << "DEBUG UPSCALE: x1f_size=" << x1f_size << " x2f_size=" << x2f_size << " x3f_size=" << x3f_size << std::endl;
        std::cout << "DEBUG UPSCALE: cis=" << cis << " cie=" << cie << " cjs=" << cjs << " cje=" << cje << " cks=" << cks << " cke=" << cke << std::endl;
    }
    
    // Process each of the 8 fine meshblocks
    for (int fmb = 0; fmb < 8; fmb++) {
        // Determine the position of this fine meshblock within the coarse meshblock
        int fmb_k = (fmb / 4);      // 0 or 1 in k direction
        int fmb_j = ((fmb / 2) % 2); // 0 or 1 in j direction  
        int fmb_i = (fmb % 2);       // 0 or 1 in i direction
        
        // Get pointers to fine face data for this meshblock
        Real* fine_mb_x1f = fine_x1f + (fine_mb_start_id + fmb) * x1f_size;
        Real* fine_mb_x2f = fine_x2f + (fine_mb_start_id + fmb) * x2f_size;
        Real* fine_mb_x3f = fine_x3f + (fine_mb_start_id + fmb) * x3f_size;
        
        // STEP 1: Prolongate shared faces between coarse and fine cells
        // Following AthenaK's logic from mesh_refinement.cpp
        
        // Prolongate X1 faces (faces perpendicular to x1 direction)
        // Loop over coarse X1 faces that map to this fine meshblock
        for (int k = cks; k <= cke; k++) {
            for (int j = cjs; j <= cje; j++) {
                // Determine which coarse X1 faces to process for this fine meshblock
                int ci_start = cis + (fmb_i == 0 ? 0 : cnx1/2);
                int ci_end = cis + (fmb_i == 0 ? cnx1/2 : cnx1);
                
                for (int i = ci_start; i <= ci_end; i++) { // Note: <= because X1 faces have nx1+1 faces
                    // Calculate fine indices 
                    int fi = cng + 2 * (i - ci_start);
                    int fj = cng + (multi_d ? 2 * ((j - cjs) % (cnx2/2)) : 0);
                    int fk = cng + (three_d ? 2 * ((k - cks) % (cnx3/2)) : 0);
                    
                    // Skip if we're beyond the active region for this fine meshblock
                    if ((multi_d && ((fmb_j == 0 && j >= cjs + cnx2/2) || (fmb_j == 1 && j < cjs + cnx2/2))) ||
                        (three_d && ((fmb_k == 0 && k >= cks + cnx3/2) || (fmb_k == 1 && k < cks + cnx3/2)))) {
                        continue;
                    }
                    
                    // Apply prolongation for shared X1 face
                    ProlongFCSharedX1FaceStandalone(coarse_x1f, fine_mb_x1f,
                                                    k, j, i, fk, fj, fi,
                                                    cnx1_tot, cnx2_tot, cnx3_tot,
                                                    cnx1_tot, cnx2_tot, cnx3_tot,
                                                    multi_d, three_d);
                }
            }
        }
        
        // Prolongate X2 faces (faces perpendicular to x2 direction)
        if (multi_d) {
            for (int k = cks; k <= cke; k++) {
                // Determine which coarse X2 faces to process for this fine meshblock
                int cj_start = cjs + (fmb_j == 0 ? 0 : cnx2/2);
                int cj_end = cjs + (fmb_j == 0 ? cnx2/2 : cnx2);
                
                for (int j = cj_start; j <= cj_end; j++) { // Note: <= because X2 faces have nx2+1 faces
                    for (int i = cis; i <= cie; i++) {
                        // Calculate fine indices
                        int fi = cng + 2 * ((i - cis) % (cnx1/2));
                        int fj = cng + 2 * (j - cj_start);
                        int fk = cng + (three_d ? 2 * ((k - cks) % (cnx3/2)) : 0);
                        
                        // Skip if we're beyond the active region for this fine meshblock
                        if ((fmb_i == 0 && i >= cis + cnx1/2) || (fmb_i == 1 && i < cis + cnx1/2) ||
                            (three_d && ((fmb_k == 0 && k >= cks + cnx3/2) || (fmb_k == 1 && k < cks + cnx3/2)))) {
                            continue;
                        }
                        
                        // Apply prolongation for shared X2 face
                        ProlongFCSharedX2FaceStandalone(coarse_x2f, fine_mb_x2f,
                                                        k, j, i, fk, fj, fi,
                                                        cnx1_tot, cnx2_tot, cnx3_tot,
                                                        cnx1_tot, cnx2_tot, cnx3_tot,
                                                        three_d);
                    }
                }
            }
        }
        
        // Prolongate X3 faces (faces perpendicular to x3 direction)
        if (three_d) {
            // Determine which coarse X3 faces to process for this fine meshblock
            int ck_start = cks + (fmb_k == 0 ? 0 : cnx3/2);
            int ck_end = cks + (fmb_k == 0 ? cnx3/2 : cnx3);
            
            for (int k = ck_start; k <= ck_end; k++) { // Note: <= because X3 faces have nx3+1 faces
                for (int j = cjs; j <= cje; j++) {
                    for (int i = cis; i <= cie; i++) {
                        // Calculate fine indices
                        int fi = cng + 2 * ((i - cis) % (cnx1/2));
                        int fj = cng + (multi_d ? 2 * ((j - cjs) % (cnx2/2)) : 0);
                        int fk = cng + 2 * (k - ck_start);
                        
                        // Skip if we're beyond the active region for this fine meshblock
                        if ((fmb_i == 0 && i >= cis + cnx1/2) || (fmb_i == 1 && i < cis + cnx1/2) ||
                            (multi_d && ((fmb_j == 0 && j >= cjs + cnx2/2) || (fmb_j == 1 && j < cjs + cnx2/2)))) {
                            continue;
                        }
                        
                        // Apply prolongation for shared X3 face
                        ProlongFCSharedX3FaceStandalone(coarse_x3f, fine_mb_x3f,
                                                        k, j, i, fk, fj, fi,
                                                        cnx1_tot, cnx2_tot, cnx3_tot,
                                                        cnx1_tot, cnx2_tot, cnx3_tot,
                                                        multi_d);
                    }
                }
            }
        }
        
        // STEP 2: Prolongate internal faces using divergence-preserving interpolation
        // This fills the faces that are internal to the fine cells
        // Loop over the cells in this fine meshblock and fill internal faces
        for (int k = cks; k <= cke; k++) {
            for (int j = cjs; j <= cje; j++) {
                for (int i = cis; i <= cie; i++) {
                    // Check if this coarse cell maps to this fine meshblock
                    bool in_i_range = (fmb_i == 0 && i < cis + cnx1/2) || (fmb_i == 1 && i >= cis + cnx1/2);
                    bool in_j_range = !multi_d || (fmb_j == 0 && j < cjs + cnx2/2) || (fmb_j == 1 && j >= cjs + cnx2/2);
                    bool in_k_range = !three_d || (fmb_k == 0 && k < cks + cnx3/2) || (fmb_k == 1 && k >= cks + cnx3/2);
                    
                    if (in_i_range && in_j_range && in_k_range) {
                        // Calculate fine cell indices for the lower-left corner of the 2x2x2 fine cells
                        int fi = cng + 2 * ((i - cis) % (cnx1/2));
                        int fj = cng + (multi_d ? 2 * ((j - cjs) % (cnx2/2)) : 0);
                        int fk = cng + (three_d ? 2 * ((k - cks) % (cnx3/2)) : 0);
                        
                        // Apply divergence-preserving prolongation for internal faces
                        bool one_d = !multi_d && !three_d;
                        if (!one_d) {  // In 1D, internal faces are trivial
                            ProlongFCInternalStandalone(fine_mb_x1f, fine_mb_x2f, fine_mb_x3f,
                                                       fk, fj, fi,
                                                       cnx1_tot, cnx2_tot, cnx3_tot,
                                                       three_d);
                        }
                    }
                }
            }
        }
    }
    
    return true;
}

bool Upscaler::UpscaleRestartFile(const std::string& output_filename) {
    const PhysicsConfig& phys_config = reader_.GetPhysicsConfig();
    const RegionIndcs& mb_indcs = reader_.GetMBIndcs();
    PhysicsReader* phys_reader = reader_.GetPhysicsReader();
    
    if (!phys_reader) {
        std::cerr << "ERROR: No physics reader available" << std::endl;
        return false;
    }
    
    int cnx1 = mb_indcs.nx1;
    int cnx2 = mb_indcs.nx2;
    int cnx3 = mb_indcs.nx3;
    int cng = mb_indcs.ng;
    
    // Total cells including ghost zones
    int cnx1_tot = cnx1 + 2*cng;
    int cnx2_tot = (cnx2 > 1) ? cnx2 + 2*cng : 1;
    int cnx3_tot = (cnx3 > 1) ? cnx3 + 2*cng : 1;
    int cells_per_mb = cnx3_tot * cnx2_tot * cnx1_tot;
    
    // Add overflow checks for memory allocation
    if (reader_.GetMyRank() == 0) {
        std::cout << "\nMemory allocation check:" << std::endl;
        std::cout << "  fine_nmb_total_ = " << fine_nmb_total_ << std::endl;
        std::cout << "  cells_per_mb = " << cells_per_mb << std::endl;
        
        // Check for potential overflow in size calculations
        size_t max_size_t = std::numeric_limits<size_t>::max();
        
        // Check hydro allocation
        if (phys_config.has_hydro) {
            size_t hydro_elements = fine_nmb_total_ * phys_config.nhydro * cells_per_mb;
            if (fine_nmb_total_ > max_size_t / (phys_config.nhydro * cells_per_mb)) {
                std::cerr << "ERROR: Integer overflow detected in hydro data allocation!" << std::endl;
                std::cerr << "  Required elements would exceed size_t maximum" << std::endl;
                return false;
            }
            std::cout << "  Hydro data: " << hydro_elements << " elements = " 
                      << (hydro_elements * sizeof(Real) / (1024.0*1024.0*1024.0)) << " GB" << std::endl;
        }
        
        // Check MHD allocation
        if (phys_config.has_mhd) {
            size_t mhd_elements = fine_nmb_total_ * phys_config.nmhd * cells_per_mb;
            if (fine_nmb_total_ > max_size_t / (phys_config.nmhd * cells_per_mb)) {
                std::cerr << "ERROR: Integer overflow detected in MHD data allocation!" << std::endl;
                return false;
            }
            
            // Check face-centered arrays
            size_t x1f_check = fine_nmb_total_ * cnx3_tot * cnx2_tot * (cnx1_tot + 1);
            size_t x2f_check = fine_nmb_total_ * cnx3_tot * (cnx2_tot + 1) * cnx1_tot;
            size_t x3f_check = fine_nmb_total_ * (cnx3_tot + 1) * cnx2_tot * cnx1_tot;
            
            std::cout << "  MHD cell-centered: " << mhd_elements << " elements = " 
                      << (mhd_elements * sizeof(Real) / (1024.0*1024.0*1024.0)) << " GB" << std::endl;
            std::cout << "  MHD face x1f: " << x1f_check << " elements = " 
                      << (x1f_check * sizeof(Real) / (1024.0*1024.0*1024.0)) << " GB" << std::endl;
            std::cout << "  MHD face x2f: " << x2f_check << " elements = " 
                      << (x2f_check * sizeof(Real) / (1024.0*1024.0*1024.0)) << " GB" << std::endl;
            std::cout << "  MHD face x3f: " << x3f_check << " elements = " 
                      << (x3f_check * sizeof(Real) / (1024.0*1024.0*1024.0)) << " GB" << std::endl;
        }
    }
    
    // Allocate storage for fine data
    if (phys_config.has_hydro) {
        fine_hydro_data_.resize(fine_nmb_total_ * phys_config.nhydro * cells_per_mb);
    }
    if (phys_config.has_mhd) {
        fine_mhd_data_.resize(fine_nmb_total_ * phys_config.nmhd * cells_per_mb);
        
        // Initialize face-centered arrays and explicitly set to zero
        size_t x1f_total_size = fine_nmb_total_ * cnx3_tot * cnx2_tot * (cnx1_tot + 1);
        size_t x2f_total_size = fine_nmb_total_ * cnx3_tot * (cnx2_tot + 1) * cnx1_tot;
        size_t x3f_total_size = fine_nmb_total_ * (cnx3_tot + 1) * cnx2_tot * cnx1_tot;
        
        fine_x1f_data_.assign(x1f_total_size, 0.0);
        fine_x2f_data_.assign(x2f_total_size, 0.0);
        fine_x3f_data_.assign(x3f_total_size, 0.0);
        
        if (reader_.GetMyRank() == 0) {
            std::cout << "DEBUG: Allocated face arrays - x1f_size=" << x1f_total_size 
                      << " x2f_size=" << x2f_total_size << " x3f_size=" << x3f_total_size << std::endl;
        }
    }
    if (phys_config.has_turbulence) {
        fine_turb_data_.resize(fine_nmb_total_ * phys_config.nforce * cells_per_mb);
    }
    // Add other physics as needed...
    
    if (reader_.GetMyRank() == 0) {
        std::cout << "\nStarting upscaling process..." << std::endl;
    }
    
    // Process each coarse meshblock
    int nmb_local = reader_.GetNMBThisRank();
    int nfine_per_coarse = 8;  // For 3D case
    for (int local_mb = 0; local_mb < nmb_local; local_mb++) {
        // For now, assume global_mb = local_mb (no MPI distribution)
        int global_mb = local_mb;
        int fine_mb_start = global_mb * nfine_per_coarse;
        
        if (reader_.GetMyRank() == 0 && local_mb % 10 == 0) {
            std::cout << "  Processing meshblock " << local_mb << "/" << nmb_local << std::endl;
        }
        
        // Upscale hydro data
        if (phys_config.has_hydro && phys_reader->GetHydroData().size() > local_mb) {
            const Real* coarse_hydro = phys_reader->GetHydroData()[local_mb].data();
            Real* fine_hydro = fine_hydro_data_.data();
            UpscaleMeshBlock(global_mb, fine_mb_start, coarse_hydro, fine_hydro, phys_config.nhydro);
        }
        
        // Upscale MHD data
        if (phys_config.has_mhd && phys_reader->GetMHDData().size() > local_mb) {
            const Real* coarse_mhd = phys_reader->GetMHDData()[local_mb].data();
            Real* fine_mhd = fine_mhd_data_.data();
            
            // DEBUG: Print coarse MHD data before upscaling
            if (reader_.GetMyRank() == 0 && local_mb == 0) {
                std::cout << "\n=== DEBUG: COARSE MHD DATA (MeshBlock " << global_mb << ") ===" << std::endl;
                const RegionIndcs& mb_indcs = reader_.GetMBIndcs();
                int cells_per_mb = (mb_indcs.nx1 + 2*mb_indcs.ng) * (mb_indcs.nx2 + 2*mb_indcs.ng) * (mb_indcs.nx3 + 2*mb_indcs.ng);
                
                std::cout << "Coarse MHD data (first 5 cells):" << std::endl;
                for (int cell = 0; cell < std::min(5, cells_per_mb); cell++) {
                    std::cout << "  Cell " << cell << ": ";
                    for (int var = 0; var < phys_config.nmhd; var++) {
                        Real value = coarse_mhd[var * cells_per_mb + cell];
                        std::cout << "var[" << var << "]=" << std::scientific << value << " ";
                    }
                    std::cout << std::endl;
                }
            }
            
            UpscaleMeshBlock(global_mb, fine_mb_start, coarse_mhd, fine_mhd, phys_config.nmhd);
            
            // DEBUG: Print fine MHD data after upscaling (only for first meshblock)
            if (reader_.GetMyRank() == 0 && local_mb == 0) {
                std::cout << "Fine MHD data after upscaling (first 5 cells of fine meshblock 0):" << std::endl;
                const RegionIndcs& mb_indcs = reader_.GetMBIndcs();
                int fine_cells_per_mb = (mb_indcs.nx1 + 2*mb_indcs.ng) * (mb_indcs.nx2 + 2*mb_indcs.ng) * (mb_indcs.nx3 + 2*mb_indcs.ng);
                
                for (int cell = 0; cell < std::min(5, fine_cells_per_mb); cell++) {
                    std::cout << "  Cell " << cell << ": ";
                    for (int var = 0; var < phys_config.nmhd; var++) {
                        Real value = fine_mhd[var * fine_cells_per_mb + cell];
                        std::cout << "var[" << var << "]=" << std::scientific << value << " ";
                    }
                    std::cout << std::endl;
                }
            }
            
            // Upscale face-centered magnetic fields
            const Real* coarse_x1f = phys_reader->GetMHDB1f()[local_mb].data();
            const Real* coarse_x2f = phys_reader->GetMHDB2f()[local_mb].data();
            const Real* coarse_x3f = phys_reader->GetMHDB3f()[local_mb].data();
            Real* fine_x1f = fine_x1f_data_.data();
            Real* fine_x2f = fine_x2f_data_.data();
            Real* fine_x3f = fine_x3f_data_.data();
            
            // DEBUG: Print coarse face-centered data
            if (reader_.GetMyRank() == 0) {
                std::cout << "Coarse face-centered B-field data (first 5 faces):" << std::endl;
                std::cout << "  Bx (x1f): ";
                for (int i = 0; i < 5; i++) {
                    std::cout << coarse_x1f[i] << " ";
                }
                std::cout << std::endl;
                std::cout << "  By (x2f): ";
                for (int i = 0; i < 5; i++) {
                    std::cout << coarse_x2f[i] << " ";
                }
                std::cout << std::endl;
                std::cout << "  Bz (x3f): ";
                for (int i = 0; i < 5; i++) {
                    std::cout << coarse_x3f[i] << " ";
                }
                std::cout << std::endl;
            }
            
            UpscaleFaceCenteredData(global_mb, fine_mb_start,
                                    coarse_x1f, coarse_x2f, coarse_x3f,
                                    fine_x1f, fine_x2f, fine_x3f);
                                   
            // DEBUG: Print fine face-centered data after upscaling
            if (reader_.GetMyRank() == 0 && local_mb == 0) {
                std::cout << "Fine face-centered B-field data after upscaling:" << std::endl;
                
                // Print data for each of the 8 fine meshblocks
                for (int fmb = 0; fmb < 8; fmb++) {
                    std::cout << "  Fine meshblock " << fmb << " (first 5 faces):" << std::endl;
                    
                    // Calculate offsets for this fine meshblock
                    // Calculate face sizes based on mesh dimensions
                    int x1f_size_local = cnx3_tot * cnx2_tot * (cnx1_tot + 1);
                    int x2f_size_local = cnx3_tot * (cnx2_tot + ((cnx2 > 1) ? 1 : 0)) * cnx1_tot;
                    int x3f_size_local = (cnx3_tot + ((cnx3 > 1) ? 1 : 0)) * cnx2_tot * cnx1_tot;
                    
                    int x1f_offset = (fine_mb_start + fmb) * x1f_size_local;
                    int x2f_offset = (fine_mb_start + fmb) * x2f_size_local;
                    int x3f_offset = (fine_mb_start + fmb) * x3f_size_local;
                    
                    std::cout << "    Bx (x1f): ";
                    for (int i = 0; i < 5; i++) {
                        std::cout << fine_x1f[x1f_offset + i] << " ";
                    }
                    std::cout << std::endl;
                    
                    std::cout << "    By (x2f): ";
                    for (int i = 0; i < 5; i++) {
                        std::cout << fine_x2f[x2f_offset + i] << " ";
                    }
                    std::cout << std::endl;
                    
                    std::cout << "    Bz (x3f): ";
                    for (int i = 0; i < 5; i++) {
                        std::cout << fine_x3f[x3f_offset + i] << " ";
                    }
                    std::cout << std::endl;
                }
                std::cout << "=== END UPSCALER DEBUG ===" << std::endl << std::endl;
            }
        }
        
        // Upscale turbulence force data
        if (phys_config.has_turbulence && phys_reader->GetTurbData().size() > local_mb) {
            const Real* coarse_turb = phys_reader->GetTurbData()[local_mb].data();
            Real* fine_turb = fine_turb_data_.data();
            UpscaleMeshBlock(global_mb, fine_mb_start, coarse_turb, fine_turb, phys_config.nforce);
        }
    }
    
    // Now write the upscaled data to output file
    if (reader_.GetMyRank() == 0) {
        std::cout << "\nWriting upscaled restart file: " << output_filename << std::endl;
    }
    
    // Create restart writer and write the file
    RestartWriter writer(reader_, *this);
    bool success = writer.WriteRestartFile(output_filename);
    
    if (success && reader_.GetMyRank() == 0) {
        std::cout << "Successfully wrote upscaled restart file!" << std::endl;
    }
    
    return success;
}