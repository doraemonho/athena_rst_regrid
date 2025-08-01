#include "physics_reader.hpp"
#include <iostream>
#include <iomanip>
#include <algorithm>

#if MPI_PARALLEL_ENABLED
#include <mpi.h>
#endif

PhysicsReader::PhysicsReader(IOWrapper& file, const PhysicsConfig& config, 
                            const MPIDistribution& mpi_dist, int my_rank,
                            int nout1, int nout2, int nout3)
    : file_(file), config_(config), mpi_dist_(mpi_dist), my_rank_(my_rank),
      nout1_(nout1), nout2_(nout2), nout3_(nout3) {
    
    nmb_thisrank_ = mpi_dist_.GetNumMeshBlocks(my_rank_);
    
    // Initialize data storage (only for MeshBlocks owned by this rank)
    if (config_.has_hydro) {
        hydro_data_.resize(nmb_thisrank_);
        for (int mb = 0; mb < nmb_thisrank_; mb++) {
            hydro_data_[mb].resize(config_.nhydro * nout1_ * nout2_ * nout3_);
        }
    }
    
    if (config_.has_mhd) {
        mhd_data_.resize(nmb_thisrank_);
        mhd_b1f_.resize(nmb_thisrank_);
        mhd_b2f_.resize(nmb_thisrank_);
        mhd_b3f_.resize(nmb_thisrank_);
        
        for (int mb = 0; mb < nmb_thisrank_; mb++) {
            mhd_data_[mb].resize(config_.nmhd * nout1_ * nout2_ * nout3_);
            mhd_b1f_[mb].resize((nout1_ + 1) * nout2_ * nout3_);
            mhd_b2f_[mb].resize(nout1_ * (nout2_ + 1) * nout3_);
            mhd_b3f_[mb].resize(nout1_ * nout2_ * (nout3_ + 1));
        }
    }
    
    if (config_.has_turbulence) {
        turb_data_.resize(nmb_thisrank_);
        for (int mb = 0; mb < nmb_thisrank_; mb++) {
            turb_data_[mb].resize(config_.nforce * nout1_ * nout2_ * nout3_);
        }
    }
}

bool PhysicsReader::ReadPhysicsData(IOWrapperSizeT physics_start, IOWrapperSizeT data_size) {
    std::cout << "Reading physics data..." << std::endl;
    // Validate and potentially adjust data size
    IOWrapperSizeT calculated_size = CalculateDataSize();
    
    std::cout << "Data size per MeshBlock: " << data_size << " bytes" << std::endl;
    
    // Get MPI distribution info
    int noutmbs_min, noutmbs_max;
    mpi_dist_.GetRankMinMax(&noutmbs_min, &noutmbs_max);
    
    // Following AthenaK pgen.cpp:289-291 exactly
    int mygids = mpi_dist_.GetStartingMeshBlockID(my_rank_);
    IOWrapperSizeT offset_myrank = physics_start + data_size * mygids;
    IOWrapperSizeT myoffset = offset_myrank;
    
    // Debug output
    std::cout << "Rank " << my_rank_ << ": mygids=" << mygids 
              << ", offset_myrank=" << offset_myrank 
              << ", nmb_thisrank=" << nmb_thisrank_ << std::endl;
    
    std::cout << "Physics data reading with PHYSICS-MODULE-INTERLEAVED layout:" << std::endl;
    std::cout << "  My rank: " << my_rank_ << std::endl;
    std::cout << "  Starting MeshBlock ID: " << mygids << std::endl;
    std::cout << "  MeshBlocks for this rank: " << nmb_thisrank_ << std::endl;
    
    // Read physics modules in AthenaK's exact order
    if (config_.has_hydro) {
        if (!ReadHydroData(myoffset, offset_myrank, data_size, noutmbs_min, noutmbs_max)) {
            return false;
        }
    }
    if (config_.has_mhd) {
        // print the offset
        std::cout << "MHD offset: " << myoffset << std::endl;
        std::cout << "MHD offset_myrank: " << offset_myrank << std::endl;
        if (!ReadMHDData(myoffset, offset_myrank, data_size, noutmbs_min, noutmbs_max)) {
            return false;
        }
    }
    
    if (config_.has_radiation) {
        SkipRadiationData(offset_myrank);
    }
    
    if (config_.has_turbulence) {
        if (!ReadTurbulenceData(myoffset, offset_myrank, data_size, noutmbs_min, noutmbs_max)) {
            return false;
        }
    }
    
    if (config_.has_z4c || config_.has_adm) {
        SkipZ4cADMData(offset_myrank);
    }
    
    std::cout << "Physics data read successfully" << std::endl;
    return true;
}

IOWrapperSizeT PhysicsReader::CalculateDataSize() const {
    IOWrapperSizeT data_size = 0;
    
    if (config_.has_hydro) {
        data_size += nout1_ * nout2_ * nout3_ * config_.nhydro * sizeof(Real);
    }
    if (config_.has_mhd) {
        data_size += nout1_ * nout2_ * nout3_ * config_.nmhd * sizeof(Real);     // MHD CC
        data_size += (nout1_ + 1) * nout2_ * nout3_ * sizeof(Real);             // MHD B1f
        data_size += nout1_ * (nout2_ + 1) * nout3_ * sizeof(Real);             // MHD B2f
        data_size += nout1_ * nout2_ * (nout3_ + 1) * sizeof(Real);             // MHD B3f
    }
    if (config_.has_turbulence) {
        data_size += nout1_ * nout2_ * nout3_ * config_.nforce * sizeof(Real);
    }
    
    return data_size;
}

bool PhysicsReader::ReadHydroData(IOWrapperSizeT& myoffset, IOWrapperSizeT& offset_myrank, 
                                 IOWrapperSizeT data_size, int noutmbs_min, int noutmbs_max) {
    std::cout << "Reading Hydro data for MeshBlocks owned by this rank..." << std::endl;
    size_t hydro_size = config_.nhydro * nout1_ * nout2_ * nout3_;
    
    // Following AthenaK's MPI pattern: read only MeshBlocks for this rank
    for (int m = 0; m < noutmbs_max; ++m) {
        if (m < noutmbs_min) {
            if (m < nmb_thisrank_) {
                if (file_.Read_Reals_at_all(hydro_data_[m].data(), hydro_size, myoffset) != hydro_size) {
                    std::cerr << "ERROR: Failed to read hydro data for local MeshBlock " << m << std::endl;
                    return false;
                }
            }
            myoffset += data_size;
        } else if (m < nmb_thisrank_) {
            if (file_.Read_Reals_at(hydro_data_[m].data(), hydro_size, myoffset) != hydro_size) {
                std::cerr << "ERROR: Failed to read hydro data for local MeshBlock " << m << std::endl;
                return false;
            }
            myoffset += data_size;
        }
    }
    
    // Move to next physics module - following pgen.cpp line 338-339 EXACTLY
    offset_myrank += nout1_ * nout2_ * nout3_ * config_.nhydro * sizeof(Real);
    myoffset = offset_myrank;  // ✓ CORRECT: No additional rank offset
    
    return true;
}

bool PhysicsReader::ReadMHDData(IOWrapperSizeT& myoffset, IOWrapperSizeT& offset_myrank, 
                               IOWrapperSizeT data_size, int noutmbs_min, int noutmbs_max) {
    // Read MHD cell-centered data
    std::cout << "Reading MHD cell-centered data..." << std::endl;
    size_t mhd_cc_size = config_.nmhd * nout1_ * nout2_ * nout3_;
    
    for (int m = 0; m < noutmbs_max; ++m) {
        if (m < noutmbs_min) {
            if (m < nmb_thisrank_) {
                if (file_.Read_Reals_at_all(mhd_data_[m].data(), mhd_cc_size, myoffset) != mhd_cc_size) {
                    std::cerr << "ERROR: Failed to read MHD CC data for local MeshBlock " << m << std::endl;
                    return false;
                }
            }
            myoffset += data_size;
        } else if (m < nmb_thisrank_) {
            if (file_.Read_Reals_at(mhd_data_[m].data(), mhd_cc_size, myoffset) != mhd_cc_size) {
                std::cerr << "ERROR: Failed to read MHD CC data for local MeshBlock " << m << std::endl;
                return false;
            }
            myoffset += data_size;
        }
    }
    
    // Move to face fields section
    offset_myrank += nout1_ * nout2_ * nout3_ * config_.nmhd * sizeof(Real);
    myoffset = offset_myrank;
    
    // Read MHD face-centered fields
    std::cout << "Reading MHD face-centered fields..." << std::endl;
    size_t b1f_size = (nout1_ + 1) * nout2_ * nout3_;
    size_t b2f_size = nout1_ * (nout2_ + 1) * nout3_;
    size_t b3f_size = nout1_ * nout2_ * (nout3_ + 1);
    
    // Following AthenaK's exact sequential reading pattern from pgen.cpp:382-460
    for (int m = 0; m < noutmbs_max; ++m) {
        // every rank has a MB to read, so read collectively
        if (m < noutmbs_min) {
            if (m < nmb_thisrank_) {
                // Sequential reading: B1f -> B2f -> B3f exactly like AthenaK
                if (file_.Read_Reals_at_all(mhd_b1f_[m].data(), b1f_size, myoffset) != b1f_size) {
                    std::cerr << "ERROR: Failed to read B1f data for local MeshBlock " << m << std::endl;
                    return false;
                }
            }
            myoffset += b1f_size * sizeof(Real);
            
            if (m < nmb_thisrank_) {
                if (file_.Read_Reals_at_all(mhd_b2f_[m].data(), b2f_size, myoffset) != b2f_size) {
                    std::cerr << "ERROR: Failed to read B2f data for local MeshBlock " << m << std::endl;
                    return false;
                }
            }
            myoffset += b2f_size * sizeof(Real);
            
            if (m < nmb_thisrank_) {
                if (file_.Read_Reals_at_all(mhd_b3f_[m].data(), b3f_size, myoffset) != b3f_size) {
                    std::cerr << "ERROR: Failed to read B3f data for local MeshBlock " << m << std::endl;
                    return false;
                }
            }
            myoffset += b3f_size * sizeof(Real);
            
            // Skip remaining data in this MeshBlock - matching pgen.cpp:421
            myoffset += data_size - (b1f_size + b2f_size + b3f_size) * sizeof(Real);
            
        } else if (m < nmb_thisrank_) {
            // Sequential reading for remaining ranks
            if (file_.Read_Reals_at(mhd_b1f_[m].data(), b1f_size, myoffset) != b1f_size) {
                std::cerr << "ERROR: Failed to read B1f data for local MeshBlock " << m << std::endl;
                return false;
            }
            myoffset += b1f_size * sizeof(Real);
            
            if (file_.Read_Reals_at(mhd_b2f_[m].data(), b2f_size, myoffset) != b2f_size) {
                std::cerr << "ERROR: Failed to read B2f data for local MeshBlock " << m << std::endl;
                return false;
            }
            myoffset += b2f_size * sizeof(Real);
            
            if (file_.Read_Reals_at(mhd_b3f_[m].data(), b3f_size, myoffset) != b3f_size) {
                std::cerr << "ERROR: Failed to read B3f data for local MeshBlock " << m << std::endl;
                return false;
            }
            myoffset += b3f_size * sizeof(Real);
            
            myoffset += data_size - (b1f_size + b2f_size + b3f_size) * sizeof(Real);
        }
    }
    
    // Move past all face field sections
    offset_myrank += (nout1_ + 1) * nout2_ * nout3_ * sizeof(Real);    // B1f
    offset_myrank += nout1_ * (nout2_ + 1) * nout3_ * sizeof(Real);    // B2f
    offset_myrank += nout1_ * nout2_ * (nout3_ + 1) * sizeof(Real);    // B3f
    myoffset = offset_myrank;
    
    return true;
}

bool PhysicsReader::ReadTurbulenceData(IOWrapperSizeT& myoffset, IOWrapperSizeT& offset_myrank, 
                                      IOWrapperSizeT data_size, int noutmbs_min, int noutmbs_max) {
    std::cout << "Reading Turbulence data for MeshBlocks owned by this rank..." << std::endl;
    size_t turb_size = config_.nforce * nout1_ * nout2_ * nout3_;
    
    for (int m = 0; m < noutmbs_max; ++m) {
        if (m < noutmbs_min) {
            if (m < nmb_thisrank_) {
                if (file_.Read_Reals_at_all(turb_data_[m].data(), turb_size, myoffset) != turb_size) {
                    std::cerr << "ERROR: Failed to read turbulence data for local MeshBlock " << m << std::endl;
                    return false;
                }
            }
            myoffset += data_size;
        } else if (m < nmb_thisrank_) {
            if (file_.Read_Reals_at(turb_data_[m].data(), turb_size, myoffset) != turb_size) {
                std::cerr << "ERROR: Failed to read turbulence data for local MeshBlock " << m << std::endl;
                return false;
            }
            myoffset += data_size;
        }
    }
    
    offset_myrank += nout1_ * nout2_ * nout3_ * config_.nforce * sizeof(Real);
    myoffset = offset_myrank;
    
    return true;
}

void PhysicsReader::SkipRadiationData(IOWrapperSizeT& offset_myrank) const {
    std::cout << "Skipping Radiation data..." << std::endl;
    offset_myrank += nout1_ * nout2_ * nout3_ * config_.nrad * sizeof(Real);
}

void PhysicsReader::SkipZ4cADMData(IOWrapperSizeT& offset_myrank) const {
    if (config_.has_z4c) {
        std::cout << "Skipping Z4c data..." << std::endl;
        offset_myrank += nout1_ * nout2_ * nout3_ * config_.nz4c * sizeof(Real);
    }
    
    if (config_.has_adm) {
        std::cout << "Skipping ADM data..." << std::endl;
        offset_myrank += nout1_ * nout2_ * nout3_ * config_.nadm * sizeof(Real);
    }
}

void PhysicsReader::DisplayPhysicsData(int local_meshblock_id) const {
    if (local_meshblock_id >= nmb_thisrank_) {
        std::cerr << "ERROR: Local MeshBlock ID " << local_meshblock_id 
                  << " out of range (0-" << nmb_thisrank_-1 << ")" << std::endl;
        return;
    }
    
    std::cout << "\n" << std::string(60, '=') << std::endl;
    std::cout << "PHYSICS DATA - Local MeshBlock " << local_meshblock_id << std::endl;
    std::cout << std::string(60, '=') << std::endl;
    
    // Display sample values from the first few cells
    int sample_cells = std::min(5, nout1_ * nout2_ * nout3_);
    
    if (config_.has_mhd && !mhd_data_.empty()) {
        std::cout << "\nMHD Cell-Centered Data (first " << sample_cells << " cells):" << std::endl;
        std::cout << "  Variables: [density, mom_x, mom_y, mom_z, energy, Bx, By, Bz]" << std::endl;
        for (int cell = 0; cell < sample_cells; cell++) {
            std::cout << "  Cell " << cell << ": ";
            for (int var = 0; var < config_.nmhd; var++) {
                int idx = var * nout1_ * nout2_ * nout3_ + cell;
                std::cout << std::scientific << std::setprecision(6) << mhd_data_[local_meshblock_id][idx] << " ";
            }
            std::cout << std::endl;
        }
        
        // Show face-centered magnetic field values
        std::cout << "\nFace-centered magnetic fields (first " << sample_cells << " cells):" << std::endl;
        if (!mhd_b1f_.empty()) {
            std::cout << "B1f: ";
            for (int cell = 0; cell < sample_cells; cell++) {
                std::cout << std::scientific << std::setprecision(6) << mhd_b1f_[local_meshblock_id][cell] << " ";
            }
            std::cout << std::endl;
        }
        if (!mhd_b2f_.empty()) {
            std::cout << "B2f: ";
            for (int cell = 0; cell < sample_cells; cell++) {
                std::cout << std::scientific << std::setprecision(6) << mhd_b2f_[local_meshblock_id][cell] << " ";
            }
            std::cout << std::endl;
        }
        if (!mhd_b3f_.empty()) {
            std::cout << "B3f: ";
            for (int cell = 0; cell < sample_cells; cell++) {
                std::cout << std::scientific << std::setprecision(6) << mhd_b3f_[local_meshblock_id][cell] << " ";
            }
            std::cout << std::endl;
        }
    }
    
    if (config_.has_turbulence && !turb_data_.empty()) {
        std::cout << "\nTurbulence Force Data (first " << sample_cells << " cells):" << std::endl;
        for (int cell = 0; cell < sample_cells; cell++) {
            std::cout << "  Cell " << cell << ": ";
            for (int var = 0; var < config_.nforce; var++) {
                int idx = var * nout1_ * nout2_ * nout3_ + cell;
                std::cout << std::scientific << std::setprecision(6) << turb_data_[local_meshblock_id][idx] << " ";
            }
            std::cout << std::endl;
        }
    }
}

void PhysicsReader::ValidateUniformValues() const {
    if (!config_.has_mhd || mhd_data_.empty()) {
        std::cout << "No MHD data available for validation" << std::endl;
        return;
    }
    
    std::cout << "\n" << std::string(60, '=') << std::endl;
    std::cout << "DATA SUMMARY FOR FIRST LOCAL MESHBLOCK" << std::endl;
    std::cout << std::string(60, '=') << std::endl;
    
    int mb = 0;
    int cell = 0;
    
    // Show first cell values
    Real density = mhd_data_[mb][0 * nout1_ * nout2_ * nout3_ + cell];
    Real energy = mhd_data_[mb][4 * nout1_ * nout2_ * nout3_ + cell];
    
    std::cout << std::fixed << std::setprecision(6);
    std::cout << "Cell-centered values (first cell):" << std::endl;
    std::cout << "  ρ (density) = " << density << std::endl;
    std::cout << "  E (energy) = " << energy << std::endl;
    
    if (!mhd_b1f_.empty()) {
        std::cout << "  B1f (x1-face) = " << mhd_b1f_[mb][cell] << std::endl;
    }
    if (!mhd_b2f_.empty()) {
        std::cout << "  B2f (x2-face) = " << mhd_b2f_[mb][cell] << std::endl;
    }
    if (!mhd_b3f_.empty()) {
        std::cout << "  B3f (x3-face) = " << mhd_b3f_[mb][cell] << std::endl;
    }
    
    // Check uniformity
    bool is_uniform = true;
    int check_cells = std::min(5, nout1_ * nout2_ * nout3_);
    
    for (int c = 1; c < check_cells; c++) {
        Real density_c = mhd_data_[mb][0 * nout1_ * nout2_ * nout3_ + c];
        Real energy_c = mhd_data_[mb][4 * nout1_ * nout2_ * nout3_ + c];
        
        if (std::abs(density_c - density) > TOLERANCE || std::abs(energy_c - energy) > TOLERANCE) {
            is_uniform = false;
            break;
        }
    }
    
    std::cout << "\nUniformity check: Data appears " << (is_uniform ? "uniform" : "non-uniform") 
              << " across first " << check_cells << " cells" << std::endl;
}