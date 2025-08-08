#include "restart_writer.hpp"
#include "upscaler.hpp"
#include "parameter_parser.hpp"
#include "simple_parameter_input.hpp"
#include <iostream>
#include <sstream>
#include <cstring>

#if MPI_PARALLEL_ENABLED
#include <mpi.h>
#endif

// RestartWriter implementation
RestartWriter::RestartWriter(RestartReader& reader, Upscaler& upscaler) 
    : reader_(reader), upscaler_(upscaler), header_offset_(0), 
      my_rank_(0), nranks_(1), nmb_thisrank_(0), 
      noutmbs_min_(0), noutmbs_max_(0) {
    
    // Get MPI info from reader
    my_rank_ = reader_.GetMyRank();
    nranks_ = reader_.GetNRanks();
    
    // Use the upscaler's fine mesh distribution
    // Each rank has already upscaled its assigned coarse meshblocks
    size_t fine_nmb_total = GetFineNMBTotal();
    nmb_thisrank_ = upscaler_.fine_nmb_thisrank_;
    
    // Setup MPI distribution that matches what upscaler created
    // Each rank gets nfine_per_coarse_ fine MBs for each of its coarse MBs
    mpi_dist_ = std::make_unique<MPIDistribution>(nranks_, fine_nmb_total);
    
    // Calculate gids and nmb for each rank based on upscaler's distribution
    std::vector<int> gids_eachrank(nranks_);
    std::vector<int> nmb_eachrank(nranks_);
    
#if MPI_PARALLEL_ENABLED
    // Gather the fine mesh distribution from all ranks
    int my_fine_start = upscaler_.fine_gids_start_;
    int my_fine_count = upscaler_.fine_nmb_thisrank_;
    
    MPI_Allgather(&my_fine_start, 1, MPI_INT, gids_eachrank.data(), 1, MPI_INT, MPI_COMM_WORLD);
    MPI_Allgather(&my_fine_count, 1, MPI_INT, nmb_eachrank.data(), 1, MPI_INT, MPI_COMM_WORLD);
#else
    gids_eachrank[0] = 0;
    nmb_eachrank[0] = fine_nmb_total;
#endif
    
    mpi_dist_->SetupDistribution(&gids_eachrank, &nmb_eachrank);
    mpi_dist_->GetRankMinMax(&noutmbs_min_, &noutmbs_max_);
    
    if (my_rank_ == 0) {
        std::cout << "RestartWriter initialized with MPI support:" << std::endl;
        std::cout << "  Number of ranks: " << nranks_ << std::endl;
        std::cout << "  Fine mesh total MeshBlocks: " << fine_nmb_total << std::endl;
        mpi_dist_->PrintDistribution();
    }
}

RegionSize RestartWriter::GetFineMeshSize() const { 
    return upscaler_.fine_mesh_size_; 
}

RegionIndcs RestartWriter::GetFineMeshIndcs() const { 
    return upscaler_.fine_mesh_indcs_; 
}

RegionIndcs RestartWriter::GetFineMBIndcs() const { 
    return upscaler_.fine_mb_indcs_; 
}

size_t RestartWriter::GetFineNMBTotal() const { 
    return upscaler_.fine_nmb_total_; 
}

const std::vector<LogicalLocation>& RestartWriter::GetFineLlocEachMB() const { 
    return upscaler_.fine_lloc_eachmb_; 
}

const std::vector<float>& RestartWriter::GetFineCostEachMB() const { 
    return upscaler_.fine_cost_eachmb_; 
}

const Real* RestartWriter::GetFineHydroData() const { 
    return upscaler_.fine_hydro_data_.data(); 
}

const Real* RestartWriter::GetFineMHDData() const { 
    return upscaler_.fine_mhd_data_.data(); 
}

const Real* RestartWriter::GetFineX1FData() const { 
    return upscaler_.fine_x1f_data_.data(); 
}

const Real* RestartWriter::GetFineX2FData() const { 
    return upscaler_.fine_x2f_data_.data(); 
}

const Real* RestartWriter::GetFineX3FData() const { 
    return upscaler_.fine_x3f_data_.data(); 
}

const Real* RestartWriter::GetFineRadData() const { 
    return upscaler_.fine_rad_data_.data(); 
}

const Real* RestartWriter::GetFineTurbData() const { 
    return upscaler_.fine_turb_data_.data(); 
}

const Real* RestartWriter::GetFineZ4cData() const { 
    return upscaler_.fine_z4c_data_.data(); 
}

const Real* RestartWriter::GetFineADMData() const { 
    return upscaler_.fine_adm_data_.data(); 
}

const std::vector<int>& RestartWriter::GetFineGidsEachRank() const { 
    return mpi_dist_->GetGidsEachRank(); 
}

const std::vector<int>& RestartWriter::GetFineNmbEachRank() const { 
    return mpi_dist_->GetNmbEachRank(); 
}

bool RestartWriter::WriteRestartFile(const std::string& filename) {
    // All ranks must call Open() for MPI_File_open() to work (collective operation)
    if (!file_.Open(filename.c_str(), IOWrapper::FileMode::write)) {
        std::cerr << "ERROR: Cannot create restart file: " << filename << std::endl;
        return false;
    }
    
    // Follow AthenaK's exact writing sequence
    if (!WriteParameterData()) {
        std::cerr << "ERROR: Failed to write parameter data" << std::endl;
        return false;
    }
    
    if (!WriteHeaderData()) {
        std::cerr << "ERROR: Failed to write header data" << std::endl;
        return false;
    }
    
    if (!WriteMeshBlockLists()) {
        std::cerr << "ERROR: Failed to write MeshBlock lists" << std::endl;
        return false;
    }
    
    if (!WriteInternalState()) {
        std::cerr << "ERROR: Failed to write internal state data" << std::endl;
        return false;
    }
    
    if (!WritePhysicsData()) {
        std::cerr << "ERROR: Failed to write physics data" << std::endl;
        return false;
    }
    
    file_.Close();
    return true;
}

bool RestartWriter::WriteParameterData() {
    if (reader_.GetMyRank() == 0) {
        std::cout << "WriteParameterData: Writing upscaled parameter data" << std::endl;
        
        std::string original_params = reader_.GetParameterString();
        std::cout << "  Original parameter string size: " << original_params.size() << " bytes" << std::endl;
        
        // Update mesh parameters for upscaling (double the resolution) using string replacement
        // This preserves the original formatting, trailing whitespace, and PAR_DUMP line
        RegionIndcs fine_indcs = GetFineMeshIndcs();
        
        // Create modified parameter string by direct replacement
        std::string param_string = original_params;
        
        // Find the <mesh> section for targeted replacement
        size_t mesh_start = param_string.find("<mesh>");
        size_t mesh_end = param_string.find("<meshblock>");
        
        if (mesh_start == std::string::npos) {
            std::cerr << "ERROR: Could not find <mesh> section in parameters" << std::endl;
            return false;
        }
        if (mesh_end == std::string::npos) {
            mesh_end = param_string.find("<", mesh_start + 6); // Find next section
        }
        
        
        // Helper lambda to replace nx values while preserving exact formatting
        auto replace_nx = [&](const std::string& param_name, int new_value) -> bool {
            // Look for pattern "nx1    = 64" in mesh section only
            std::string search_pattern = param_name + "    = ";
            size_t pos = param_string.find(search_pattern, mesh_start);
            
            
            if (pos != std::string::npos && pos < mesh_end) {
                size_t value_start = pos + search_pattern.length();
                size_t value_end = param_string.find_first_of(" \t", value_start);
                if (value_end == std::string::npos || value_end > mesh_end) {
                    value_end = param_string.find('\n', value_start);
                }
                
                std::string old_value = param_string.substr(value_start, value_end - value_start);
                param_string.replace(value_start, value_end - value_start, std::to_string(new_value));
                return true;
            }
            std::cerr << "    WARNING: Could not find " << param_name << " in mesh section" << std::endl;
            return false;
        };
        
        // Update the mesh parameters (only in <mesh> section, not <meshblock>)
        replace_nx("nx1", fine_indcs.nx1);
        replace_nx("nx2", fine_indcs.nx2);
        replace_nx("nx3", fine_indcs.nx3);
        
        // Note: We do NOT update parameters in <meshblock> as the meshblock size remains constant
        // The upscaling maintains the same meshblock size but creates more meshblocks
        std::cout << "  Modified parameter string size: " << param_string.size() << " bytes" << std::endl;
        
        
        // Write the updated parameter string (preserves original format including <par_end>\n)
        IOWrapperSizeT bytes_written = file_.Write_any_type(param_string.c_str(), 
                                                           param_string.size(), "byte");
        if (bytes_written != param_string.size()) {
            std::cerr << "ERROR: Failed to write parameter data (wrote " << bytes_written 
                      << " bytes, expected " << param_string.size() << ")" << std::endl;
            return false;
        }
        
        // Store the header offset
        header_offset_ = file_.GetPosition();
        std::cout << "  Header offset set to: " << header_offset_ << std::endl;
    }
    
#if MPI_PARALLEL_ENABLED
    // Ensure all ranks wait for rank 0 to finish writing parameters
    MPI_Barrier(MPI_COMM_WORLD);
    // Broadcast header offset to all ranks
    MPI_Bcast(&header_offset_, sizeof(IOWrapperSizeT), MPI_BYTE, 0, MPI_COMM_WORLD);
#endif
    
    return true;
}

bool RestartWriter::WriteHeaderData() {
    if (reader_.GetMyRank() == 0) {
        // Write header data following AthenaK format
        int fine_nmb_total = static_cast<int>(GetFineNMBTotal());
        int root_level = reader_.GetRootLevel() + 1;  // Root level increases by 1 when doubling resolution
        RegionSize fine_mesh_size = GetFineMeshSize();
        RegionIndcs fine_mesh_indcs = GetFineMeshIndcs();
        RegionIndcs fine_mb_indcs = GetFineMBIndcs();
        Real time = reader_.GetTime();
        Real dt = reader_.GetDt();  // Keep original timestep - AthenaK will adjust if needed
        int ncycle = reader_.GetNCycle();
        
        // Debug output
        std::cout << "Writing header data:" << std::endl;
        std::cout << "  fine_nmb_total = " << fine_nmb_total << std::endl;
        std::cout << "  root_level = " << root_level << " (original was " << reader_.GetRootLevel() << ")" << std::endl;
        std::cout << "  fine_mesh_indcs.nx1 = " << fine_mesh_indcs.nx1 << std::endl;
        
        // Debug: Check file position before writing
        IOWrapperSizeT pos_before = file_.GetPosition();
        std::cout << "  File position before header write: " << pos_before << std::endl;
        
        IOWrapperSizeT bytes_written;
        bytes_written = file_.Write_any_type(&fine_nmb_total, sizeof(int), "byte");
        std::cout << "  Wrote nmb_total: " << bytes_written << " bytes" << std::endl;
        
        bytes_written = file_.Write_any_type(&root_level, sizeof(int), "byte");
        std::cout << "  Wrote root_level: " << bytes_written << " bytes" << std::endl;
        file_.Write_any_type(&fine_mesh_size, sizeof(RegionSize), "byte");
        file_.Write_any_type(&fine_mesh_indcs, sizeof(RegionIndcs), "byte");
        file_.Write_any_type(&fine_mb_indcs, sizeof(RegionIndcs), "byte");
        file_.Write_any_type(&time, sizeof(Real), "byte");
        file_.Write_any_type(&dt, sizeof(Real), "byte");
        file_.Write_any_type(&ncycle, sizeof(int), "byte");
    }
    
#if MPI_PARALLEL_ENABLED
    // Ensure all ranks wait for rank 0 to finish writing header
    MPI_Barrier(MPI_COMM_WORLD);
#endif
    
    return true;
}

bool RestartWriter::WriteMeshBlockLists() {
    if (reader_.GetMyRank() == 0) {
        const std::vector<LogicalLocation>& fine_llocs = GetFineLlocEachMB();
        const std::vector<float>& fine_costs = GetFineCostEachMB();
        size_t fine_nmb_total = GetFineNMBTotal();
        
        file_.Write_any_type(fine_llocs.data(), fine_nmb_total * sizeof(LogicalLocation), "byte");
        file_.Write_any_type(fine_costs.data(), fine_nmb_total * sizeof(float), "byte");
    }
    
    return true;
}

bool RestartWriter::WriteInternalState() {
    if (reader_.GetMyRank() == 0) {
        const PhysicsConfig& phys_config = reader_.GetPhysicsConfig();
        
        // Write internal state data following AthenaK format (src/outputs/restart.cpp:219-235)
        
        // Write Z4c last output time if present
        if (phys_config.has_z4c) {
            Real last_output_time = reader_.GetZ4cLastOutputTime();
            file_.Write_any_type(&last_output_time, sizeof(Real), "byte");
        }
        
        // Write compact object tracker data if present
        if (phys_config.nco > 0) {
            const std::vector<Real>& tracker_data = reader_.GetTrackerData();
            file_.Write_any_type(tracker_data.data(), tracker_data.size() * sizeof(Real), "byte");
        }
        
        // Write turbulence RNG state if present
        if (phys_config.has_turbulence) {
            const RNG_State& rng_state = reader_.GetRNGState();
            file_.Write_any_type(&rng_state, sizeof(RNG_State), "byte");
        }
    }
    
#if MPI_PARALLEL_ENABLED
    // Ensure all ranks wait for rank 0 to finish writing internal state
    MPI_Barrier(MPI_COMM_WORLD);
#endif
    
    return true;
}

bool RestartWriter::WritePhysicsData() {
    const PhysicsConfig& phys_config = reader_.GetPhysicsConfig();
    RegionIndcs mb_indcs = GetFineMBIndcs();
    size_t fine_nmb_total = GetFineNMBTotal();
    
    int nx1 = mb_indcs.nx1;
    int nx2 = mb_indcs.nx2;
    int nx3 = mb_indcs.nx3;
    int ng = mb_indcs.ng;
    
    // Total cells including ghost zones
    int nout1 = nx1 + 2*ng;
    int nout2 = (nx2 > 1) ? nx2 + 2*ng : 1;
    int nout3 = (nx3 > 1) ? nx3 + 2*ng : 1;
    
    // Calculate data size per meshblock
    IOWrapperSizeT data_size = 0;
    if (phys_config.has_hydro) {
        data_size += nout1*nout2*nout3*phys_config.nhydro*sizeof(Real);
    }
    if (phys_config.has_mhd) {
        data_size += nout1*nout2*nout3*phys_config.nmhd*sizeof(Real);   // mhd u0
        data_size += (nout1+1)*nout2*nout3*sizeof(Real);    // mhd b0.x1f
        data_size += nout1*(nout2+1)*nout3*sizeof(Real);    // mhd b0.x2f
        data_size += nout1*nout2*(nout3+1)*sizeof(Real);    // mhd b0.x3f
    }
    if (phys_config.has_radiation) {
        data_size += nout1*nout2*nout3*phys_config.nrad*sizeof(Real);
    }
    if (phys_config.has_turbulence) {
        data_size += nout1*nout2*nout3*phys_config.nforce*sizeof(Real);
    }
    if (phys_config.has_z4c) {
        data_size += nout1*nout2*nout3*phys_config.nz4c*sizeof(Real);
    } else if (phys_config.has_adm) {
        data_size += nout1*nout2*nout3*phys_config.nadm*sizeof(Real);
    }
    
    // Write data size
    if (my_rank_ == 0) {
        file_.Write_any_type(&data_size, sizeof(IOWrapperSizeT), "byte");
    }
    
#if MPI_PARALLEL_ENABLED
    // Ensure all ranks wait for rank 0 to write data size
    MPI_Barrier(MPI_COMM_WORLD);
#endif
    
    // Calculate offset sizes following AthenaK's pattern (src/outputs/restart.cpp:267-278)
    IOWrapperSizeT step1size = header_offset_;  // Parameter data size
    IOWrapperSizeT step2size = fine_nmb_total * (sizeof(LogicalLocation) + sizeof(float));  // MB lists
    IOWrapperSizeT step3size = 0;  // Internal state data size
    if (phys_config.has_z4c) step3size += sizeof(Real);
    if (phys_config.nco > 0) step3size += 3 * phys_config.nco * sizeof(Real);
    if (phys_config.has_turbulence) step3size += sizeof(RNG_State);
    
    // Add header data size (3*sizeof(int) + 2*sizeof(Real) + sizeof(RegionSize) + 2*sizeof(RegionIndcs))
    step1size += 3*sizeof(int) + 2*sizeof(Real) + sizeof(RegionSize) + 2*sizeof(RegionIndcs);
    
    // Calculate offset following AthenaK's pattern (src/outputs/restart.cpp:277-279)
    // Each rank calculates its starting offset based on its starting global ID
    IOWrapperSizeT offset_myrank = step1size + step2size + step3size + sizeof(IOWrapperSizeT) 
                                  + data_size * mpi_dist_->GetStartingMeshBlockID(my_rank_);
    IOWrapperSizeT myoffset = offset_myrank;
    
    // Write physics data following AthenaK's exact pattern with MPI
    int cells_per_mb = nout1 * nout2 * nout3;
    
    // Write hydro data following AthenaK's pattern (src/outputs/restart.cpp:284-317)
    if (phys_config.has_hydro) {
        const Real* hydro_data = GetFineHydroData();
        // Each rank writes its assigned MeshBlocks
        for (int m = 0; m < noutmbs_max_; ++m) {
            // Every rank has a MB to write, so write collectively
            if (m < noutmbs_min_) {
                if (m < nmb_thisrank_) {
                    const Real* mb_data = hydro_data + m * phys_config.nhydro * cells_per_mb;
                    IOWrapperSizeT mbcnt = phys_config.nhydro * cells_per_mb;
                    if (file_.Write_any_type_at_all(mb_data, mbcnt, myoffset, "Real") != mbcnt) {
                        std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
                                  << std::endl << "cell-centered hydro data not written correctly to rst file, "
                                  << "restart file is broken." << std::endl;
                        return false;
                    }
                }
                myoffset += data_size;
            // Some ranks are finished writing, so use non-collective write
            } else if (m < nmb_thisrank_) {
                const Real* mb_data = hydro_data + m * phys_config.nhydro * cells_per_mb;
                IOWrapperSizeT mbcnt = phys_config.nhydro * cells_per_mb;
                if (file_.Write_any_type_at(mb_data, mbcnt, myoffset, "Real") != mbcnt) {
                    std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
                              << std::endl << "cell-centered hydro data not written correctly to rst file, "
                              << "restart file is broken." << std::endl;
                    return false;
                }
                myoffset += data_size;
            }
        }
        offset_myrank += nout1*nout2*nout3*phys_config.nhydro*sizeof(Real); // hydro u0
        myoffset = offset_myrank;
        
    }
    
    // Write MHD data following AthenaK's pattern (src/outputs/restart.cpp:318-350)
    // IMPORTANT: AthenaK writes cell-centered data for ALL meshblocks first,
    // THEN face-centered data for ALL meshblocks
    if (phys_config.has_mhd) {
        // First write ALL cell-centered MHD data
        const Real* mhd_data = GetFineMHDData();
        for (int m = 0; m < noutmbs_max_; ++m) {
            // Every rank has a MB to write, so write collectively
            if (m < noutmbs_min_) {
                if (m < nmb_thisrank_) {
                    const Real* mb_data = mhd_data + m * phys_config.nmhd * cells_per_mb;
                    IOWrapperSizeT mbcnt = phys_config.nmhd * cells_per_mb;
                    if (file_.Write_any_type_at_all(mb_data, mbcnt, myoffset, "Real") != mbcnt) {
                        std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
                                  << std::endl << "cell-centered mhd data not written correctly to rst file, "
                                  << "restart file is broken." << std::endl;
                        return false;
                    }
                }
                myoffset += data_size;
            // Some ranks are finished writing, so use non-collective write
            } else if (m < nmb_thisrank_) {
                const Real* mb_data = mhd_data + m * phys_config.nmhd * cells_per_mb;
                IOWrapperSizeT mbcnt = phys_config.nmhd * cells_per_mb;
                if (file_.Write_any_type_at(mb_data, mbcnt, myoffset, "Real") != mbcnt) {
                    std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
                              << std::endl << "cell-centered mhd data not written correctly to rst file, "
                              << "restart file is broken." << std::endl;
                    return false;
                }
                myoffset += data_size;
            }
        }
        offset_myrank += nout1*nout2*nout3*phys_config.nmhd*sizeof(Real);   // mhd u0
        myoffset = offset_myrank;
        
        // Write MHD face-centered data following AthenaK's pattern (src/outputs/restart.cpp:352-432)
        // Then write ALL face-centered data
        const Real* x1f_data = GetFineX1FData();
        const Real* x2f_data = GetFineX2FData();
        const Real* x3f_data = GetFineX3FData();
        
        int x1f_size = nout3 * nout2 * (nout1 + 1);
        int x2f_size = nout3 * (nout2 + 1) * nout1;
        int x3f_size = (nout3 + 1) * nout2 * nout1;
        
        for (int m = 0; m < noutmbs_max_; ++m) {
            // Every rank has a MB to write, so write collectively
            if (m < noutmbs_min_) {
                if (m < nmb_thisrank_) {
                    // Write x1f, x2f, x3f for this meshblock in sequence following AthenaK pattern
                    const Real* x1f_mb = x1f_data + m * x1f_size;
                    if (file_.Write_any_type_at_all(x1f_mb, x1f_size, myoffset, "Real") != x1f_size) {
                        std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
                                  << std::endl << "b0.x1f data not written correctly to rst file, "
                                  << "restart file is broken." << std::endl;
                        return false;
                    }
                    myoffset += x1f_size*sizeof(Real);
                    
                    const Real* x2f_mb = x2f_data + m * x2f_size;
                    if (file_.Write_any_type_at_all(x2f_mb, x2f_size, myoffset, "Real") != x2f_size) {
                        std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
                                  << std::endl << "b0.x2f data not written correctly to rst file, "
                                  << "restart file is broken." << std::endl;
                        return false;
                    }
                    myoffset += x2f_size*sizeof(Real);
                    
                    const Real* x3f_mb = x3f_data + m * x3f_size;
                    if (file_.Write_any_type_at_all(x3f_mb, x3f_size, myoffset, "Real") != x3f_size) {
                        std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
                                  << std::endl << "b0.x3f data not written correctly to rst file, "
                                  << "restart file is broken." << std::endl;
                        return false;
                    }
                    myoffset += x3f_size*sizeof(Real);
                    
                    // Adjust offset to account for other physics data in this meshblock
                    myoffset += data_size-(x1f_size+x2f_size+x3f_size)*sizeof(Real);
                } else {
                    // Rank doesn't have this MB but needs to advance offset for collective write
                    myoffset += data_size;
                }
            // Some ranks are finished writing, so use non-collective write
            } else if (m < nmb_thisrank_) {
                // Write x1f, x2f, x3f for this meshblock in sequence following AthenaK pattern
                const Real* x1f_mb = x1f_data + m * x1f_size;
                if (file_.Write_any_type_at(x1f_mb, x1f_size, myoffset, "Real") != x1f_size) {
                    std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
                              << std::endl << "b0.x1f data not written correctly to rst file, "
                              << "restart file is broken." << std::endl;
                    return false;
                }
                myoffset += x1f_size*sizeof(Real);
                
                const Real* x2f_mb = x2f_data + m * x2f_size;
                if (file_.Write_any_type_at(x2f_mb, x2f_size, myoffset, "Real") != x2f_size) {
                    std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
                              << std::endl << "b0.x2f data not written correctly to rst file, "
                              << "restart file is broken." << std::endl;
                    return false;
                }
                myoffset += x2f_size*sizeof(Real);
                
                const Real* x3f_mb = x3f_data + m * x3f_size;
                if (file_.Write_any_type_at(x3f_mb, x3f_size, myoffset, "Real") != x3f_size) {
                    std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
                              << std::endl << "b0.x3f data not written correctly to rst file, "
                              << "restart file is broken." << std::endl;
                    return false;
                }
                myoffset += x3f_size*sizeof(Real);
                
                // Adjust offset to account for other physics data in this meshblock
                myoffset += data_size-(x1f_size+x2f_size+x3f_size)*sizeof(Real);
            }
        }
        offset_myrank += (nout1+1)*nout2*nout3*sizeof(Real);    // mhd b0.x1f
        offset_myrank += nout1*(nout2+1)*nout3*sizeof(Real);    // mhd b0.x2f
        offset_myrank += nout1*nout2*(nout3+1)*sizeof(Real);    // mhd b0.x3f
        myoffset = offset_myrank;
    }
        
    // Write turbulence force data following AthenaK's pattern (src/outputs/restart.cpp:469-502)
    if (phys_config.has_turbulence) {
        const Real* turb_data = GetFineTurbData();
        for (int m = 0; m < noutmbs_max_; ++m) {
            // Every rank has a MB to write, so write collectively
            if (m < noutmbs_min_) {
                if (m < nmb_thisrank_) {
                    const Real* mb_data = turb_data + m * phys_config.nforce * cells_per_mb;
                    IOWrapperSizeT mbcnt = phys_config.nforce * cells_per_mb;
                    if (file_.Write_any_type_at_all(mb_data, mbcnt, myoffset, "Real") != mbcnt) {
                        std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
                                  << std::endl << "cell-centered turb data not written correctly to rst file, "
                                  << "restart file is broken." << std::endl;
                        return false;
                    }
                }
                myoffset += data_size;
            // Some ranks are finished writing, so use non-collective write
            } else if (m < nmb_thisrank_) {
                const Real* mb_data = turb_data + m * phys_config.nforce * cells_per_mb;
                IOWrapperSizeT mbcnt = phys_config.nforce * cells_per_mb;
                if (file_.Write_any_type_at(mb_data, mbcnt, myoffset, "Real") != mbcnt) {
                    std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
                              << std::endl << "cell-centered turb data not written correctly to rst file, "
                              << "restart file is broken." << std::endl;
                    return false;
                }
                myoffset += data_size;
            }
        }
        offset_myrank += nout1*nout2*nout3*phys_config.nforce*sizeof(Real); // forcing
        myoffset = offset_myrank;
    }
        
    // Write Z4c data following AthenaK's pattern (src/outputs/restart.cpp:504-536)
    if (phys_config.has_z4c) {
        const Real* z4c_data = GetFineZ4cData();
        for (int m = 0; m < noutmbs_max_; ++m) {
            // Every rank has a MB to write, so write collectively
            if (m < noutmbs_min_) {
                if (m < nmb_thisrank_) {
                    const Real* mb_data = z4c_data + m * phys_config.nz4c * cells_per_mb;
                    IOWrapperSizeT mbcnt = phys_config.nz4c * cells_per_mb;
                    if (file_.Write_any_type_at_all(mb_data, mbcnt, myoffset, "Real") != mbcnt) {
                        std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
                                  << std::endl << "cell-centered z4c data not written correctly to rst file, "
                                  << "restart file is broken." << std::endl;
                        return false;
                    }
                }
                myoffset += data_size;
            // Some ranks are finished writing, so use non-collective write
            } else if (m < nmb_thisrank_) {
                const Real* mb_data = z4c_data + m * phys_config.nz4c * cells_per_mb;
                IOWrapperSizeT mbcnt = phys_config.nz4c * cells_per_mb;
                if (file_.Write_any_type_at(mb_data, mbcnt, myoffset, "Real") != mbcnt) {
                    std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
                              << std::endl << "cell-centered z4c data not written correctly to rst file, "
                              << "restart file is broken." << std::endl;
                    return false;
                }
                myoffset += data_size;
            }
        }
        offset_myrank += nout1*nout2*nout3*phys_config.nz4c*sizeof(Real); // z4c u0
        myoffset = offset_myrank;
    } else if (phys_config.has_adm) {
        const Real* adm_data = GetFineADMData();
        for (int m = 0; m < noutmbs_max_; ++m) {
            // Every rank has a MB to write, so write collectively
            if (m < noutmbs_min_) {
                if (m < nmb_thisrank_) {
                    const Real* mb_data = adm_data + m * phys_config.nadm * cells_per_mb;
                    IOWrapperSizeT mbcnt = phys_config.nadm * cells_per_mb;
                    if (file_.Write_any_type_at_all(mb_data, mbcnt, myoffset, "Real") != mbcnt) {
                        std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
                                  << std::endl << "cell-centered adm data not written correctly to rst file, "
                                  << "restart file is broken." << std::endl;
                        return false;
                    }
                }
                myoffset += data_size;
            // Some ranks are finished writing, so use non-collective write
            } else if (m < nmb_thisrank_) {
                const Real* mb_data = adm_data + m * phys_config.nadm * cells_per_mb;
                IOWrapperSizeT mbcnt = phys_config.nadm * cells_per_mb;
                if (file_.Write_any_type_at(mb_data, mbcnt, myoffset, "Real") != mbcnt) {
                    std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
                              << std::endl << "cell-centered adm data not written correctly to rst file, "
                              << "restart file is broken." << std::endl;
                    return false;
                }
                myoffset += data_size;
            }
        }
        offset_myrank += nout1*nout2*nout3*phys_config.nadm*sizeof(Real); // adm u_adm
        myoffset = offset_myrank;
    }
    
    // Write radiation data following AthenaK's pattern (src/outputs/restart.cpp:434-467)
    if (phys_config.has_radiation) {
        const Real* rad_data = GetFineRadData();
        for (int m = 0; m < noutmbs_max_; ++m) {
            // Every rank has a MB to write, so write collectively
            if (m < noutmbs_min_) {
                if (m < nmb_thisrank_) {
                    const Real* mb_data = rad_data + m * phys_config.nrad * cells_per_mb;
                    IOWrapperSizeT mbcnt = phys_config.nrad * cells_per_mb;
                    if (file_.Write_any_type_at_all(mb_data, mbcnt, myoffset, "Real") != mbcnt) {
                        std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
                                  << std::endl << "cell-centered rad data not written correctly to rst file, "
                                  << "restart file is broken." << std::endl;
                        return false;
                    }
                }
                myoffset += data_size;
            // Some ranks are finished writing, so use non-collective write
            } else if (m < nmb_thisrank_) {
                const Real* mb_data = rad_data + m * phys_config.nrad * cells_per_mb;
                IOWrapperSizeT mbcnt = phys_config.nrad * cells_per_mb;
                if (file_.Write_any_type_at(mb_data, mbcnt, myoffset, "Real") != mbcnt) {
                    std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
                              << std::endl << "cell-centered rad data not written correctly to rst file, "
                              << "restart file is broken." << std::endl;
                    return false;
                }
                myoffset += data_size;
            }
        }
        offset_myrank += nout1*nout2*nout3*phys_config.nrad*sizeof(Real);   // radiation i0
        myoffset = offset_myrank;
    }
    
#if MPI_PARALLEL_ENABLED
    // Ensure all ranks finish writing
    MPI_Barrier(MPI_COMM_WORLD);
#endif
    
    return true;
}