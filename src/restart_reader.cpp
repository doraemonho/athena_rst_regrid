#include "restart_reader.hpp"
#include <iomanip>
#include <sstream>
#include <cstring>

#if MPI_PARALLEL_ENABLED
#include <mpi.h>
#endif

RestartReader::RestartReader() 
    : header_offset_(0), my_rank_(0), nranks_(1), nmb_total_(0), root_level_(0), 
      time_(0.0), dt_(0.0), ncycle_(0), nout1_(0), nout2_(0), nout3_(0),
      z4c_last_output_time_(0.0) {
    
    // Initialize structures
    std::memset(&mesh_size_, 0, sizeof(RegionSize));
    std::memset(&mesh_indcs_, 0, sizeof(RegionIndcs));
    std::memset(&mb_indcs_, 0, sizeof(RegionIndcs));
    std::memset(&rng_state_, 0, sizeof(RNG_State));
}

RestartReader::~RestartReader() {
    // All ranks opened the file (collective operation), so all ranks must close it
    file_.Close();
}

bool RestartReader::ReadRestartFile(const std::string& filename, int my_rank, int nranks, 
                                   const std::vector<int>* gids_eachrank,
                                   const std::vector<int>* nmb_eachrank) {
    filename_ = filename;
    my_rank_ = my_rank;
    nranks_ = nranks;
    
    if (my_rank_ == 0) {
        std::cout << "Reading AthenaK restart file: " << filename << std::endl;
    }
    
    // All ranks must call Open() for MPI_File_open() to work (collective operation)
    if (!file_.Open(filename.c_str(), IOWrapper::FileMode::read)) {
        std::cerr << "ERROR: Cannot open restart file: " << filename << std::endl;
        return false;
    }
    
    // Follow AthenaK's exact reading sequence
    if (!ReadParameterData()) {
        std::cerr << "ERROR: Failed to read parameter data" << std::endl;
        return false;
    }
    
    if (!ReadHeaderData()) {
        std::cerr << "ERROR: Failed to read header data" << std::endl;
        return false;
    }
    
    if (!ReadMeshBlockLists()) {
        std::cerr << "ERROR: Failed to read MeshBlock lists" << std::endl;
        return false;
    }
    
    if (!ReadInternalState()) {
        std::cerr << "ERROR: Failed to read internal state data" << std::endl;
        return false;
    }
    
    // Setup MPI distribution after reading mesh structure
    mpi_dist_ = std::make_unique<MPIDistribution>(nranks_, nmb_total_);
    mpi_dist_->SetupDistribution(gids_eachrank, nmb_eachrank);

    // read the data size and advance the file pointer - MPI-aware
    if (my_rank_ == 0) {
        if (file_.Read_bytes(&data_size, sizeof(IOWrapperSizeT), 1) != 1) {
            std::cerr << "ERROR: Failed to read data size" << std::endl;
            return false;
        }
    }
#if MPI_PARALLEL_ENABLED
    MPI_Bcast(&data_size, sizeof(IOWrapperSizeT), MPI_BYTE, 0, MPI_COMM_WORLD);
#endif
    
    // Debug output
    std::cout << "Rank " << my_rank_ << ": data_size=" << data_size << std::endl;
    
    if (!SetupPhysicsReader()) {
        std::cerr << "ERROR: Failed to setup physics reader" << std::endl;
        return false;
    }
    
    if (my_rank_ == 0) {
        std::cout << "Successfully read restart file!" << std::endl;
    }
    return true;
}

bool RestartReader::ReadParameterData() {
    if (my_rank_ == 0) {
        std::cout << "Reading parameter data..." << std::endl;
    }
    
    // Read parameter data until <par_end> marker is found - MPI-aware
    std::stringstream param_stream;
    constexpr int kBufSize = 4096;
    char buf[kBufSize];
    IOWrapperSizeT header = 0;
    
    // Only rank 0 reads the parameter data
    if (my_rank_ == 0) {
        // Search for <par_end> marker - following AthenaK parameter_input.cpp:179-204 exactly
        IOWrapperSizeT ret;
        do {
            ret = file_.Read_bytes(buf, sizeof(char), kBufSize);
            if (ret == 0) {
                std::cerr << "ERROR: Reached end of file without finding <par_end>" << std::endl;
                return false;
            }
            
            param_stream.write(buf, ret);  // add the buffer into the stream
            header += ret;
            std::string sbuf = param_stream.str();  // create string for search
            size_t loc = sbuf.find("<par_end>", 0); // search from the top of the stream
            if (loc != std::string::npos) {  // found <par_end>
                header = loc + 10;  // store the header length (+10 for "<par_end>" length)
                break;
            }
            if (header > kBufSize*10) {  // AthenaK's 40KB limit check
                std::cerr << "ERROR: <par_end> is not found in the first 40KBytes. "
                          << "Probably the file is broken or the wrong file is specified" << std::endl;
                return false;
            }
        } while (ret == kBufSize);  // till EOF (or par_end is found) - AthenaK's exact condition
    }
    
    // Broadcast header offset and parameter string to all ranks
    if (my_rank_ == 0) {
        std::cout << "Rank 0: About to broadcast header offset: " << header << std::endl;
    } else {
        std::cout << "Rank " << my_rank_ << ": Waiting for header broadcast" << std::endl;
    }
#if MPI_PARALLEL_ENABLED
    MPI_Bcast(&header, sizeof(IOWrapperSizeT), MPI_BYTE, 0, MPI_COMM_WORLD);
#endif
    if (my_rank_ == 0) {
        std::cout << "Rank 0: Header broadcast complete" << std::endl;
    } else {
        std::cout << "Rank " << my_rank_ << ": Received header: " << header << std::endl;
    }
    
    // Store header offset (position after <par_end>)
    header_offset_ = header;
    if (my_rank_ == 0) {
        std::cout << "Header offset: " << header_offset_ << std::endl;
    }
    
    // Broadcast parameter string to all ranks
    if (my_rank_ == 0) {
        // Only get the parameter portion (up to header_offset_), not the entire buffer
        std::string full_buffer = param_stream.str();
        parameter_string_ = full_buffer.substr(0, header_offset_);
    }
    
#if MPI_PARALLEL_ENABLED
    // Broadcast string length first
    size_t param_length = parameter_string_.length();
    MPI_Bcast(&param_length, sizeof(size_t), MPI_BYTE, 0, MPI_COMM_WORLD);
    
    // Resize string on other ranks
    if (my_rank_ != 0) {
        parameter_string_.resize(param_length);
    }
    
    // Broadcast string content
    MPI_Bcast(const_cast<char*>(parameter_string_.data()), param_length, MPI_CHAR, 0, MPI_COMM_WORLD);
#endif
    
    ShowParameterSample(parameter_string_);
    
    physics_config_ = ParameterParser::ParseParameters(parameter_string_);
    
    if (my_rank_ == 0) {
        std::cout << "Parameter data size: " << header << " bytes" << std::endl;
    }
    return true;
}

bool RestartReader::ReadHeaderData() {
    std::cout << "Reading header data..." << std::endl;
    
    // Read fixed-size header data following build_tree.cpp:BuildTreeFromRestart() - MPI-aware
    IOWrapperSizeT headersize = 3*sizeof(int) + 2*sizeof(Real) + 
                               sizeof(RegionSize) + 2*sizeof(RegionIndcs);
    std::vector<char> headerdata(headersize);
    
    // Only rank 0 reads the header data
    if (my_rank_ == 0) {
        // Position file pointer after parameter data
        file_.Seek(header_offset_);
        
        if (file_.Read_bytes(headerdata.data(), 1, headersize) != headersize) {
            std::cerr << "ERROR: Header size read incorrectly" << std::endl;
            return false;
        }
    }
    
    // Broadcast header data to all ranks
#if MPI_PARALLEL_ENABLED
    MPI_Bcast(headerdata.data(), headersize, MPI_BYTE, 0, MPI_COMM_WORLD);
#endif
    
    // Extract header data in the exact order from AthenaK
    IOWrapperSizeT hdos = 0;
    std::memcpy(&nmb_total_, &(headerdata[hdos]), sizeof(int));
    hdos += sizeof(int);
    std::memcpy(&root_level_, &(headerdata[hdos]), sizeof(int));
    hdos += sizeof(int);
    std::memcpy(&mesh_size_, &(headerdata[hdos]), sizeof(RegionSize));
    hdos += sizeof(RegionSize);
    std::memcpy(&mesh_indcs_, &(headerdata[hdos]), sizeof(RegionIndcs));
    hdos += sizeof(RegionIndcs);
    std::memcpy(&mb_indcs_, &(headerdata[hdos]), sizeof(RegionIndcs));
    hdos += sizeof(RegionIndcs);
    std::memcpy(&time_, &(headerdata[hdos]), sizeof(Real));
    hdos += sizeof(Real);
    std::memcpy(&dt_, &(headerdata[hdos]), sizeof(Real));
    hdos += sizeof(Real);
    std::memcpy(&ncycle_, &(headerdata[hdos]), sizeof(int));
    
    // Calculate grid dimensions including ghost zones
    nout1_ = mb_indcs_.nx1 + 2 * mb_indcs_.ng;
    nout2_ = (mb_indcs_.nx2 > 1) ? (mb_indcs_.nx2 + 2 * mb_indcs_.ng) : 1;
    nout3_ = (mb_indcs_.nx3 > 1) ? (mb_indcs_.nx3 + 2 * mb_indcs_.ng) : 1;
    
    std::cout << "Header data read successfully:" << std::endl;
    std::cout << "  Total MeshBlocks: " << nmb_total_ << std::endl;
    std::cout << "  Root level: " << root_level_ << std::endl;
    std::cout << "  Time: " << time_ << std::endl;
    std::cout << "  Grid dimensions (with ghosts): " << nout1_ << " x " << nout2_ << " x " << nout3_ << std::endl;
    
    return true;
}

bool RestartReader::ReadMeshBlockLists() {
    std::cout << "Reading MeshBlock lists..." << std::endl;
    
    // Read logical locations and costs - MPI-aware
    IOWrapperSizeT listsize = sizeof(LogicalLocation) + sizeof(float);
    std::vector<char> idlist(listsize * nmb_total_);
    
    // Only rank 0 reads the mesh block lists
    if (my_rank_ == 0) {
        if (file_.Read_bytes(idlist.data(), listsize, nmb_total_) != static_cast<size_t>(nmb_total_)) {
            std::cerr << "ERROR: Incorrect number of MeshBlocks read" << std::endl;
            return false;
        }
    }
    
    // Broadcast mesh block lists to all ranks
#if MPI_PARALLEL_ENABLED
    MPI_Bcast(idlist.data(), listsize * nmb_total_, MPI_BYTE, 0, MPI_COMM_WORLD);
#endif
    // Parse logical locations and costs
    lloc_eachmb_.resize(nmb_total_);
    cost_eachmb_.resize(nmb_total_);
    
    int os = 0;
    for (int i = 0; i < nmb_total_; i++) {
        std::memcpy(&(lloc_eachmb_[i]), &(idlist[os]), sizeof(LogicalLocation));
        os += sizeof(LogicalLocation);
    }
    for (int i = 0; i < nmb_total_; i++) {
        std::memcpy(&(cost_eachmb_[i]), &(idlist[os]), sizeof(float));
        os += sizeof(float);
    }
    
    std::cout << "MeshBlock lists read successfully" << std::endl;
    return true;
}

bool RestartReader::ReadInternalState() {
    std::cout << "Reading internal state data..." << std::endl;
    std::cout << "Current position: " << file_.GetPosition() << std::endl;

    // Read and store Z4c data if present - MPI-aware
    if (physics_config_.has_z4c) {
        if (my_rank_ == 0) {
            if (file_.Read_Reals(&z4c_last_output_time_, 1) != 1) {
                std::cerr << "ERROR: Failed to read Z4c last_output_time" << std::endl;
                return false;
            }
        }
#if MPI_PARALLEL_ENABLED
        MPI_Bcast(&z4c_last_output_time_, sizeof(Real), MPI_BYTE, 0, MPI_COMM_WORLD);
#endif
    }
    
    // Read and store compact object tracker data - MPI-aware
    if (physics_config_.nco > 0) {
        tracker_data_.resize(3 * physics_config_.nco);
        if (my_rank_ == 0) {
            if (file_.Read_Reals(tracker_data_.data(), 3 * physics_config_.nco) != 
                static_cast<size_t>(3 * physics_config_.nco)) {
                std::cerr << "ERROR: Failed to read tracker data" << std::endl;
                return false;
            }
        }
#if MPI_PARALLEL_ENABLED
        MPI_Bcast(tracker_data_.data(), 3 * physics_config_.nco * sizeof(Real), MPI_BYTE, 0, MPI_COMM_WORLD);
#endif
    }
    
    // Read and store turbulence RNG state - MPI-aware
    if (physics_config_.has_turbulence) {
        if (my_rank_ == 0) {
            if (file_.Read_bytes(&rng_state_, sizeof(RNG_State), 1) != 1) {
                std::cerr << "ERROR: Failed to read RNG state" << std::endl;
                return false;
            }
        }
#if MPI_PARALLEL_ENABLED
        MPI_Bcast(&rng_state_, sizeof(RNG_State), MPI_BYTE, 0, MPI_COMM_WORLD);
#endif
    }
    
    std::cout << "Internal state data processed" << std::endl;
    std::cout << "Current position: " << file_.GetPosition() << std::endl;
    return true;
}

bool RestartReader::SetupPhysicsReader() {
    std::cout << "Setting up physics reader..." << std::endl;
    
    physics_reader_ = std::make_unique<PhysicsReader>(
        file_, physics_config_, *mpi_dist_, my_rank_, nout1_, nout2_, nout3_);
    
    // Get physics start position - MPI-aware
    IOWrapperSizeT physics_start = 0;
    if (my_rank_ == 0) {
        physics_start = file_.GetPosition();
    }
#if MPI_PARALLEL_ENABLED
    MPI_Bcast(&physics_start, sizeof(IOWrapperSizeT), MPI_BYTE, 0, MPI_COMM_WORLD);
#endif
    std::cout << "Rank " << my_rank_ << ": physics_start=" << physics_start << std::endl;
    return physics_reader_->ReadPhysicsData(physics_start, data_size);
}

void RestartReader::ShowParameterSample(const std::string& param_string) const {
    std::cout << "Parameter string sample (first 5000 chars):" << std::endl;
    std::string sample = param_string.substr(0, 5000);
    std::cout << sample << std::endl;
    if (param_string.length() > 5000) std::cout << "..." << std::endl;
}

void RestartReader::DisplayMeshInfo() const {
    std::cout << "\n" << std::string(60, '=') << std::endl;
    std::cout << "RESTART FILE INFORMATION" << std::endl;
    std::cout << std::string(60, '=') << std::endl;
    
    std::cout << "File: " << filename_ << std::endl;
    std::cout << "Time: " << std::scientific << std::setprecision(6) << time_ << std::endl;
    std::cout << "Timestep: " << std::scientific << std::setprecision(6) << dt_ << std::endl;
    std::cout << "Cycle: " << ncycle_ << std::endl;
    std::cout << "Total MeshBlocks: " << nmb_total_ << std::endl;
    std::cout << "Root level: " << root_level_ << std::endl;
    
    std::cout << "\nMesh dimensions:" << std::endl;
    std::cout << "  Physical domain: [" << mesh_size_.x1min << ", " << mesh_size_.x1max << "] x ["
              << mesh_size_.x2min << ", " << mesh_size_.x2max << "] x ["
              << mesh_size_.x3min << ", " << mesh_size_.x3max << "]" << std::endl;
    std::cout << "  Grid spacing: " << mesh_size_.dx1 << " x " << mesh_size_.dx2 << " x " << mesh_size_.dx3 << std::endl;
    std::cout << "  Mesh cells: " << mesh_indcs_.nx1 << " x " << mesh_indcs_.nx2 << " x " << mesh_indcs_.nx3 << std::endl;
    std::cout << "  MeshBlock cells: " << mb_indcs_.nx1 << " x " << mb_indcs_.nx2 << " x " << mb_indcs_.nx3 << std::endl;
    std::cout << "  Ghost cells: " << mb_indcs_.ng << std::endl;
    std::cout << "  Output dimensions: " << nout1_ << " x " << nout2_ << " x " << nout3_ << std::endl;
    
    std::cout << "\nPhysics modules:" << std::endl;
    std::cout << "  Hydro: " << (physics_config_.has_hydro ? "Yes" : "No") << " (nhydro=" << physics_config_.nhydro << ")" << std::endl;
    std::cout << "  MHD: " << (physics_config_.has_mhd ? "Yes" : "No") << " (nmhd=" << physics_config_.nmhd << ")" << std::endl;
    std::cout << "  Turbulence: " << (physics_config_.has_turbulence ? "Yes" : "No") << " (nforce=" << physics_config_.nforce << ")" << std::endl;
    std::cout << "  Radiation: " << (physics_config_.has_radiation ? "Yes" : "No") << std::endl;
    std::cout << "  Z4c: " << (physics_config_.has_z4c ? "Yes" : "No") << std::endl;
    
    if (mpi_dist_) {
        std::cout << "\nMPI Information:" << std::endl; 
        std::cout << "  My rank: " << my_rank_ << " / " << nranks_ << std::endl;
        std::cout << "  MeshBlocks for this rank: " << GetNMBThisRank() << std::endl;
        mpi_dist_->PrintDistribution();
    }
}

void RestartReader::DisplayPhysicsData(int local_meshblock_id) const {
    if (physics_reader_) {
        physics_reader_->DisplayPhysicsData(local_meshblock_id);
    } else {
        std::cout << "No physics data available" << std::endl;
    }
}

void RestartReader::ValidateUniformValues() const {
    if (physics_reader_) {
        physics_reader_->ValidateUniformValues();
    } else {
        std::cout << "No physics data available for validation" << std::endl;
    }
}