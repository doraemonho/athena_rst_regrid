#ifndef RESTART_READER_HPP_
#define RESTART_READER_HPP_

#include <iostream>
#include <string>
#include <vector>
#include <memory>

#include "common.hpp"
#include "io_wrapper.hpp"
#include "parameter_parser.hpp"
#include "mpi_distribution.hpp"
#include "physics_reader.hpp"

//----------------------------------------------------------------------------------------
// Main restart reader class (cleaned up and modular)
//----------------------------------------------------------------------------------------
class RestartReader {
public:
    RestartReader();
    ~RestartReader();
    
    bool ReadRestartFile(const std::string& filename, int my_rank = 0, int nranks = 1, 
                         const std::vector<int>* gids_eachrank = nullptr,
                         const std::vector<int>* nmb_eachrank = nullptr);
    
    void DisplayMeshInfo() const;
    void DisplayPhysicsData(int local_meshblock_id = 0) const;
    void ValidateUniformValues() const;
    
    // Accessors for MPI information
    int GetNMBThisRank() const { return mpi_dist_ ? mpi_dist_->GetNumMeshBlocks(my_rank_) : 0; }
    int GetMyRank() const { return my_rank_; }
    int GetNRanks() const { return nranks_; }
    
    // Accessors for mesh information
    const RegionIndcs& GetMeshIndcs() const { return mesh_indcs_; }
    const RegionIndcs& GetMBIndcs() const { return mb_indcs_; }
    const RegionSize& GetMeshSize() const { return mesh_size_; }
    int GetNMBTotal() const { return nmb_total_; }
    int GetRootLevel() const { return root_level_; }
    Real GetTime() const { return time_; }
    Real GetDt() const { return dt_; }
    int GetNCycle() const { return ncycle_; }
    const std::vector<LogicalLocation>& GetLlocEachMB() const { return lloc_eachmb_; }
    const std::vector<float>& GetCostEachMB() const { return cost_eachmb_; }
    const PhysicsConfig& GetPhysicsConfig() const { return physics_config_; }
    PhysicsReader* GetPhysicsReader() const { return physics_reader_.get(); }
    const std::string& GetParameterString() const { return parameter_string_; }
    
    // Accessors for internal state data
    Real GetZ4cLastOutputTime() const { return z4c_last_output_time_; }
    const std::vector<Real>& GetTrackerData() const { return tracker_data_; }
    const RNG_State& GetRNGState() const { return rng_state_; }

private:
    // Core data
    IOWrapper file_;
    std::string filename_;
    IOWrapperSizeT header_offset_;
    IOWrapperSizeT data_size;

    
    // MPI support
    int my_rank_;
    int nranks_;
    std::unique_ptr<MPIDistribution> mpi_dist_;
    
    // Mesh information
    int nmb_total_;
    int root_level_;
    RegionSize mesh_size_;
    RegionIndcs mesh_indcs_;
    RegionIndcs mb_indcs_;
    Real time_;
    Real dt_;
    int ncycle_;
    
    // MeshBlock information
    std::vector<LogicalLocation> lloc_eachmb_;
    std::vector<float> cost_eachmb_;
    
    // Physics configuration
    PhysicsConfig physics_config_;
    
    // Grid dimensions (including ghost zones)
    int nout1_, nout2_, nout3_;
    
    // Physics data reader
    std::unique_ptr<PhysicsReader> physics_reader_;
    
    // Parameter string
    std::string parameter_string_;
    
    // Internal state data
    Real z4c_last_output_time_;
    std::vector<Real> tracker_data_;
    RNG_State rng_state_;
    
    // Private reading functions
    bool ReadParameterData();
    bool ReadHeaderData();
    bool ReadMeshBlockLists();
    bool ReadInternalState();
    bool SetupPhysicsReader();
    
    // Helper function
    void ShowParameterSample(const std::string& param_string) const;
};

#endif // RESTART_READER_HPP_