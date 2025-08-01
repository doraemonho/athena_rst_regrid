#ifndef PHYSICS_READER_HPP_
#define PHYSICS_READER_HPP_

#include <vector>
#include "common.hpp"
#include "io_wrapper.hpp"
#include "mpi_distribution.hpp"

//----------------------------------------------------------------------------------------
// Physics data reading functionality
//----------------------------------------------------------------------------------------
class PhysicsReader {
public:
    PhysicsReader(IOWrapper& file, const PhysicsConfig& config, 
                  const MPIDistribution& mpi_dist, int my_rank,
                  int nout1, int nout2, int nout3);
    
    bool ReadPhysicsData(IOWrapperSizeT physics_start, IOWrapperSizeT data_size);
    
    // Data accessors
    const std::vector<std::vector<Real>>& GetHydroData() const { return hydro_data_; }
    const std::vector<std::vector<Real>>& GetMHDData() const { return mhd_data_; }
    const std::vector<std::vector<Real>>& GetMHDB1f() const { return mhd_b1f_; }
    const std::vector<std::vector<Real>>& GetMHDB2f() const { return mhd_b2f_; }
    const std::vector<std::vector<Real>>& GetMHDB3f() const { return mhd_b3f_; }
    const std::vector<std::vector<Real>>& GetTurbData() const { return turb_data_; }
    
    void DisplayPhysicsData(int local_meshblock_id) const;
    void ValidateUniformValues() const;
    
private:
    IOWrapper& file_;
    const PhysicsConfig& config_;
    const MPIDistribution& mpi_dist_;
    int my_rank_;
    int nout1_, nout2_, nout3_;
    int nmb_thisrank_;
    
    // Physics data storage (only for MeshBlocks owned by this rank)
    std::vector<std::vector<Real>> hydro_data_;
    std::vector<std::vector<Real>> mhd_data_;
    std::vector<std::vector<Real>> mhd_b1f_;
    std::vector<std::vector<Real>> mhd_b2f_;
    std::vector<std::vector<Real>> mhd_b3f_;
    std::vector<std::vector<Real>> turb_data_;
    
    // Helper functions
    IOWrapperSizeT CalculateDataSize() const;
    bool ReadHydroData(IOWrapperSizeT& myoffset, IOWrapperSizeT& offset_myrank, 
                       IOWrapperSizeT data_size, int noutmbs_min, int noutmbs_max);
    bool ReadMHDData(IOWrapperSizeT& myoffset, IOWrapperSizeT& offset_myrank, 
                     IOWrapperSizeT data_size, int noutmbs_min, int noutmbs_max);
    bool ReadTurbulenceData(IOWrapperSizeT& myoffset, IOWrapperSizeT& offset_myrank, 
                           IOWrapperSizeT data_size, int noutmbs_min, int noutmbs_max);
    void SkipRadiationData(IOWrapperSizeT& offset_myrank) const;
    void SkipZ4cADMData(IOWrapperSizeT& offset_myrank) const;
};

#endif // PHYSICS_READER_HPP_