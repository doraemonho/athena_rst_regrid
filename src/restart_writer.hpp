#ifndef RESTART_WRITER_HPP_
#define RESTART_WRITER_HPP_

#include <string>
#include <memory>
#include "common.hpp"
#include "io_wrapper.hpp"
#include "restart_reader.hpp"
#include "mpi_distribution.hpp"

// Forward declaration
class Upscaler;

//----------------------------------------------------------------------------------------
// RestartWriter class - writes restart files in AthenaK format
//----------------------------------------------------------------------------------------
class RestartWriter {
public:
    RestartWriter(RestartReader& reader, Upscaler& upscaler);
    ~RestartWriter() = default;
    
    bool WriteRestartFile(const std::string& filename);
    
private:
    RestartReader& reader_;
    Upscaler& upscaler_;
    IOWrapper file_;
    IOWrapperSizeT header_offset_;
    
    // MPI-related members
    int my_rank_;
    int nranks_;
    std::unique_ptr<MPIDistribution> mpi_dist_;
    int nmb_thisrank_;  // number of MeshBlocks assigned to this rank
    int noutmbs_min_;   // min number of MeshBlocks across all ranks
    int noutmbs_max_;   // max number of MeshBlocks across all ranks
    
    // Private methods
    bool WriteParameterData();
    bool WriteHeaderData();
    bool WriteMeshBlockLists();
    bool WriteInternalState();
    bool WritePhysicsData();
    
    // Get upscaled mesh configuration from upscaler
    RegionSize GetFineMeshSize() const;
    RegionIndcs GetFineMeshIndcs() const;
    RegionIndcs GetFineMBIndcs() const;
    size_t GetFineNMBTotal() const;
    const std::vector<LogicalLocation>& GetFineLlocEachMB() const;
    const std::vector<float>& GetFineCostEachMB() const;
    
    // Get upscaled data arrays from upscaler
    const Real* GetFineHydroData() const;
    const Real* GetFineMHDData() const;
    const Real* GetFineX1FData() const;
    const Real* GetFineX2FData() const;
    const Real* GetFineX3FData() const;
    const Real* GetFineRadData() const;
    const Real* GetFineTurbData() const;
    const Real* GetFineZ4cData() const;
    const Real* GetFineADMData() const;
    
    // Get MPI distribution info for fine mesh
    const std::vector<int>& GetFineGidsEachRank() const;
    const std::vector<int>& GetFineNmbEachRank() const;
};

#endif // RESTART_WRITER_HPP_