#ifndef RESTART_WRITER_HPP_
#define RESTART_WRITER_HPP_

#include <string>
#include "common.hpp"
#include "io_wrapper.hpp"
#include "restart_reader.hpp"

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
    int GetFineNMBTotal() const;
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
};

#endif // RESTART_WRITER_HPP_