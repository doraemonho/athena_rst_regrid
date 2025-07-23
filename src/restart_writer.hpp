#ifndef RESTART_WRITER_HPP_
#define RESTART_WRITER_HPP_
//========================================================================================
// AthenaK Regridding Tool - RestartWriter Class
// Licensed under the 3-clause BSD License (the "LICENSE")
//========================================================================================
//! \file restart_writer.hpp
//  \brief Class for writing regridded data to AthenaK restart files

#include <string>
#include <sstream>

// Standalone includes
#include "restart_reader.hpp"
#include "mesh_regrid.hpp"
#include "data_interpolator.hpp"
#include "io_wrapper.hpp"

//----------------------------------------------------------------------------------------
//! \class RestartWriter
//! \brief Handles writing regridded data to AthenaK restart files
class RestartWriter {
 public:
  RestartWriter() = default;
  ~RestartWriter() = default;
  
  // Main interface functions
  bool WriteRestartFile(const std::string& output_filename,
                       const RestartData& input_data,
                       const NewMeshData& mesh_data,
                       const InterpolatedData& physics_data);
  
  // Configuration options
  void SetOutputFileNumber(int file_num) { output_file_number_ = file_num; }
  void SetUpdateTime(bool update) { update_time_ = update; }
  void SetNewTime(Real new_time) { new_time_ = new_time; }
  
  // Error handling
  const std::string& GetLastError() const { return last_error_; }
  
 private:
  // Write file sections (matching AthenaK format)
  bool WriteHeader(IOWrapper& file, const RestartData& input_data, 
                  const NewMeshData& mesh_data);
  bool WriteMeshStructure(IOWrapper& file, const NewMeshData& mesh_data);
  bool WriteInternalState(IOWrapper& file, const RestartData& input_data);
  bool WritePhysicsData(IOWrapper& file, const NewMeshData& mesh_data,
                       const InterpolatedData& physics_data);
  
  // Header generation
  std::string GenerateUpdatedInputParameters(const RestartData& input_data,
                                            const NewMeshData& mesh_data);
  bool UpdateParameterValue(std::string& params, const std::string& section,
                           const std::string& key, const std::string& value);
  
  // Data size calculations
  IOWrapperSizeT CalculateDataSize(const NewMeshData& mesh_data,
                                  const InterpolatedData& physics_data) const;
  IOWrapperSizeT CalculateOffsetSize(const RestartData& input_data) const;
  
  // Validation
  bool ValidateOutputData(const NewMeshData& mesh_data,
                         const InterpolatedData& physics_data) const;
  bool ValidateFileStructure(const std::string& filename) const;
  
  // Error reporting
  void SetError(const std::string& error) const { last_error_ = error; }
  mutable std::string last_error_;
  
  // Configuration
  int output_file_number_ = 0;
  bool update_time_ = false;
  Real new_time_ = 0.0;
};

#endif // RESTART_WRITER_HPP_