#ifndef REGRID_DRIVER_HPP_
#define REGRID_DRIVER_HPP_
//========================================================================================
// AthenaK Regridding Tool - RegridDriver Class
// Licensed under the 3-clause BSD License (the "LICENSE")
//========================================================================================
//! \file regrid_driver.hpp
//  \brief Main driver class for the AthenaK regridding tool

#include <string>
#include <chrono>

// Tool components
#include "restart_reader.hpp"
#include "mesh_regrid.hpp"
#include "data_interpolator.hpp"
#include "restart_writer.hpp"
#include "verification.hpp"

//----------------------------------------------------------------------------------------
//! \struct RegridOptions
//! \brief Configuration options for regridding process
struct RegridOptions {
  // Input/Output files
  std::string input_file;
  std::string output_file;
  
  // Regridding options
  int refinement_factor = 2;        // Resolution multiplier (default 2x)
  bool preserve_time = true;        // Keep original time
  
  // Verification options
  bool verify_conservation = true;   // Check conservation laws
  bool verify_divergence_b = true;   // Check ∇·B = 0
  bool verbose = false;              // Detailed output
  
  // Output options
  int output_file_number = 0;       // New file number (0 = auto)
  Real new_time = 0.0;              // New time value (if !preserve_time)
  
  // Performance options
  bool show_progress = true;        // Show progress information
  bool benchmark = false;           // Show timing information
};

//----------------------------------------------------------------------------------------
//! \class RegridDriver
//! \brief Main driver class that orchestrates the regridding process
class RegridDriver {
 public:
  RegridDriver() = default;
  ~RegridDriver() = default;
  
  // Main interface
  bool RunRegridding(const RegridOptions& options);
  
  // Individual steps (for testing/debugging)
  bool ReadInputFile(const std::string& filename);
  bool PerformRegridding();
  bool WriteOutputFile(const std::string& filename);
  bool RunVerification();
  
  // Information functions
  void PrintSummary() const;
  void PrintTimingInfo() const;
  void PrintMemoryUsage() const;
  
  // Error handling
  const std::string& GetLastError() const { return last_error_; }
  bool HasError() const { return !last_error_.empty(); }
  
 private:
  // Core components
  RestartReader reader_;
  MeshRegridder mesh_regridder_;
  DataInterpolator data_interpolator_;
  RestartWriter writer_;
  Verification verifier_;
  
  // Data containers
  RestartData input_data_;
  NewMeshData mesh_data_;
  InterpolatedData physics_data_;
  
  // Configuration
  RegridOptions options_;
  
  // Timing information
  std::chrono::high_resolution_clock::time_point start_time_;
  double read_time_ = 0.0;
  double mesh_time_ = 0.0;
  double interpolation_time_ = 0.0;
  double write_time_ = 0.0;
  double verification_time_ = 0.0;
  double total_time_ = 0.0;
  
  // Memory usage tracking
  size_t peak_memory_usage_ = 0;
  size_t input_memory_usage_ = 0;
  size_t output_memory_usage_ = 0;
  
  // Helper functions
  void StartTimer() { start_time_ = std::chrono::high_resolution_clock::now(); }
  double StopTimer() const {
    auto end_time = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(end_time - start_time_);
    return duration.count() / 1000.0;
  }
  
  void UpdateMemoryUsage();
  size_t EstimateMemoryUsage(const RestartData& data) const;
  size_t EstimateMemoryUsage(const InterpolatedData& data) const;
  
  bool ValidateOptions() const;
  void ConfigureComponents();
  
  // Error reporting
  void SetError(const std::string& error) const { last_error_ = error; }
  void AppendError(const std::string& component, const std::string& error) {
    if (!last_error_.empty()) last_error_ += "; ";
    last_error_ += component + ": " + error;
  }
  mutable std::string last_error_;
};

#endif // REGRID_DRIVER_HPP_