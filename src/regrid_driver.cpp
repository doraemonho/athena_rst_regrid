//========================================================================================
// AthenaK Regridding Tool - RegridDriver Implementation  
// Licensed under the 3-clause BSD License (the "LICENSE")
//========================================================================================
//! \file regrid_driver.cpp
//  \brief Implementation of main driver for AthenaK regridding tool

#include "regrid_driver.hpp"
#include <iostream>
#include <iomanip>
#include <cmath>

//----------------------------------------------------------------------------------------
//! \fn bool RegridDriver::RunRegridding()
//! \brief Main function that runs the complete regridding process
bool RegridDriver::RunRegridding(const RegridOptions& options) {
  SetError("");
  options_ = options;
  
  // Validate options
  if (!ValidateOptions()) {
    return false;
  }
  
  // Configure components
  ConfigureComponents();
  
  std::cout << "=== AthenaK MHD Restart File Regridding Tool ===" << std::endl;
  std::cout << "Input file: " << options_.input_file << std::endl;
  std::cout << "Output file: " << options_.output_file << std::endl;
  std::cout << "Refinement factor: " << options_.refinement_factor << "x" << std::endl;
  std::cout << "=================================================" << std::endl;
  
  auto total_start = std::chrono::high_resolution_clock::now();
  
  // Step 1: Read input file
  if (options_.show_progress) std::cout << "\n[1/5] Reading input restart file..." << std::endl;
  if (!ReadInputFile(options_.input_file)) {
    return false;
  }
  
  // Step 2: Regrid mesh structure
  if (options_.show_progress) std::cout << "\n[2/5] Regridding mesh structure..." << std::endl;
  if (!PerformRegridding()) {
    return false;
  }
  
  // Step 3: Write output file
  if (options_.show_progress) std::cout << "\n[3/5] Writing regridded restart file..." << std::endl;
  if (!WriteOutputFile(options_.output_file)) {
    return false;
  }
  
  // Step 4: Run verification
  if (options_.verify_conservation || options_.verify_divergence_b) {
    if (options_.show_progress) std::cout << "\n[4/5] Running verification checks..." << std::endl;
    if (!RunVerification()) {
      std::cout << "Warning: Verification failed, but output file was created" << std::endl;
    }
  }
  
  // Step 5: Print summary
  if (options_.show_progress) std::cout << "\n[5/5] Generating summary..." << std::endl;
  
  auto total_end = std::chrono::high_resolution_clock::now();
  total_time_ = std::chrono::duration_cast<std::chrono::milliseconds>(total_end - total_start).count() / 1000.0;
  
  PrintSummary();
  
  if (options_.benchmark) {
    PrintTimingInfo();
  }
  
  std::cout << "\n=== Regridding completed successfully! ===" << std::endl;
  return true;
}

//----------------------------------------------------------------------------------------
//! \fn bool RegridDriver::ReadInputFile()
//! \brief Read and validate input restart file
bool RegridDriver::ReadInputFile(const std::string& filename) {
  StartTimer();
  
  if (!reader_.ReadRestartFile(filename, input_data_)) {
    AppendError("Reader", reader_.GetLastError());
    return false;
  }
  
  // Validate this is an MHD file
  if (!reader_.ValidateMHDFile(input_data_)) {
    SetError("Input file is not a valid MHD restart file");
    return false;
  }
  
  // Print file information if verbose
  if (options_.verbose) {
    reader_.PrintFileInfo(input_data_);
  }
  
  // Update memory usage
  input_memory_usage_ = EstimateMemoryUsage(input_data_);
  UpdateMemoryUsage();
  
  read_time_ = StopTimer();
  
  if (options_.show_progress) {
    std::cout << "Successfully read " << input_data_.nmb_total 
              << " MeshBlocks (" << std::fixed << std::setprecision(2) 
              << read_time_ << "s)" << std::endl;
  }
  
  return true;
}

//----------------------------------------------------------------------------------------
//! \fn bool RegridDriver::PerformRegridding()
//! \brief Perform mesh regridding and data interpolation
bool RegridDriver::PerformRegridding() {
  // Step 2a: Regrid mesh structure
  StartTimer();
  
  if (!mesh_regridder_.RegridMesh(input_data_, mesh_data_, options_.refinement_factor)) {
    AppendError("MeshRegridder", mesh_regridder_.GetLastError());
    return false;
  }
  
  if (options_.verbose) {
    mesh_regridder_.PrintRegridInfo(input_data_, mesh_data_);
  }
  
  mesh_time_ = StopTimer();
  
  if (options_.show_progress) {
    std::cout << "Mesh regridding: " << input_data_.nmb_total << " → " 
              << mesh_data_.nmb_total_new << " MeshBlocks (" 
              << std::fixed << std::setprecision(2) << mesh_time_ << "s)" << std::endl;
  }
  
  // Step 2b: Interpolate physics data
  StartTimer();
  
  if (!data_interpolator_.InterpolateAllData(input_data_, mesh_data_, physics_data_)) {
    AppendError("DataInterpolator", data_interpolator_.GetLastError());
    return false;
  }
  
  // Update memory usage
  output_memory_usage_ = EstimateMemoryUsage(physics_data_);
  UpdateMemoryUsage();
  
  interpolation_time_ = StopTimer();
  
  if (options_.show_progress) {
    std::cout << "Data interpolation completed (" 
              << std::fixed << std::setprecision(2) << interpolation_time_ << "s)" << std::endl;
  }
  
  return true;
}

//----------------------------------------------------------------------------------------
//! \fn bool RegridDriver::WriteOutputFile()
//! \brief Write regridded data to output file
bool RegridDriver::WriteOutputFile(const std::string& filename) {
  StartTimer();
  
  // Configure writer
  if (options_.output_file_number > 0) {
    writer_.SetOutputFileNumber(options_.output_file_number);
  }
  
  if (!options_.preserve_time) {
    writer_.SetUpdateTime(true);
    writer_.SetNewTime(options_.new_time);
  }
  
  if (!writer_.WriteRestartFile(filename, input_data_, mesh_data_, physics_data_)) {
    AppendError("RestartWriter", writer_.GetLastError());
    return false;
  }
  
  write_time_ = StopTimer();
  
  if (options_.show_progress) {
    std::cout << "Output file written (" 
              << std::fixed << std::setprecision(2) << write_time_ << "s)" << std::endl;
  }
  
  return true;
}

//----------------------------------------------------------------------------------------
//! \fn bool RegridDriver::RunVerification()
//! \brief Run verification checks on the regridded data
bool RegridDriver::RunVerification() {
  StartTimer();
  
  bool verification_passed = true;
  
  // Check conservation laws
  if (options_.verify_conservation) {
    if (!data_interpolator_.VerifyConservation(input_data_, physics_data_)) {
      AppendError("Conservation", data_interpolator_.GetLastError());
      verification_passed = false;
    }
  }
  
  // Check divergence-free constraint for magnetic field
  if (options_.verify_divergence_b) {
    if (!data_interpolator_.VerifyDivergenceB(physics_data_)) {
      AppendError("Divergence", data_interpolator_.GetLastError());
      verification_passed = false;
    }
  }
  
  verification_time_ = StopTimer();
  
  if (options_.show_progress) {
    std::cout << "Verification " << (verification_passed ? "passed" : "failed") 
              << " (" << std::fixed << std::setprecision(2) << verification_time_ << "s)" 
              << std::endl;
  }
  
  return verification_passed;
}

//----------------------------------------------------------------------------------------
//! \fn void RegridDriver::PrintSummary()
//! \brief Print summary of regridding operation
void RegridDriver::PrintSummary() const {
  std::cout << "\n=== Regridding Summary ===" << std::endl;
  std::cout << "Original resolution: " << input_data_.mesh_indcs.nx1 << "³" << std::endl;
  std::cout << "New resolution: " << mesh_data_.mesh_indcs_new.nx1 << "³" << std::endl;
  std::cout << "Refinement factor: " << options_.refinement_factor << "x per dimension" << std::endl;
  std::cout << "Original MeshBlocks: " << input_data_.nmb_total << std::endl;
  std::cout << "New MeshBlocks: " << mesh_data_.nmb_total_new << std::endl;
  std::cout << "MeshBlock size: " << mesh_data_.mb_indcs_new.nx1 << "³ (unchanged)" << std::endl;
  
  // Calculate resolution increase
  long long old_cells = static_cast<long long>(input_data_.mesh_indcs.nx1) * 
                       input_data_.mesh_indcs.nx2 * input_data_.mesh_indcs.nx3;
  long long new_cells = static_cast<long long>(mesh_data_.mesh_indcs_new.nx1) * 
                       mesh_data_.mesh_indcs_new.nx2 * mesh_data_.mesh_indcs_new.nx3;
  
  std::cout << "Total cells: " << old_cells << " → " << new_cells 
            << " (" << new_cells/old_cells << "x increase)" << std::endl;
  
  // Memory usage
  std::cout << "Peak memory usage: " << peak_memory_usage_ / (1024*1024) << " MB" << std::endl;
  
  std::cout << "Total time: " << std::fixed << std::setprecision(2) 
            << total_time_ << "s" << std::endl;
  std::cout << "=========================" << std::endl;
}

//----------------------------------------------------------------------------------------
//! \fn void RegridDriver::PrintTimingInfo()
//! \brief Print detailed timing information
void RegridDriver::PrintTimingInfo() const {
  std::cout << "\n=== Timing Breakdown ===" << std::endl;
  std::cout << std::fixed << std::setprecision(3);
  std::cout << "Reading input file:    " << std::setw(8) << read_time_ << "s ("
            << std::setw(5) << std::setprecision(1) << (read_time_/total_time_*100) << "%)" << std::endl;
  std::cout << "Mesh regridding:       " << std::setw(8) << mesh_time_ << "s ("
            << std::setw(5) << (mesh_time_/total_time_*100) << "%)" << std::endl;
  std::cout << "Data interpolation:    " << std::setw(8) << interpolation_time_ << "s ("
            << std::setw(5) << (interpolation_time_/total_time_*100) << "%)" << std::endl;
  std::cout << "Writing output file:   " << std::setw(8) << write_time_ << "s ("
            << std::setw(5) << (write_time_/total_time_*100) << "%)" << std::endl;
  std::cout << "Verification:          " << std::setw(8) << verification_time_ << "s ("
            << std::setw(5) << (verification_time_/total_time_*100) << "%)" << std::endl;
  std::cout << "------------------------" << std::endl;
  std::cout << "Total:                 " << std::setw(8) << total_time_ << "s" << std::endl;
  std::cout << "========================" << std::endl;
}

//----------------------------------------------------------------------------------------
//! \fn bool RegridDriver::ValidateOptions()
//! \brief Validate regridding options
bool RegridDriver::ValidateOptions() const {
  if (options_.input_file.empty()) {
    SetError("Input file not specified");
    return false;
  }
  
  if (options_.output_file.empty()) {
    SetError("Output file not specified");
    return false;
  }
  
  if (options_.refinement_factor < 1 || options_.refinement_factor > 4) {
    SetError("Refinement factor must be between 1 and 4");
    return false;
  }
  
  if (options_.input_file == options_.output_file) {
    SetError("Input and output files cannot be the same");
    return false;
  }
  
  return true;
}

//----------------------------------------------------------------------------------------
//! \fn void RegridDriver::ConfigureComponents()
//! \brief Configure tool components based on options
void RegridDriver::ConfigureComponents() {
  // Configure data interpolator
}

//----------------------------------------------------------------------------------------
//! \fn void RegridDriver::UpdateMemoryUsage()
//! \brief Update peak memory usage tracking
void RegridDriver::UpdateMemoryUsage() {
  size_t current_usage = input_memory_usage_ + output_memory_usage_;
  peak_memory_usage_ = std::max(peak_memory_usage_, current_usage);
}

//----------------------------------------------------------------------------------------
//! \fn size_t RegridDriver::EstimateMemoryUsage()
//! \brief Estimate memory usage for RestartData
size_t RegridDriver::EstimateMemoryUsage(const RestartData& data) const {
  size_t usage = 0;
  
  usage += data.hydro_data.size() * sizeof(Real);
  usage += data.mhd_data.size() * sizeof(Real);
  usage += data.mhd_b1f_data.size() * sizeof(Real);
  usage += data.mhd_b2f_data.size() * sizeof(Real);
  usage += data.mhd_b3f_data.size() * sizeof(Real);
  usage += data.rad_data.size() * sizeof(Real);
  usage += data.force_data.size() * sizeof(Real);
  usage += data.z4c_data.size() * sizeof(Real);
  usage += data.adm_data.size() * sizeof(Real);
  
  // Add overhead for other data structures
  usage += data.lloc_eachmb.size() * sizeof(LogicalLocation);
  usage += data.cost_eachmb.size() * sizeof(float);
  usage += data.input_params.size();
  
  return usage;
}

//----------------------------------------------------------------------------------------
//! \fn size_t RegridDriver::EstimateMemoryUsage()
//! \brief Estimate memory usage for InterpolatedData  
size_t RegridDriver::EstimateMemoryUsage(const InterpolatedData& data) const {
  size_t usage = 0;
  
  usage += data.hydro_data.size() * sizeof(Real);
  usage += data.mhd_data.size() * sizeof(Real);
  usage += data.mhd_b1f_data.size() * sizeof(Real);
  usage += data.mhd_b2f_data.size() * sizeof(Real);
  usage += data.mhd_b3f_data.size() * sizeof(Real);
  usage += data.rad_data.size() * sizeof(Real);
  usage += data.force_data.size() * sizeof(Real);
  usage += data.z4c_data.size() * sizeof(Real);
  usage += data.adm_data.size() * sizeof(Real);
  
  return usage;
}