//========================================================================================
// AthenaK Regridding Tool - RestartWriter Implementation
// Licensed under the 3-clause BSD License (the "LICENSE")
//========================================================================================
//! \file restart_writer.cpp
//  \brief Implementation of RestartWriter class

#include "restart_writer.hpp"
#include <iostream>
#include <iomanip>
#include <sstream>
#include <algorithm>
#include <cstdio>
#include <cstring>
#include <sys/stat.h>

//----------------------------------------------------------------------------------------
//! \fn bool RestartWriter::WriteRestartFile()
//! \brief Main function to write regridded restart file
bool RestartWriter::WriteRestartFile(const std::string& output_filename,
                                    const RestartData& input_data,
                                    const NewMeshData& mesh_data,
                                    const InterpolatedData& physics_data) {
  SetError("");
  
  std::cout << "Writing regridded restart file: " << output_filename << std::endl;
  
  // Validate input data
  if (!ValidateOutputData(mesh_data, physics_data)) {
    return false;
  }
  
  // Create output directory if it doesn't exist
  size_t last_slash = output_filename.find_last_of("/");
  if (last_slash != std::string::npos) {
    std::string dir_name = output_filename.substr(0, last_slash);
    if (!dir_name.empty()) {
      mkdir(dir_name.c_str(), 0775);
    }
  }
  
  // Open the output file
  IOWrapper file;
  if (file.Open(output_filename.c_str(), IOWrapper::FileMode::write) != 0) {
    SetError("Failed to create output file: " + output_filename + " (" + file.GetLastError() + ")");
    return false;
  }
  
  // Write file sections in the same order as AthenaK
  if (!WriteHeader(file, input_data, mesh_data)) {
    file.Close();
    return false;
  }
  
  if (!WriteMeshStructure(file, mesh_data)) {
    file.Close();
    return false;
  }
  
  if (!WriteInternalState(file, input_data)) {
    file.Close();
    return false;
  }
  
  if (!WritePhysicsData(file, mesh_data, physics_data)) {
    file.Close();
    return false;
  }
  
  file.Close();
  
  // Validate the written file
  if (!ValidateFileStructure(output_filename)) {
    return false;
  }
  
  std::cout << "Successfully wrote regridded restart file with " 
            << mesh_data.nmb_total_new << " MeshBlocks" << std::endl;
  
  return true;
}

//----------------------------------------------------------------------------------------
//! \fn bool RestartWriter::WriteHeader()
//! \brief Write restart file header with updated parameters
bool RestartWriter::WriteHeader(IOWrapper& file, const RestartData& input_data, 
                               const NewMeshData& mesh_data) {
  // Generate updated input parameters
  std::string updated_params = GenerateUpdatedInputParameters(input_data, mesh_data);
  
  // Write parameter string
  if (file.Write_bytes(updated_params.c_str(), 1, updated_params.size()) != 
      updated_params.size()) {
    SetError("Failed to write input parameters");
    return false;
  }
  
  // Write mesh information
  
  // Write exactly like AthenaK does
  if (file.Write_any_type(&mesh_data.nmb_total_new, sizeof(int), "byte") != sizeof(int)) {
    SetError("Failed to write nmb_total");
    return false;
  }
  
  if (file.Write_any_type(&mesh_data.root_level_new, sizeof(int), "byte") != sizeof(int)) {
    SetError("Failed to write root_level");
    return false;
  }
  
  if (file.Write_any_type(&mesh_data.mesh_size_new, sizeof(RegionSize), "byte") != sizeof(RegionSize)) {
    SetError("Failed to write mesh_size");
    return false;
  }
  
  if (file.Write_any_type(&mesh_data.mesh_indcs_new, sizeof(RegionIndcs), "byte") != 
      sizeof(RegionIndcs)) {
    SetError("Failed to write mesh_indcs");
    return false;
  }
  
  if (file.Write_any_type(&mesh_data.mb_indcs_new, sizeof(RegionIndcs), "byte") != 
      sizeof(RegionIndcs)) {
    SetError("Failed to write mb_indcs");
    return false;
  }
  
  // Write time information (optionally updated)
  Real time_to_write = update_time_ ? new_time_ : input_data.time;
  if (file.Write_any_type(&time_to_write, sizeof(Real), "byte") != sizeof(Real)) {
    SetError("Failed to write time");
    return false;
  }
  
  if (file.Write_any_type(&mesh_data.dt_new, sizeof(Real), "byte") != sizeof(Real)) {
    SetError("Failed to write dt");
    return false;
  }
  
  if (file.Write_any_type(&input_data.ncycle, sizeof(int), "byte") != sizeof(int)) {
    SetError("Failed to write ncycle");
    return false;
  }
  
  return true;
}

//----------------------------------------------------------------------------------------
//! \fn bool RestartWriter::WriteMeshStructure()
//! \brief Write MeshBlock structure information
bool RestartWriter::WriteMeshStructure(IOWrapper& file, const NewMeshData& mesh_data) {
  // Write logical locations
  size_t lloc_bytes = mesh_data.nmb_total_new * sizeof(LogicalLocation);
            
  if (file.Write_any_type(&mesh_data.lloc_eachmb_new[0], lloc_bytes, "byte") != lloc_bytes) {
    SetError("Failed to write logical locations");
    return false;
  }
  
  // Write costs
  size_t cost_bytes = mesh_data.nmb_total_new * sizeof(float);
  if (file.Write_any_type(&mesh_data.cost_eachmb_new[0], cost_bytes, "byte") != cost_bytes) {
    SetError("Failed to write MeshBlock costs");
    return false;
  }
  
  return true;
}

//----------------------------------------------------------------------------------------
//! \fn bool RestartWriter::WriteInternalState()
//! \brief Write internal state of physics modules
bool RestartWriter::WriteInternalState(IOWrapper& file, const RestartData& input_data) {
  // Write Z4c information if present
  if (input_data.has_z4c) {
    if (file.Write_any_type(&input_data.z4c_last_output_time, sizeof(Real), "byte") != 
        sizeof(Real)) {
      SetError("Failed to write Z4c last output time");
      return false;
    }
  }
  
  // Write puncture positions if present
  if (input_data.has_punctures) {
    size_t punct_bytes = input_data.puncture_positions.size() * sizeof(Real);
    if (file.Write_any_type(input_data.puncture_positions.data(), punct_bytes, "byte") != 
        punct_bytes) {
      SetError("Failed to write puncture positions");
      return false;
    }
  }
  
  // Write turbulence RNG state if present (exactly like AthenaK does)
  if (input_data.has_turb) {
    if (file.Write_any_type(&input_data.turb_rng_state, sizeof(RNG_State), "byte") != 
        sizeof(RNG_State)) {
      SetError("Failed to write turbulence RNG state");
      return false;
    }
  }
  
  return true;
}

//----------------------------------------------------------------------------------------
//! \fn bool RestartWriter::WritePhysicsData()
//! \brief Write physics data exactly like AthenaK: one MeshBlock at a time with interleaved format
bool RestartWriter::WritePhysicsData(IOWrapper& file, const NewMeshData& mesh_data,
                                    const InterpolatedData& physics_data) {
  // Calculate total data size per MeshBlock (exactly like AthenaK does it)
  IOWrapperSizeT data_size = 0;
  int nout1 = physics_data.nout1;
  int nout2 = physics_data.nout2; 
  int nout3 = physics_data.nout3;
  
  if (physics_data.nhydro > 0) {
    data_size += nout1*nout2*nout3*physics_data.nhydro*sizeof(Real); // hydro u0
  }
  if (physics_data.nmhd > 0) {
    data_size += nout1*nout2*nout3*physics_data.nmhd*sizeof(Real);   // mhd u0
    data_size += (nout1+1)*nout2*nout3*sizeof(Real);                // mhd b0.x1f
    data_size += nout1*(nout2+1)*nout3*sizeof(Real);                // mhd b0.x2f  
    data_size += nout1*nout2*(nout3+1)*sizeof(Real);                // mhd b0.x3f
  }
  if (physics_data.nrad > 0) {
    data_size += nout1*nout2*nout3*physics_data.nrad*sizeof(Real);   // radiation i0
  }
  if (!physics_data.force_data.empty()) {
    data_size += nout1*nout2*nout3*3*sizeof(Real);                  // forcing
  }
  if (physics_data.nz4c > 0) {
    data_size += nout1*nout2*nout3*physics_data.nz4c*sizeof(Real);  // z4c u0
  }
  if (physics_data.nadm > 0) {
    data_size += nout1*nout2*nout3*physics_data.nadm*sizeof(Real);  // adm u_adm
  }
  
  // Write the data size (this is written once, like AthenaK does)
  if (file.Write_any_type(&data_size, sizeof(IOWrapperSizeT), "byte") != sizeof(IOWrapperSizeT)) {
    SetError("Failed to write data size");
    return false;
  }
  
  // Write data in AthenaK's EXACT MeshBlock-interleaved format matching restart.cpp:284-570
  // Each MeshBlock gets ALL its physics data written before moving to next MeshBlock
  
  // MeshBlock-interleaved format: MB0[hydro+mhd+b1f+b2f+b3f+turb], MB1[...], MB2[...]
  for (int mb = 0; mb < physics_data.nmb_total; mb++) {
    
    // 1. Write hydro data for this MeshBlock (if present)
    if (physics_data.nhydro > 0) {
      size_t mb_start = mb * physics_data.nhydro * nout1 * nout2 * nout3;
      size_t mb_bytes = physics_data.nhydro * nout1 * nout2 * nout3 * sizeof(Real);
      
      if (file.Write_bytes(physics_data.hydro_data.data() + mb_start, 1, mb_bytes) != mb_bytes) {
        SetError("Failed to write hydro data for MeshBlock " + std::to_string(mb));
        return false;
      }
    }
    
    // 2. Write MHD conserved data for this MeshBlock (if present)
    if (physics_data.nmhd > 0) {
      size_t mb_start = mb * physics_data.nmhd * nout1 * nout2 * nout3;
      size_t mb_bytes = physics_data.nmhd * nout1 * nout2 * nout3 * sizeof(Real);
      
      if (file.Write_bytes(physics_data.mhd_data.data() + mb_start, 1, mb_bytes) != mb_bytes) {
        SetError("Failed to write MHD data for MeshBlock " + std::to_string(mb));
        return false;
      }
      
      // 3. Write B-field components for this MeshBlock: B1f, then B2f, then B3f
      
      // B1f for this MeshBlock
      size_t b1f_start = mb * (nout1+1) * nout2 * nout3;
      size_t b1f_bytes = (nout1+1) * nout2 * nout3 * sizeof(Real);
      
      if (file.Write_bytes(physics_data.mhd_b1f_data.data() + b1f_start, 1, b1f_bytes) != b1f_bytes) {
        SetError("Failed to write B1f data for MeshBlock " + std::to_string(mb));
        return false;
      }
      
      // B2f for this MeshBlock
      size_t b2f_start = mb * nout1 * (nout2+1) * nout3;
      size_t b2f_bytes = nout1 * (nout2+1) * nout3 * sizeof(Real);
      
      if (file.Write_bytes(physics_data.mhd_b2f_data.data() + b2f_start, 1, b2f_bytes) != b2f_bytes) {
        SetError("Failed to write B2f data for MeshBlock " + std::to_string(mb));
        return false;
      }
      
      // B3f for this MeshBlock
      size_t b3f_start = mb * nout1 * nout2 * (nout3+1);
      size_t b3f_bytes = nout1 * nout2 * (nout3+1) * sizeof(Real);
      
      if (file.Write_bytes(physics_data.mhd_b3f_data.data() + b3f_start, 1, b3f_bytes) != b3f_bytes) {
        SetError("Failed to write B3f data for MeshBlock " + std::to_string(mb));
        return false;
      }
    }
    
    // 4. Write radiation data for this MeshBlock (if present)
    if (physics_data.nrad > 0) {
      size_t mb_start = mb * physics_data.nrad * nout1 * nout2 * nout3;
      size_t mb_bytes = physics_data.nrad * nout1 * nout2 * nout3 * sizeof(Real);
      
      if (file.Write_bytes(physics_data.rad_data.data() + mb_start, 1, mb_bytes) != mb_bytes) {
        SetError("Failed to write radiation data for MeshBlock " + std::to_string(mb));
        return false;
      }
    }
    
    // 5. Write turbulence forcing data for this MeshBlock (if present)
    if (!physics_data.force_data.empty()) {
      size_t mb_start = mb * 3 * nout1 * nout2 * nout3;
      size_t mb_bytes = 3 * nout1 * nout2 * nout3 * sizeof(Real);
      
      if (file.Write_bytes(physics_data.force_data.data() + mb_start, 1, mb_bytes) != mb_bytes) {
        SetError("Failed to write turbulence data for MeshBlock " + std::to_string(mb));
        return false;
      }
    }
    
    // 6. Write Z4c data for this MeshBlock (if present)
    if (physics_data.nz4c > 0) {
      size_t mb_start = mb * physics_data.nz4c * nout1 * nout2 * nout3;
      size_t mb_bytes = physics_data.nz4c * nout1 * nout2 * nout3 * sizeof(Real);
      
      if (file.Write_bytes(physics_data.z4c_data.data() + mb_start, 1, mb_bytes) != mb_bytes) {
        SetError("Failed to write Z4c data for MeshBlock " + std::to_string(mb));
        return false;
      }
    }
    
    // 7. Write ADM data for this MeshBlock (if present)
    if (physics_data.nadm > 0) {
      size_t mb_start = mb * physics_data.nadm * nout1 * nout2 * nout3;
      size_t mb_bytes = physics_data.nadm * nout1 * nout2 * nout3 * sizeof(Real);
      
      if (file.Write_bytes(physics_data.adm_data.data() + mb_start, 1, mb_bytes) != mb_bytes) {
        SetError("Failed to write ADM data for MeshBlock " + std::to_string(mb));
        return false;
      }
    }
  }
  
  return true;
}

//----------------------------------------------------------------------------------------
//! \fn std::string RestartWriter::GenerateUpdatedInputParameters()
//! \brief Generate updated input parameter string with new resolution
std::string RestartWriter::GenerateUpdatedInputParameters(const RestartData& input_data,
                                                          const NewMeshData& mesh_data) {
  std::string updated_params = input_data.input_params;
  
  // Update mesh resolution
  UpdateParameterValue(updated_params, "mesh", "nx1", 
                      std::to_string(mesh_data.mesh_indcs_new.nx1));
  UpdateParameterValue(updated_params, "mesh", "nx2", 
                      std::to_string(mesh_data.mesh_indcs_new.nx2));
  UpdateParameterValue(updated_params, "mesh", "nx3", 
                      std::to_string(mesh_data.mesh_indcs_new.nx3));
  
  // Update file number if specified
  if (output_file_number_ > 0) {
    UpdateParameterValue(updated_params, "output1", "file_number", 
                        std::to_string(output_file_number_));
  }
  
  // Restore the PAR_DUMP footer that was stripped during reading
  updated_params += "#------------------------- PAR_DUMP -------------------------\n";
  updated_params += "<par_end>\n";
  
  return updated_params;
}

//----------------------------------------------------------------------------------------
//! \fn bool RestartWriter::UpdateParameterValue()
//! \brief Update a specific parameter value in the parameter string
bool RestartWriter::UpdateParameterValue(std::string& params, const std::string& section,
                                        const std::string& key, const std::string& value) {
  std::istringstream iss(params);
  std::ostringstream oss;
  std::string line;
  bool in_section = false;
  bool found_key = false;
  
  while (std::getline(iss, line)) {
    // Check if we're entering the target section
    if (line.find("<" + section + ">") != std::string::npos) {
      in_section = true;
      oss << line << "\n";
      continue;
    }
    
    // Check if we're leaving the section (explicit closing tag or new section)
    if (in_section && (line.find("</") != std::string::npos || 
                       (line.find("<") != std::string::npos && 
                        line.find("<" + section + ">") == std::string::npos))) {
      // If we didn't find the key, add it before leaving the section
      if (!found_key) {
        oss << key << " = " << value << "\n";
      }
      in_section = false;
      oss << line << "\n";
      continue;
    }
    
    // Check if this is the key we want to update (must be at start of line)
    if (in_section) {
      // Find the parameter name at the start of the line (ignoring whitespace)
      size_t start = line.find_first_not_of(" \t");
      if (start != std::string::npos && 
          line.substr(start, key.length()) == key &&
          (start + key.length() >= line.length() || 
           !std::isalnum(line[start + key.length()]))) {
        // Check if followed by = or whitespace
        size_t after_key = start + key.length();
        if (after_key < line.length() && 
            (line[after_key] == '=' || std::isspace(line[after_key]))) {
          // Find the equals sign and value position
          size_t eq_pos = line.find('=', after_key);
          if (eq_pos != std::string::npos) {
            // Find the end of the value (before any comment)
            size_t value_start = line.find_first_not_of(" \t", eq_pos + 1);
            size_t comment_pos = line.find('#', value_start);
            
            // Preserve original formatting: keep everything before value and after value
            std::string prefix = line.substr(0, value_start);
            std::string suffix = (comment_pos != std::string::npos) ? 
                               line.substr(comment_pos) : "";
            
            // Reconstruct line with new value and proper spacing
            if (!suffix.empty()) {
              oss << prefix << value << "       " << suffix << "\n";
            } else {
              oss << prefix << value << "\n";
            }
            found_key = true;
            continue;
          }
        }
      }
    }
    
    // Copy line as-is
    oss << line << "\n";
  }
  
  params = oss.str();
  return found_key || !in_section; // Success if found or if section didn't exist
}

//----------------------------------------------------------------------------------------
//! \fn IOWrapperSizeT RestartWriter::CalculateDataSize()
//! \brief Calculate total size of physics data to be written
IOWrapperSizeT RestartWriter::CalculateDataSize(const NewMeshData& mesh_data,
                                               const InterpolatedData& physics_data) const {
  int nout1 = physics_data.nout1;
  int nout2 = physics_data.nout2;
  int nout3 = physics_data.nout3;
  
  // Per MeshBlock data size (AthenaK format)
  IOWrapperSizeT mb_data_size = 0;
  
  if (physics_data.nhydro > 0) {
    IOWrapperSizeT hydro_size = nout1 * nout2 * nout3 * physics_data.nhydro * sizeof(Real);
    mb_data_size += hydro_size;
  }
  
  if (physics_data.nmhd > 0) {
    IOWrapperSizeT mhd_cons = nout1 * nout2 * nout3 * physics_data.nmhd * sizeof(Real);
    IOWrapperSizeT b1f_size = (nout1 + 1) * nout2 * nout3 * sizeof(Real);
    IOWrapperSizeT b2f_size = nout1 * (nout2 + 1) * nout3 * sizeof(Real);
    IOWrapperSizeT b3f_size = nout1 * nout2 * (nout3 + 1) * sizeof(Real);
    
    mb_data_size += mhd_cons + b1f_size + b2f_size + b3f_size;
  }
  
  if (physics_data.nrad > 0) {
    mb_data_size += nout1 * nout2 * nout3 * physics_data.nrad * sizeof(Real);
  }
  
  // Only include turbulence data if it was actually present in the original file
  if (!physics_data.force_data.empty()) {
    IOWrapperSizeT force_size = nout1 * nout2 * nout3 * physics_data.nforce * sizeof(Real);
    mb_data_size += force_size;
  }
  
  if (physics_data.nz4c > 0) {
    mb_data_size += nout1 * nout2 * nout3 * physics_data.nz4c * sizeof(Real);
  }
  
  if (physics_data.nadm > 0) {
    mb_data_size += nout1 * nout2 * nout3 * physics_data.nadm * sizeof(Real);
  }
  
  return mb_data_size;
}

//----------------------------------------------------------------------------------------
//! \fn bool RestartWriter::ValidateOutputData()
//! \brief Validate output data before writing
bool RestartWriter::ValidateOutputData(const NewMeshData& mesh_data,
                                      const InterpolatedData& physics_data) const {
  // Check mesh data consistency
  if (mesh_data.nmb_total_new <= 0) {
    SetError("Invalid number of output MeshBlocks");
    return false;
  }
  
  if (mesh_data.lloc_eachmb_new.size() != static_cast<size_t>(mesh_data.nmb_total_new)) {
    SetError("Inconsistent logical location array size");
    return false;
  }
  
  if (mesh_data.cost_eachmb_new.size() != static_cast<size_t>(mesh_data.nmb_total_new)) {
    SetError("Inconsistent cost array size");
    return false;
  }
  
  // Check physics data consistency
  if (physics_data.nmb_total != mesh_data.nmb_total_new) {
    SetError("Physics data MeshBlock count doesn't match mesh data");
    return false;
  }
  
  // Check array sizes for each physics module
  if (physics_data.nmhd > 0) {
    size_t expected_mhd = physics_data.nmb_total * physics_data.nmhd * 
                         physics_data.nout1 * physics_data.nout2 * physics_data.nout3;
    if (physics_data.mhd_data.size() != expected_mhd) {
      SetError("Inconsistent MHD data array size");
      return false;
    }
    
    size_t expected_b1f = physics_data.nmb_total * physics_data.nout3 * 
                         physics_data.nout2 * (physics_data.nout1 + 1);
    if (physics_data.mhd_b1f_data.size() != expected_b1f) {
      SetError("Inconsistent B1f data array size");
      return false;
    }
    
    size_t expected_b2f = physics_data.nmb_total * physics_data.nout3 * 
                         (physics_data.nout2 + 1) * physics_data.nout1;
    if (physics_data.mhd_b2f_data.size() != expected_b2f) {
      SetError("Inconsistent B2f data array size");
      return false;
    }
    
    size_t expected_b3f = physics_data.nmb_total * (physics_data.nout3 + 1) * 
                         physics_data.nout2 * physics_data.nout1;
    if (physics_data.mhd_b3f_data.size() != expected_b3f) {
      SetError("Inconsistent B3f data array size");
      return false;
    }
  }
  
  return true;
}

//----------------------------------------------------------------------------------------
//! \fn bool RestartWriter::ValidateFileStructure()
//! \brief Validate the written file structure
bool RestartWriter::ValidateFileStructure(const std::string& filename) const {
  // Try to reopen the file to check it was written correctly
  IOWrapper test_file;
  if (test_file.Open(filename.c_str(), IOWrapper::FileMode::read) != 0) {
    SetError("Cannot reopen written file for validation");
    return false;
  }
  
  // Read a few bytes to check file is accessible
  char buffer[1024];
  size_t bytes_read = test_file.Read_bytes(buffer, 1, sizeof(buffer));
  test_file.Close();
  
  if (bytes_read == 0) {
    SetError("Written file appears to be empty or corrupted");
    return false;
  }
  
  return true;
}