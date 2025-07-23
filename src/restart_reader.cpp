//========================================================================================
// AthenaK Regridding Tool - RestartReader Implementation
// Licensed under the 3-clause BSD License (the "LICENSE")
//========================================================================================
//! \file restart_reader.cpp
//  \brief Implementation of RestartReader class for reading AthenaK restart files

#include "restart_reader.hpp"
#include <iostream>
#include <fstream>
#include <sstream>
#include <cstring>
#include <algorithm>

//----------------------------------------------------------------------------------------
//! \fn bool RestartReader::ReadRestartFile()
//! \brief Main function to read an AthenaK restart file
bool RestartReader::ReadRestartFile(const std::string& filename, RestartData& data) {
  SetError("");
  
  // Open the restart file
  IOWrapper file;
  if (file.Open(filename.c_str(), IOWrapper::FileMode::read) != 0) {
    SetError("Failed to open restart file: " + filename);
    return false;
  }
  
  std::cout << "Reading restart file: " << filename << std::endl;
  
  // Read file sections in AthenaK format order
  if (!ReadHeader(file, data)) {
    file.Close();
    return false;
  }
  
  if (!ReadMeshStructure(file, data)) {
    file.Close();
    return false;
  }
  
  if (!ReadInternalState(file, data)) {
    file.Close();
    return false;
  }
  
  if (!ReadPhysicsData(file, data)) {
    file.Close();
    return false;
  }
  
  file.Close();
  
  // Validate the data we read
  if (!ValidateDataConsistency(data)) {
    return false;
  }
  
  std::cout << "Successfully read restart file with " << data.nmb_total 
            << " MeshBlocks" << std::endl;
  return true;
}

//----------------------------------------------------------------------------------------
//! \fn bool RestartReader::ReadHeader()
//! \brief Read the restart file header following AthenaK format
bool RestartReader::ReadHeader(IOWrapper& file, RestartData& data) {
  // Read input parameters as variable-length string until <par_end> marker
  std::string params_str;
  char buffer[8192];
  bool found_end = false;
  
  while (!found_end) {
    size_t bytes_read = file.Read_bytes(buffer, 1, sizeof(buffer)-1);
    if (bytes_read == 0) {
      SetError("Unexpected end of file while reading input parameters");
      return false;
    }
    
    buffer[bytes_read] = '\0';
    params_str += buffer;
    
    // Look for <par_end> marker
    size_t pos = params_str.find("<par_end>");
    if (pos != std::string::npos) {
      // Trim to end at <par_end> marker
      params_str = params_str.substr(0, pos);
      found_end = true;
      
      // Seek to position after <par_end> + newline
      IOWrapperSizeT seek_pos = pos + 10; // +9 for "<par_end>" +1 for newline
      file.Seek(seek_pos);
    }
  }
  
  data.input_params = params_str;
  
  // Parse input parameters to extract physics module information
  if (!ParseInputParameters(params_str, data)) {
    return false;
  }
  
  // Read mesh metadata in exact AthenaK order
  if (file.Read_bytes(&data.nmb_total, sizeof(int), 1) != sizeof(int)) {
    SetError("Failed to read nmb_total");
    return false;
  }
  
  if (file.Read_bytes(&data.root_level, sizeof(int), 1) != sizeof(int)) {
    SetError("Failed to read root_level");
    return false;
  }
  
  // Read RegionSize structure (domain bounds + grid spacing)
  if (file.Read_bytes(&data.mesh_size, sizeof(RegionSize), 1) != sizeof(RegionSize)) {
    SetError("Failed to read mesh_size");
    return false;
  }
  
  // Read mesh-level indices
  if (file.Read_bytes(&data.mesh_indcs, sizeof(RegionIndcs), 1) != sizeof(RegionIndcs)) {
    SetError("Failed to read mesh_indcs");
    return false;
  }
  
  // Read MeshBlock-level indices  
  if (file.Read_bytes(&data.mb_indcs, sizeof(RegionIndcs), 1) != sizeof(RegionIndcs)) {
    SetError("Failed to read mb_indcs");
    return false;
  }
  
  std::cout << "Debug: Read mesh dimensions: " << data.mesh_indcs.nx1 << " x " 
            << data.mesh_indcs.nx2 << " x " << data.mesh_indcs.nx3 << std::endl;
  std::cout << "Debug: Read MeshBlock dimensions: " << data.mb_indcs.nx1 << " x " 
            << data.mb_indcs.nx2 << " x " << data.mb_indcs.nx3 << std::endl;
  std::cout << "Debug: Total MeshBlocks: " << data.nmb_total << std::endl;
  
  // Read simulation time
  if (file.Read_bytes(&data.time, sizeof(Real), 1) != sizeof(Real)) {
    SetError("Failed to read time");
    return false;
  }
  
  // Read timestep
  if (file.Read_bytes(&data.dt, sizeof(Real), 1) != sizeof(Real)) {
    SetError("Failed to read dt");
    return false;
  }
  
  // Read cycle number
  if (file.Read_bytes(&data.ncycle, sizeof(int), 1) != sizeof(int)) {
    SetError("Failed to read ncycle");
    return false;
  }
  
  // Calculate output dimensions (including ghost zones)
  data.nout1 = data.mb_indcs.nx1 + 2*data.mb_indcs.ng;
  data.nout2 = (data.mb_indcs.nx2 > 1) ? (data.mb_indcs.nx2 + 2*data.mb_indcs.ng) : 1;
  data.nout3 = (data.mb_indcs.nx3 > 1) ? (data.mb_indcs.nx3 + 2*data.mb_indcs.ng) : 1;
  
  return ValidateHeader(data);
}

//----------------------------------------------------------------------------------------
//! \fn bool RestartReader::ReadMeshStructure()
//! \brief Read MeshBlock structure information
bool RestartReader::ReadMeshStructure(IOWrapper& file, RestartData& data) {
  // Read logical locations array
  data.lloc_eachmb.resize(data.nmb_total);
  size_t lloc_bytes = data.nmb_total * sizeof(LogicalLocation);
  if (file.Read_bytes(data.lloc_eachmb.data(), 1, lloc_bytes) != lloc_bytes) {
    SetError("Failed to read logical locations");
    return false;
  }
  
  // Read costs array
  data.cost_eachmb.resize(data.nmb_total);
  size_t cost_bytes = data.nmb_total * sizeof(float);
  if (file.Read_bytes(data.cost_eachmb.data(), 1, cost_bytes) != cost_bytes) {
    SetError("Failed to read MeshBlock costs");
    return false;
  }
  
  return ValidateMeshStructure(data);
}

//----------------------------------------------------------------------------------------
//! \fn bool RestartReader::ReadInternalState()
//! \brief Read internal state of physics modules
bool RestartReader::ReadInternalState(IOWrapper& file, RestartData& data) {
  // Read Z4c last output time if Z4c physics is enabled
  if (data.has_z4c) {
    if (file.Read_bytes(&data.z4c_last_output_time, sizeof(Real), 1) != sizeof(Real)) {
      SetError("Failed to read Z4c last output time");
      return false;
    }
  }
  
  // Read puncture tracker positions if present
  if (data.has_punctures) {
    size_t nco = data.puncture_positions.size() / 3;
    size_t punct_bytes = 3 * nco * sizeof(Real);
    if (file.Read_bytes(data.puncture_positions.data(), 1, punct_bytes) != punct_bytes) {
      SetError("Failed to read puncture positions");
      return false;
    }
  }
  
  // Read turbulence RNG state if turbulence driving is enabled
  if (data.has_turb) {
    if (file.Read_bytes(&data.turb_rng_state, sizeof(RNG_State), 1) != sizeof(RNG_State)) {
      SetError("Failed to read turbulence RNG state");
      return false;
    }
  }
  
  return true;
}

//----------------------------------------------------------------------------------------
//! \fn bool RestartReader::ReadPhysicsData()
//! \brief Read physics variable arrays in MeshBlock-interleaved format
bool RestartReader::ReadPhysicsData(IOWrapper& file, RestartData& data) {
  // Read data size header for validation
  IOWrapperSizeT data_size;
  if (file.Read_bytes(&data_size, sizeof(IOWrapperSizeT), 1) != sizeof(IOWrapperSizeT)) {
    SetError("Failed to read data size");
    return false;
  }
  
  // Calculate array sizes
  size_t total_cells = data.nmb_total * data.nout1 * data.nout2 * data.nout3;
  
  // Allocate arrays for physics data
  if (data.nhydro > 0) {
    data.hydro_data.resize(total_cells * data.nhydro);
  }
  
  if (data.nmhd > 0) {
    data.mhd_data.resize(total_cells * data.nmhd);
    data.mhd_b1f_data.resize(data.nmb_total * (data.nout1+1) * data.nout2 * data.nout3);
    data.mhd_b2f_data.resize(data.nmb_total * data.nout1 * (data.nout2+1) * data.nout3);
    data.mhd_b3f_data.resize(data.nmb_total * data.nout1 * data.nout2 * (data.nout3+1));
  }
  
  if (data.nrad > 0) {
    data.rad_data.resize(total_cells * data.nrad);
  }
  
  if (data.has_turb) {
    data.force_data.resize(total_cells * data.nforce);
  }
  
  if (data.nz4c > 0) {
    data.z4c_data.resize(total_cells * data.nz4c);
  }
  
  if (data.nadm > 0) {
    data.adm_data.resize(total_cells * data.nadm);
  }
  
  // Read physics data in MeshBlock-interleaved format
  // Each MeshBlock contains ALL its physics data before moving to next MeshBlock
  
  for (int mb = 0; mb < data.nmb_total; mb++) {
    size_t mb_cells = data.nout1 * data.nout2 * data.nout3;
    
    // Read hydro data for this MeshBlock
    if (data.nhydro > 0) {
      size_t mb_start = mb * data.nhydro * mb_cells;
      size_t mb_bytes = data.nhydro * mb_cells * sizeof(Real);
      
      if (file.Read_bytes(data.hydro_data.data() + mb_start, 1, mb_bytes) != mb_bytes) {
        SetError("Failed to read hydro data for MeshBlock " + std::to_string(mb));
        return false;
      }
    }
    
    // Read MHD conserved variables for this MeshBlock
    if (data.nmhd > 0) {
      size_t mb_start = mb * data.nmhd * mb_cells;
      size_t mb_bytes = data.nmhd * mb_cells * sizeof(Real);
      
      if (file.Read_bytes(data.mhd_data.data() + mb_start, 1, mb_bytes) != mb_bytes) {
        SetError("Failed to read MHD data for MeshBlock " + std::to_string(mb));
        return false;
      }
      
      // Read face-centered B-fields for this MeshBlock
      
      // B1f (x-faces): (nx1+1) x ny x nz
      size_t b1f_start = mb * (data.nout1+1) * data.nout2 * data.nout3;
      size_t b1f_bytes = (data.nout1+1) * data.nout2 * data.nout3 * sizeof(Real);
      
      if (file.Read_bytes(data.mhd_b1f_data.data() + b1f_start, 1, b1f_bytes) != b1f_bytes) {
        SetError("Failed to read B1f data for MeshBlock " + std::to_string(mb));
        return false;
      }
      
      // B2f (y-faces): nx x (ny+1) x nz  
      size_t b2f_start = mb * data.nout1 * (data.nout2+1) * data.nout3;
      size_t b2f_bytes = data.nout1 * (data.nout2+1) * data.nout3 * sizeof(Real);
      
      if (file.Read_bytes(data.mhd_b2f_data.data() + b2f_start, 1, b2f_bytes) != b2f_bytes) {
        SetError("Failed to read B2f data for MeshBlock " + std::to_string(mb));
        return false;
      }
      
      // B3f (z-faces): nx x ny x (nz+1)
      size_t b3f_start = mb * data.nout1 * data.nout2 * (data.nout3+1);
      size_t b3f_bytes = data.nout1 * data.nout2 * (data.nout3+1) * sizeof(Real);
      
      if (file.Read_bytes(data.mhd_b3f_data.data() + b3f_start, 1, b3f_bytes) != b3f_bytes) {
        SetError("Failed to read B3f data for MeshBlock " + std::to_string(mb));
        return false;
      }
    }
    
    // Read radiation data for this MeshBlock
    if (data.nrad > 0) {
      size_t mb_start = mb * data.nrad * mb_cells;
      size_t mb_bytes = data.nrad * mb_cells * sizeof(Real);
      
      if (file.Read_bytes(data.rad_data.data() + mb_start, 1, mb_bytes) != mb_bytes) {
        SetError("Failed to read radiation data for MeshBlock " + std::to_string(mb));
        return false;
      }
    }
    
    // Read turbulence forcing data for this MeshBlock
    if (data.has_turb) {
      size_t mb_start = mb * data.nforce * mb_cells;
      size_t mb_bytes = data.nforce * mb_cells * sizeof(Real);
      
      if (file.Read_bytes(data.force_data.data() + mb_start, 1, mb_bytes) != mb_bytes) {
        SetError("Failed to read turbulence data for MeshBlock " + std::to_string(mb));
        return false;
      }
    }
    
    // Read Z4c data for this MeshBlock
    if (data.nz4c > 0) {
      size_t mb_start = mb * data.nz4c * mb_cells;
      size_t mb_bytes = data.nz4c * mb_cells * sizeof(Real);
      
      if (file.Read_bytes(data.z4c_data.data() + mb_start, 1, mb_bytes) != mb_bytes) {
        SetError("Failed to read Z4c data for MeshBlock " + std::to_string(mb));
        return false;
      }
    }
    
    // Read ADM data for this MeshBlock  
    if (data.nadm > 0) {
      size_t mb_start = mb * data.nadm * mb_cells;
      size_t mb_bytes = data.nadm * mb_cells * sizeof(Real);
      
      if (file.Read_bytes(data.adm_data.data() + mb_start, 1, mb_bytes) != mb_bytes) {
        SetError("Failed to read ADM data for MeshBlock " + std::to_string(mb));
        return false;
      }
    }
  }
  
  return true;
}

//----------------------------------------------------------------------------------------
//! \fn bool RestartReader::ParseInputParameters()
//! \brief Parse input parameter string to extract physics module info
bool RestartReader::ParseInputParameters(const std::string& params, RestartData& data) {
  std::istringstream iss(params);
  std::string line;
  
  // Initialize physics module settings
  data.nhydro = 0;
  data.nmhd = 0;
  data.nrad = 0;
  data.nforce = 3;  // turbulence always has 3 force components
  data.nz4c = 0;
  data.nadm = 0;
  data.has_z4c = false;
  data.has_punctures = false;
  data.has_turb = false;
  
  while (std::getline(iss, line)) {
    // Trim whitespace
    line.erase(0, line.find_first_not_of(" \t"));
    line.erase(line.find_last_not_of(" \t") + 1);
    
    // Skip empty lines and comments
    if (line.empty() || line[0] == '#') continue;
    
    // Check for physics module section headers
    if (line == "<hydro>") {
      data.nhydro = 5;  // Standard ideal gas hydro: density, momentum(3), energy
    } else if (line == "<mhd>") {
      data.nmhd = 5;   // Standard ideal gas MHD: same 5 conserved variables
    } else if (line == "<radiation>") {
      data.nrad = 1;   // Default radiation transport
    } else if (line == "<z4c>") {
      data.has_z4c = true;
      data.nz4c = 25;  // Standard Z4c spacetime variables
    } else if (line == "<adm>") {
      data.nadm = 10;  // Standard ADM variables
    } else if (line == "<turb_driving>") {
      data.has_turb = true;
    }
  }
  
  return true;
}

//----------------------------------------------------------------------------------------
//! \fn void RestartReader::CalculateDataSizes()
//! \brief Calculate expected data sizes for validation
void RestartReader::CalculateDataSizes(RestartData& data) {
  // Calculate expected total data size for validation
  IOWrapperSizeT expected_size = 0;
  size_t mb_cells = data.nout1 * data.nout2 * data.nout3;
  
  if (data.nhydro > 0) {
    expected_size += data.nmb_total * data.nhydro * mb_cells * sizeof(Real);
  }
  
  if (data.nmhd > 0) {
    expected_size += data.nmb_total * data.nmhd * mb_cells * sizeof(Real);
    expected_size += data.nmb_total * (data.nout1+1) * data.nout2 * data.nout3 * sizeof(Real);
    expected_size += data.nmb_total * data.nout1 * (data.nout2+1) * data.nout3 * sizeof(Real);
    expected_size += data.nmb_total * data.nout1 * data.nout2 * (data.nout3+1) * sizeof(Real);
  }
  
  if (data.nrad > 0) {
    expected_size += data.nmb_total * data.nrad * mb_cells * sizeof(Real);
  }
  
  if (data.has_turb) {
    expected_size += data.nmb_total * data.nforce * mb_cells * sizeof(Real);
  }
  
  if (data.nz4c > 0) {
    expected_size += data.nmb_total * data.nz4c * mb_cells * sizeof(Real);
  }
  
  if (data.nadm > 0) {
    expected_size += data.nmb_total * data.nadm * mb_cells * sizeof(Real);
  }
}

//----------------------------------------------------------------------------------------
//! \fn bool RestartReader::ValidateMHDFile()
//! \brief Validate that this is a proper MHD restart file
bool RestartReader::ValidateMHDFile(const RestartData& data) const {
  if (data.nmhd == 0) {
    return false;  // Not an MHD file
  }
  
  if (data.mhd_data.empty() || data.mhd_b1f_data.empty() || 
      data.mhd_b2f_data.empty() || data.mhd_b3f_data.empty()) {
    return false;  // Missing MHD data
  }
  
  return true;
}

//----------------------------------------------------------------------------------------
//! \fn bool RestartReader::ValidateHeader()
//! \brief Validate header data consistency
bool RestartReader::ValidateHeader(const RestartData& data) const {
  if (data.nmb_total <= 0) {
    SetError("Invalid number of MeshBlocks");
    return false;
  }
  
  if (data.mb_indcs.nx1 <= 0 || data.mb_indcs.ng < 2) {
    SetError("Invalid MeshBlock dimensions");
    return false;
  }
  
  if (data.nout1 <= 0 || data.nout2 <= 0 || data.nout3 <= 0) {
    SetError("Invalid output dimensions");
    return false;
  }
  
  return true;
}

//----------------------------------------------------------------------------------------
//! \fn bool RestartReader::ValidateMeshStructure()
//! \brief Validate mesh structure data
bool RestartReader::ValidateMeshStructure(const RestartData& data) const {
  if (data.lloc_eachmb.size() != static_cast<size_t>(data.nmb_total)) {
    SetError("Inconsistent logical location array size");
    return false;
  }
  
  if (data.cost_eachmb.size() != static_cast<size_t>(data.nmb_total)) {
    SetError("Inconsistent cost array size");
    return false;
  }
  
  return true;
}

//----------------------------------------------------------------------------------------
//! \fn bool RestartReader::ValidateDataConsistency()
//! \brief Validate overall data consistency
bool RestartReader::ValidateDataConsistency(const RestartData& data) const {
  // Require MHD data for this tool
  if (data.nmhd == 0) {
    SetError("This tool requires MHD data");
    return false;
  }
  
  // Check array sizes are consistent with expected dimensions
  if (data.nmhd > 0) {
    size_t expected_mhd = data.nmb_total * data.nmhd * data.nout1 * data.nout2 * data.nout3;
    if (data.mhd_data.size() != expected_mhd) {
      SetError("Inconsistent MHD data array size");
      return false;
    }
    
    size_t expected_b1f = data.nmb_total * (data.nout1+1) * data.nout2 * data.nout3;
    if (data.mhd_b1f_data.size() != expected_b1f) {
      SetError("Inconsistent B1f data array size");
      return false;
    }
    
    size_t expected_b2f = data.nmb_total * data.nout1 * (data.nout2+1) * data.nout3;
    if (data.mhd_b2f_data.size() != expected_b2f) {
      SetError("Inconsistent B2f data array size");
      return false;
    }
    
    size_t expected_b3f = data.nmb_total * data.nout1 * data.nout2 * (data.nout3+1);
    if (data.mhd_b3f_data.size() != expected_b3f) {
      SetError("Inconsistent B3f data array size");
      return false;
    }
  }
  
  return true;
}

//----------------------------------------------------------------------------------------
//! \fn void RestartReader::PrintFileInfo()
//! \brief Print information about the restart file
void RestartReader::PrintFileInfo(const RestartData& data) const {
  std::cout << "\n=== Restart File Information ===" << std::endl;
  std::cout << "Number of MeshBlocks: " << data.nmb_total << std::endl;
  std::cout << "Root level: " << data.root_level << std::endl;
  std::cout << "Time: " << data.time << std::endl;
  std::cout << "Time step: " << data.dt << std::endl;
  std::cout << "Cycle: " << data.ncycle << std::endl;
  
  std::cout << "\nMeshBlock dimensions (including ghost zones):" << std::endl;
  std::cout << "  nx1 = " << data.nout1 << " (active: " << data.mb_indcs.nx1 << ")" << std::endl;
  std::cout << "  nx2 = " << data.nout2 << " (active: " << data.mb_indcs.nx2 << ")" << std::endl;
  std::cout << "  nx3 = " << data.nout3 << " (active: " << data.mb_indcs.nx3 << ")" << std::endl;
  std::cout << "  ng = " << data.mb_indcs.ng << std::endl;
  
  std::cout << "\nMesh size:" << std::endl;
  std::cout << "  x1: [" << data.mesh_size.x1min << ", " << data.mesh_size.x1max << "]" << std::endl;
  std::cout << "  x2: [" << data.mesh_size.x2min << ", " << data.mesh_size.x2max << "]" << std::endl;
  std::cout << "  x3: [" << data.mesh_size.x3min << ", " << data.mesh_size.x3max << "]" << std::endl;
  
  std::cout << "\nPhysics modules:" << std::endl;
  if (data.nhydro > 0) std::cout << "  Hydro: " << data.nhydro << " variables" << std::endl;
  if (data.nmhd > 0) std::cout << "  MHD: " << data.nmhd << " variables" << std::endl;
  if (data.nrad > 0) std::cout << "  Radiation: " << data.nrad << " angles" << std::endl;
  if (data.has_turb) std::cout << "  Turbulence driver: active" << std::endl;
  if (data.nz4c > 0) std::cout << "  Z4c: " << data.nz4c << " variables" << std::endl;
  if (data.nadm > 0) std::cout << "  ADM: " << data.nadm << " variables" << std::endl;
  
  std::cout << "================================\n" << std::endl;
}