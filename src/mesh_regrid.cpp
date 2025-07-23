//========================================================================================
// AthenaK Regridding Tool - MeshRegridder Implementation
// Licensed under the 3-clause BSD License (the "LICENSE")
//========================================================================================
//! \file mesh_regrid.cpp
//  \brief Implementation of MeshRegridder class

#include "mesh_regrid.hpp"
#include <iostream>
#include <algorithm>
#include <cmath>
#include <cassert>
#include <numeric>

//----------------------------------------------------------------------------------------
//! \fn bool MeshRegridder::RegridMesh()
//! \brief Main function to regrid mesh from N³ to (2N)³ resolution
bool MeshRegridder::RegridMesh(const RestartData& input_data, NewMeshData& output_data, int refinement_factor) {
  SetError("");
  
  std::cout << "Regridding mesh from " << input_data.mesh_indcs.nx1 << "×" 
            << input_data.mesh_indcs.nx2 << "×" << input_data.mesh_indcs.nx3 << " to " 
            << (refinement_factor*input_data.mesh_indcs.nx1) << "×"
            << (refinement_factor*input_data.mesh_indcs.nx2) << "×"
            << (refinement_factor*input_data.mesh_indcs.nx3) << " resolution..." << std::endl;
  
  // Validate input mesh
  if (!ValidateInputMesh(input_data)) {
    return false;
  }
  
  // Update mesh parameters for doubled resolution
  if (!UpdateMeshParameters(input_data, output_data, refinement_factor)) {
    return false;
  }
  
  // Create new MeshBlock tree structure
  if (!CreateNewMeshBlockTree(input_data, output_data)) {
    return false;
  }
  
  // Distribute new MeshBlocks across MPI ranks
  if (!DistributeMeshBlocks(output_data)) {
    return false;
  }
  
  // Validate the regridded mesh
  if (!ValidateOutputMesh(input_data, output_data)) {
    return false;
  }
  
  std::cout << "Successfully regridded mesh: " << input_data.nmb_total 
            << " → " << output_data.nmb_total_new << " MeshBlocks" << std::endl;
  
  return true;
}

//----------------------------------------------------------------------------------------
//! \fn bool MeshRegridder::UpdateMeshParameters()
//! \brief Update mesh parameters for doubled resolution
bool MeshRegridder::UpdateMeshParameters(const RestartData& input_data, NewMeshData& output_data, int refinement_factor) {
  // For N³→(rN)³: expand root grid by refinement_factor³ MeshBlocks
  int children_per_parent = refinement_factor * refinement_factor * refinement_factor;
  output_data.nmb_total_new = input_data.nmb_total * children_per_parent;
  
  std::cout << "Debug: input nmb_total = " << input_data.nmb_total << std::endl;
  std::cout << "Debug: children_per_parent = " << children_per_parent << std::endl;
  std::cout << "Debug: output nmb_total_new = " << output_data.nmb_total_new << std::endl;
  
  // New root level is one higher than the current maximum level
  int max_level = input_data.root_level;
  for (const auto& loc : input_data.lloc_eachmb) {
    max_level = std::max(max_level, loc.level);
  }
  output_data.root_level_new = max_level + 1;
  
  // Update mesh indices for refined resolution
  // Use mesh_indcs to represent the total mesh dimensions
  output_data.mesh_indcs_new = input_data.mesh_indcs;
  output_data.mesh_indcs_new.nx1 *= refinement_factor;  
  output_data.mesh_indcs_new.nx2 *= refinement_factor;  
  output_data.mesh_indcs_new.nx3 *= refinement_factor;
  
  // Preserve original domain bounds from input
  output_data.mesh_size_new.x1min = input_data.mesh_size.x1min;  
  output_data.mesh_size_new.x1max = input_data.mesh_size.x1max;   
  output_data.mesh_size_new.x2min = input_data.mesh_size.x2min;  
  output_data.mesh_size_new.x2max = input_data.mesh_size.x2max;   
  output_data.mesh_size_new.x3min = input_data.mesh_size.x3min;  
  output_data.mesh_size_new.x3max = input_data.mesh_size.x3max;
  
  
  // Recalculate mesh grid spacing correctly for doubled resolution
  Real domain_x1 = output_data.mesh_size_new.x1max - output_data.mesh_size_new.x1min;  // 200
  Real domain_x2 = output_data.mesh_size_new.x2max - output_data.mesh_size_new.x2min;  // 200  
  Real domain_x3 = output_data.mesh_size_new.x3max - output_data.mesh_size_new.x3min;  // 200
  
  output_data.mesh_size_new.dx1 = domain_x1 / output_data.mesh_indcs_new.nx1;  // 200/256 = 0.78125
  output_data.mesh_size_new.dx2 = domain_x2 / output_data.mesh_indcs_new.nx2;  // 200/256 = 0.78125
  output_data.mesh_size_new.dx3 = domain_x3 / output_data.mesh_indcs_new.nx3;  // 200/256 = 0.78125
  
  // MeshBlock size stays the same (64³ + ghost zones), we just have 8x more of them
  output_data.mb_indcs_new = input_data.mb_indcs;
  
  // Output dimensions stay the same (each MeshBlock is still 72³ including ghost zones)
  output_data.nout1 = input_data.nout1;
  output_data.nout2 = input_data.nout2; 
  output_data.nout3 = input_data.nout3;
  
  // Update time step for finer resolution (smaller dt for stability)
  output_data.dt_new = input_data.dt / refinement_factor;
  
  
  return true;
}

//----------------------------------------------------------------------------------------
//! \fn bool MeshRegridder::CreateNewMeshBlockTree()
//! \brief Create new MeshBlock tree - expand 1×1×1 to 2×2×2 root grid (8 MeshBlocks)
bool MeshRegridder::CreateNewMeshBlockTree(const RestartData& input_data, NewMeshData& output_data) {
  // Reserve space for new arrays (8x more blocks)
  std::cout << "Debug: Allocating space for " << output_data.nmb_total_new << " MeshBlocks" << std::endl;
  std::cout << "Debug: sizeof(LogicalLocation) = " << sizeof(LogicalLocation) << std::endl;
  std::cout << "Debug: Total memory for lloc_eachmb_new = " << output_data.nmb_total_new * sizeof(LogicalLocation) << " bytes" << std::endl;
  
  try {
    output_data.lloc_eachmb_new.reserve(output_data.nmb_total_new);
    output_data.cost_eachmb_new.reserve(output_data.nmb_total_new);
    output_data.old_to_new_map.resize(input_data.nmb_total);
  } catch (const std::exception& e) {
    SetError("Failed to allocate memory for new MeshBlock arrays: " + std::string(e.what()));
    return false;
  }
  
  // Create 8 MeshBlocks from each original MeshBlock (2×2×2 expansion)
  int new_gid = 0;
  for (int old_gid = 0; old_gid < input_data.nmb_total; ++old_gid) {
    const LogicalLocation& parent_loc = input_data.lloc_eachmb[old_gid];
    float parent_cost = input_data.cost_eachmb[old_gid];
    
    // The array is already sized to 8 elements
    
    // Create 8 child MeshBlocks in 2×2×2 arrangement
    for (int child_id = 0; child_id < 8; ++child_id) {
      // Create logical location for child block
      LogicalLocation child_loc = RefineLogicalLocation(parent_loc, child_id);
      output_data.lloc_eachmb_new.push_back(child_loc);
      
      // Each child MeshBlock gets the same cost as parent
      output_data.cost_eachmb_new.push_back(parent_cost);
      
      // Store mapping
      output_data.old_to_new_map[old_gid][child_id] = new_gid;
      
      new_gid++;
    }
  }
  
  // Verify we created the expected number of MeshBlocks
  if (new_gid != output_data.nmb_total_new) {
    SetError("Mismatch in created MeshBlock count");
    return false;
  }
  
  
  return ValidateTreeStructure(output_data.lloc_eachmb_new);
}

//----------------------------------------------------------------------------------------
//! \fn LogicalLocation MeshRegridder::RefineLogicalLocation()
//! \brief Create refined logical location for a child MeshBlock
//! For regridding N³→(2N)³: create 8 child blocks at root level that tile the domain
LogicalLocation MeshRegridder::RefineLogicalLocation(const LogicalLocation& parent_loc, 
                                                   int child_id) const {
  LogicalLocation child_loc;
  
  // Get 3D child indices (0 or 1 in each direction)
  int ic, jc, kc;
  ChildIDToIndex3D(child_id, ic, jc, kc);
  
  // Properly refine the parent's logical location
  // Child level is parent level + 1
  child_loc.level = parent_loc.level + 1;
  
  // Child coordinates are parent coordinates * 2 + child offset
  child_loc.lx1 = parent_loc.lx1 * 2 + ic;
  child_loc.lx2 = parent_loc.lx2 * 2 + jc;  
  child_loc.lx3 = parent_loc.lx3 * 2 + kc;
  
  return child_loc;
}

//----------------------------------------------------------------------------------------
//! \fn bool MeshRegridder::DistributeMeshBlocks()
//! \brief Distribute MeshBlocks across MPI ranks for optimal load balancing
bool MeshRegridder::DistributeMeshBlocks(NewMeshData& output_data, int nranks) {
  // For simplicity, assume single rank for now
  // In a full implementation, this would distribute based on cost
  
  output_data.nmb_eachrank_new.resize(nranks, 0);
  output_data.gids_eachrank_new.resize(nranks + 1, 0);
  
  if (nranks == 1) {
    output_data.nmb_eachrank_new[0] = output_data.nmb_total_new;
    output_data.gids_eachrank_new[0] = 0;
    output_data.gids_eachrank_new[1] = output_data.nmb_total_new;
  } else {
    // Distribute MeshBlocks evenly across ranks
    int mbs_per_rank = output_data.nmb_total_new / nranks;
    int remainder = output_data.nmb_total_new % nranks;
    
    int current_gid = 0;
    for (int rank = 0; rank < nranks; ++rank) {
      output_data.gids_eachrank_new[rank] = current_gid;
      output_data.nmb_eachrank_new[rank] = mbs_per_rank;
      if (rank < remainder) {
        output_data.nmb_eachrank_new[rank]++;
      }
      current_gid += output_data.nmb_eachrank_new[rank];
    }
    output_data.gids_eachrank_new[nranks] = output_data.nmb_total_new;
  }
  
  return true;
}

//----------------------------------------------------------------------------------------
//! \fn void MeshRegridder::SortByZOrder()
//! \brief Sort MeshBlocks by Z-order (Morton order) for better cache locality
void MeshRegridder::SortByZOrder(std::vector<LogicalLocation>& locs, 
                                std::vector<int>& indices) const {
  std::sort(indices.begin(), indices.end(), [&](int a, int b) {
    return ComputeZOrder(locs[a]) < ComputeZOrder(locs[b]);
  });
  
  // Apply the sorting to the logical locations
  std::vector<LogicalLocation> temp_locs = locs;
  for (size_t i = 0; i < locs.size(); ++i) {
    locs[i] = temp_locs[indices[i]];
  }
}

//----------------------------------------------------------------------------------------
//! \fn std::uint64_t MeshRegridder::ComputeZOrder()
//! \brief Compute Z-order (Morton order) value for a logical location
std::uint64_t MeshRegridder::ComputeZOrder(const LogicalLocation& loc) const {
  // Simple Z-order computation for 3D coordinates
  // Interleave bits of lx1, lx2, lx3
  std::uint64_t x = static_cast<std::uint64_t>(loc.lx1);
  std::uint64_t y = static_cast<std::uint64_t>(loc.lx2);
  std::uint64_t z = static_cast<std::uint64_t>(loc.lx3);
  
  std::uint64_t result = 0;
  for (int i = 0; i < 21; ++i) { // 21 bits per dimension for 64-bit result
    result |= ((x & (1ULL << i)) << (2*i)) |
              ((y & (1ULL << i)) << (2*i + 1)) |
              ((z & (1ULL << i)) << (2*i + 2));
  }
  
  // Include level in high-order bits
  result |= (static_cast<std::uint64_t>(loc.level) << 63);
  
  return result;
}

//----------------------------------------------------------------------------------------
//! \fn bool MeshRegridder::ValidateTreeStructure()
//! \brief Validate the structure of the MeshBlock tree
bool MeshRegridder::ValidateTreeStructure(const std::vector<LogicalLocation>& locs) const {
  for (const auto& loc : locs) {
    if (!IsValidLogicalLocation(loc)) {
      SetError("Invalid logical location found in tree");
      return false;
    }
  }
  
  // Check for duplicates
  std::vector<LogicalLocation> sorted_locs = locs;
  std::sort(sorted_locs.begin(), sorted_locs.end(), 
    [](const LogicalLocation& a, const LogicalLocation& b) {
      if (a.level != b.level) return a.level < b.level;
      if (a.lx1 != b.lx1) return a.lx1 < b.lx1;
      if (a.lx2 != b.lx2) return a.lx2 < b.lx2;
      return a.lx3 < b.lx3;
    });
  
  for (size_t i = 1; i < sorted_locs.size(); ++i) {
    const LogicalLocation& prev = sorted_locs[i-1];
    const LogicalLocation& curr = sorted_locs[i];
    if (prev.level == curr.level && prev.lx1 == curr.lx1 && 
        prev.lx2 == curr.lx2 && prev.lx3 == curr.lx3) {
      SetError("Duplicate logical location found in tree");
      return false;
    }
  }
  
  return true;
}

//----------------------------------------------------------------------------------------
//! \fn bool MeshRegridder::ValidateInputMesh()
//! \brief Validate input mesh for regridding compatibility
bool MeshRegridder::ValidateInputMesh(const RestartData& input_data) const {
  // Check that mesh dimensions are even (can be divided by 2)
  // Use mb_indcs (MeshBlock indices) which has the correct active zone dimensions
  if (input_data.mb_indcs.nx1 % 2 != 0 || 
      input_data.mb_indcs.nx2 % 2 != 0 || 
      input_data.mb_indcs.nx3 % 2 != 0) {
    SetError("Input mesh dimensions must be even for 2x regridding");
    return false;
  }
  
  // Check that we have MHD data
  if (input_data.nmhd == 0) {
    SetError("Input file must contain MHD data");
    return false;
  }
  
  // Check MeshBlock structure is valid
  if (input_data.nmb_total <= 0) {
    SetError("Invalid number of MeshBlocks");
    return false;
  }
  
  // Validate mesh dimensions match MeshBlock count
  int expected_mb_x = input_data.mesh_indcs.nx1 / input_data.mb_indcs.nx1;
  int expected_mb_y = input_data.mesh_indcs.nx2 / input_data.mb_indcs.nx2;
  int expected_mb_z = input_data.mesh_indcs.nx3 / input_data.mb_indcs.nx3;
  int expected_nmb = expected_mb_x * expected_mb_y * expected_mb_z;
  
  if (expected_nmb != input_data.nmb_total) {
    SetError("Mesh dimensions inconsistent with MeshBlock count. Expected " + 
             std::to_string(expected_nmb) + " MeshBlocks for " +
             std::to_string(input_data.mesh_indcs.nx1) + "×" +
             std::to_string(input_data.mesh_indcs.nx2) + "×" +
             std::to_string(input_data.mesh_indcs.nx3) + " mesh with " +
             std::to_string(input_data.mb_indcs.nx1) + "×" +
             std::to_string(input_data.mb_indcs.nx2) + "×" +
             std::to_string(input_data.mb_indcs.nx3) + " MeshBlocks, but found " +
             std::to_string(input_data.nmb_total));
    return false;
  }
  
  // Check ghost zones
  if (input_data.mb_indcs.ng < 2) {
    SetError("Insufficient ghost zones for prolongation operations");
    return false;
  }
  
  return true;
}

//----------------------------------------------------------------------------------------
//! \fn bool MeshRegridder::ValidateOutputMesh()
//! \brief Validate output mesh structure
bool MeshRegridder::ValidateOutputMesh(const RestartData& input_data, const NewMeshData& output_data) const {
  // Check total number of MeshBlocks
  if (output_data.nmb_total_new <= 0) {
    SetError("Invalid number of output MeshBlocks");
    return false;
  }
  
  // Check array sizes
  if (output_data.lloc_eachmb_new.size() != static_cast<size_t>(output_data.nmb_total_new)) {
    SetError("Inconsistent logical location array size");
    return false;
  }
  
  if (output_data.cost_eachmb_new.size() != static_cast<size_t>(output_data.nmb_total_new)) {
    SetError("Inconsistent cost array size");
    return false;
  }
  
  // Check mesh parameters are reasonable
  if (output_data.mesh_indcs_new.nx1 <= 0 || 
      output_data.mesh_indcs_new.nx2 <= 0 || 
      output_data.mesh_indcs_new.nx3 <= 0) {
    SetError("Invalid output mesh dimensions");
    return false;
  }
  
  return true;
}

//----------------------------------------------------------------------------------------
//! \fn int MeshRegridder::CalculateRefinementFactor()
//! \brief Calculate the refinement factor for the mesh
int MeshRegridder::CalculateRefinementFactor(const RestartData& input_data) {
  // For this implementation, we always refine by factor of 2 in each dimension
  return 2;
}

//----------------------------------------------------------------------------------------
//! \fn bool MeshRegridder::IsValidForRegridding()
//! \brief Check if input data is valid for regridding
bool MeshRegridder::IsValidForRegridding(const RestartData& input_data) {
  MeshRegridder regridder;
  return regridder.ValidateInputMesh(input_data);
}

//----------------------------------------------------------------------------------------
//! \fn bool MeshRegridder::ValidateRegridding()
//! \brief Validate the regridding was performed correctly
bool MeshRegridder::ValidateRegridding(const RestartData& input_data, 
                                     const NewMeshData& output_data) const {
  // Check that we have 8x more MeshBlocks
  if (output_data.nmb_total_new != input_data.nmb_total * 8) {
    SetError("Incorrect number of MeshBlocks after regridding");
    return false;
  }
  
  // Check that mesh resolution doubled
  if (output_data.mesh_indcs_new.nx1 != input_data.mesh_indcs.nx1 * 2 ||
      output_data.mesh_indcs_new.nx2 != input_data.mesh_indcs.nx2 * 2 ||
      output_data.mesh_indcs_new.nx3 != input_data.mesh_indcs.nx3 * 2) {
    SetError("Mesh resolution not correctly doubled");
    return false;
  }
  
  // Check that MeshBlock size stayed the same
  if (output_data.mb_indcs_new.nx1 != input_data.mb_indcs.nx1 ||
      output_data.mb_indcs_new.nx2 != input_data.mb_indcs.nx2 ||
      output_data.mb_indcs_new.nx3 != input_data.mb_indcs.nx3) {
    SetError("MeshBlock size should remain unchanged");
    return false;
  }
  
  // Check physical domain size is preserved
  Real dx_ratio = output_data.mesh_size_new.dx1 / input_data.mesh_size.dx1;
  if (std::abs(dx_ratio - 0.5) > 1e-12) {
    SetError("Grid spacing not correctly halved");
    return false;
  }
  
  return true;
}

//----------------------------------------------------------------------------------------
//! \fn void MeshRegridder::PrintRegridInfo()
//! \brief Print information about the regridding operation
void MeshRegridder::PrintRegridInfo(const RestartData& input_data, 
                                   const NewMeshData& output_data) const {
  std::cout << "\n=== Mesh Regridding Summary ===" << std::endl;
  std::cout << "Original mesh:" << std::endl;
  std::cout << "  Resolution: " << input_data.mesh_indcs.nx1 << " × " 
            << input_data.mesh_indcs.nx2 << " × " << input_data.mesh_indcs.nx3 << std::endl;
  std::cout << "  MeshBlocks: " << input_data.nmb_total << std::endl;
  std::cout << "  Grid spacing: " << input_data.mesh_size.dx1 << " × " 
            << input_data.mesh_size.dx2 << " × " << input_data.mesh_size.dx3 << std::endl;
  
  std::cout << "\nRegridded mesh:" << std::endl;
  std::cout << "  Resolution: " << output_data.mesh_indcs_new.nx1 << " × " 
            << output_data.mesh_indcs_new.nx2 << " × " << output_data.mesh_indcs_new.nx3 << std::endl;
  std::cout << "  MeshBlocks: " << output_data.nmb_total_new << std::endl;
  std::cout << "  Grid spacing: " << output_data.mesh_size_new.dx1 << " × " 
            << output_data.mesh_size_new.dx2 << " × " << output_data.mesh_size_new.dx3 << std::endl;
  
  std::cout << "\nMeshBlock details:" << std::endl;
  std::cout << "  Size: " << output_data.mb_indcs_new.nx1 << " × " 
            << output_data.mb_indcs_new.nx2 << " × " << output_data.mb_indcs_new.nx3 << std::endl;
  std::cout << "  Ghost zones: " << output_data.mb_indcs_new.ng << std::endl;
  std::cout << "  Refinement factor: 2× in each dimension" << std::endl;
  std::cout << "  Total cells: " << output_data.mesh_indcs_new.nx1 * 
                                     output_data.mesh_indcs_new.nx2 * 
                                     output_data.mesh_indcs_new.nx3 << std::endl;
  std::cout << "==============================\n" << std::endl;
}