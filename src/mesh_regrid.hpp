#ifndef MESH_REGRID_HPP_
#define MESH_REGRID_HPP_
//========================================================================================
// AthenaK Regridding Tool - MeshRegridder Class
// Licensed under the 3-clause BSD License (the "LICENSE")
//========================================================================================
//! \file mesh_regrid.hpp
//  \brief Class for regridding mesh from N³ to (2N)³ resolution

#include <vector>
#include <array>

// Standalone includes
#include "restart_reader.hpp"

//----------------------------------------------------------------------------------------
//! \struct NewMeshData
//! \brief Container for regridded mesh information
struct NewMeshData {
  // Updated mesh parameters
  int nmb_total_new;
  int root_level_new;
  RegionSize mesh_size_new;
  RegionIndcs mesh_indcs_new;
  RegionIndcs mb_indcs_new;  // stays the same
  Real dt_new;  // updated time step for finer resolution
  
  // New MeshBlock structure
  std::vector<LogicalLocation> lloc_eachmb_new;
  std::vector<float> cost_eachmb_new;
  
  // MPI distribution for new mesh (for output)
  std::vector<int> nmb_eachrank_new;
  std::vector<int> gids_eachrank_new;
  
  // Mapping from old to new MeshBlocks
  std::vector<std::array<int, 8>> old_to_new_map;  // Each old MB maps to 8 new MBs
  
  // Output dimensions (same as input since meshblock size unchanged)
  int nout1, nout2, nout3;
};

//----------------------------------------------------------------------------------------
//! \class MeshRegridder
//! \brief Handles mesh regridding from N³ to (2N)³ resolution
class MeshRegridder {
 public:
  MeshRegridder() = default;
  ~MeshRegridder() = default;
  
  // Main interface functions  
  bool RegridMesh(const RestartData& input_data, NewMeshData& output_data, int refinement_factor = 2);
  bool ValidateRegridding(const RestartData& input_data, const NewMeshData& output_data) const;
  void PrintRegridInfo(const RestartData& input_data, const NewMeshData& output_data) const;
  
  // Utility functions
  static int CalculateRefinementFactor(const RestartData& input_data);
  static bool IsValidForRegridding(const RestartData& input_data);
  
  // Error handling
  const std::string& GetLastError() const { return last_error_; }
  
 private:
  // Core regridding functions
  bool UpdateMeshParameters(const RestartData& input_data, NewMeshData& output_data, int refinement_factor = 2);
  bool CreateNewMeshBlockTree(const RestartData& input_data, NewMeshData& output_data);
  bool DistributeMeshBlocks(NewMeshData& output_data, int nranks = 1);
  
  // MeshBlock creation helpers
  LogicalLocation RefineLogicalLocation(const LogicalLocation& old_loc, int child_id) const;
  void CreateChildMeshBlocks(const LogicalLocation& parent_loc, int parent_gid,
                           std::vector<LogicalLocation>& new_locs,
                           std::vector<int>& parent_map) const;
  
  // Tree structure helpers
  bool ValidateTreeStructure(const std::vector<LogicalLocation>& locs) const;
  void SortByZOrder(std::vector<LogicalLocation>& locs, std::vector<int>& indices) const;
  std::uint64_t ComputeZOrder(const LogicalLocation& loc) const;
  
  // Validation helpers
  bool ValidateInputMesh(const RestartData& input_data) const;
  bool ValidateOutputMesh(const RestartData& input_data, const NewMeshData& output_data) const;
  bool CheckNeighborConsistency(const NewMeshData& output_data) const;
  
  // Error reporting
  void SetError(const std::string& error) const { last_error_ = error; }
  mutable std::string last_error_;
};

//----------------------------------------------------------------------------------------
//! Helper functions for logical location manipulation

//! Convert 3D child index to linear child ID (0-7)
inline int ChildIndex3DToID(int ic, int jc, int kc) {
  return ic + 2*jc + 4*kc;
}

//! Convert linear child ID to 3D child index
inline void ChildIDToIndex3D(int child_id, int& ic, int& jc, int& kc) {
  ic = child_id & 1;
  jc = (child_id >> 1) & 1;
  kc = (child_id >> 2) & 1;
}

//! Check if logical location is valid
inline bool IsValidLogicalLocation(const LogicalLocation& loc) {
  return (loc.level >= 0) && (loc.lx1 >= 0) && (loc.lx2 >= 0) && (loc.lx3 >= 0);
}

//! Get maximum logical location at a given level
inline LogicalLocation GetMaxLogicalLocation(int level, const RegionIndcs& mesh_indcs) {
  LogicalLocation max_loc;
  max_loc.level = level;
  max_loc.lx1 = mesh_indcs.nx1 * (1 << level) / mesh_indcs.nx1 - 1;
  max_loc.lx2 = mesh_indcs.nx2 * (1 << level) / mesh_indcs.nx2 - 1;
  max_loc.lx3 = mesh_indcs.nx3 * (1 << level) / mesh_indcs.nx3 - 1;
  return max_loc;
}

#endif // MESH_REGRID_HPP_