#ifndef DATA_INTERPOLATOR_HPP_
#define DATA_INTERPOLATOR_HPP_
//========================================================================================
// AthenaK Regridding Tool - DataInterpolator Class
// Licensed under the 3-clause BSD License (the "LICENSE")
//========================================================================================
//! \file data_interpolator.hpp
//  \brief Class for interpolating physics data using prolongation operators

#include <vector>
#include <array>

// Standalone includes
#include "restart_reader.hpp"
#include "mesh_regrid.hpp"

//----------------------------------------------------------------------------------------
//! \struct InterpolatedData
//! \brief Container for interpolated physics data at new resolution
struct InterpolatedData {
  // Interpolated physics data arrays (same structure as RestartData but for new mesh)
  std::vector<Real> hydro_data;      // [nmb_new][nhydro][nout3][nout2][nout1]
  std::vector<Real> mhd_data;        // [nmb_new][nmhd][nout3][nout2][nout1] 
  std::vector<Real> mhd_b1f_data;    // [nmb_new][nout3][nout2][nout1+1]
  std::vector<Real> mhd_b2f_data;    // [nmb_new][nout3][nout2+1][nout1]
  std::vector<Real> mhd_b3f_data;    // [nmb_new][nout3+1][nout2][nout1]
  std::vector<Real> rad_data;        // [nmb_new][nrad][nout3][nout2][nout1]
  std::vector<Real> force_data;      // [nmb_new][nforce][nout3][nout2][nout1]
  std::vector<Real> z4c_data;        // [nmb_new][nz4c][nout3][nout2][nout1]
  std::vector<Real> adm_data;        // [nmb_new][nadm][nout3][nout2][nout1]
  
  // Data dimensions
  int nmb_total, nhydro, nmhd, nrad, nforce, nz4c, nadm;
  int nout1, nout2, nout3;
  
  // Helper functions for array indexing (same as RestartData)
  size_t GetCCIndex(int mb, int var, int k, int j, int i, int nvar) const {
    return mb*(size_t)nvar*nout3*nout2*nout1 + (size_t)var*nout3*nout2*nout1 + 
           (size_t)k*nout2*nout1 + (size_t)j*nout1 + (size_t)i;
  }
  
  size_t GetFCX1Index(int mb, int k, int j, int i) const {
    return mb*(size_t)nout3*nout2*(nout1+1) + (size_t)k*nout2*(nout1+1) + 
           (size_t)j*(nout1+1) + (size_t)i;
  }
  
  size_t GetFCX2Index(int mb, int k, int j, int i) const {
    return mb*(size_t)nout3*(nout2+1)*nout1 + (size_t)k*(nout2+1)*nout1 + 
           (size_t)j*nout1 + (size_t)i;
  }
  
  size_t GetFCX3Index(int mb, int k, int j, int i) const {
    return mb*(size_t)(nout3+1)*nout2*nout1 + (size_t)k*nout2*nout1 + 
           (size_t)j*nout1 + (size_t)i;
  }
};

//----------------------------------------------------------------------------------------
//! \class DataInterpolator
//! \brief Handles physics data interpolation using AthenaK prolongation operators
class DataInterpolator {
 public:
  DataInterpolator() = default;
  ~DataInterpolator() = default;
  
  // Main interface functions
  bool InterpolateAllData(const RestartData& input_data, const NewMeshData& mesh_data,
                         InterpolatedData& output_data);
  
  // Individual interpolation functions
  bool InterpolateCellCenteredData(const RestartData& input_data, const NewMeshData& mesh_data,
                                  InterpolatedData& output_data);
  bool InterpolateFaceCenteredData(const RestartData& input_data, const NewMeshData& mesh_data,
                                  InterpolatedData& output_data);
  
  // Verification functions
  bool VerifyConservation(const RestartData& input_data, const InterpolatedData& output_data) const;
  bool VerifyDivergenceB(const InterpolatedData& output_data) const;
  
  // Configuration
  
  // Error handling
  const std::string& GetLastError() const { return last_error_; }
  
 private:
  // Core prolongation functions (based on AthenaK algorithms)
  void ProlongCC(const Real* coarse_data, Real* fine_data, 
                int nvar, int nx1, int nx2, int nx3, int ng,
                int coarse_i, int coarse_j, int coarse_k) const;
  
  // Face-centered field prolongation functions
  void ProlongFC1(const Real* coarse_b1f, Real* fine_b1f,
                  int nx1, int nx2, int nx3, int ng,
                  int coarse_i, int coarse_j, int coarse_k) const;
  
  void ProlongFC2(const Real* coarse_b2f, Real* fine_b2f,
                  int nx1, int nx2, int nx3, int ng,
                  int coarse_i, int coarse_j, int coarse_k) const;
  
  void ProlongFC3(const Real* coarse_b3f, Real* fine_b3f,
                  int nx1, int nx2, int nx3, int ng,
                  int coarse_i, int coarse_j, int coarse_k) const;
  
  void ProlongFCSharedX1Face(const Real* coarse_data, Real* fine_data,
                           int nx1, int nx2, int nx3, int ng,
                           int coarse_i, int coarse_j, int coarse_k) const;
  
  void ProlongFCSharedX2Face(const Real* coarse_data, Real* fine_data,
                           int nx1, int nx2, int nx3, int ng,
                           int coarse_i, int coarse_j, int coarse_k) const;
  
  void ProlongFCSharedX3Face(const Real* coarse_data, Real* fine_data,
                           int nx1, int nx2, int nx3, int ng,
                           int coarse_i, int coarse_j, int coarse_k) const;
  
  void ProlongFC(const Real* coarse_data, Real* fine_data,
                int nvar, int nx1, int nx2, int nx3, int ng,
                int child_i, int child_j, int child_k,
                int child_gid, int parent_gid, 
                const RestartData& input_data, int face_dir) const;
  
  void ProlongFCInternal(Real* b1f_data, Real* b2f_data, Real* b3f_data,
                        int nx1, int nx2, int nx3, int ng,
                        int fine_i, int fine_j, int fine_k) const;
  
  // Helper functions
  void AllocateOutputArrays(const RestartData& input_data, const NewMeshData& mesh_data,
                           InterpolatedData& output_data) const;
  
  void CopyPhysicsParameters(const RestartData& input_data, InterpolatedData& output_data) const;
  
  // Index mapping functions
  void GetCoarseIndices(int child_id, int fine_i, int fine_j, int fine_k,
                       int& coarse_i, int& coarse_j, int& coarse_k) const;
  
  void GetFineIndices(int child_id, int coarse_i, int coarse_j, int coarse_k,
                     int& fine_i_start, int& fine_j_start, int& fine_k_start) const;
  
  // Primitive/conserved variable conversion
  bool ConvertConsToProms(Real* u_data, Real* w_data, const Real* b_data,
                         int nvar, int nx1, int nx2, int nx3, int ng) const;
  
  bool ConvertPrimsTocons(const Real* w_data, const Real* b_data, Real* u_data,
                         int nvar, int nx1, int nx2, int nx3, int ng) const;
  
  
  // Array indexing helpers
  size_t GetIndex4D(int n, int k, int j, int i, int nx1, int nx2, int nx3) const {
    return n*(size_t)nx3*nx2*nx1 + (size_t)k*nx2*nx1 + (size_t)j*nx1 + (size_t)i;
  }
  
  size_t GetIndex5D(int m, int n, int k, int j, int i, 
                   int nvar, int nx1, int nx2, int nx3) const {
    return m*(size_t)nvar*nx3*nx2*nx1 + (size_t)n*nx3*nx2*nx1 + 
           (size_t)k*nx2*nx1 + (size_t)j*nx1 + (size_t)i;
  }
  
  // Helper functions matching AthenaK exactly
  Real Sign(Real x) const { return (x < 0.0) ? -1.0 : 1.0; }  // AthenaK's SIGN macro
  
  // Configuration
  
  // Error reporting
  void SetError(const std::string& error) const { last_error_ = error; }
  mutable std::string last_error_;
};

#endif // DATA_INTERPOLATOR_HPP_