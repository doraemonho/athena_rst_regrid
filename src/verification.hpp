#ifndef VERIFICATION_HPP_
#define VERIFICATION_HPP_
//========================================================================================
// AthenaK Regridding Tool - Verification Module
// Licensed under the 3-clause BSD License (the "LICENSE")  
//========================================================================================
//! \file verification.hpp
//  \brief Verification and validation functions for regridded data

#include <string>
#include <vector>

// Standalone includes
#include "restart_reader.hpp"
#include "data_interpolator.hpp"

//----------------------------------------------------------------------------------------
//! \struct VerificationResults
//! \brief Results from verification checks
struct VerificationResults {
  // Conservation checks
  bool mass_conserved = true;
  bool momentum_conserved = true;
  bool energy_conserved = true;
  bool magnetic_flux_conserved = true;
  
  Real mass_error = 0.0;
  Real momentum_error = 0.0;
  Real energy_error = 0.0;
  Real magnetic_flux_error = 0.0;
  
  // Divergence checks
  bool divergence_free = true;
  Real max_div_b = 0.0;
  Real avg_div_b = 0.0;
  
  // Accuracy checks  
  Real interpolation_error = 0.0;
  int failed_cells = 0;
  
  // Overall result
  bool all_passed = true;
  std::string summary;
};

//----------------------------------------------------------------------------------------
//! \class Verification
//! \brief Handles verification and validation of regridded data
class Verification {
 public:
  Verification() = default;
  ~Verification() = default;
  
  // Main verification functions
  bool RunAllChecks(const RestartData& input_data, const InterpolatedData& output_data,
                   VerificationResults& results);
  
  // Individual verification functions
  bool CheckConservation(const RestartData& input_data, const InterpolatedData& output_data,
                        VerificationResults& results);
  bool CheckDivergenceB(const InterpolatedData& output_data, VerificationResults& results);
  bool CheckInterpolationAccuracy(const RestartData& input_data, 
                                 const InterpolatedData& output_data,
                                 VerificationResults& results);
  
  // Reporting functions
  void PrintResults(const VerificationResults& results) const;
  void WriteReport(const std::string& filename, const VerificationResults& results) const;
  
  // Configuration
  void SetTolerance(Real tol) { tolerance_ = tol; }
  void SetVerbose(bool verbose) { verbose_ = verbose; }
  
  // Error handling
  const std::string& GetLastError() const { return last_error_; }
  
 private:
  // Conservation check helpers
  Real CalculateTotalMass(const std::vector<Real>& data, int nmb, int nout1, int nout2, int nout3);
  Real CalculateTotalMomentum(const std::vector<Real>& data, int nmb, int nout1, int nout2, int nout3, int component);
  Real CalculateTotalEnergy(const std::vector<Real>& data, int nmb, int nout1, int nout2, int nout3);
  Real CalculateMagneticFlux(const std::vector<Real>& b_data, int nmb, int nout1, int nout2, int nout3);
  
  // Divergence check helpers
  Real CalculateDivergenceB(const InterpolatedData& data, int mb, int i, int j, int k) const;
  void ComputeDivergenceStatistics(const InterpolatedData& data, Real& max_div, Real& avg_div) const;
  
  // Accuracy check helpers
  Real ComputeInterpolationError(const RestartData& input_data, 
                                const InterpolatedData& output_data) const;
  
  // Utility functions
  size_t GetCCIndex(int mb, int var, int k, int j, int i, int nvar, int nx1, int nx2, int nx3) const;
  size_t GetFCIndex(int mb, int k, int j, int i, int nx1, int nx2, int nx3) const;
  
  // Error reporting
  void SetError(const std::string& error) const { last_error_ = error; }
  mutable std::string last_error_;
  
  // Configuration
  Real tolerance_ = 1e-12;
  bool verbose_ = false;
};

#endif // VERIFICATION_HPP_