//========================================================================================
// AthenaK Regridding Tool - Verification Implementation
// Licensed under the 3-clause BSD License (the "LICENSE")
//========================================================================================
//! \file verification.cpp
//  \brief Implementation of verification functions for regridded data

#include "verification.hpp"
#include <iostream>
#include <iomanip>
#include <fstream>
#include <cmath>
#include <algorithm>

//----------------------------------------------------------------------------------------
//! \fn bool Verification::RunAllChecks()
//! \brief Run all verification checks
bool Verification::RunAllChecks(const RestartData& input_data, 
                                const InterpolatedData& output_data,
                                VerificationResults& results) {
  SetError("");
  
  if (verbose_) {
    std::cout << "Running verification checks..." << std::endl;
  }
  
  // Initialize results
  results = VerificationResults{};
  
  // Check conservation laws
  if (!CheckConservation(input_data, output_data, results)) {
    results.all_passed = false;
  }
  
  // Check divergence-free constraint for B
  if (output_data.nmhd > 0) {
    if (!CheckDivergenceB(output_data, results)) {
      results.all_passed = false;
    }
  }
  
  // Check interpolation accuracy
  if (!CheckInterpolationAccuracy(input_data, output_data, results)) {
    results.all_passed = false;
  }
  
  // Generate summary
  std::ostringstream summary;
  summary << "Verification " << (results.all_passed ? "PASSED" : "FAILED");
  if (!results.all_passed) {
    summary << " - ";
    if (!results.mass_conserved) summary << "mass ";
    if (!results.momentum_conserved) summary << "momentum ";
    if (!results.energy_conserved) summary << "energy ";
    if (!results.divergence_free) summary << "divB ";
    summary << "check(s) failed";
  }
  results.summary = summary.str();
  
  if (verbose_) {
    PrintResults(results);
  }
  
  return results.all_passed;
}

//----------------------------------------------------------------------------------------
//! \fn bool Verification::CheckConservation()
//! \brief Check conservation of mass, momentum, and energy
bool Verification::CheckConservation(const RestartData& input_data, 
                                     const InterpolatedData& output_data,
                                     VerificationResults& results) {
  if (output_data.nmhd == 0) return true; // No MHD data to check
  
  bool all_conserved = true;
  
  // Calculate total quantities from input data
  Real input_mass = CalculateTotalMass(input_data.mhd_data, input_data.nmb_total,
                                      input_data.nout1, input_data.nout2, input_data.nout3);
  
  Real input_momentum_x = CalculateTotalMomentum(input_data.mhd_data, input_data.nmb_total,
                                                input_data.nout1, input_data.nout2, input_data.nout3, IM1);
  
  Real input_energy = CalculateTotalEnergy(input_data.mhd_data, input_data.nmb_total,
                                          input_data.nout1, input_data.nout2, input_data.nout3);
  
  // Calculate total quantities from output data
  Real output_mass = CalculateTotalMass(output_data.mhd_data, output_data.nmb_total,
                                       output_data.nout1, output_data.nout2, output_data.nout3);
  
  Real output_momentum_x = CalculateTotalMomentum(output_data.mhd_data, output_data.nmb_total,
                                                 output_data.nout1, output_data.nout2, output_data.nout3, IM1);
  
  Real output_energy = CalculateTotalEnergy(output_data.mhd_data, output_data.nmb_total,
                                           output_data.nout1, output_data.nout2, output_data.nout3);
  
  // Calculate relative errors
  results.mass_error = std::abs(output_mass - input_mass) / std::max(std::abs(input_mass), 1e-15);
  results.momentum_error = std::abs(output_momentum_x - input_momentum_x) / 
                          std::max(std::abs(input_momentum_x), 1e-15);
  results.energy_error = std::abs(output_energy - input_energy) / 
                        std::max(std::abs(input_energy), 1e-15);
  
  // Check against tolerance
  results.mass_conserved = (results.mass_error < tolerance_);
  results.momentum_conserved = (results.momentum_error < tolerance_);
  results.energy_conserved = (results.energy_error < tolerance_);
  
  if (!results.mass_conserved || !results.momentum_conserved || !results.energy_conserved) {
    all_conserved = false;
  }
  
  if (verbose_) {
    std::cout << "Conservation check:" << std::endl;
    std::cout << "  Mass error: " << std::scientific << results.mass_error 
              << " (" << (results.mass_conserved ? "PASS" : "FAIL") << ")" << std::endl;
    std::cout << "  Momentum error: " << results.momentum_error
              << " (" << (results.momentum_conserved ? "PASS" : "FAIL") << ")" << std::endl;
    std::cout << "  Energy error: " << results.energy_error
              << " (" << (results.energy_conserved ? "PASS" : "FAIL") << ")" << std::endl;
  }
  
  return all_conserved;
}

//----------------------------------------------------------------------------------------
//! \fn bool Verification::CheckDivergenceB()
//! \brief Check that ∇·B = 0 is maintained
bool Verification::CheckDivergenceB(const InterpolatedData& output_data,
                                   VerificationResults& results) {
  if (output_data.nmhd == 0) return true;
  
  ComputeDivergenceStatistics(output_data, results.max_div_b, results.avg_div_b);
  
  results.divergence_free = (results.max_div_b < tolerance_);
  
  if (verbose_) {
    std::cout << "Divergence check:" << std::endl;
    std::cout << "  Max |∇·B|: " << std::scientific << results.max_div_b 
              << " (" << (results.divergence_free ? "PASS" : "FAIL") << ")" << std::endl;
    std::cout << "  Avg |∇·B|: " << results.avg_div_b << std::endl;
  }
  
  return results.divergence_free;
}

//----------------------------------------------------------------------------------------
//! \fn bool Verification::CheckInterpolationAccuracy()
//! \brief Check interpolation accuracy (placeholder)
bool Verification::CheckInterpolationAccuracy(const RestartData& input_data,
                                              const InterpolatedData& output_data,
                                              VerificationResults& results) {
  // This is a placeholder implementation
  // In a full implementation, this would check against analytical solutions
  // or perform other accuracy tests
  
  results.interpolation_error = 0.0;
  results.failed_cells = 0;
  
  if (verbose_) {
    std::cout << "Interpolation accuracy check: SKIPPED (not implemented)" << std::endl;
  }
  
  return true;
}

//----------------------------------------------------------------------------------------
//! \fn Real Verification::CalculateTotalMass()
//! \brief Calculate total mass in the domain
Real Verification::CalculateTotalMass(const std::vector<Real>& data, int nmb, 
                                     int nout1, int nout2, int nout3) {
  Real total_mass = 0.0;
  
  for (int mb = 0; mb < nmb; ++mb) {
    for (int k = 0; k < nout3; ++k) {
      for (int j = 0; j < nout2; ++j) {
        for (int i = 0; i < nout1; ++i) {
          size_t idx = GetCCIndex(mb, IDN, k, j, i, 5, nout1, nout2, nout3);
          total_mass += data[idx];
        }
      }
    }
  }
  
  return total_mass;
}

//----------------------------------------------------------------------------------------
//! \fn Real Verification::CalculateTotalMomentum()
//! \brief Calculate total momentum component in the domain
Real Verification::CalculateTotalMomentum(const std::vector<Real>& data, int nmb,
                                         int nout1, int nout2, int nout3, int component) {
  Real total_momentum = 0.0;
  
  for (int mb = 0; mb < nmb; ++mb) {
    for (int k = 0; k < nout3; ++k) {
      for (int j = 0; j < nout2; ++j) {
        for (int i = 0; i < nout1; ++i) {
          size_t idx = GetCCIndex(mb, component, k, j, i, 5, nout1, nout2, nout3);
          total_momentum += data[idx];
        }
      }
    }
  }
  
  return total_momentum;
}

//----------------------------------------------------------------------------------------
//! \fn Real Verification::CalculateTotalEnergy()
//! \brief Calculate total energy in the domain
Real Verification::CalculateTotalEnergy(const std::vector<Real>& data, int nmb,
                                       int nout1, int nout2, int nout3) {
  Real total_energy = 0.0;
  
  for (int mb = 0; mb < nmb; ++mb) {
    for (int k = 0; k < nout3; ++k) {
      for (int j = 0; j < nout2; ++j) {
        for (int i = 0; i < nout1; ++i) {
          size_t idx = GetCCIndex(mb, IEN, k, j, i, 5, nout1, nout2, nout3);
          total_energy += data[idx];
        }
      }
    }
  }
  
  return total_energy;
}

//----------------------------------------------------------------------------------------
//! \fn void Verification::ComputeDivergenceStatistics()
//! \brief Compute statistics for ∇·B
void Verification::ComputeDivergenceStatistics(const InterpolatedData& data, 
                                              Real& max_div, Real& avg_div) const {
  max_div = 0.0;
  avg_div = 0.0;
  int cell_count = 0;
  
  int ng = 2; // Assuming ghost zones = 2
  int nx1 = data.nout1 - 2*ng;
  int nx2 = (data.nout2 > 1) ? data.nout2 - 2*ng : 1;
  int nx3 = (data.nout3 > 1) ? data.nout3 - 2*ng : 1;
  
  for (int mb = 0; mb < data.nmb_total; ++mb) {
    for (int k = ng; k < ng + nx3; ++k) {
      for (int j = ng; j < ng + nx2; ++j) {
        for (int i = ng; i < ng + nx1; ++i) {
          Real div_b = CalculateDivergenceB(data, mb, i, j, k);
          Real abs_div_b = std::abs(div_b);
          
          max_div = std::max(max_div, abs_div_b);
          avg_div += abs_div_b;
          cell_count++;
        }
      }
    }
  }
  
  if (cell_count > 0) {
    avg_div /= cell_count;
  }
}

//----------------------------------------------------------------------------------------
//! \fn Real Verification::CalculateDivergenceB()
//! \brief Calculate ∇·B at a specific cell
Real Verification::CalculateDivergenceB(const InterpolatedData& data, int mb, 
                                       int i, int j, int k) const {
  int nout1 = data.nout1;
  int nout2 = data.nout2;
  int nout3 = data.nout3;
  
  // Calculate divergence using face-centered fields
  size_t idx_b1f_right = mb * nout3 * nout2 * (nout1+1) + k * nout2 * (nout1+1) + j * (nout1+1) + (i+1);
  size_t idx_b1f_left = mb * nout3 * nout2 * (nout1+1) + k * nout2 * (nout1+1) + j * (nout1+1) + i;
  Real db1_dx = data.mhd_b1f_data[idx_b1f_right] - data.mhd_b1f_data[idx_b1f_left];
  
  Real db2_dy = 0.0;
  if (nout2 > 1) {
    size_t idx_b2f_up = mb * nout3 * (nout2+1) * nout1 + k * (nout2+1) * nout1 + (j+1) * nout1 + i;
    size_t idx_b2f_down = mb * nout3 * (nout2+1) * nout1 + k * (nout2+1) * nout1 + j * nout1 + i;
    db2_dy = data.mhd_b2f_data[idx_b2f_up] - data.mhd_b2f_data[idx_b2f_down];
  }
  
  Real db3_dz = 0.0;
  if (nout3 > 1) {
    size_t idx_b3f_front = mb * (nout3+1) * nout2 * nout1 + (k+1) * nout2 * nout1 + j * nout1 + i;
    size_t idx_b3f_back = mb * (nout3+1) * nout2 * nout1 + k * nout2 * nout1 + j * nout1 + i;
    db3_dz = data.mhd_b3f_data[idx_b3f_front] - data.mhd_b3f_data[idx_b3f_back];
  }
  
  return db1_dx + db2_dy + db3_dz;
}

//----------------------------------------------------------------------------------------
//! \fn size_t Verification::GetCCIndex()
//! \brief Get index for cell-centered variable
size_t Verification::GetCCIndex(int mb, int var, int k, int j, int i, int nvar,
                                int nx1, int nx2, int nx3) const {
  return mb * (size_t)nvar * nx3 * nx2 * nx1 + var * (size_t)nx3 * nx2 * nx1 + 
         (size_t)k * nx2 * nx1 + (size_t)j * nx1 + (size_t)i;
}

//----------------------------------------------------------------------------------------
//! \fn void Verification::PrintResults()
//! \brief Print verification results
void Verification::PrintResults(const VerificationResults& results) const {
  std::cout << "\n=== Verification Results ===" << std::endl;
  std::cout << "Overall status: " << (results.all_passed ? "PASSED" : "FAILED") << std::endl;
  
  std::cout << "\nConservation checks:" << std::endl;
  std::cout << "  Mass:     " << (results.mass_conserved ? "PASS" : "FAIL")
            << " (error: " << std::scientific << results.mass_error << ")" << std::endl;
  std::cout << "  Momentum: " << (results.momentum_conserved ? "PASS" : "FAIL")
            << " (error: " << results.momentum_error << ")" << std::endl;
  std::cout << "  Energy:   " << (results.energy_conserved ? "PASS" : "FAIL")
            << " (error: " << results.energy_error << ")" << std::endl;
  
  std::cout << "\nMagnetic field checks:" << std::endl;
  std::cout << "  ∇·B = 0:  " << (results.divergence_free ? "PASS" : "FAIL")
            << " (max |∇·B|: " << results.max_div_b << ")" << std::endl;
  
  std::cout << "============================" << std::endl;
}