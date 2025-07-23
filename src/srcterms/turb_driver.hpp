#ifndef SRCTERMS_TURB_DRIVER_HPP_
#define SRCTERMS_TURB_DRIVER_HPP_
//========================================================================================
// AthenaXXX astrophysical plasma code
// Copyright(C) 2020 James M. Stone <jmstone@ias.edu> and the Athena code team
// Licensed under the 3-clause BSD License (the "LICENSE")
//========================================================================================
//! \file turb_driver.hpp
//  \brief defines turbulence driver class, which implements data and functions for
//  randomly forced turbulence which evolves via an Ornstein-Uhlenbeck stochastic process

#include <memory>

#include "athena.hpp"
#include "mesh/mesh.hpp"
#include "parameter_input.hpp"
#include "utils/random.hpp"

//----------------------------------------------------------------------------------------
//! \class TurbulenceDriver

class TurbulenceDriver {
 public:
  TurbulenceDriver(MeshBlockPack *pp, ParameterInput *pin);
  ~TurbulenceDriver();

  DvceArray5D<Real> force;  // arrays used for turb forcing
  RNG_State rstate;                    // random state

  DualArray1D<Real> kx_mode, ky_mode, kz_mode;
  DualArray1D<Real> OUPhase0_x, OUPhase0_y, OUPhase0_z, OUPhase1_x, OUPhase1_y, OUPhase1_z; 
  DualArray1D<Real> Phase0_x, Phase0_y, Phase0_z, Phase1_x, Phase1_y, Phase1_z;
  DualArray1D<Real> Ampl;

  // parameters of driving
  int nlow, nhigh;
  Real nmid;
  int mode_count;
  Real tcorr, dedt;
  Real dtdrive;
  Real drive_time = 0.0;
  Real SolWeight;
  Real WeigthNorm;
  Real OUVar;
  int drive_type;
  bool flag_drive;
  int rseed;

  // functions
  void IncludeInitializeModesTask(std::shared_ptr<TaskList> tl, TaskID start);
  void IncludeAddForcingTask(std::shared_ptr<TaskList> tl, TaskID start);
  TaskStatus InitializeModes(Driver *pdrive, int stage);
  TaskStatus AddForcing(Driver *pdrive, int stage);
  void UpdateBlockCount();
  void Initialize();
 private:
  bool first_time = true;   // flag to enable initialization on first call
  MeshBlockPack *pmy_pack;  // ptr to MeshBlockPack containing this TurbulenceDriver


//----------------------------------------------------------------------------------------
//! \fn  static Real RandomPhase
// \brief Generates a normally distributed random number using the
//        specific Box-Muller form from the Fortran RandomPhase function.
//
// input: state Pointer to the RNG state (must be initialized).
// return: double A random number from a standard normal distribution.

KOKKOS_INLINE_FUNCTION
static Real RandomPhase(RNG_State *state) {
    Real r0, r1;
    // Generate two uniform random numbers in (0, 1)
    // Need to ensure r0 is not exactly 0 for log(r0).
    // RanSt returns values up to RNMX = 1.0 - DBL_EPSILON.
    // It might theoretically return 0, though highly unlikely.
    
    r0 = RanSt(state);
    if (r0 == 0.0) {
      r0 = 1e-10;
    }

    r1 = RanSt(state);

    // Apply the exact formula from the Fortran code
    return sqrt(-2.0 * log(r0)) * cos(2.0*M_PI * r1);
}
};

#endif  // SRCTERMS_TURB_DRIVER_HPP_
