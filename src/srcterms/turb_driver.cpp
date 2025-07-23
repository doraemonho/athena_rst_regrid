//========================================================================================
// AthenaXXX astrophysical plasma code
// Copyright(C) 2020 James M. Stone <jmstone@ias.edu> and the Athena code team
// Licensed under the 3-clause BSD License (the "LICENSE")
//========================================================================================
//! \file turb_driver.cpp
//  \brief implementation of functions in TurbulenceDriver

#include <algorithm>
#include <iostream>
#include <limits>
#include <memory>

#include "athena.hpp"
#include "parameter_input.hpp"
#include "coordinates/cell_locations.hpp"
#include "mesh/mesh.hpp"
#include "hydro/hydro.hpp"
#include "mhd/mhd.hpp"
#include "ion-neutral/ion-neutral.hpp"
#include "driver/driver.hpp"
#include "utils/random.hpp"
#include "eos/eos.hpp"
#include "eos/ideal_c2p_hyd.hpp"
#include "eos/ideal_c2p_mhd.hpp"
#include "turb_driver.hpp"

//----------------------------------------------------------------------------------------
// constructor, initializes data structures and parameters

TurbulenceDriver::TurbulenceDriver(MeshBlockPack *pp, ParameterInput *pin) :
  pmy_pack(pp),
  force("force",1,1,1,1,1),
  Ampl("Ampl",1),
  kx_mode("kx_mode",1),ky_mode("ky_mode",1),kz_mode("kz_mode",1),
  OUPhase0_x("OUPhase0_x",1),OUPhase0_y("OUPhase0_y",1),OUPhase0_z("OUPhase0_z",1),
  OUPhase1_x("OUPhase1_x",1),OUPhase1_y("OUPhase1_y",1),OUPhase1_z("OUPhase1_z",1),
  Phase0_x("Phase0_x",1),Phase0_y("Phase0_y",1),Phase0_z("Phase0_z",1),
  Phase1_x("Phase1_x",1),Phase1_y("Phase1_y",1),Phase1_z("Phase1_z",1) {
  // allocate memory for force registers
  Mesh *pm = pmy_pack->pmesh;
  int nmb = pmy_pack->nmb_thispack;
  auto &indcs = pmy_pack->pmesh->mb_indcs;
  auto &one_d = pmy_pack->pmesh->one_d;
  auto &three_d = pmy_pack->pmesh->three_d;
  int ncells1 = indcs.nx1 + 2*(indcs.ng);
  int ncells2 = (indcs.nx2 > 1)? (indcs.nx2 + 2*(indcs.ng)) : 1;
  int ncells3 = (indcs.nx3 > 1)? (indcs.nx3 + 2*(indcs.ng)) : 1;
  Real lx = pm->mesh_size.x1max - pm->mesh_size.x1min;
  Real ly = pm->mesh_size.x2max - pm->mesh_size.x2min;
  Real lz = pm->mesh_size.x3max - pm->mesh_size.x3min;

  Kokkos::realloc(force, nmb, 3, ncells3, ncells2, ncells1);

  // range of modes including, corresponding to kmin and kmax 
  nlow  = pin->GetOrAddInteger("turb_driving", "nlow", 1);
  nhigh = pin->GetOrAddInteger("turb_driving", "nhigh", 2);
  nmid  = pin->GetOrAddReal("turb_driving", "nmid", 0.5*(nlow + nhigh));
  // energy injection rate
  dedt = pin->GetOrAddReal("turb_driving", "dedt", 0.0);
  // correlation time
  tcorr = pin->GetOrAddReal("turb_driving", "tcorr", 0.0);
  // driving interval
  dtdrive = pin->GetOrAddReal("turb_driving", "dtdrive", 0.0);
  // driving ratio (compressive to solenoidal)
  SolWeight = pin->GetOrAddReal("turb_driving", "SolWeight", 1.0);
  // driving type
  drive_type = pin->GetOrAddInteger("turb_driving", "drive_type", 0);
  // rseed
  rseed = pin->GetOrAddInteger("turb_driving", "rseed", -1);
  if (dtdrive > 0.0) {
    drive_time = pmy_pack->pmesh->time;
  }  

  if (one_d) {
    std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__ << std::endl
        << "One-dimensional driving is not supported" << std::endl;
    std::exit(EXIT_FAILURE);
  }

  Real dkx = 2.0*M_PI/lx;
  Real dky = 2.0*M_PI/ly;
  Real dkz = 2.0*M_PI/lz;
  Real dkl = std::max(dkx, std::max(dky, dkz));
  
  // Initialize OUVar  
  Real rho0 = pin->GetReal("turb_driving", "rho0");
  Real P0 = pin->GetReal("turb_driving", "P0");
  EquationOfState *peos;
  if (pmy_pack->phydro != nullptr) peos = (pmy_pack->phydro->peos);
  if (pmy_pack->pmhd != nullptr) peos = pmy_pack->pmhd->peos;
  Real gamma = peos->eos_data.gamma;

  Real c_s = sqrt(gamma*P0/rho0);
  Real kc = 0.5*(nlow + nhigh)*dkl;
  Real Mach = pow(2.0*M_PI * dedt /(rho0 * kc), 1.0/3.0) / c_s;
  Real TDecay = (2.0*M_PI)/(Mach*c_s*kc);
  OUVar = (Mach*c_s)/TDecay;

  Real nlow_sqr = nlow*nlow;
  Real nhigh_sqr = nhigh*nhigh;

  mode_count = 0;

  int nkx, nky, nkz;
  Real nsqr;
  // count the number of driving modes
  int modesPerK = 2;
  if (three_d)
    modesPerK = 4;

  for (nkx = 0; nkx <= nhigh; nkx++) {
    for (nky = 0; nky <= nhigh; nky++) {
      for (nkz = 0; nkz <= nhigh; nkz++) {
        if (nkx == 0 && nky == 0 && nkz == 0) continue;
        nsqr = SQR(nkx) + SQR(nky) + SQR(nkz);
        if (nsqr >= nlow_sqr && nsqr <= nhigh_sqr) {
          mode_count += modesPerK;
        }
      }
    }
  }
  Kokkos::realloc(Ampl, mode_count);

  // Mode between Fortran and C++
  // Mode(1,:) == kx_mode
  // Mode(2,:) == ky_mode
  // Mode(3,:) == kz_mode
  Kokkos::realloc(kx_mode, mode_count);
  Kokkos::realloc(ky_mode, mode_count);
  Kokkos::realloc(kz_mode, mode_count);

  Kokkos::realloc(OUPhase0_x, mode_count);
  Kokkos::realloc(OUPhase0_y, mode_count);
  Kokkos::realloc(OUPhase0_z, mode_count);
  Kokkos::realloc(OUPhase1_x, mode_count);
  Kokkos::realloc(OUPhase1_y, mode_count);
  Kokkos::realloc(OUPhase1_z, mode_count);

  Kokkos::realloc(Phase0_x, mode_count);
  Kokkos::realloc(Phase0_y, mode_count);
  Kokkos::realloc(Phase0_z, mode_count);
  Kokkos::realloc(Phase1_x, mode_count);
  Kokkos::realloc(Phase1_y, mode_count);
  Kokkos::realloc(Phase1_z, mode_count);

  Initialize();
}

//----------------------------------------------------------------------------------------
// destructor

TurbulenceDriver::~TurbulenceDriver() {
}

//----------------------------------------------------------------------------------------
//! \fn  noid Initialize
//  \brief Function to initialize the driver

void TurbulenceDriver::Initialize() {
  Mesh *pm = pmy_pack->pmesh;
  int nmb = pmy_pack->nmb_thispack;
  auto &indcs = pmy_pack->pmesh->mb_indcs;
  auto &one_d = pmy_pack->pmesh->one_d;
  auto &three_d = pmy_pack->pmesh->three_d;
  int ncells1 = indcs.nx1 + 2*(indcs.ng);
  int ncells2 = (indcs.nx2 > 1)? (indcs.nx2 + 2*(indcs.ng)) : 1;
  int ncells3 = (indcs.nx3 > 1)? (indcs.nx3 + 2*(indcs.ng)) : 1;
  //TODO: sync all the phases arrays to the host
  auto force_ = force;
  par_for("force_init_pgen",DevExeSpace(),
          0,nmb-1,0,2,0,ncells3-1,0,ncells2-1,0,ncells1-1,
  KOKKOS_LAMBDA(int m, int n, int k, int j, int i) {
    force_(m,n,k,j,i) = 0.0;
  });

  rstate.idum = rseed;

  auto kx_mode_ = kx_mode;
  auto ky_mode_ = ky_mode;
  auto kz_mode_ = kz_mode;
  auto Ampl_ = Ampl;

  Real dkx, dky, dkz, kx, ky, kz;
  Real lx = pm->mesh_size.x1max - pm->mesh_size.x1min;
  Real ly = pm->mesh_size.x2max - pm->mesh_size.x2min;
  Real lz = pm->mesh_size.x3max - pm->mesh_size.x3min;
  dkx = 2.0*M_PI/lx;
  dky = 2.0*M_PI/ly;
  dkz = 2.0*M_PI/lz;

  int nmode = 0;
  int nkx, nky, nkz;
  Real nsqr;
  Real nlow_sqr = nlow*nlow;
  Real nhigh_sqr = nhigh*nhigh;
  Real norm, kiso;
  auto mode_count_ = mode_count;
  for (nkx = 0; nkx <= nhigh; nkx++) {
    for (nky = 0; nky <= nhigh; nky++) {
      for (nkz = 0; nkz <= nhigh; nkz++) {
        if (nkx == 0 && nky == 0 && nkz == 0) continue;
        nsqr = SQR(nkx) + SQR(nky) + SQR(nkz);
        if (nsqr >= nlow_sqr && nsqr <= nhigh_sqr) {
          kx = dkx*nkx;
          ky = dky*nky;
          kz = dkz*nkz;
          Real dkl = std::max(dkx, std::max(dky, dkz));
          Real kc = nmid*dkl;
          Real km = (nhigh-nlow)*dkl;
          // Generate Fourier amplitudes
          kiso = sqrt(SQR(kx) + SQR(ky) + SQR(kz)); 
          // parabola profile Schmidt ea. (2006) Computers & Fluids 35 (2006) 353–371.
          
          if (kiso > 1e-16) {
            if (drive_type == 0) {
              norm = 1.0 - 4.0*SQR(kiso-kc)/SQR(km);
            } else if (drive_type == 1) {
              norm = 1.0 - SQR(kiso-nlow*dkl)/SQR(dkl);
            } else if (drive_type == 2) {
              norm = 1.0 - 0.25*SQR((kiso-dkl)/dkl);
            }
          } else {
            norm = 0.0;
          }
          Ampl_.h_view(nmode) = norm;
          kx_mode_.h_view(nmode) = kx;
          ky_mode_.h_view(nmode) = ky;
          kz_mode_.h_view(nmode) = kz;
          nmode++;

          Ampl_.h_view(nmode) = norm;
          kx_mode_.h_view(nmode) =  kx;
          ky_mode_.h_view(nmode) = -ky;
          kz_mode_.h_view(nmode) =  kz;
          nmode++;
          
          if (three_d) {
            Ampl_.h_view(nmode) = norm;
            kx_mode_.h_view(nmode) =  kx;
            ky_mode_.h_view(nmode) =  ky;
            kz_mode_.h_view(nmode) = -kz;
            nmode++;

            Ampl_.h_view(nmode) = norm;
            kx_mode_.h_view(nmode) = -kx;
            ky_mode_.h_view(nmode) =  ky;
            kz_mode_.h_view(nmode) = -kz;
            nmode++;
          }
        }
      }
    }
  }

  kx_mode_.template modify<HostMemSpace>();
  kx_mode_.template sync<DevExeSpace>();
  ky_mode_.template modify<HostMemSpace>();
  ky_mode_.template sync<DevExeSpace>();
  kz_mode_.template modify<HostMemSpace>();
  kz_mode_.template sync<DevExeSpace>();
  Ampl_.template modify<HostMemSpace>();
  Ampl_.template sync<DevExeSpace>();

  // Initialize OU process
  auto OUPhase0_x_ = OUPhase0_x;
  auto OUPhase0_y_ = OUPhase0_y;
  auto OUPhase0_z_ = OUPhase0_z;
  auto OUPhase1_x_ = OUPhase1_x;
  auto OUPhase1_y_ = OUPhase1_y;
  auto OUPhase1_z_ = OUPhase1_z;
  
  for (int n=0; n<mode_count_; n++) {
    //TODO: init the phase second time and compute the OU process
    //TODO: RandomPhase and OUVar <- RandomPhase done and OUVar wait to be done
    // First component in x,y,z
    OUPhase0_x_.h_view(n) = OUVar*RandomPhase(&rstate);
    OUPhase0_y_.h_view(n) = OUVar*RandomPhase(&rstate);
    OUPhase0_z_.h_view(n) = OUVar*RandomPhase(&rstate);
    // Second component in x,y,z
    OUPhase1_x_.h_view(n) = OUVar*RandomPhase(&rstate);
    OUPhase1_y_.h_view(n) = OUVar*RandomPhase(&rstate);
    OUPhase1_z_.h_view(n) = OUVar*RandomPhase(&rstate);
  }

  OUPhase0_x_.template modify<HostMemSpace>();
  OUPhase0_x_.template sync<DevExeSpace>();
  OUPhase0_y_.template modify<HostMemSpace>();
  OUPhase0_y_.template sync<DevExeSpace>();
  OUPhase0_z_.template modify<HostMemSpace>();
  OUPhase0_z_.template sync<DevExeSpace>(); 
  OUPhase1_x_.template modify<HostMemSpace>();
  OUPhase1_x_.template sync<DevExeSpace>();
  OUPhase1_y_.template modify<HostMemSpace>();
  OUPhase1_y_.template sync<DevExeSpace>();
  OUPhase1_z_.template modify<HostMemSpace>();
  OUPhase1_z_.template sync<DevExeSpace>();

  return;
}

//----------------------------------------------------------------------------------------
//! \fn  void IncludeModeEvolutionTasks
//  \brief Includes task in the operator split task list that constructs new modes with
//  random amplitudes and phases that can be used to evolve the force via an O-U process
//  Called by MeshBlockPack::AddPhysics() function

void TurbulenceDriver::IncludeInitializeModesTask(std::shared_ptr<TaskList> tl,
                                                  TaskID start) {
  auto id_init = tl->AddTask(&TurbulenceDriver::InitializeModes, this, start);
  auto id_add = tl->AddTask(&TurbulenceDriver::AddForcing, this, id_init);
  return;
}

//----------------------------------------------------------------------------------------
//! \fn  void IncludeAddForcingTask
//  \brief includes task in the stage_run task list for adding random forcing to fluid
//  as an explicit source terms in each stage of integrator
//  Called by MeshBlockPack::AddPhysics() function

void TurbulenceDriver::IncludeAddForcingTask(std::shared_ptr<TaskList> tl, TaskID start) {
  // These must be inserted after update task, but before send_u
  if (pmy_pack->pionn == nullptr) {
    if (pmy_pack->phydro != nullptr) {
      auto id = tl->InsertTask(&TurbulenceDriver::AddForcing, this,
                              pmy_pack->phydro->id.flux, pmy_pack->phydro->id.rkupdt);
    }
    if (pmy_pack->pmhd != nullptr) {
      auto id = tl->InsertTask(&TurbulenceDriver::AddForcing, this,
                              pmy_pack->pmhd->id.flux, pmy_pack->pmhd->id.rkupdt);
    }
  } else {
    auto id = tl->InsertTask(&TurbulenceDriver::AddForcing, this,
                            pmy_pack->pionn->id.n_flux, pmy_pack->pionn->id.n_rkupdt);
  }

  return;
}

//----------------------------------------------------------------------------------------
//! \fn InitializeModes()
// \brief Initializes driving, and so is only executed once at start of calc.
// Cannot be included in constructor since (it seems) Kokkos::par_for not allowed in cons.

TaskStatus TurbulenceDriver::InitializeModes(Driver *pdrive, int stage) {
  Mesh *pm = pmy_pack->pmesh;
  auto &indcs = pmy_pack->pmesh->mb_indcs;
  int is = indcs.is, ie = indcs.ie;
  int js = indcs.js, je = indcs.je;
  int ks = indcs.ks, ke = indcs.ke;
  int &nx1 = indcs.nx1;
  int &nx2 = indcs.nx2;
  int &nx3 = indcs.nx3;
  auto &gindcs = pm->mesh_indcs;
  
  Real lx = pm->mesh_size.x1max - pm->mesh_size.x1min;
  Real ly = pm->mesh_size.x2max - pm->mesh_size.x2min;
  Real lz = pm->mesh_size.x3max - pm->mesh_size.x3min;
  int gnx1 = gindcs.nx1;
  int gnx2 = gindcs.nx2;
  int gnx3 = gindcs.nx3;

  int ncells1 = indcs.nx1 + 2*(indcs.ng);
  int ncells2 = (indcs.nx2 > 1)? (indcs.nx2 + 2*(indcs.ng)) : 1;
  int ncells3 = (indcs.nx3 > 1)? (indcs.nx3 + 2*(indcs.ng)) : 1;

  flag_drive = true;
  // implusive driving check
  if (dtdrive > 0.0 &&  pm->time <= drive_time){
    flag_drive = false;
  }else if (dtdrive > 0.0 && pm->time > drive_time){
    drive_time += dtdrive;
  }
  Real dt = dtdrive > 0.0 ? dtdrive : pm->dt; // Use dtdrive directly instead of max(dtdrive, pm->dt)

  if (flag_drive){    
    // Zero out new force array
    auto force_ = force;
    int &nmb = pmy_pack->nmb_thispack;
    par_for("force_init", DevExeSpace(),0,nmb-1,0,2,ks,ke,js,je,is,ie,
    KOKKOS_LAMBDA(int m, int n, int k, int j, int i) {
      force_(m,n,k,j,i) = 0.0;
    });

    Real nlow_sqr = SQR(nlow);
    Real nhigh_sqr = SQR(nhigh);
    auto mode_count_ = mode_count;

    // Update the phases
    Real fcorr, gcorr;
    if (tcorr <= 1e-6) {  // use whitenoise
      fcorr = 0.0;
      gcorr = 1.0;
    } else {
      fcorr = std::exp(-dt/tcorr);
      gcorr = std::sqrt(1.0 - fcorr*fcorr);
    }

    // OU process
    auto OUPhase0_x_ = OUPhase0_x;
    auto OUPhase0_y_ = OUPhase0_y;
    auto OUPhase0_z_ = OUPhase0_z;
    auto OUPhase1_x_ = OUPhase1_x;
    auto OUPhase1_y_ = OUPhase1_y;
    auto OUPhase1_z_ = OUPhase1_z;
    for (int n=0; n<mode_count_; n++) {
      OUPhase0_x_.h_view(n) = fcorr*OUPhase0_x_.h_view(n) + gcorr*OUVar*RandomPhase(&rstate);
      OUPhase0_y_.h_view(n) = fcorr*OUPhase0_y_.h_view(n) + gcorr*OUVar*RandomPhase(&rstate);
      OUPhase0_z_.h_view(n) = fcorr*OUPhase0_z_.h_view(n) + gcorr*OUVar*RandomPhase(&rstate);
      OUPhase1_x_.h_view(n) = fcorr*OUPhase1_x_.h_view(n) + gcorr*OUVar*RandomPhase(&rstate);
      OUPhase1_y_.h_view(n) = fcorr*OUPhase1_y_.h_view(n) + gcorr*OUVar*RandomPhase(&rstate);
      OUPhase1_z_.h_view(n) = fcorr*OUPhase1_z_.h_view(n) + gcorr*OUVar*RandomPhase(&rstate);
    }

    OUPhase0_x_.template modify<HostMemSpace>();
    OUPhase0_x_.template sync<DevExeSpace>();
    OUPhase0_y_.template modify<HostMemSpace>();
    OUPhase0_y_.template sync<DevExeSpace>();
    OUPhase0_z_.template modify<HostMemSpace>();
    OUPhase0_z_.template sync<DevExeSpace>();

    OUPhase1_x_.template modify<HostMemSpace>();
    OUPhase1_x_.template sync<DevExeSpace>();
    OUPhase1_y_.template modify<HostMemSpace>();
    OUPhase1_y_.template sync<DevExeSpace>();
    OUPhase1_z_.template modify<HostMemSpace>();
    OUPhase1_z_.template sync<DevExeSpace>();

    // Calculate phases and compressible/incompressible modes
    auto Phase0_x_ = Phase0_x;
    auto Phase0_y_ = Phase0_y;
    auto Phase0_z_ = Phase0_z;
    auto Phase1_x_ = Phase1_x;
    auto Phase1_y_ = Phase1_y;
    auto Phase1_z_ = Phase1_z;
    auto kx_mode_ = kx_mode;
    auto ky_mode_ = ky_mode;
    auto kz_mode_ = kz_mode;
    for (int n=0; n<mode_count_; n++) {
      Real ka = 0.0;
      Real kb = 0.0;
      Real kk = 0.0;
      
      // kk = k \cdot k
      kk = kx_mode_.h_view(n)*kx_mode_.h_view(n) + 
           ky_mode_.h_view(n)*ky_mode_.h_view(n) + 
           kz_mode_.h_view(n)*kz_mode_.h_view(n);

      // ka = k \cdot Phase0
      ka = kx_mode_.h_view(n)*OUPhase0_x_.h_view(n) + 
           ky_mode_.h_view(n)*OUPhase0_y_.h_view(n) + 
           kz_mode_.h_view(n)*OUPhase0_z_.h_view(n);

      // kb = k \cdot Phase1
      kb = kx_mode_.h_view(n)*OUPhase1_x_.h_view(n) + 
           ky_mode_.h_view(n)*OUPhase1_y_.h_view(n) + 
           kz_mode_.h_view(n)*OUPhase1_z_.h_view(n);

      Real kk1 = 1.0/kk;

      Real diva_x  = kx_mode_.h_view(n)*ka*kk1;
      Real divb_x  = kx_mode_.h_view(n)*kb*kk1;
      Real curla_x = OUPhase0_x_.h_view(n) - divb_x;
      Real curlb_x = OUPhase1_x_.h_view(n) - diva_x;
      Phase0_x_.h_view(n) = SolWeight*curla_x + (1.0 - SolWeight)*divb_x;
      Phase1_x_.h_view(n) = SolWeight*curlb_x + (1.0 - SolWeight)*diva_x;

      Real diva_y  = ky_mode_.h_view(n)*ka*kk1;
      Real divb_y  = ky_mode_.h_view(n)*kb*kk1;
      Real curla_y = OUPhase0_y_.h_view(n) - divb_y;
      Real curlb_y = OUPhase1_y_.h_view(n) - diva_y;
      Phase0_y_.h_view(n) = SolWeight*curla_y + (1.0 - SolWeight)*divb_y;
      Phase1_y_.h_view(n) = SolWeight*curlb_y + (1.0 - SolWeight)*diva_y;

      Real diva_z  = kz_mode_.h_view(n)*ka*kk1;
      Real divb_z  = kz_mode_.h_view(n)*kb*kk1;
      Real curla_z = OUPhase0_z_.h_view(n) - divb_z;
      Real curlb_z = OUPhase1_z_.h_view(n) - diva_z;
      Phase0_z_.h_view(n) = SolWeight*curla_z + (1.0 - SolWeight)*divb_z;
      Phase1_z_.h_view(n) = SolWeight*curlb_z + (1.0 - SolWeight)*diva_z;
    }

    Phase0_x_.template modify<HostMemSpace>();
    Phase0_x_.template sync<DevExeSpace>();
    Phase0_y_.template modify<HostMemSpace>();
    Phase0_y_.template sync<DevExeSpace>();
    Phase0_z_.template modify<HostMemSpace>();
    Phase0_z_.template sync<DevExeSpace>();

    Phase1_x_.template modify<HostMemSpace>();
    Phase1_x_.template sync<DevExeSpace>();
    Phase1_y_.template modify<HostMemSpace>();
    Phase1_y_.template sync<DevExeSpace>();
    Phase1_z_.template modify<HostMemSpace>();
    Phase1_z_.template sync<DevExeSpace>();

    // Assemble Forcing Array from Phases array
    auto &size = pmy_pack->pmb->mb_size;
    auto Ampl_ = Ampl;
    for (int n=0; n<mode_count_; n++) {
      if (Ampl_.h_view(n) > 1e-16) { // only compute for non-zero modes for efficiency
        par_for("force_compute", DevExeSpace(),0,nmb-1,ks,ke,js,je,is,ie,
        KOKKOS_LAMBDA(int m, int k, int j, int i) {
          // kr = k \cdot r = kx*x + ky*y + kz*z
          Real kr = 0.0;

          if (ncells1-1 > 0) {
            Real &x1min = size.d_view(m).x1min;
            Real &x1max = size.d_view(m).x1max;
            Real x1v = CellCenterX(i-is, nx1, x1min, x1max);
            Real k1v = kx_mode_.d_view(n);
            kr += k1v*x1v;
          }

          if (ncells2-1 > 0) {
            Real &x2min = size.d_view(m).x2min;
            Real &x2max = size.d_view(m).x2max;
            Real x2v = CellCenterX(j-js, nx2, x2min, x2max);
            Real k2v = ky_mode_.d_view(n);
            kr += k2v*x2v;
          }

          if (ncells3-1 > 0) {
            Real &x3min = size.d_view(m).x3min;
            Real &x3max = size.d_view(m).x3max;
            Real x3v = CellCenterX(k-ks, nx3, x3min, x3max);
            Real k3v = kz_mode_.d_view(n);
            kr += k3v*x3v;
          }

          Real coskr = cos(kr);
          Real sinkr = sin(kr);
          force_(m,0,k,j,i)  += Ampl_.d_view(n)*( Phase0_x_.d_view(n)*coskr -
                                                Phase1_x_.d_view(n)*sinkr );
          force_(m,1,k,j,i)  += Ampl_.d_view(n)*( Phase0_y_.d_view(n)*coskr -
                                                Phase1_y_.d_view(n)*sinkr );
          force_(m,2,k,j,i)  += Ampl_.d_view(n)*( Phase0_z_.d_view(n)*coskr -
                                                Phase1_z_.d_view(n)*sinkr );
        });
      }
    }
  }
  return TaskStatus::complete;
}

//----------------------------------------------------------------------------------------
//! \fn apply forcing

TaskStatus TurbulenceDriver::AddForcing(Driver *pdrive, int stage) {
  Mesh *pm = pmy_pack->pmesh;
  auto &indcs = pmy_pack->pmesh->mb_indcs;
  int is = indcs.is, ie = indcs.ie;
  int js = indcs.js, je = indcs.je;
  int ks = indcs.ks, ke = indcs.ke;
  int &nmb = pmy_pack->nmb_thispack;
  int &nx1 = indcs.nx1;
  int &nx2 = indcs.nx2;
  int &nx3 = indcs.nx3;
  auto &gindcs = pm->mesh_indcs;
  
  Real lx = pm->mesh_size.x1max - pm->mesh_size.x1min;
  Real ly = pm->mesh_size.x2max - pm->mesh_size.x2min;
  Real lz = pm->mesh_size.x3max - pm->mesh_size.x3min;
  int gnx1 = gindcs.nx1;
  int gnx2 = gindcs.nx2;
  int gnx3 = gindcs.nx3;


  Real dt = dtdrive > 0.0 ? std::max(dtdrive, pm->dt) : pm->dt; // checking if dtdrive is set

  EquationOfState *peos;

  DvceArray5D<Real> u0, u0_;
  DvceArray5D<Real> w0;
  if (pmy_pack->phydro != nullptr) u0 = (pmy_pack->phydro->u0);
  if (pmy_pack->phydro != nullptr) peos = (pmy_pack->phydro->peos);
  if (pmy_pack->pmhd != nullptr) u0 = (pmy_pack->pmhd->u0);
  if (pmy_pack->pmhd != nullptr) peos = pmy_pack->pmhd->peos;

  bool flag_ideal = false;

  if (pmy_pack->phydro != nullptr){
    if (pmy_pack->phydro->peos->eos_data.is_ideal)
      flag_ideal = true;
  }
  if (pmy_pack->pmhd != nullptr){
    if (pmy_pack->pmhd->peos->eos_data.is_ideal)
      flag_ideal = true;
  }
  if (flag_drive){

    DvceArray5D<Real> u0, u0_;
    if (pmy_pack->phydro != nullptr) u0 = (pmy_pack->phydro->u0);
    if (pmy_pack->pmhd != nullptr) u0 = (pmy_pack->pmhd->u0);

    // compute the net momentum
    auto force_ = force;
    const int nmkji = nmb*nx3*nx2*nx1;
    const int nkji = nx3*nx2*nx1;
    const int nji  = nx2*nx1;
    Real t0 = 0.0, t1 = 0.0, t2 = 0.0, t3 = 0.0;

    Kokkos::parallel_reduce("net_mom_1", Kokkos::RangePolicy<>(DevExeSpace(),0,nmkji),
    KOKKOS_LAMBDA(const int &idx, Real &sum_t0, Real &sum_t1,
                                  Real &sum_t2, Real &sum_t3) {
      // compute n,k,j,i indices of thread
      int m = (idx)/nkji;
      int k = (idx - m*nkji)/nji;
      int j = (idx - m*nkji - k*nji)/nx1;
      int i = (idx - m*nkji - k*nji - j*nx1) + is;
      k += ks;
      j += js;
      Real den = u0(m,IDN,k,j,i);

      sum_t0 += den;
      sum_t1 += den*force_(m,0,k,j,i);
      sum_t2 += den*force_(m,1,k,j,i);
      sum_t3 += den*force_(m,2,k,j,i);
    }, Kokkos::Sum<Real>(t0), Kokkos::Sum<Real>(t1),
      Kokkos::Sum<Real>(t2), Kokkos::Sum<Real>(t3));


  #if MPI_PARALLEL_ENABLED
    Real m[4], gm[4];
    m[0] = t0; m[1] = t1; m[2] = t2; m[3] = t3;
    MPI_Allreduce(m, gm, 4, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    t0 = gm[0]; t1 = gm[1]; t2 = gm[2]; t3 = gm[3];
  #endif

    par_for("force_remove_net_mom", DevExeSpace(),0,nmb-1,ks,ke,js,je,is,ie,
    KOKKOS_LAMBDA(int m, int k, int j, int i) {
      force_(m,0,k,j,i) -= t1/t0;
      force_(m,1,k,j,i) -= t2/t0;
      force_(m,2,k,j,i) -= t3/t0;
    });

   // adjust the innjeciton rate
   //TODO: check if the injection rate is correct
    t0 = 0.0;
    t1 = 0.0;
    auto &size = pmy_pack->pmb->mb_size;
    Kokkos::parallel_reduce("net_mom_2", Kokkos::RangePolicy<>(DevExeSpace(),0,nmkji),
    KOKKOS_LAMBDA(const int &idx, Real &sum_t0, Real &sum_t1) {
      // compute n,k,j,i indices of thread
      int m = (idx)/nkji;
      int k = (idx - m*nkji)/nji;
      int j = (idx - m*nkji - k*nji)/nx1;
      int i = (idx - m*nkji - k*nji - j*nx1) + is;
      k += ks;
      j += js;

      Real den  = u0(m,IDN,k,j,i);
      Real mom1 = u0(m,IM1,k,j,i);
      Real mom2 = u0(m,IM2,k,j,i);
      Real mom3 = u0(m,IM3,k,j,i);

      Real dv1 = force_(m,0,k,j,i);
      Real dv2 = force_(m,1,k,j,i);
      Real dv3 = force_(m,2,k,j,i);
      
      Real dvol = size.d_view(m).dx1*size.d_view(m).dx2*size.d_view(m).dx3;

      sum_t0 += 0.5*den*(dv1*dv1 + dv2*dv2 + dv3*dv3)*dvol;
      sum_t1 += (mom1*dv1 + mom2*dv2 + mom3*dv3)*dvol;
    }, Kokkos::Sum<Real>(t0), Kokkos::Sum<Real>(t1));

  #if MPI_PARALLEL_ENABLED
    m[0] = t0; m[1] = t1;
    MPI_Allreduce(m, gm, 2, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    t0 = gm[0]; t1 = gm[1];
  #endif

    t0 = std::max(t0, 1.0e-20);
    t1 = std::max(t1, 1.0e-20);

    Real Vol = lx*ly*lz;
    Real a = t0*dt;
    Real b = t1; 
    Real c = -dedt*Vol; // TODO: check if this is correct
    Real sd = sqrt(SQR(b) - 4.0*a*c);
    Real Acorr = (sd - b) / (2.0*a);
    // TODO: Normalize the forcing, dt added in fortran but c++ moved to AddForcing
    par_for("force_norm", DevExeSpace(),0,nmb-1,ks,ke,js,je,is,ie,
    KOKKOS_LAMBDA(int m, int k, int j, int i) {
      force_(m,0,k,j,i) *= Acorr;
      force_(m,1,k,j,i) *= Acorr;
      force_(m,2,k,j,i) *= Acorr;
    });

    // push the force
    par_for("push",DevExeSpace(),0,nmb-1,ks,ke,js,je,is,ie,
    KOKKOS_LAMBDA(int m, int k, int j, int i) {
      Real v1 = force_(m,0,k,j,i);
      Real v2 = force_(m,1,k,j,i);
      Real v3 = force_(m,2,k,j,i);

      Real den = u0(m,IDN,k,j,i);

      Real M1   = u0(m,IM1,k,j,i); 
      Real M2   = u0(m,IM2,k,j,i);
      Real M3   = u0(m,IM3,k,j,i);
      Real dv1  = v1*dt;
      Real dv2  = v2*dt;
      Real dv3  = v3*dt;

      if (flag_ideal) {
        u0(m,IEN,k,j,i) +=     (dv1*M1  +  dv2*M2 +  dv3*M3) + 
                          0.5*(dv1*dv1 + dv2*dv2 + dv3*dv3)*den;
      }

      u0(m,IM1,k,j,i) += den*v1*dt;
      u0(m,IM2,k,j,i) += den*v2*dt;
      u0(m,IM3,k,j,i) += den*v3*dt;

    });

/*
    // compute the net momentum
    t0 = 0.0, t1 = 0.0, t2 = 0.0, t3 = 0.0;
    Kokkos::parallel_reduce("net_mom_2", Kokkos::RangePolicy<>(DevExeSpace(),0,nmkji),
    KOKKOS_LAMBDA(const int &idx, Real &sum_t0, Real &sum_t1,
                                  Real &sum_t2, Real &sum_t3) {
      // compute n,k,j,i indices of thread
      int m = (idx)/nkji;
      int k = (idx - m*nkji)/nji;
      int j = (idx - m*nkji - k*nji)/nx1;
      int i = (idx - m*nkji - k*nji - j*nx1) + is;
      k += ks;
      j += js;
      Real den = u0(m,IDN,k,j,i);

      sum_t0 += den;
      sum_t1 += u0(m,IM1,k,j,i);
      sum_t2 += u0(m,IM2,k,j,i);
      sum_t3 += u0(m,IM3,k,j,i);
    }, Kokkos::Sum<Real>(t0), Kokkos::Sum<Real>(t1),
      Kokkos::Sum<Real>(t2), Kokkos::Sum<Real>(t3));

    #if MPI_PARALLEL_ENABLED
      m[0] = t0; m[1] = t1; m[2] = t2; m[3] = t3;
      MPI_Allreduce(m, gm, 4, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
      t0 = gm[0]; t1 = gm[1]; t2 = gm[2]; t3 = gm[3];
    #endif
    std::cout << "t0: " << t0 << " net mom1: " << t1 << " net mom2: " << t2 << " net mom3: " << t3 << std::endl;


    t0 = 0.0, t1 = 0.0, t2 = 0.0, t3 = 0.0;

    Kokkos::parallel_reduce("net_mom_1", Kokkos::RangePolicy<>(DevExeSpace(),0,nmkji),
    KOKKOS_LAMBDA(const int &idx, Real &sum_t0, Real &sum_t1,
                                  Real &sum_t2, Real &sum_t3) {
      // compute n,k,j,i indices of thread
      int m = (idx)/nkji;
      int k = (idx - m*nkji)/nji;
      int j = (idx - m*nkji - k*nji)/nx1;
      int i = (idx - m*nkji - k*nji - j*nx1) + is;
      k += ks;
      j += js;
      Real den = u0(m,IDN,k,j,i);

      sum_t0 += den;
      sum_t1 += den*force_(m,0,k,j,i);
      sum_t2 += den*force_(m,1,k,j,i);
      sum_t3 += den*force_(m,2,k,j,i);
    }, Kokkos::Sum<Real>(t0), Kokkos::Sum<Real>(t1),
      Kokkos::Sum<Real>(t2), Kokkos::Sum<Real>(t3));


    #if MPI_PARALLEL_ENABLED
      m[0] = t0; m[1] = t1; m[2] = t2; m[3] = t3;
      MPI_Allreduce(m, gm, 4, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
      t0 = gm[0]; t1 = gm[1]; t2 = gm[2]; t3 = gm[3];
    #endif
    std::cout << "t0: " << t0 << " net force1: " << t1 << " net force2: " << t2 << " net force3: " << t3 << std::endl;
*/
  }
  return TaskStatus::complete;
}