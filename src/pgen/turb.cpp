//========================================================================================
// AthenaXXX astrophysical plasma code
// Copyright(C) 2020 James M. Stone <jmstone@ias.edu> and the Athena code team
// Licensed under the 3-clause BSD License (the "LICENSE")
//========================================================================================
//! \file turb.cpp
//  \brief Problem generator for turbulence
#include <iostream> // cout

#include "athena.hpp"
#include "globals.hpp"
#include "parameter_input.hpp"
#include "coordinates/cell_locations.hpp"
#include "mesh/mesh.hpp"
#include "eos/eos.hpp"
#include "hydro/hydro.hpp"
#include "mhd/mhd.hpp"
#include "pgen.hpp"
#include "srcterms/turb_driver.hpp"

// User-defined history functions
void TurbulentHistory(HistoryData *pdata, Mesh *pm);
void RefinementCondition(MeshBlockPack* pmbp);
bool RefineFlag1 = false;
bool RefineFlag2 = false;
auto AMR_t1 = 0.0;
auto AMR_t2 = 0.0;
//----------------------------------------------------------------------------------------
//! \fn void MeshBlock::Turb_()
//  \brief Problem Generator for turbulence

void ProblemGenerator::UserProblem(ParameterInput *pin, const bool restart) {
  user_ref_func  = RefinementCondition;

  MeshBlockPack *pmbp = pmy_mesh_->pmb_pack;
  auto &size = pmbp->pmb->mb_size;
  auto &indcs = pmy_mesh_->mb_indcs;

  if (pmbp->phydro == nullptr && pmbp->pmhd == nullptr) {
    std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__ << std::endl
       << "Turbulence problem generator can only be run with Hydro and/or MHD, but no "
       << "<hydro> or <mhd> block in input file" << std::endl;
    exit(EXIT_FAILURE);
  }

  // enroll user history function
  user_hist_func = TurbulentHistory;
  user_hist = pin->GetOrAddBoolean("problem","user_hist",false);

  // capture variables for kernel
  int &is = indcs.is; int &ie = indcs.ie;
  int &js = indcs.js; int &je = indcs.je;
  int &ks = indcs.ks; int &ke = indcs.ke;

  if (restart) return;
  
  Real cs = pin->GetOrAddReal("eos","iso_sound_speed",1.0);
  //Real beta = pin->GetOrAddReal("problem","beta",1.0);
  Real B0   = pin->GetReal("problem","B0");
  Real t0   = pin->GetOrAddReal("problem","Btheta",0.0); // 0 = z-axis, pi/2 = y-axis
  if (abs(t0) > M_PI/2.0) {
    std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__ << std::endl
              << "Btheta must be between 0 and pi/2" << std::endl;
    exit(EXIT_FAILURE);
  }
  
  AMR_t1 = pin->GetOrAddReal("problem","AMR_t1",0.0);
  AMR_t2 = pin->GetOrAddReal("problem","AMR_t2",0.0);

  // Initialize Hydro variables -------------------------------
  if (pmbp->phydro != nullptr) {
    Real d_i = pin->GetOrAddReal("problem","d_i",1.0);
    Real d_n = pin->GetOrAddReal("problem","d_n",1.0);
    auto &u0 = pmbp->phydro->u0;
    EOS_Data &eos = pmbp->phydro->peos->eos_data;
    Real gm1 = eos.gamma - 1.0;
    Real p0 = d_n*cs*cs/eos.gamma;
    // Set initial conditions
    par_for("pgen_turb", DevExeSpace(),0,(pmbp->nmb_thispack-1),ks,ke,js,je,is,ie,
    KOKKOS_LAMBDA(int m, int k, int j, int i) {
      u0(m,IDN,k,j,i) = d_n;
      u0(m,IM1,k,j,i) = 0.0;
      u0(m,IM2,k,j,i) = 0.0;
      u0(m,IM3,k,j,i) = 0.0;
      if (eos.is_ideal) {
        u0(m,IEN,k,j,i) = p0/gm1 +
           0.5*(SQR(u0(m,IM1,k,j,i)) + SQR(u0(m,IM2,k,j,i)) +
           SQR(u0(m,IM3,k,j,i)))/u0(m,IDN,k,j,i);
      }
    });
  }

  // Initialize MHD variables ---------------------------------
  if (pmbp->pmhd != nullptr) {
    Real d_n = pin->GetOrAddReal("problem","d_n",1.0);
    //Real B0 = cs*std::sqrt(2.0*d_n/beta);
    auto &u0 = pmbp->pmhd->u0;
    auto &b0 = pmbp->pmhd->b0;
    EOS_Data &eos = pmbp->pmhd->peos->eos_data;
    Real gm1 = eos.gamma - 1.0;
    Real p0 = d_n*cs*cs/eos.gamma;

    bool TGtest = pin->GetOrAddBoolean("problem","TGtest",false);
    if (TGtest) {
      std::cout << "TGtest is true" << std::endl;
    }
    Real lx, ly, lz, U0;
    lx = pmy_mesh_->mesh_size.x1max - pmy_mesh_->mesh_size.x1min;
    ly = pmy_mesh_->mesh_size.x2max - pmy_mesh_->mesh_size.x2min;
    lz = pmy_mesh_->mesh_size.x3max - pmy_mesh_->mesh_size.x3min;
    U0 = pin->GetOrAddReal("problem","U0", 1.0);
    int nx1 = indcs.nx1;
    int nx2 = indcs.nx2;
    int nx3 = indcs.nx3;
    Real pfactor = 4.0/3.0*std::sqrt(2.0/3.0);
    Real sqrt3 = std::sqrt(3.0);
    Real theta1 = std::asin(-sqrt3/2.0);
    Real phi1 = std::asin(sqrt3/2.0);
    Real phi2 = M_PI/2.0;
    if (TGtest) {
      // print the parameters if TGtest is true
      std::cout << "pfactor: " << pfactor << std::endl;
      std::cout << "theta1: " << theta1 << " phi1: " << phi1 << " phi2: " << phi2 << std::endl;
    }
    // Set initial conditions
    par_for("pgen_turb", DevExeSpace(),0,(pmbp->nmb_thispack-1),ks,ke,js,je,is,ie,
    KOKKOS_LAMBDA(int m, int k, int j, int i) {
      u0(m,IDN,k,j,i) = d_n;

      if (TGtest) {
        Real &x1min = size.d_view(m).x1min;
        Real &x1max = size.d_view(m).x1max;
        Real x1v = CellCenterX(i-is, nx1, x1min, x1max);

        Real &x2min = size.d_view(m).x2min;
        Real &x2max = size.d_view(m).x2max;
        Real x2v = CellCenterX(j-js, nx2, x2min, x2max);

        Real &x3min = size.d_view(m).x3min;
        Real &x3max = size.d_view(m).x3max;
        Real x3v = CellCenterX(k-ks, nx3, x3min, x3max);

        Real kx = 2.0*M_PI*x1v/lx;
        Real ky = 2.0*M_PI*x2v/ly;
        Real kz = 2.0*M_PI*x3v/lz;
        u0(m,IM1,k,j,i) =  d_n*U0*pfactor*(sin(kx+theta1)*cos(ky+phi1)*sin(kz+phi2) - 
                                          cos(kz+theta1)*sin(kx+phi1)*sin(ky+phi2));
        u0(m,IM2,k,j,i) =  d_n*U0*pfactor*(sin(ky+theta1)*cos(kz+phi1)*sin(kx+phi2) - 
                                          cos(kx+theta1)*sin(ky+phi1)*sin(kz+phi2));
        u0(m,IM3,k,j,i) =  d_n*U0*pfactor*(sin(kz+theta1)*cos(kx+phi1)*sin(ky+phi2) - 
                                          cos(ky+theta1)*sin(kz+phi1)*sin(kx+phi2));
      }else{
        u0(m,IM1,k,j,i) = 0.0;
        u0(m,IM2,k,j,i) = 0.0;
        u0(m,IM3,k,j,i) = 0.0;
      }

      // initialize B
      b0.x1f(m,k,j,i) = 0.0;
      b0.x2f(m,k,j,i) = B0*sin(t0);
      b0.x3f(m,k,j,i) = B0*cos(t0);
      if (i==ie) {b0.x1f(m,k,j,i+1) = 0.0;}
      if (j==je) {b0.x2f(m,k,j+1,i) = B0*sin(t0);}
      if (k==ke) {b0.x3f(m,k+1,j,i) = B0*cos(t0);}

      if (eos.is_ideal) {
        u0(m,IEN,k,j,i) = p0/gm1 + 0.5*B0*B0 + // fix contribution from dB
           0.5*(SQR(u0(m,IM1,k,j,i)) + SQR(u0(m,IM2,k,j,i)) +
           SQR(u0(m,IM3,k,j,i)))/u0(m,IDN,k,j,i);
      }
    });
  }

  // Initialize ion-neutral variables -------------------------
  if (pmbp->pionn != nullptr) {
    Real d_i = pin->GetOrAddReal("problem","d_i",1.0);
    Real d_n = pin->GetOrAddReal("problem","d_n",1.0);
    //Real B0 = cs*std::sqrt(2.0*(d_i+d_n)/beta);

    // MHD
    auto &u0 = pmbp->pmhd->u0;
    auto &b0 = pmbp->pmhd->b0;
    EOS_Data &eos = pmbp->pmhd->peos->eos_data;
    Real gm1 = eos.gamma - 1.0;
    Real p0 = d_i/eos.gamma; // TODO(@user): multiply by ionized density

    // Set initial conditions
    par_for("pgen_turb_mhd", DevExeSpace(),0,(pmbp->nmb_thispack-1),ks,ke,js,je,is,ie,
    KOKKOS_LAMBDA(int m, int k, int j, int i) {
      u0(m,IDN,k,j,i) = d_i;
      u0(m,IM1,k,j,i) = 0.0;
      u0(m,IM2,k,j,i) = 0.0;
      u0(m,IM3,k,j,i) = 0.0;

      // initialize B
      b0.x1f(m,k,j,i) = 0.0;
      b0.x2f(m,k,j,i) = 0.0;
      b0.x3f(m,k,j,i) = B0;
      if (i==ie) {b0.x1f(m,k,j,i+1) = 0.0;}
      if (j==je) {b0.x2f(m,k,j+1,i) = 0.0;}
      if (k==ke) {b0.x3f(m,k+1,j,i) = B0;}

      if (eos.is_ideal) {
        u0(m,IEN,k,j,i) = p0/gm1 + 0.5*B0*B0 + // fix contribution from dB
           0.5*(SQR(u0(m,IM1,k,j,i)) + SQR(u0(m,IM2,k,j,i)) +
           SQR(u0(m,IM3,k,j,i)))/u0(m,IDN,k,j,i);
      }
    });
    // Hydro
    auto &u0_ = pmbp->phydro->u0;
    EOS_Data &eos_ = pmbp->phydro->peos->eos_data;
    Real gm1_ = eos_.gamma - 1.0;
    Real p0_ = d_n/eos_.gamma; // TODO(@user): multiply by neutral density

    // Set initial conditions
    par_for("pgen_turb_hydro", DevExeSpace(),0,(pmbp->nmb_thispack-1),ks,ke,js,je,is,ie,
    KOKKOS_LAMBDA(int m, int k, int j, int i) {
      u0_(m,IDN,k,j,i) = d_n;
      u0_(m,IM1,k,j,i) = 0.0;
      u0_(m,IM2,k,j,i) = 0.0;
      u0_(m,IM3,k,j,i) = 0.0;
      if (eos_.is_ideal) {
        u0_(m,IEN,k,j,i) = p0_/gm1_ +
            0.5*(SQR(u0_(m,IM1,k,j,i)) + SQR(u0_(m,IM2,k,j,i)) +
            SQR(u0_(m,IM3,k,j,i)))/u0_(m,IDN,k,j,i);
      }
    });
  }

  return;
}


//----------------------------------------------------------------------------------------
// Function for computing history variables
void TurbulentHistory(HistoryData *pdata, Mesh *pm) {

  pdata->nhist = 13;
  // Volume-averaged quantities
  pdata->label[0] = "|∇⋅B|";      // absdivB
  pdata->label[1] = "|B2/ρ|";     // B2overRho
  pdata->label[2] = "|ρ2|";       // absrho2
  pdata->label[3] = "|∇⋅V|";      // absdivV
  pdata->label[4] = "|∇⋅V|2";     // absdivV2
  pdata->label[5] = "|p∇⋅V|";     // abspdivV
  pdata->label[6] = "|∇XV|2";     // curlU2
  pdata->label[7] = "|ρu⋅a|";     // rhoudota_OU
  pdata->label[8] = "Ms^2";       // Ms2
  pdata->label[9] = "Ma^2";       // Ma2
  pdata->label[10] = "u^2";       // u2
  pdata->label[11] = "eint";      // eint
  pdata->label[12] = "|∇xB2|";     // curlB2

  // Capture class variables for kernel
  auto &bcc =  pm->pmb_pack->pmhd->bcc0;
  auto &b = pm->pmb_pack->pmhd->b0;
  auto &w0_ = pm->pmb_pack->pmhd->w0;
  auto &u0 = pm->pmb_pack->pmhd->u0;
  auto &force = pm->pmb_pack->pturb->force;
  auto &size = pm->pmb_pack->pmb->mb_size;
  auto &meshsize = pm->pmb_pack->pmesh->mesh_size;
  int &nhist_ = pdata->nhist;
      
  Real gamma = pm->pmb_pack->pmhd->peos->eos_data.gamma;
  bool is_ideal = pm->pmb_pack->pmhd->peos->eos_data.is_ideal;
  Real cs0 = 0.0;
  if (!is_ideal) {
    cs0 = pm->pmb_pack->pmhd->peos->eos_data.iso_cs;
  }
  Real gm1 = gamma - 1.0;
  // Loop over all MeshBlocks in this pack
  auto &indcs = pm->pmb_pack->pmesh->mb_indcs;
  int is = indcs.is; int nx1 = indcs.nx1;
  int js = indcs.js; int nx2 = indcs.nx2;
  int ks = indcs.ks; int nx3 = indcs.nx3;
  const int nmkji = (pm->pmb_pack->nmb_thispack)*nx3*nx2*nx1;
  const int nkji = nx3*nx2*nx1;
  const int nji  = nx2*nx1;

  // Calculate the total volume of the mesh
  Real lx = pm->mesh_size.x1max - pm->mesh_size.x1min;
  Real ly = pm->mesh_size.x2max - pm->mesh_size.x2min;
  Real lz = pm->mesh_size.x3max - pm->mesh_size.x3min;
  Real TotalVol = lx*ly*lz;

  array_sum::GlobalSum sum_this_mb;
  Kokkos::parallel_reduce("HistSums",Kokkos::RangePolicy<>(DevExeSpace(), 0, nmkji),
  KOKKOS_LAMBDA(const int &idx, array_sum::GlobalSum &mb_sum) {
    // Compute n,k,j,i indices of thread
    int m = (idx)/nkji;
    int k = (idx - m*nkji)/nji;
    int j = (idx - m*nkji - k*nji)/nx1;
    int i = (idx - m*nkji - k*nji - j*nx1) + is;
    k += ks;
    j += js;

    Real dv = size.d_view(m).dx1*size.d_view(m).dx2*size.d_view(m).dx3;
    Real dx1 = size.d_view(m).dx1;
    Real dx2 = size.d_view(m).dx2;
    Real dx3 = size.d_view(m).dx3;
    Real dx_squared = dx1 * dx1;


    // MHD conserved variables:
    array_sum::GlobalSum hvars;

    // Common calculations for both MHD and Hydro
    
    // Calculate velocity magnitude squared
    Real u2 = (w0_(m,IVX,k,j,i)*w0_(m,IVX,k,j,i))
            + (w0_(m,IVY,k,j,i)*w0_(m,IVY,k,j,i))
            + (w0_(m,IVZ,k,j,i)*w0_(m,IVZ,k,j,i));
    
    // Calculate density squared
    Real rho = u0(m,IDN,k,j,i);
    Real rho2 = rho * rho;

    // Calculate Ms^2 (Mach number squared)
    Real p = 0.0;
    Real cs2 = 0.0;
    if (is_ideal) {
      p = w0_(m,IEN,k,j,i)*gm1;
      cs2 = gamma * p / rho;
    }else{
      p = rho*cs0*cs0;
      cs2 = cs0*cs0;
    }
    Real Ms2 = u2 / cs2;

    // Calculate divergence of velocity
    // Approximating face areas with cell dimensions for simplicity
    Real vel_divx = (w0_(m,IVX,k,j,i+1) - w0_(m,IVX,k,j,i-1)) / (2.0 * dx1);
    Real vel_divy = (w0_(m,IVY,k,j+1,i) - w0_(m,IVY,k,j-1,i)) / (2.0 * dx2);
    Real vel_divz = (w0_(m,IVZ,k+1,j,i) - w0_(m,IVZ,k-1,j,i)) / (2.0 * dx3);
    Real div_v = vel_divx + vel_divy + vel_divz;
    
    // Calculate ρu⋅a (density times velocity dot turbulent acceleration)
    Real rhoudota = rho * (w0_(m,IVX,k,j,i)*force(m,0,k,j,i) + 
                           w0_(m,IVY,k,j,i)*force(m,1,k,j,i) + 
                           w0_(m,IVZ,k,j,i)*force(m,2,k,j,i));

    // Calculate curl of velocity (vorticity)
    Real curl_vx = ((w0_(m,IVZ,k,j+1,i) - w0_(m,IVZ,k,j,i)) / (dx2)) - 
                   ((w0_(m,IVY,k+1,j,i) - w0_(m,IVY,k,j,i)) / (dx3));
    Real curl_vy = ((w0_(m,IVX,k+1,j,i) - w0_(m,IVX,k,j,i)) / (dx3)) -
                   ((w0_(m,IVZ,k,j,i+1) - w0_(m,IVZ,k,j,i)) / (dx1));
    Real curl_vz = ((w0_(m,IVY,k,j,i+1) - w0_(m,IVY,k,j,i)) / (dx1)) -
                   ((w0_(m,IVX,k,j+1,i) - w0_(m,IVX,k,j,i)) / (dx2));
    Real curl_v2 = curl_vx*curl_vx + curl_vy*curl_vy + curl_vz*curl_vz;

    // Calculate |∇xB^2|
    Real div_bx = (b.x1f(m,k,j,i+1) - b.x1f(m,k,j,i)) / dx1;
    Real div_by = (b.x2f(m,k,j+1,i) - b.x2f(m,k,j,i)) / dx2;
    Real div_bz = (b.x3f(m,k+1,j,i) - b.x3f(m,k,j,i)) / dx3;
    Real div_b = div_bx + div_by + div_bz;

    Real Jx = (b.x3f(m,k,j,i) - b.x3f(m,k,j-1,i)) / dx2 - 
              (b.x2f(m,k,j,i) - b.x2f(m,k-1,j,i)) / dx3;
    Real Jy = (b.x1f(m,k,j,i) - b.x1f(m,k-1,j,i)) / dx3 - 
              (b.x3f(m,k,j,i) - b.x3f(m,k,j,i-1)) / dx1;
    Real Jz = (b.x2f(m,k,j,i) - b.x2f(m,k,j,i-1)) / dx1 -
              (b.x1f(m,k,j,i) - b.x1f(m,k,j-1,i)) / dx2;
    Real J2 = Jx*Jx + Jy*Jy + Jz*Jz;

    // Calculate B^2
    Real b_sq = 0.25*( SQR(b.x1f(m,k,j,i) + b.x1f(m,k,j,i+1)) +
                        SQR(b.x2f(m,k,j,i) + b.x2f(m,k,j+1,i)) +
                        SQR(b.x3f(m,k,j,i) + b.x3f(m,k+1,j,i)) );
    
    // Calculate |B2/ρ|
    Real b2_over_rho = b_sq / rho;
    
    // Calculate Ma^2 (Alfvenic Mach number squared)
    Real va2 = b_sq / rho;
    Real Ma2 = u2 / va2;
    
    // Assign values to history variables for MHD
    hvars.the_array[0] = dv/TotalVol*abs(div_b);
    hvars.the_array[1] = dv/TotalVol*b2_over_rho;
    hvars.the_array[2] = dv/TotalVol*rho2;
    hvars.the_array[3] = dv/TotalVol*abs(div_v);
    hvars.the_array[4] = dv/TotalVol*div_v * div_v;
    hvars.the_array[5] = dv/TotalVol*p * div_v;
    hvars.the_array[6] = dv/TotalVol*curl_v2;
    hvars.the_array[7] = dv/TotalVol*rhoudota;
    hvars.the_array[8] = dv/TotalVol*Ms2;
    hvars.the_array[9] = dv/TotalVol*Ma2;
    hvars.the_array[10] = dv/TotalVol*u2;
    hvars.the_array[11] = dv/TotalVol*p/gm1;
    hvars.the_array[12] = dv/TotalVol*J2;

    // Sum into parallel reduce
    mb_sum += hvars;
  }, Kokkos::Sum<array_sum::GlobalSum>(sum_this_mb));
  Kokkos::fence();

  // Store data into hdata array
  for (int n=0; n<pdata->nhist; ++n) {
    pdata->hdata[n] = sum_this_mb.the_array[n];
  }
  return;
}

// how decide the refinement
void RefinementCondition(MeshBlockPack* pmbp) {
  Real time = pmbp->pmesh->time;
  auto& refine_flag = pmbp->pmesh->pmr->refine_flag;  // This is a DualView
  int mbs = pmbp->pmesh->gids_eachrank[global_variable::my_rank];
  int nmb = pmbp->nmb_thispack;

  // For DualView, we can directly access the host view
  // Mark the host view as modified
  refine_flag.modify_host();

  // Refine each meshblock
  if (time > AMR_t1 && RefineFlag1 == false) {
    RefineFlag1 = true;
    std::cout << "=---------RefinementCondition Entered-------------=" << std::endl;
    for (int m = 0; m < nmb; ++m) {
      refine_flag.h_view(m+mbs) = 1;
    }
  }
  if (time > AMR_t2 && RefineFlag2 == false) {
    RefineFlag2 = true;
    for (int m = 0; m < nmb; ++m) {
      refine_flag.h_view(m+mbs) = 1;
    }
  }

  // Sync the changes to device
  refine_flag.sync_device();

  return;
}