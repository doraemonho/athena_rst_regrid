//========================================================================================
// AthenaXXX astrophysical plasma code
// Copyright(C) 2020 James M. Stone <jmstone@ias.edu> and the Athena code team
// Licensed under the 3-clause BSD License (the "LICENSE")
//========================================================================================
//! \file kink.cpp
//  \brief Problem generator for non-realisitic kink instability

#include <iostream>
#include <sstream>

#include "athena.hpp"
#include "parameter_input.hpp"
#include "mesh/mesh.hpp"
#include "coordinates/coordinates.hpp"
#include "coordinates/cartesian_ks.hpp"
#include "coordinates/cell_locations.hpp"
#include "eos/eos.hpp"
#include "hydro/hydro.hpp"
#include "mhd/mhd.hpp"

#include <Kokkos_Random.hpp>

void DivB(HistoryData *pdata, Mesh *pm);
void RefinementCondition(MeshBlockPack* pmbp);
bool RefineFlag1 = false;
bool RefineFlag2 = false;
auto AMR_t1 = 0.0;
auto AMR_t2 = 0.0;

//----------------------------------------------------------------------------------------
//! \fn void MeshBlock::Kink_()
//  \brief Problem Generator for turbulence

void ProblemGenerator::UserProblem(ParameterInput *pin, const bool restart) {
  if (restart) return;
  user_ref_func  = RefinementCondition;

  MeshBlockPack *pmbp = pmy_mesh_->pmb_pack;
  auto &indcs = pmy_mesh_->mb_indcs;
  auto &size = pmbp->pmb->mb_size;

  if (pmbp->phydro == nullptr && pmbp->pmhd == nullptr) {
    std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__ << std::endl
      << "kink problem generator can only be run with Hydro and/or MHD, but no " 
      << "<hydro> or <mhd> block in input file" << std::endl;
    exit(EXIT_FAILURE);
  }
  
  AMR_t1 = pin->GetOrAddReal("problem","AMR_t1",0.0);
  AMR_t2 = pin->GetOrAddReal("problem","AMR_t2",0.0);

  // capture variables for kernel
  int &is = indcs.is; int &ie = indcs.ie;
  int &js = indcs.js; int &je = indcs.je;
  int &ks = indcs.ks; int &ke = indcs.ke;

  // Set up history output function 
  //user_hist_func = DivB;

  // Initialize MHD variables ---------------------------------
  if (pmbp->pmhd != nullptr) {

    Real  d0 = pin->GetReal("problem", "d0");
    Real  p0 = pin->GetReal("problem", "p0");
    Real  B0 = pin->GetReal("problem", "b0");
    Real  r0 = pin->GetReal("problem", "r0");
    Real   L = pin->GetReal("mesh", "x3max") - pin->GetReal("mesh", "x3min");
    
    // kink modes
    Real    amp = pin->GetReal("problem", "amp");
    int  nmodes = pin->GetOrAddInteger("problem", "nmodes",1);

    // create kokkos array from it
    Kokkos::View<int*> nmode_arr("nmode_arr", nmodes);
    auto nmode_arr_host = Kokkos::create_mirror_view(nmode_arr);  // Create a host mirror of the Kokkos View
    for (int i = 0; i < nmodes; ++i) {
        nmode_arr_host(i) = pin->GetInteger("problem", "n_" + std::to_string(i+1));
    }

    // Deep copy from host to device
    Kokkos::deep_copy(nmode_arr, nmode_arr_host);

    auto &u0 = pmbp->pmhd->u0;
    auto &b0 = pmbp->pmhd->b0;
    Real pi = 3.14159265;
    EOS_Data &eos = pmbp->pmhd->peos->eos_data;
    Real gm1 = eos.gamma - 1.0;

    // Set up vector potential
    int ncells1 = indcs.nx1 + 2*(indcs.ng);
    int ncells2 = (indcs.nx2 > 1)? (indcs.nx2 + 2*(indcs.ng)) : 1;
    int ncells3 = (indcs.nx3 > 1)? (indcs.nx3 + 2*(indcs.ng)) : 1;
    int nmb     = pmbp->nmb_thispack;
    DvceArray4D<Real> az;
    Kokkos::realloc(az, nmb, ncells3, ncells2, ncells1);

    // random perturbation
    Kokkos::Random_XorShift64_Pool<> rand_pool64(pmbp->gids);

    par_for("pgen_vector_potential", DevExeSpace(), 0,nmb-1,ks-1,ke+2,js-1,je+2,is-1,ie+2,
    KOKKOS_LAMBDA(int m, int k, int j, int i) {
      int nx1 = indcs.nx1;
      int nx2 = indcs.nx2;
      Real x1f     = LeftEdgeX(i-is, nx1, size.d_view(m).x1min, size.d_view(m).x1max);
      Real x2f     = LeftEdgeX(j-js, nx2, size.d_view(m).x2min, size.d_view(m).x2max);
      Real  r      = sqrt(x1f*x1f  +  x2f*x2f);
      az(m,k,j,i)  = -B0*r0/2.0*log(1.0 + SQR(r/r0));
    });

    // Set initial conditions
    par_for("pgen_kink", DevExeSpace(),0,(pmbp->nmb_thispack-1),ks,ke,js,je,is,ie,
    KOKKOS_LAMBDA(int m, int k, int j, int i) {

      // Compute face-centered fields & cell-centered coordinates
      int nx1 = indcs.nx1;
      int nx2 = indcs.nx2;
      int nx3 = indcs.nx3;
      Real x1f   = CellCenterX(i -is, nx1, size.d_view(m).x1min, size.d_view(m).x1max);
      Real x2f   = CellCenterX(j -js, nx2, size.d_view(m).x2min, size.d_view(m).x2max);
      Real x3f   = CellCenterX(k -ks, nx3, size.d_view(m).x3min, size.d_view(m).x3max);

      // compute variables in cylindrical coordinates
      Real r    = sqrt(x1f  *x1f   + x2f  *x2f);
      Real phi  = atan2(x2f, x1f);

      auto rand_gen = rand_pool64.get_state();  // get random number state this thread
      Real vr   = 0.0;
      for (int n_i=0; n_i<nmodes; ++n_i){
        Real n = nmode_arr(n_i);
        Real rval = 1.0 + 0.5*(rand_gen.frand() - 0.5);
        vr += rval*exp(-r/r0)*cos(phi)*sin(2.0*pi*n*x3f/L);
      }
      vr*=amp;
      rand_pool64.free_state(rand_gen);  // free state for use by other threads

      u0(m,IDN,k,j,i) = d0;
      u0(m,IM1,k,j,i) = d0*amp*vr*-sin(phi);
      u0(m,IM2,k,j,i) = d0*amp*vr*+cos(phi);
      u0(m,IM3,k,j,i) = 0;      
    });

    // initialize B from vector potenetial
    // B_z   = B0/(1 + (r/r0)^2)
    // B_phi = B0/(1 + (r/r0)^2)*(r/r0)
    // x,y,z -> r,theta,phi
    // We use vector potential approach (0,0,az) for Bphi to ensure ∇⋅B = 0
    // B_phi = B0/(1 + (r/r0)^2)*(r/r0) -> az = -B0*r0/2 *ln(1 + (r/r0)^2)
    par_for("pgen_kink_b1", DevExeSpace(),0,(pmbp->nmb_thispack-1),ks,ke,js,je,is,ie+1,
    KOKKOS_LAMBDA(int m, int k, int j, int i) {
      Real dx2   = size.d_view(m).dx2;
      b0.x1f(m,k,j,i) = (az(m,k,j+1,i) - az(m,k,j,i  ))/dx2;
    });

    par_for("pgen_kink_b2", DevExeSpace(),0,(pmbp->nmb_thispack-1),ks,ke,js,je+1,is,ie,
    KOKKOS_LAMBDA(int m, int k, int j, int i) {
      Real dx1   = size.d_view(m).dx1;
      b0.x2f(m,k,j,i) = (az(m,k,j  ,i) - az(m,k,j,i+1))/dx1;
    });

    par_for("pgen_kink_b3", DevExeSpace(),0,(pmbp->nmb_thispack-1),ks,ke+1,js,je,is,ie,
    KOKKOS_LAMBDA(int m, int k, int j, int i) {
      Real x1f   = LeftEdgeX(i -is, indcs.nx1, size.d_view(m).x1min, size.d_view(m).x1max);
      Real x2f   = LeftEdgeX(j -js, indcs.nx2, size.d_view(m).x2min, size.d_view(m).x2max);
      Real r    = sqrt(x1f  *x1f   + x2f  *x2f);
      // no r dependence in x3 for z component of B
      b0.x3f(m,k,j,i) = B0/(1 + SQR(r/r0));
    });

    par_for("pgen_kink", DevExeSpace(),0,(pmbp->nmb_thispack-1),ks,ke,js,je,is,ie,
    KOKKOS_LAMBDA(int m, int k, int j, int i) {
      if (eos.is_ideal) {
        u0(m,IEN,k,j,i) = p0/gm1 + 
                    0.5*(SQR(0.5*(b0.x1f(m,k,j,i) + b0.x1f(m,k,j,i+1))) +
                          SQR(0.5*(b0.x2f(m,k,j,i) + b0.x2f(m,k,j+1,i))) +
                          SQR(0.5*(b0.x3f(m,k,j,i) + b0.x3f(m,k+1,j,i)))) +
                    0.5*(SQR(u0(m,IM1,k,j,i)) + SQR(u0(m,IM2,k,j,i)) +
                          SQR(u0(m,IM3,k,j,i)))/u0(m,IDN,k,j,i);
      }
      if (fabs(b0.x1f(m,k,j,i)) > 5.0){
        printf("b.x1f = %f\n", b0.x1f(m,k,j,i));}
      if (fabs(b0.x2f(m,k,j,i)) > 5.0){
        printf("b.x2f = %f\n", b0.x2f(m,k,j,i));}

      if ( (i==ie) && (fabs(b0.x1f(m,k,j,i+1)) > 5.0)){
        printf("b.x1f = %f\n", b0.x1f(m,k,j,i+1));}
      if ( (j==je) && (fabs(b0.x2f(m,k,j+1,i)) > 5.0)){
        printf("b.x2f = %f\n", b0.x2f(m,k,j+1,i));}
    });
  }
  return;
}

void DivB(HistoryData *pdata, Mesh *pm) {

  pdata->nhist = 1;
  // div B in math labal
  pdata->label[0] = "<|∇⋅B|>";

  // capture class variabels for kernel
  auto &bcc = pm->pmb_pack->pmhd->bcc0;
  auto &b = pm->pmb_pack->pmhd->b0;
  auto &w0_ = pm->pmb_pack->pmhd->w0;
  auto &size = pm->pmb_pack->pmb->mb_size;
  int &nhist_ = pdata->nhist;

  // loop over all MeshBlocks in this pack
  auto &indcs = pm->pmb_pack->pmesh->mb_indcs;
  int is = indcs.is; int nx1 = indcs.nx1;
  int js = indcs.js; int nx2 = indcs.nx2;
  int ks = indcs.ks; int nx3 = indcs.nx3;
  const int nmkji = (pm->pmb_pack->nmb_thispack)*nx3*nx2*nx1;
  const int nkji = nx3*nx2*nx1;
  const int nji  = nx2*nx1;
  const int Nx = indcs.nx1;
  const int Ny = (indcs.nx2 > 1)? indcs.nx2 : 1;
  const int Nz = (indcs.nx3 > 1)? indcs.nx3 : 1;

  array_sum::GlobalSum sum_this_mb;
  Kokkos::parallel_reduce("HistSums",Kokkos::RangePolicy<>(DevExeSpace(), 0, nmkji),
  KOKKOS_LAMBDA(const int &idx, array_sum::GlobalSum &mb_sum) {
    // compute n,k,j,i indices of thread
    int m = (idx)/nkji;
    int k = (idx - m*nkji)/nji;
    int j = (idx - m*nkji - k*nji)/nx1;
    int i = (idx - m*nkji - k*nji - j*nx1) + is;
    k += ks;
    j += js;

    Real dx1 = size.d_view(m).dx1;
    Real dx2 = size.d_view(m).dx2;
    Real dx3 = size.d_view(m).dx3;

    Real face1 = dx2*dx3;
    Real face2 = dx1*dx3;
    Real face3 = dx1*dx2;
    
    array_sum::GlobalSum hvars;

    hvars.the_array[0] += fabs(  face1*(b.x1f(m,k,j,i+1) - b.x1f(m,k,j,i))
                              +  face2*(b.x2f(m,k,j+1,i) - b.x2f(m,k,j,i))
                              +  face3*(b.x3f(m,k+1,j,i) - b.x3f(m,k,j,i)))/Nx/Ny/Nz;

    // fill rest of the_array with zeros, if nhist < NHISTORY_VARIABLES
    for (int n=nhist_; n<NHISTORY_VARIABLES; ++n) {
      hvars.the_array[n] = 0.0;
    }

    // sum into parallel reduce
    mb_sum += hvars;

  }, Kokkos::Sum<array_sum::GlobalSum>(sum_this_mb));
  Kokkos::fence();

  // store data into hdata array
  for (int n=0; n<pdata->nhist; ++n) {
    pdata->hdata[n] = sum_this_mb.the_array[n];
  }

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
    //std::cout << "=---------RefinementCondition Entered-------------=" << std::endl;
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