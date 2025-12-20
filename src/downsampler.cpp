#include "downsampler.hpp"

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <cstring>
#include <iostream>
#include <vector>

#include "binary_writer.hpp"
#include "mpi_distribution.hpp"
#include "parameter_parser.hpp"

#if MPI_PARALLEL_ENABLED
#include <mpi.h>
#endif

namespace {
struct ChunkHeader {
  std::int32_t coarse_gid;
  std::int32_t child_index;
};

constexpr int kVarDens = 0;
constexpr int kVarVelx = 1;
constexpr int kVarVely = 2;
constexpr int kVarVelz = 3;
constexpr int kVarEint = 4;  // present only for ideal EOS (AthenaK mhd_w/mhd_w_bcc)

int FindRankForGid(const std::vector<int>& gids_eachrank,
                   const std::vector<int>& nmb_eachrank, int gid) {
  auto it = std::upper_bound(gids_eachrank.begin(), gids_eachrank.end(), gid);
  int rank = static_cast<int>(it - gids_eachrank.begin()) - 1;
  rank = std::max(rank, 0);
  if (rank >= static_cast<int>(gids_eachrank.size())) {
    return static_cast<int>(gids_eachrank.size()) - 1;
  }
  const int start = gids_eachrank[rank];
  const int end = start + nmb_eachrank[rank];
  if (gid < start || gid >= end) {
    for (int r = 0; r < static_cast<int>(gids_eachrank.size()); ++r) {
      const int s = gids_eachrank[r];
      const int e = s + nmb_eachrank[r];
      if (gid >= s && gid < e) {
        return r;
      }
    }
  }
  return rank;
}
}  // namespace

Downsampler::Downsampler(RestartReader& reader, int downsample_factor)
    : reader_(reader), downsample_factor_(downsample_factor) {
  const auto& mesh_indcs = reader_.GetMeshIndcs();
  multi_d_ = (mesh_indcs.nx2 > 1);
  three_d_ = (mesh_indcs.nx3 > 1);

  if (downsample_factor_ == 2) {
    downsample_levels_ = 1;
  } else if (downsample_factor_ == 4) {
    downsample_levels_ = 2;
  } else {
    downsample_levels_ = 0;
    nfine_per_coarse_ = 0;
    return;
  }

  nfine_per_coarse_ = downsample_factor_;
  if (multi_d_) {
    nfine_per_coarse_ *= downsample_factor_;
  }
  if (three_d_) {
    nfine_per_coarse_ *= downsample_factor_;
  }
}

int Downsampler::ChildOffsetX1(int child_index) const {
  const int dims = three_d_ ? 3 : (multi_d_ ? 2 : 1);
  int ox = 0;
  for (int bit = 0; bit < downsample_levels_; ++bit) {
    const int pos = bit * dims;
    ox |= ((child_index >> pos) & 1) << bit;
  }
  return ox;
}

int Downsampler::ChildOffsetX2(int child_index) const {
  const int dims = three_d_ ? 3 : (multi_d_ ? 2 : 1);
  int ox = 0;
  for (int bit = 0; bit < downsample_levels_; ++bit) {
    const int pos = bit * dims + 1;
    ox |= ((child_index >> pos) & 1) << bit;
  }
  return ox;
}

int Downsampler::ChildOffsetX3(int child_index) const {
  const int dims = three_d_ ? 3 : (multi_d_ ? 2 : 1);
  int ox = 0;
  for (int bit = 0; bit < downsample_levels_; ++bit) {
    const int pos = bit * dims + 2;
    ox |= ((child_index >> pos) & 1) << bit;
  }
  return ox;
}

LogicalLocation Downsampler::CoarseLocation(const LogicalLocation& fine) const {
  LogicalLocation coarse;
  coarse.level = fine.level - downsample_levels_;
  coarse.lx1 = fine.lx1 / downsample_factor_;
  coarse.lx2 = fine.lx2 / downsample_factor_;
  coarse.lx3 = fine.lx3 / downsample_factor_;
  return coarse;
}

bool Downsampler::SameLocation(const LogicalLocation& a, const LogicalLocation& b) {
  return a.level == b.level && a.lx1 == b.lx1 && a.lx2 == b.lx2 && a.lx3 == b.lx3;
}

std::uint64_t Downsampler::CellIndex(int i, int j, int k, int nx1, int nx2) {
  return static_cast<std::uint64_t>((k * nx2 + j) * nx1 + i);
}

std::uint64_t Downsampler::FaceIndexX1(int i, int j, int k, int nx1p1, int nx2) {
  return static_cast<std::uint64_t>((k * nx2 + j) * nx1p1 + i);
}

std::uint64_t Downsampler::FaceIndexX2(int i, int j, int k, int nx1, int nx2p1) {
  return static_cast<std::uint64_t>((k * nx2p1 + j) * nx1 + i);
}

std::uint64_t Downsampler::FaceIndexX3(int i, int j, int k, int nx1, int nx2) {
  return static_cast<std::uint64_t>((k * nx2 + j) * nx1 + i);
}

bool Downsampler::ParseEOS() {
  eos_ = "ideal";
  dfloor_ = static_cast<Real>(FLT_MIN);
  pfloor_ = static_cast<Real>(FLT_MIN);
  tfloor_ = static_cast<Real>(FLT_MIN);
  sfloor_ = static_cast<Real>(FLT_MIN);
  tmax_ = static_cast<Real>(FLT_MAX);

  const std::string& params = reader_.GetParameterString();
  auto eos_opt = ParameterParser::ExtractStringInBlock(params, "mhd", "eos");
  if (eos_opt.has_value()) {
    eos_ = *eos_opt;
  }

  auto gamma_opt = ParameterParser::ExtractRealInBlock(params, "mhd", "gamma");
  if (!gamma_opt.has_value()) {
    gamma_opt = ParameterParser::ExtractRealInBlock(params, "hydro", "gamma");
  }
  if (gamma_opt.has_value()) {
    gamma_ = *gamma_opt;
  }

  auto cs_opt = ParameterParser::ExtractRealInBlock(params, "mhd", "iso_sound_speed");
  if (cs_opt.has_value()) {
    iso_sound_speed_ = *cs_opt;
  }

  auto dfloor_opt = ParameterParser::ExtractRealInBlock(params, "mhd", "dfloor");
  if (!dfloor_opt.has_value()) {
    dfloor_opt = ParameterParser::ExtractRealInBlock(params, "hydro", "dfloor");
  }
  if (dfloor_opt.has_value()) {
    dfloor_ = *dfloor_opt;
  }

  auto pfloor_opt = ParameterParser::ExtractRealInBlock(params, "mhd", "pfloor");
  if (!pfloor_opt.has_value()) {
    pfloor_opt = ParameterParser::ExtractRealInBlock(params, "hydro", "pfloor");
  }
  if (pfloor_opt.has_value()) {
    pfloor_ = *pfloor_opt;
  }

  auto tfloor_opt = ParameterParser::ExtractRealInBlock(params, "mhd", "tfloor");
  if (!tfloor_opt.has_value()) {
    tfloor_opt = ParameterParser::ExtractRealInBlock(params, "hydro", "tfloor");
  }
  if (tfloor_opt.has_value()) {
    tfloor_ = *tfloor_opt;
  }

  auto sfloor_opt = ParameterParser::ExtractRealInBlock(params, "mhd", "sfloor");
  if (!sfloor_opt.has_value()) {
    sfloor_opt = ParameterParser::ExtractRealInBlock(params, "hydro", "sfloor");
  }
  if (sfloor_opt.has_value()) {
    sfloor_ = *sfloor_opt;
  }

  auto tmax_opt = ParameterParser::ExtractRealInBlock(params, "mhd", "tmax");
  if (!tmax_opt.has_value()) {
    tmax_opt = ParameterParser::ExtractRealInBlock(params, "hydro", "tmax");
  }
  if (tmax_opt.has_value()) {
    tmax_ = *tmax_opt;
  }
  return true;
}

bool Downsampler::EnsureCoarseSetup() {
  if (coarse_ready_) {
    return true;
  }
  ComputeCoarseMeshConfig();
  if (!BuildCoarseLogicalLocations()) {
    return false;
  }
  coarse_ready_ = true;
  return true;
}

void Downsampler::ComputeCoarseMeshConfig() {
  const RegionSize& fine_mesh_size = reader_.GetMeshSize();
  const RegionIndcs& fine_mesh_indcs = reader_.GetMeshIndcs();
  const RegionIndcs& fine_mb_indcs = reader_.GetMBIndcs();

  coarse_mesh_size_ = fine_mesh_size;
  coarse_mesh_size_.dx1 = fine_mesh_size.dx1 * static_cast<Real>(downsample_factor_);
  coarse_mesh_size_.dx2 =
      fine_mesh_size.dx2 * static_cast<Real>(multi_d_ ? downsample_factor_ : 1);
  coarse_mesh_size_.dx3 =
      fine_mesh_size.dx3 * static_cast<Real>(three_d_ ? downsample_factor_ : 1);

  coarse_mesh_indcs_ = fine_mesh_indcs;
  coarse_mesh_indcs_.nx1 = fine_mesh_indcs.nx1 / downsample_factor_;
  coarse_mesh_indcs_.nx2 = multi_d_ ? (fine_mesh_indcs.nx2 / downsample_factor_) : 1;
  coarse_mesh_indcs_.nx3 = three_d_ ? (fine_mesh_indcs.nx3 / downsample_factor_) : 1;

  coarse_mesh_indcs_.is = fine_mesh_indcs.ng;
  coarse_mesh_indcs_.ie = coarse_mesh_indcs_.is + coarse_mesh_indcs_.nx1 - 1;
  coarse_mesh_indcs_.js = fine_mesh_indcs.ng;
  coarse_mesh_indcs_.je = coarse_mesh_indcs_.js + coarse_mesh_indcs_.nx2 - 1;
  coarse_mesh_indcs_.ks = fine_mesh_indcs.ng;
  coarse_mesh_indcs_.ke = coarse_mesh_indcs_.ks + coarse_mesh_indcs_.nx3 - 1;

  coarse_mb_indcs_ = fine_mb_indcs;

  coarse_nmb_total_ = reader_.GetNMBTotal() / nfine_per_coarse_;
  coarse_root_level_ = reader_.GetRootLevel() - downsample_levels_;
}

bool Downsampler::BuildCoarseLogicalLocations() {
  const auto& fine_llocs = reader_.GetLlocEachMB();
  if (fine_llocs.empty()) {
    return false;
  }

  const int fine_nmb_total = reader_.GetNMBTotal();
  if (fine_nmb_total % nfine_per_coarse_ != 0) {
    if (reader_.GetMyRank() == 0) {
      std::cerr << "ERROR: nmb_total=" << fine_nmb_total
                << " not divisible by nfine_per_coarse=" << nfine_per_coarse_
                << std::endl;
    }
    return false;
  }

  if (reader_.GetRootLevel() < downsample_levels_) {
    if (reader_.GetMyRank() == 0) {
      std::cerr << "ERROR: root_level=" << reader_.GetRootLevel()
                << " too small for " << downsample_factor_ << "x downsampling"
                << std::endl;
    }
    return false;
  }

  coarse_lloc_eachmb_.resize(coarse_nmb_total_);

  for (int cg = 0; cg < coarse_nmb_total_; ++cg) {
    const int fg0 = cg * nfine_per_coarse_;
    const LogicalLocation parent = CoarseLocation(fine_llocs[fg0]);
    coarse_lloc_eachmb_[cg] = parent;

    for (int c = 0; c < nfine_per_coarse_; ++c) {
      const int fg = fg0 + c;
      const LogicalLocation& child = fine_llocs[fg];
      if (child.level != fine_llocs[fg0].level) {
        if (reader_.GetMyRank() == 0) {
          std::cerr << "ERROR: Mixed levels within child group at fine gid " << fg0
                    << std::endl;
        }
        return false;
      }

      const LogicalLocation computed_parent = CoarseLocation(child);
      if (!SameLocation(computed_parent, parent)) {
        if (reader_.GetMyRank() == 0) {
          std::cerr << "ERROR: Non-parent grouping at coarse gid " << cg
                    << " fine gid " << fg << std::endl;
        }
        return false;
      }

      const int ox1 = ChildOffsetX1(c);
      const int ox2 = multi_d_ ? ChildOffsetX2(c) : 0;
      const int ox3 = three_d_ ? ChildOffsetX3(c) : 0;

      if (child.lx1 != downsample_factor_ * parent.lx1 + ox1
          || child.lx2 != downsample_factor_ * parent.lx2 + ox2
          || child.lx3 != downsample_factor_ * parent.lx3 + ox3) {
        if (reader_.GetMyRank() == 0) {
          std::cerr << "ERROR: Fine gids not in expected child order at coarse gid " << cg
                    << std::endl;
        }
        return false;
      }
    }
  }

  return true;
}

bool Downsampler::ComputeDownsampledChunk(int local_fine_mb,
                                         std::vector<float>* chunk) const {
  const auto* phys = reader_.GetPhysicsReader();
  if (phys == nullptr) {
    return false;
  }

  const auto& config = reader_.GetPhysicsConfig();
  if (!config.has_mhd) {
    return false;
  }

  const bool has_energy = (config.nmhd >= 5);
  const bool ideal_eos = (eos_ == "ideal");
  const bool isothermal_eos = (eos_ == "isothermal");

  if (ideal_eos && !has_energy) {
    return false;
  }
  if (!ideal_eos && !isothermal_eos) {
    return false;
  }

  // Match AthenaK `mhd_w_bcc` bin output ordering:
  // - ideal EOS: dens, velx, vely, velz, eint, bcc1, bcc2, bcc3  (8 vars)
  // - isothermal EOS: dens, velx, vely, velz, bcc1, bcc2, bcc3   (7 vars)
  const int bcc_index_base = ideal_eos ? 5 : 4;
  const int out_vars = ideal_eos ? 8 : 7;

  const auto& mb_indcs = reader_.GetMBIndcs();
  const int ng = mb_indcs.ng;
  const int nx1 = mb_indcs.nx1;
  const int nx2 = (mb_indcs.nx2 > 1) ? mb_indcs.nx2 : 1;
  const int nx3 = (mb_indcs.nx3 > 1) ? mb_indcs.nx3 : 1;

  if ((nx1 % downsample_factor_) != 0
      || (multi_d_ && (nx2 % downsample_factor_) != 0)
      || (three_d_ && (nx3 % downsample_factor_) != 0)) {
    return false;
  }

  const int nx1_sub = nx1 / downsample_factor_;
  const int nx2_sub = multi_d_ ? (nx2 / downsample_factor_) : 1;
  const int nx3_sub = three_d_ ? (nx3 / downsample_factor_) : 1;

  const int nout1 = nx1 + 2 * ng;
  const int nout2 = multi_d_ ? (nx2 + 2 * ng) : 1;
  const int nout3 = three_d_ ? (nx3 + 2 * ng) : 1;

  const auto& mhd = phys->GetMHDData().at(local_fine_mb);
  const auto& b1f = phys->GetMHDB1f().at(local_fine_mb);
  const auto& b2f = phys->GetMHDB2f().at(local_fine_mb);
  const auto& b3f = phys->GetMHDB3f().at(local_fine_mb);

  const std::uint64_t sub_cells =
      static_cast<std::uint64_t>(nx1_sub) * nx2_sub * nx3_sub;
  chunk->assign(static_cast<std::size_t>(out_vars) * sub_cells, 0.0F);

  std::vector<Real> b1c(static_cast<std::size_t>(nx3_sub) * nx2_sub * (nx1_sub + 1), 0.0);
  std::vector<Real> b2c(static_cast<std::size_t>(nx3_sub) * (nx2_sub + 1) * nx1_sub, 0.0);
  std::vector<Real> b3c(static_cast<std::size_t>(nx3_sub + 1) * nx2_sub * nx1_sub, 0.0);

  const int nout1p1 = nout1 + 1;
  const int nout2p1 = nout2 + 1;
  const int nout3p1 = nout3 + 1;

  for (int k = 0; k < nx3_sub; ++k) {
    for (int j = 0; j < nx2_sub; ++j) {
      for (int i = 0; i <= nx1_sub; ++i) {
        const int fi = ng + downsample_factor_ * i;
        const int fj = multi_d_ ? (ng + downsample_factor_ * j) : 0;
        const int fk = three_d_ ? (ng + downsample_factor_ * k) : 0;

        Real sum = 0.0;
        if (three_d_) {
          for (int kk = 0; kk < downsample_factor_; ++kk) {
            for (int jj = 0; jj < downsample_factor_; ++jj) {
              const std::uint64_t idx =
                  FaceIndexX1(fi, fj + jj, fk + kk, nout1p1, nout2);
              sum += b1f[idx];
            }
          }
          sum *= 1.0
                 / static_cast<Real>(downsample_factor_ * downsample_factor_);
        } else if (multi_d_) {
          for (int jj = 0; jj < downsample_factor_; ++jj) {
            sum += b1f[FaceIndexX1(fi, fj + jj, 0, nout1p1, nout2)];
          }
          sum *= 1.0 / static_cast<Real>(downsample_factor_);
        } else {
          sum = b1f[FaceIndexX1(fi, 0, 0, nout1p1, 1)];
        }

        const std::uint64_t cidx = FaceIndexX1(i, j, k, nx1_sub + 1, nx2_sub);
        b1c[cidx] = sum;
      }
    }
  }

  for (int k = 0; k < nx3_sub; ++k) {
    for (int j = 0; j <= nx2_sub; ++j) {
      for (int i = 0; i < nx1_sub; ++i) {
        const int fi = ng + downsample_factor_ * i;
        const int fj = multi_d_ ? (ng + downsample_factor_ * j) : 0;
        const int fk = three_d_ ? (ng + downsample_factor_ * k) : 0;

        Real sum = 0.0;
        if (three_d_) {
          for (int kk = 0; kk < downsample_factor_; ++kk) {
            for (int ii = 0; ii < downsample_factor_; ++ii) {
              const std::uint64_t idx =
                  FaceIndexX2(fi + ii, fj, fk + kk, nout1, nout2p1);
              sum += b2f[idx];
            }
          }
          sum *= 1.0
                 / static_cast<Real>(downsample_factor_ * downsample_factor_);
        } else if (multi_d_) {
          for (int ii = 0; ii < downsample_factor_; ++ii) {
            sum += b2f[FaceIndexX2(fi + ii, fj, 0, nout1, nout2p1)];
          }
          sum *= 1.0 / static_cast<Real>(downsample_factor_);
        } else {
          for (int ii = 0; ii < downsample_factor_; ++ii) {
            sum += b2f[FaceIndexX2(fi + ii, 0, 0, nout1, 2)];
          }
          sum *= 1.0 / static_cast<Real>(downsample_factor_);
        }

        const std::uint64_t cidx = FaceIndexX2(i, j, k, nx1_sub, nx2_sub + 1);
        b2c[cidx] = sum;
      }
    }
  }

  for (int k = 0; k <= nx3_sub; ++k) {
    for (int j = 0; j < nx2_sub; ++j) {
      for (int i = 0; i < nx1_sub; ++i) {
        const int fi = ng + downsample_factor_ * i;
        const int fj = multi_d_ ? (ng + downsample_factor_ * j) : 0;
        const int fk = three_d_ ? (ng + downsample_factor_ * k) : 0;

        Real sum = 0.0;
        if (three_d_) {
          for (int jj = 0; jj < downsample_factor_; ++jj) {
            for (int ii = 0; ii < downsample_factor_; ++ii) {
              const std::uint64_t idx =
                  FaceIndexX3(fi + ii, fj + jj, fk, nout1, nout2);
              sum += b3f[idx];
            }
          }
          sum *= 1.0
                 / static_cast<Real>(downsample_factor_ * downsample_factor_);
        } else if (multi_d_) {
          for (int jj = 0; jj < downsample_factor_; ++jj) {
            for (int ii = 0; ii < downsample_factor_; ++ii) {
              sum += b3f[FaceIndexX3(fi + ii, fj + jj, 0, nout1, nout2)];
            }
          }
          sum *= 1.0
                 / static_cast<Real>(downsample_factor_ * downsample_factor_);
        } else {
          for (int ii = 0; ii < downsample_factor_; ++ii) {
            sum += b3f[FaceIndexX3(fi + ii, 0, 0, nout1, nout2)];
          }
          sum *= 1.0 / static_cast<Real>(downsample_factor_);
        }

        const std::uint64_t cidx = FaceIndexX3(i, j, k, nx1_sub, nx2_sub);
        b3c[cidx] = sum;
      }
    }
  }

  const int ncells = nout1 * nout2 * nout3;

  for (int k = 0; k < nx3_sub; ++k) {
    for (int j = 0; j < nx2_sub; ++j) {
      for (int i = 0; i < nx1_sub; ++i) {
        const int fi = ng + downsample_factor_ * i;
        const int fj = multi_d_ ? (ng + downsample_factor_ * j) : 0;
        const int fk = three_d_ ? (ng + downsample_factor_ * k) : 0;

        Real rho_sum = 0.0;
        Real m1_sum = 0.0;
        Real m2_sum = 0.0;
        Real m3_sum = 0.0;
        Real eng_sum = 0.0;

        auto read_var = [&](int v, std::uint64_t idx) -> Real {
          return mhd[static_cast<std::uint64_t>(v) * ncells + idx];
        };

        const int jj_max = multi_d_ ? downsample_factor_ : 1;
        const int kk_max = three_d_ ? downsample_factor_ : 1;
        for (int kk = 0; kk < kk_max; ++kk) {
          for (int jj = 0; jj < jj_max; ++jj) {
            for (int ii = 0; ii < downsample_factor_; ++ii) {
              const std::uint64_t idx =
                  CellIndex(fi + ii, fj + jj, fk + kk, nout1, nout2);
              rho_sum += read_var(0, idx);
              m1_sum += read_var(1, idx);
              m2_sum += read_var(2, idx);
              m3_sum += read_var(3, idx);
              if (has_energy) {
                eng_sum += read_var(4, idx);
              }
            }
          }
        }

        const Real inv_vol = 1.0 / static_cast<Real>(nfine_per_coarse_);
        const Real rho = inv_vol * rho_sum;
        const Real m1 = inv_vol * m1_sum;
        const Real m2 = inv_vol * m2_sum;
        const Real m3 = inv_vol * m3_sum;
        const Real eng = inv_vol * eng_sum;

        const std::uint64_t b1_left = FaceIndexX1(i, j, k, nx1_sub + 1, nx2_sub);
        const std::uint64_t b1_right = FaceIndexX1(i + 1, j, k, nx1_sub + 1, nx2_sub);
        const std::uint64_t b2_bot = FaceIndexX2(i, j, k, nx1_sub, nx2_sub + 1);
        const std::uint64_t b2_top = FaceIndexX2(i, j + 1, k, nx1_sub, nx2_sub + 1);
        const std::uint64_t b3_bot = FaceIndexX3(i, j, k, nx1_sub, nx2_sub);
        const std::uint64_t b3_top = FaceIndexX3(i, j, k + 1, nx1_sub, nx2_sub);

        const Real bx = 0.5 * (b1c[b1_left] + b1c[b1_right]);
        const Real by = 0.5 * (b2c[b2_bot] + b2c[b2_top]);
        const Real bz = 0.5 * (b3c[b3_bot] + b3c[b3_top]);

        Real rho_out = rho;
        if (rho_out < dfloor_) {
          rho_out = dfloor_;
        }

        Real v1 = 0.0;
        Real v2 = 0.0;
        Real v3 = 0.0;
        if (rho_out != 0.0) {
          v1 = m1 / rho_out;
          v2 = m2 / rho_out;
          v3 = m3 / rho_out;
        }

        Real eint_out = 0.0;
        if (ideal_eos) {
          const Real gm1 = gamma_ - 1.0;
          if (gm1 <= 0.0) {
            return false;
          }
          const Real di = 1.0 / rho_out;
          const Real ke = 0.5 * di * (m1 * m1 + m2 * m2 + m3 * m3);
          const Real me = 0.5 * (bx * bx + by * by + bz * bz);

          Real eint = eng - ke - me;
          const Real efloor = pfloor_ / gm1;
          if (eint < efloor) {
            eint = efloor;
          }
          if (gm1 * eint * di < tfloor_) {
            eint = rho_out * tfloor_ / gm1;
          }
          if (gm1 * eint * di > tmax_) {
            eint = rho_out * tmax_ / gm1;
          }

          const Real spe_over_eps = gm1 / std::pow(rho_out, gm1);
          const Real spe = spe_over_eps * eint * di;
          if (spe <= sfloor_) {
            eint = rho_out * sfloor_ / spe_over_eps;
          }
          eint_out = eint;
        }

        const std::uint64_t cidx = CellIndex(i, j, k, nx1_sub, nx2_sub);
        (*chunk)[static_cast<std::size_t>(kVarDens) * sub_cells + cidx] =
            static_cast<float>(rho_out);
        (*chunk)[static_cast<std::size_t>(kVarVelx) * sub_cells + cidx] =
            static_cast<float>(v1);
        (*chunk)[static_cast<std::size_t>(kVarVely) * sub_cells + cidx] =
            static_cast<float>(v2);
        (*chunk)[static_cast<std::size_t>(kVarVelz) * sub_cells + cidx] =
            static_cast<float>(v3);
        if (ideal_eos) {
          (*chunk)[static_cast<std::size_t>(kVarEint) * sub_cells + cidx] =
              static_cast<float>(eint_out);
        }
        (*chunk)[static_cast<std::size_t>(bcc_index_base + 0) * sub_cells + cidx] =
            static_cast<float>(bx);
        (*chunk)[static_cast<std::size_t>(bcc_index_base + 1) * sub_cells + cidx] =
            static_cast<float>(by);
        (*chunk)[static_cast<std::size_t>(bcc_index_base + 2) * sub_cells + cidx] =
            static_cast<float>(bz);
      }
    }
  }

  return true;
}

bool Downsampler::ComputeDownsampledTurbulenceChunk(int local_fine_mb,
                                                    std::vector<float>* chunk) const {
  const auto* phys = reader_.GetPhysicsReader();
  if (phys == nullptr) {
    return false;
  }

  const auto& config = reader_.GetPhysicsConfig();
  if (!config.has_turbulence) {
    return false;
  }
  const int out_vars = config.nforce;
  if (out_vars <= 0) {
    return false;
  }

  const auto& mb_indcs = reader_.GetMBIndcs();
  const int ng = mb_indcs.ng;
  const int nx1 = mb_indcs.nx1;
  const int nx2 = (mb_indcs.nx2 > 1) ? mb_indcs.nx2 : 1;
  const int nx3 = (mb_indcs.nx3 > 1) ? mb_indcs.nx3 : 1;

  if ((nx1 % downsample_factor_) != 0
      || (multi_d_ && (nx2 % downsample_factor_) != 0)
      || (three_d_ && (nx3 % downsample_factor_) != 0)) {
    return false;
  }

  const int nx1_sub = nx1 / downsample_factor_;
  const int nx2_sub = multi_d_ ? (nx2 / downsample_factor_) : 1;
  const int nx3_sub = three_d_ ? (nx3 / downsample_factor_) : 1;

  const int nout1 = nx1 + 2 * ng;
  const int nout2 = multi_d_ ? (nx2 + 2 * ng) : 1;
  const int nout3 = three_d_ ? (nx3 + 2 * ng) : 1;
  const int ncells = nout1 * nout2 * nout3;

  const auto& turb = phys->GetTurbData().at(local_fine_mb);
  const std::uint64_t expected =
      static_cast<std::uint64_t>(out_vars) * static_cast<std::uint64_t>(ncells);
  if (turb.size() != expected) {
    return false;
  }

  const std::uint64_t sub_cells =
      static_cast<std::uint64_t>(nx1_sub) * nx2_sub * nx3_sub;
  chunk->assign(static_cast<std::size_t>(out_vars) * sub_cells, 0.0F);

  const Real inv_vol = 1.0 / static_cast<Real>(nfine_per_coarse_);
  const int jj_max = multi_d_ ? downsample_factor_ : 1;
  const int kk_max = three_d_ ? downsample_factor_ : 1;

  auto read_var = [&](int v, std::uint64_t idx) -> Real {
    return turb[static_cast<std::uint64_t>(v) * ncells + idx];
  };

  for (int k = 0; k < nx3_sub; ++k) {
    for (int j = 0; j < nx2_sub; ++j) {
      for (int i = 0; i < nx1_sub; ++i) {
        const int fi = ng + downsample_factor_ * i;
        const int fj = multi_d_ ? (ng + downsample_factor_ * j) : 0;
        const int fk = three_d_ ? (ng + downsample_factor_ * k) : 0;

        const std::uint64_t cidx = CellIndex(i, j, k, nx1_sub, nx2_sub);
        for (int v = 0; v < out_vars; ++v) {
          Real sum = 0.0;
          for (int kk = 0; kk < kk_max; ++kk) {
            for (int jj = 0; jj < jj_max; ++jj) {
              for (int ii = 0; ii < downsample_factor_; ++ii) {
                const std::uint64_t idx =
                    CellIndex(fi + ii, fj + jj, fk + kk, nout1, nout2);
                sum += read_var(v, idx);
              }
            }
          }
          (*chunk)[static_cast<std::size_t>(v) * sub_cells + cidx] =
              static_cast<float>(inv_vol * sum);
        }
      }
    }
  }

  return true;
}

bool Downsampler::DownsampleToBinary(const std::string& output_filename) {
  if (downsample_factor_ != 2 && downsample_factor_ != 4) {
    if (reader_.GetMyRank() == 0) {
      std::cerr << "ERROR: Unsupported downsample factor " << downsample_factor_
                << " (supported: 2 or 4)" << std::endl;
    }
    return false;
  }

  if (!ParseEOS()) {
    return false;
  }
  if (!EnsureCoarseSetup()) {
    return false;
  }

  const auto& mesh_indcs = reader_.GetMeshIndcs();
  if ((mesh_indcs.nx1 % downsample_factor_) != 0
      || (multi_d_ && (mesh_indcs.nx2 % downsample_factor_) != 0)
      || (three_d_ && (mesh_indcs.nx3 % downsample_factor_) != 0)) {
    if (reader_.GetMyRank() == 0) {
      std::cerr << "ERROR: Global mesh dimensions must be divisible by "
                << downsample_factor_ << std::endl;
    }
    return false;
  }

  const RegionIndcs& mb_indcs = reader_.GetMBIndcs();
  if ((mb_indcs.nx1 % downsample_factor_) != 0
      || (multi_d_ && (mb_indcs.nx2 % downsample_factor_) != 0)
      || (three_d_ && (mb_indcs.nx3 % downsample_factor_) != 0)) {
    if (reader_.GetMyRank() == 0) {
      std::cerr << "ERROR: MeshBlock cell counts must be divisible by "
                << downsample_factor_ << std::endl;
    }
    return false;
  }

  const int my_rank = reader_.GetMyRank();
  const int nranks = reader_.GetNRanks();

  MPIDistribution coarse_dist(nranks, coarse_nmb_total_);
  coarse_dist.SetupDistribution();

  const int coarse_gid_start = coarse_dist.GetStartingMeshBlockID(my_rank);
  const int coarse_nmb_thisrank = coarse_dist.GetNumMeshBlocks(my_rank);

  const int nx1 = mb_indcs.nx1;
  const int nx2 = (mb_indcs.nx2 > 1) ? mb_indcs.nx2 : 1;
  const int nx3 = (mb_indcs.nx3 > 1) ? mb_indcs.nx3 : 1;

  const std::uint64_t cells =
      static_cast<std::uint64_t>(nx1) * nx2 * nx3;

  const auto& config = reader_.GetPhysicsConfig();
  const bool has_energy = (config.nmhd >= 5);
  const bool ideal_eos = (eos_ == "ideal");
  const bool isothermal_eos = (eos_ == "isothermal");
  if (!ideal_eos && !isothermal_eos) {
    if (my_rank == 0) {
      std::cerr << "ERROR: Unsupported EOS '" << eos_ << "'" << std::endl;
    }
    return false;
  }
  if (ideal_eos && !has_energy) {
    if (my_rank == 0) {
      std::cerr << "ERROR: Ideal EOS requires energy (nmhd >= 5)" << std::endl;
    }
    return false;
  }
  const int out_vars = ideal_eos ? 8 : 7;

  std::vector<float> local_data(
      static_cast<std::size_t>(coarse_nmb_thisrank) * out_vars * cells, 0.0F);

  std::vector<int> recv_counts(coarse_nmb_thisrank, 0);

  const int fine_gid_start =
      reader_.GetMPIDistribution()->GetStartingMeshBlockID(my_rank);
  const int fine_nmb_thisrank = reader_.GetNMBThisRank();

  const int nx1_sub = nx1 / downsample_factor_;
  const int nx2_sub = multi_d_ ? (nx2 / downsample_factor_) : 1;
  const int nx3_sub = three_d_ ? (nx3 / downsample_factor_) : 1;
  const std::uint64_t sub_cells =
      static_cast<std::uint64_t>(nx1_sub) * nx2_sub * nx3_sub;

  const auto& fine_gids_eachrank = reader_.GetMPIDistribution()->GetGidsEachRank();
  const auto& fine_nmb_eachrank = reader_.GetMPIDistribution()->GetNmbEachRank();
  const auto& coarse_gids_eachrank = coarse_dist.GetGidsEachRank();
  const auto& coarse_nmb_eachrank = coarse_dist.GetNmbEachRank();

  std::vector<std::vector<char>> send_buffers(nranks);
  std::vector<int> send_chunk_counts(nranks, 0);

  std::uint64_t local_chunks_inserted = 0;

  for (int lf = 0; lf < fine_nmb_thisrank; ++lf) {
    const int fine_gid = fine_gid_start + lf;
    const int coarse_gid = fine_gid / nfine_per_coarse_;
    const int child_index = fine_gid % nfine_per_coarse_;

    const int owner_rank =
        FindRankForGid(coarse_gids_eachrank, coarse_nmb_eachrank, coarse_gid);

    std::vector<float> chunk;
    if (!ComputeDownsampledChunk(lf, &chunk)) {
      if (my_rank == 0) {
        std::cerr << "ERROR: Failed to downsample fine MeshBlock gid " << fine_gid
                  << std::endl;
      }
      return false;
    }
    if (chunk.size() != static_cast<std::size_t>(out_vars) * sub_cells) {
      if (my_rank == 0) {
        std::cerr << "ERROR: Unexpected downsample chunk size " << chunk.size()
                  << std::endl;
      }
      return false;
    }

    if (owner_rank == my_rank) {
      const int local_cg = coarse_gid - coarse_gid_start;
      if (local_cg < 0 || local_cg >= coarse_nmb_thisrank) {
        return false;
      }

      const int ox1 = ChildOffsetX1(child_index);
      const int ox2 = multi_d_ ? ChildOffsetX2(child_index) : 0;
      const int ox3 = three_d_ ? ChildOffsetX3(child_index) : 0;

      const int i_off = ox1 * nx1_sub;
      const int j_off = ox2 * nx2_sub;
      const int k_off = ox3 * nx3_sub;

      const std::uint64_t block_offset =
          static_cast<std::uint64_t>(local_cg) * out_vars * cells;
      for (int v = 0; v < out_vars; ++v) {
        const float* chunk_var = chunk.data() + static_cast<std::uint64_t>(v) * sub_cells;
        float* out_var = local_data.data() + block_offset
                         + static_cast<std::uint64_t>(v) * cells;

        for (int k = 0; k < nx3_sub; ++k) {
          for (int j = 0; j < nx2_sub; ++j) {
            for (int i = 0; i < nx1_sub; ++i) {
              const std::uint64_t cidx = CellIndex(i, j, k, nx1_sub, nx2_sub);
              const std::uint64_t pidx =
                  CellIndex(i_off + i, j_off + j, k_off + k, nx1, nx2);
              out_var[pidx] = chunk_var[cidx];
            }
          }
        }
      }

      recv_counts[local_cg] += 1;
      local_chunks_inserted += 1;
    } else {
      if (send_chunk_counts[owner_rank] == 0) {
        send_buffers[owner_rank].resize(sizeof(std::int32_t), 0);
      }

      ChunkHeader header;
      header.coarse_gid = coarse_gid;
      header.child_index = child_index;

      auto& buf = send_buffers[owner_rank];
      const std::size_t start = buf.size();
      buf.resize(start + sizeof(header) + sizeof(float) * chunk.size());
      std::memcpy(buf.data() + start, &header, sizeof(header));
      std::memcpy(buf.data() + start + sizeof(header), chunk.data(),
                  sizeof(float) * chunk.size());
      send_chunk_counts[owner_rank] += 1;
    }
  }

  for (int r = 0; r < nranks; ++r) {
    if (send_chunk_counts[r] == 0) {
      continue;
    }
    std::memcpy(send_buffers[r].data(), &send_chunk_counts[r], sizeof(std::int32_t));
  }

#if MPI_PARALLEL_ENABLED
  std::vector<MPI_Request> send_reqs;
  send_reqs.reserve(nranks);
  for (int r = 0; r < nranks; ++r) {
    if (send_chunk_counts[r] == 0) {
      continue;
    }
    MPI_Request req;
    MPI_Isend(send_buffers[r].data(), static_cast<int>(send_buffers[r].size()), MPI_BYTE,
              r, 900, MPI_COMM_WORLD, &req);
    send_reqs.push_back(req);
  }
#else
  (void)send_buffers;
#endif

  const std::uint64_t total_needed =
      static_cast<std::uint64_t>(coarse_nmb_thisrank) * nfine_per_coarse_;
  const std::uint64_t remote_needed = total_needed - local_chunks_inserted;

  std::uint64_t remote_received = 0;

#if MPI_PARALLEL_ENABLED
  while (remote_received < remote_needed) {
    MPI_Status status;
    MPI_Probe(MPI_ANY_SOURCE, 900, MPI_COMM_WORLD, &status);
    int nbytes = 0;
    MPI_Get_count(&status, MPI_BYTE, &nbytes);
    if (nbytes <= 0) {
      break;
    }
    std::vector<char> recv_buf(static_cast<std::size_t>(nbytes));
    MPI_Recv(recv_buf.data(), nbytes, MPI_BYTE, status.MPI_SOURCE, 900, MPI_COMM_WORLD,
             MPI_STATUS_IGNORE);

    std::size_t pos = 0;
    std::int32_t n_chunks = 0;
    std::memcpy(&n_chunks, recv_buf.data(), sizeof(n_chunks));
    pos += sizeof(n_chunks);

    const std::size_t floats_per_chunk = static_cast<std::size_t>(out_vars) * sub_cells;
    const std::size_t bytes_per_chunk =
        sizeof(ChunkHeader) + sizeof(float) * floats_per_chunk;
    std::vector<float> chunk_storage(floats_per_chunk);
    for (int c = 0; c < n_chunks; ++c) {
      if (pos + bytes_per_chunk > recv_buf.size()) {
        return false;
      }
      ChunkHeader header;
      std::memcpy(&header, recv_buf.data() + pos, sizeof(header));
      pos += sizeof(header);

      std::memcpy(chunk_storage.data(), recv_buf.data() + pos,
                  sizeof(float) * floats_per_chunk);
      pos += sizeof(float) * floats_per_chunk;

      const int coarse_gid = header.coarse_gid;
      const int child_index = header.child_index;

      if (coarse_gid < coarse_gid_start
          || coarse_gid >= coarse_gid_start + coarse_nmb_thisrank) {
        return false;
      }
      const int local_cg = coarse_gid - coarse_gid_start;

      const int ox1 = ChildOffsetX1(child_index);
      const int ox2 = multi_d_ ? ChildOffsetX2(child_index) : 0;
      const int ox3 = three_d_ ? ChildOffsetX3(child_index) : 0;

      const int i_off = ox1 * nx1_sub;
      const int j_off = ox2 * nx2_sub;
      const int k_off = ox3 * nx3_sub;

      const std::uint64_t block_offset =
          static_cast<std::uint64_t>(local_cg) * out_vars * cells;
      for (int v = 0; v < out_vars; ++v) {
        const float* chunk_var =
            chunk_storage.data() + static_cast<std::uint64_t>(v) * sub_cells;
        float* out_var = local_data.data() + block_offset
                         + static_cast<std::uint64_t>(v) * cells;
        for (int k = 0; k < nx3_sub; ++k) {
          for (int j = 0; j < nx2_sub; ++j) {
            for (int i = 0; i < nx1_sub; ++i) {
              const std::uint64_t cidx = CellIndex(i, j, k, nx1_sub, nx2_sub);
              const std::uint64_t pidx =
                  CellIndex(i_off + i, j_off + j, k_off + k, nx1, nx2);
              out_var[pidx] = chunk_var[cidx];
            }
          }
        }
      }

      recv_counts[local_cg] += 1;
      remote_received += 1;
    }
  }

  if (!send_reqs.empty()) {
    MPI_Waitall(static_cast<int>(send_reqs.size()), send_reqs.data(),
                MPI_STATUSES_IGNORE);
  }
#else
  if (remote_needed != 0) {
    return false;
  }
#endif

  for (int m = 0; m < coarse_nmb_thisrank; ++m) {
    if (recv_counts[m] != nfine_per_coarse_) {
      if (my_rank == 0) {
        std::cerr << "ERROR: Coarse MeshBlock " << (coarse_gid_start + m)
                  << " received " << recv_counts[m] << " children, expected "
                  << nfine_per_coarse_ << std::endl;
      }
      return false;
    }
  }

  const int nmb_rootx1 = coarse_mesh_indcs_.nx1 / coarse_mb_indcs_.nx1;
  const int nmb_rootx2 = multi_d_ ? (coarse_mesh_indcs_.nx2 / coarse_mb_indcs_.nx2) : 1;
  const int nmb_rootx3 = three_d_ ? (coarse_mesh_indcs_.nx3 / coarse_mb_indcs_.nx3) : 1;

  std::string param_out =
      ParameterParser::ReplaceMeshNx(reader_.GetParameterString(), coarse_mesh_indcs_.nx1,
                                     coarse_mesh_indcs_.nx2, coarse_mesh_indcs_.nx3);

  std::vector<std::string> labels = {"dens", "velx", "vely", "velz"};
  if (ideal_eos) {
    labels.push_back("eint");
  }
  labels.push_back("bcc1");
  labels.push_back("bcc2");
  labels.push_back("bcc3");

  BinaryWriter writer(my_rank, nranks);
  return writer.WriteBinaryFile(output_filename, labels, param_out, reader_.GetTime(),
                               reader_.GetNCycle(), coarse_mesh_size_, coarse_root_level_,
                               nmb_rootx1, nmb_rootx2, nmb_rootx3, coarse_mb_indcs_,
                               coarse_lloc_eachmb_, coarse_dist, local_data);
}

bool Downsampler::DownsampleTurbulenceToBinary(const std::string& output_filename) {
  if (downsample_factor_ != 2 && downsample_factor_ != 4) {
    if (reader_.GetMyRank() == 0) {
      std::cerr << "ERROR: Unsupported downsample factor " << downsample_factor_
                << " (supported: 2 or 4)" << std::endl;
    }
    return false;
  }

  const auto& config = reader_.GetPhysicsConfig();
  if (!config.has_turbulence) {
    if (reader_.GetMyRank() == 0) {
      std::cerr << "ERROR: Restart file does not contain turbulence force data"
                << std::endl;
    }
    return false;
  }

  if (config.nforce <= 0) {
    if (reader_.GetMyRank() == 0) {
      std::cerr << "ERROR: Invalid turbulence force variable count nforce="
                << config.nforce << std::endl;
    }
    return false;
  }

  if (!EnsureCoarseSetup()) {
    return false;
  }

  const auto& mesh_indcs = reader_.GetMeshIndcs();
  if ((mesh_indcs.nx1 % downsample_factor_) != 0
      || (multi_d_ && (mesh_indcs.nx2 % downsample_factor_) != 0)
      || (three_d_ && (mesh_indcs.nx3 % downsample_factor_) != 0)) {
    if (reader_.GetMyRank() == 0) {
      std::cerr << "ERROR: Global mesh dimensions must be divisible by "
                << downsample_factor_ << std::endl;
    }
    return false;
  }

  const RegionIndcs& mb_indcs = reader_.GetMBIndcs();
  if ((mb_indcs.nx1 % downsample_factor_) != 0
      || (multi_d_ && (mb_indcs.nx2 % downsample_factor_) != 0)
      || (three_d_ && (mb_indcs.nx3 % downsample_factor_) != 0)) {
    if (reader_.GetMyRank() == 0) {
      std::cerr << "ERROR: MeshBlock cell counts must be divisible by "
                << downsample_factor_ << std::endl;
    }
    return false;
  }

  const int my_rank = reader_.GetMyRank();
  const int nranks = reader_.GetNRanks();

  MPIDistribution coarse_dist(nranks, coarse_nmb_total_);
  coarse_dist.SetupDistribution();

  const int coarse_gid_start = coarse_dist.GetStartingMeshBlockID(my_rank);
  const int coarse_nmb_thisrank = coarse_dist.GetNumMeshBlocks(my_rank);

  const int nx1 = mb_indcs.nx1;
  const int nx2 = (mb_indcs.nx2 > 1) ? mb_indcs.nx2 : 1;
  const int nx3 = (mb_indcs.nx3 > 1) ? mb_indcs.nx3 : 1;

  const std::uint64_t cells =
      static_cast<std::uint64_t>(nx1) * nx2 * nx3;

  const int out_vars = config.nforce;
  std::vector<float> local_data(
      static_cast<std::size_t>(coarse_nmb_thisrank) * out_vars * cells, 0.0F);

  std::vector<int> recv_counts(coarse_nmb_thisrank, 0);

  const int fine_gid_start =
      reader_.GetMPIDistribution()->GetStartingMeshBlockID(my_rank);
  const int fine_nmb_thisrank = reader_.GetNMBThisRank();

  const int nx1_sub = nx1 / downsample_factor_;
  const int nx2_sub = multi_d_ ? (nx2 / downsample_factor_) : 1;
  const int nx3_sub = three_d_ ? (nx3 / downsample_factor_) : 1;
  const std::uint64_t sub_cells =
      static_cast<std::uint64_t>(nx1_sub) * nx2_sub * nx3_sub;

  const auto& coarse_gids_eachrank = coarse_dist.GetGidsEachRank();
  const auto& coarse_nmb_eachrank = coarse_dist.GetNmbEachRank();

  std::vector<std::vector<char>> send_buffers(nranks);
  std::vector<int> send_chunk_counts(nranks, 0);

  std::uint64_t local_chunks_inserted = 0;

  for (int lf = 0; lf < fine_nmb_thisrank; ++lf) {
    const int fine_gid = fine_gid_start + lf;
    const int coarse_gid = fine_gid / nfine_per_coarse_;
    const int child_index = fine_gid % nfine_per_coarse_;

    const int owner_rank =
        FindRankForGid(coarse_gids_eachrank, coarse_nmb_eachrank, coarse_gid);

    std::vector<float> chunk;
    if (!ComputeDownsampledTurbulenceChunk(lf, &chunk)) {
      if (my_rank == 0) {
        std::cerr << "ERROR: Failed to downsample turbulence for fine MeshBlock gid "
                  << fine_gid << std::endl;
      }
      return false;
    }
    if (chunk.size() != static_cast<std::size_t>(out_vars) * sub_cells) {
      if (my_rank == 0) {
        std::cerr << "ERROR: Unexpected turbulence chunk size " << chunk.size()
                  << std::endl;
      }
      return false;
    }

    if (owner_rank == my_rank) {
      const int local_cg = coarse_gid - coarse_gid_start;
      if (local_cg < 0 || local_cg >= coarse_nmb_thisrank) {
        return false;
      }

      const int ox1 = ChildOffsetX1(child_index);
      const int ox2 = multi_d_ ? ChildOffsetX2(child_index) : 0;
      const int ox3 = three_d_ ? ChildOffsetX3(child_index) : 0;

      const int i_off = ox1 * nx1_sub;
      const int j_off = ox2 * nx2_sub;
      const int k_off = ox3 * nx3_sub;

      const std::uint64_t block_offset =
          static_cast<std::uint64_t>(local_cg) * out_vars * cells;
      for (int v = 0; v < out_vars; ++v) {
        const float* chunk_var =
            chunk.data() + static_cast<std::uint64_t>(v) * sub_cells;
        float* out_var = local_data.data() + block_offset
                         + static_cast<std::uint64_t>(v) * cells;

        for (int k = 0; k < nx3_sub; ++k) {
          for (int j = 0; j < nx2_sub; ++j) {
            for (int i = 0; i < nx1_sub; ++i) {
              const std::uint64_t cidx = CellIndex(i, j, k, nx1_sub, nx2_sub);
              const std::uint64_t pidx =
                  CellIndex(i_off + i, j_off + j, k_off + k, nx1, nx2);
              out_var[pidx] = chunk_var[cidx];
            }
          }
        }
      }

      recv_counts[local_cg] += 1;
      local_chunks_inserted += 1;
    } else {
      if (send_chunk_counts[owner_rank] == 0) {
        send_buffers[owner_rank].resize(sizeof(std::int32_t), 0);
      }

      ChunkHeader header;
      header.coarse_gid = coarse_gid;
      header.child_index = child_index;

      auto& buf = send_buffers[owner_rank];
      const std::size_t start = buf.size();
      buf.resize(start + sizeof(header) + sizeof(float) * chunk.size());
      std::memcpy(buf.data() + start, &header, sizeof(header));
      std::memcpy(buf.data() + start + sizeof(header), chunk.data(),
                  sizeof(float) * chunk.size());
      send_chunk_counts[owner_rank] += 1;
    }
  }

  for (int r = 0; r < nranks; ++r) {
    if (send_chunk_counts[r] == 0) {
      continue;
    }
    std::memcpy(send_buffers[r].data(), &send_chunk_counts[r], sizeof(std::int32_t));
  }

#if MPI_PARALLEL_ENABLED
  std::vector<MPI_Request> send_reqs;
  send_reqs.reserve(nranks);
  for (int r = 0; r < nranks; ++r) {
    if (send_chunk_counts[r] == 0) {
      continue;
    }
    MPI_Request req;
    MPI_Isend(send_buffers[r].data(), static_cast<int>(send_buffers[r].size()), MPI_BYTE,
              r, 901, MPI_COMM_WORLD, &req);
    send_reqs.push_back(req);
  }
#else
  (void)send_buffers;
#endif

  const std::uint64_t total_needed =
      static_cast<std::uint64_t>(coarse_nmb_thisrank) * nfine_per_coarse_;
  const std::uint64_t remote_needed = total_needed - local_chunks_inserted;

  std::uint64_t remote_received = 0;

#if MPI_PARALLEL_ENABLED
  while (remote_received < remote_needed) {
    MPI_Status status;
    MPI_Probe(MPI_ANY_SOURCE, 901, MPI_COMM_WORLD, &status);
    int nbytes = 0;
    MPI_Get_count(&status, MPI_BYTE, &nbytes);
    if (nbytes <= 0) {
      break;
    }
    std::vector<char> recv_buf(static_cast<std::size_t>(nbytes));
    MPI_Recv(recv_buf.data(), nbytes, MPI_BYTE, status.MPI_SOURCE, 901, MPI_COMM_WORLD,
             MPI_STATUS_IGNORE);

    std::size_t pos = 0;
    std::int32_t n_chunks = 0;
    if (recv_buf.size() < sizeof(n_chunks)) {
      return false;
    }
    std::memcpy(&n_chunks, recv_buf.data(), sizeof(n_chunks));
    pos += sizeof(n_chunks);

    const std::size_t floats_per_chunk = static_cast<std::size_t>(out_vars) * sub_cells;
    const std::size_t bytes_per_chunk =
        sizeof(ChunkHeader) + sizeof(float) * floats_per_chunk;
    std::vector<float> chunk_storage(floats_per_chunk);
    for (int c = 0; c < n_chunks; ++c) {
      if (pos + bytes_per_chunk > recv_buf.size()) {
        return false;
      }
      ChunkHeader header;
      std::memcpy(&header, recv_buf.data() + pos, sizeof(header));
      pos += sizeof(header);

      std::memcpy(chunk_storage.data(), recv_buf.data() + pos,
                  sizeof(float) * floats_per_chunk);
      pos += sizeof(float) * floats_per_chunk;

      const int coarse_gid = header.coarse_gid;
      const int child_index = header.child_index;

      if (coarse_gid < coarse_gid_start
          || coarse_gid >= coarse_gid_start + coarse_nmb_thisrank) {
        return false;
      }
      const int local_cg = coarse_gid - coarse_gid_start;

      const int ox1 = ChildOffsetX1(child_index);
      const int ox2 = multi_d_ ? ChildOffsetX2(child_index) : 0;
      const int ox3 = three_d_ ? ChildOffsetX3(child_index) : 0;

      const int i_off = ox1 * nx1_sub;
      const int j_off = ox2 * nx2_sub;
      const int k_off = ox3 * nx3_sub;

      const std::uint64_t block_offset =
          static_cast<std::uint64_t>(local_cg) * out_vars * cells;
      for (int v = 0; v < out_vars; ++v) {
        const float* chunk_var =
            chunk_storage.data() + static_cast<std::uint64_t>(v) * sub_cells;
        float* out_var = local_data.data() + block_offset
                         + static_cast<std::uint64_t>(v) * cells;
        for (int k = 0; k < nx3_sub; ++k) {
          for (int j = 0; j < nx2_sub; ++j) {
            for (int i = 0; i < nx1_sub; ++i) {
              const std::uint64_t cidx = CellIndex(i, j, k, nx1_sub, nx2_sub);
              const std::uint64_t pidx =
                  CellIndex(i_off + i, j_off + j, k_off + k, nx1, nx2);
              out_var[pidx] = chunk_var[cidx];
            }
          }
        }
      }

      recv_counts[local_cg] += 1;
      remote_received += 1;
    }
  }

  if (!send_reqs.empty()) {
    MPI_Waitall(static_cast<int>(send_reqs.size()), send_reqs.data(),
                MPI_STATUSES_IGNORE);
  }
#else
  if (remote_needed != 0) {
    return false;
  }
#endif

  for (int m = 0; m < coarse_nmb_thisrank; ++m) {
    if (recv_counts[m] != nfine_per_coarse_) {
      if (my_rank == 0) {
        std::cerr << "ERROR: Coarse MeshBlock " << (coarse_gid_start + m)
                  << " received " << recv_counts[m] << " children, expected "
                  << nfine_per_coarse_ << std::endl;
      }
      return false;
    }
  }

  const int nmb_rootx1 = coarse_mesh_indcs_.nx1 / coarse_mb_indcs_.nx1;
  const int nmb_rootx2 = multi_d_ ? (coarse_mesh_indcs_.nx2 / coarse_mb_indcs_.nx2) : 1;
  const int nmb_rootx3 = three_d_ ? (coarse_mesh_indcs_.nx3 / coarse_mb_indcs_.nx3) : 1;

  std::string param_out =
      ParameterParser::ReplaceMeshNx(reader_.GetParameterString(), coarse_mesh_indcs_.nx1,
                                     coarse_mesh_indcs_.nx2, coarse_mesh_indcs_.nx3);

  std::vector<std::string> labels;
  labels.reserve(static_cast<std::size_t>(out_vars));
  for (int v = 0; v < out_vars; ++v) {
    labels.push_back("force" + std::to_string(v + 1));
  }

  BinaryWriter writer(my_rank, nranks);
  return writer.WriteBinaryFile(output_filename, labels, param_out, reader_.GetTime(),
                               reader_.GetNCycle(), coarse_mesh_size_, coarse_root_level_,
                               nmb_rootx1, nmb_rootx2, nmb_rootx3, coarse_mb_indcs_,
                               coarse_lloc_eachmb_, coarse_dist, local_data);
}
