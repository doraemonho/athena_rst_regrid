#include "downsampler.hpp"

#include <algorithm>
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
constexpr int kVarVel1 = 1;
constexpr int kVarVel2 = 2;
constexpr int kVarVel3 = 3;
constexpr int kVarPress = 4;
constexpr int kVarBcc1 = 5;
constexpr int kVarBcc2 = 6;
constexpr int kVarBcc3 = 7;

constexpr int kOutVars = 8;

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

Downsampler::Downsampler(RestartReader& reader) : reader_(reader) {
  const auto& mesh_indcs = reader_.GetMeshIndcs();
  multi_d_ = (mesh_indcs.nx2 > 1);
  three_d_ = (mesh_indcs.nx3 > 1);
  nfine_per_coarse_ = 2;
  if (multi_d_) {
    nfine_per_coarse_ = 4;
  }
  if (three_d_) {
    nfine_per_coarse_ = 8;
  }
}

LogicalLocation Downsampler::ParentLocation(const LogicalLocation& fine) {
  LogicalLocation parent;
  parent.level = fine.level - 1;
  parent.lx1 = fine.lx1 / 2;
  parent.lx2 = fine.lx2 / 2;
  parent.lx3 = fine.lx3 / 2;
  return parent;
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
  return true;
}

void Downsampler::ComputeCoarseMeshConfig() {
  const RegionSize& fine_mesh_size = reader_.GetMeshSize();
  const RegionIndcs& fine_mesh_indcs = reader_.GetMeshIndcs();
  const RegionIndcs& fine_mb_indcs = reader_.GetMBIndcs();

  coarse_mesh_size_ = fine_mesh_size;
  coarse_mesh_size_.dx1 = fine_mesh_size.dx1 * 2.0;
  coarse_mesh_size_.dx2 = fine_mesh_size.dx2 * (multi_d_ ? 2.0 : 1.0);
  coarse_mesh_size_.dx3 = fine_mesh_size.dx3 * (three_d_ ? 2.0 : 1.0);

  coarse_mesh_indcs_ = fine_mesh_indcs;
  coarse_mesh_indcs_.nx1 = fine_mesh_indcs.nx1 / 2;
  coarse_mesh_indcs_.nx2 = multi_d_ ? (fine_mesh_indcs.nx2 / 2) : 1;
  coarse_mesh_indcs_.nx3 = three_d_ ? (fine_mesh_indcs.nx3 / 2) : 1;

  coarse_mesh_indcs_.is = fine_mesh_indcs.ng;
  coarse_mesh_indcs_.ie = coarse_mesh_indcs_.is + coarse_mesh_indcs_.nx1 - 1;
  coarse_mesh_indcs_.js = fine_mesh_indcs.ng;
  coarse_mesh_indcs_.je = coarse_mesh_indcs_.js + coarse_mesh_indcs_.nx2 - 1;
  coarse_mesh_indcs_.ks = fine_mesh_indcs.ng;
  coarse_mesh_indcs_.ke = coarse_mesh_indcs_.ks + coarse_mesh_indcs_.nx3 - 1;

  coarse_mb_indcs_ = fine_mb_indcs;

  coarse_nmb_total_ = reader_.GetNMBTotal() / nfine_per_coarse_;
  coarse_root_level_ = reader_.GetRootLevel() - 1;
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

  if (reader_.GetRootLevel() < 1) {
    if (reader_.GetMyRank() == 0) {
      std::cerr << "ERROR: root_level=" << reader_.GetRootLevel()
                << " too small for 2x downsampling" << std::endl;
    }
    return false;
  }

  coarse_lloc_eachmb_.resize(coarse_nmb_total_);

  for (int cg = 0; cg < coarse_nmb_total_; ++cg) {
    const int fg0 = cg * nfine_per_coarse_;
    const LogicalLocation parent = ParentLocation(fine_llocs[fg0]);
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

      const LogicalLocation computed_parent = ParentLocation(child);
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

      if (child.lx1 != 2 * parent.lx1 + ox1 || child.lx2 != 2 * parent.lx2 + ox2
          || child.lx3 != 2 * parent.lx3 + ox3) {
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

  const auto& mb_indcs = reader_.GetMBIndcs();
  const int ng = mb_indcs.ng;
  const int nx1 = mb_indcs.nx1;
  const int nx2 = (mb_indcs.nx2 > 1) ? mb_indcs.nx2 : 1;
  const int nx3 = (mb_indcs.nx3 > 1) ? mb_indcs.nx3 : 1;

  if ((nx1 % 2) != 0 || (multi_d_ && (nx2 % 2) != 0) || (three_d_ && (nx3 % 2) != 0)) {
    return false;
  }

  const int nx1_sub = nx1 / 2;
  const int nx2_sub = multi_d_ ? (nx2 / 2) : 1;
  const int nx3_sub = three_d_ ? (nx3 / 2) : 1;

  const int nout1 = nx1 + 2 * ng;
  const int nout2 = multi_d_ ? (nx2 + 2 * ng) : 1;
  const int nout3 = three_d_ ? (nx3 + 2 * ng) : 1;

  const auto& mhd = phys->GetMHDData().at(local_fine_mb);
  const auto& b1f = phys->GetMHDB1f().at(local_fine_mb);
  const auto& b2f = phys->GetMHDB2f().at(local_fine_mb);
  const auto& b3f = phys->GetMHDB3f().at(local_fine_mb);

  const std::uint64_t sub_cells =
      static_cast<std::uint64_t>(nx1_sub) * nx2_sub * nx3_sub;
  chunk->assign(static_cast<std::size_t>(kOutVars) * sub_cells, 0.0F);

  std::vector<Real> b1c(static_cast<std::size_t>(nx3_sub) * nx2_sub * (nx1_sub + 1), 0.0);
  std::vector<Real> b2c(static_cast<std::size_t>(nx3_sub) * (nx2_sub + 1) * nx1_sub, 0.0);
  std::vector<Real> b3c(static_cast<std::size_t>(nx3_sub + 1) * nx2_sub * nx1_sub, 0.0);

  const int nout1p1 = nout1 + 1;
  const int nout2p1 = nout2 + 1;
  const int nout3p1 = nout3 + 1;

  for (int k = 0; k < nx3_sub; ++k) {
    for (int j = 0; j < nx2_sub; ++j) {
      for (int i = 0; i <= nx1_sub; ++i) {
        const int fi = ng + 2 * i;
        const int fj = multi_d_ ? (ng + 2 * j) : 0;
        const int fk = three_d_ ? (ng + 2 * k) : 0;

        Real sum = 0.0;
        if (three_d_) {
          for (int kk = 0; kk < 2; ++kk) {
            for (int jj = 0; jj < 2; ++jj) {
              const std::uint64_t idx =
                  FaceIndexX1(fi, fj + jj, fk + kk, nout1p1, nout2);
              sum += b1f[idx];
            }
          }
          sum *= 0.25;
        } else if (multi_d_) {
          sum = 0.5 * (b1f[FaceIndexX1(fi, fj, 0, nout1p1, nout2)]
                       + b1f[FaceIndexX1(fi, fj + 1, 0, nout1p1, nout2)]);
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
        const int fi = ng + 2 * i;
        const int fj = multi_d_ ? (ng + 2 * j) : 0;
        const int fk = three_d_ ? (ng + 2 * k) : 0;

        Real sum = 0.0;
        if (three_d_) {
          for (int kk = 0; kk < 2; ++kk) {
            for (int ii = 0; ii < 2; ++ii) {
              const std::uint64_t idx =
                  FaceIndexX2(fi + ii, fj, fk + kk, nout1, nout2p1);
              sum += b2f[idx];
            }
          }
          sum *= 0.25;
        } else if (multi_d_) {
          sum = 0.5 * (b2f[FaceIndexX2(fi, fj, 0, nout1, nout2p1)]
                       + b2f[FaceIndexX2(fi + 1, fj, 0, nout1, nout2p1)]);
        } else {
          sum = 0.5 * (b2f[FaceIndexX2(fi, 0, 0, nout1, 2)]
                       + b2f[FaceIndexX2(fi + 1, 0, 0, nout1, 2)]);
        }

        const std::uint64_t cidx = FaceIndexX2(i, j, k, nx1_sub, nx2_sub + 1);
        b2c[cidx] = sum;
      }
    }
  }

  for (int k = 0; k <= nx3_sub; ++k) {
    for (int j = 0; j < nx2_sub; ++j) {
      for (int i = 0; i < nx1_sub; ++i) {
        const int fi = ng + 2 * i;
        const int fj = multi_d_ ? (ng + 2 * j) : 0;
        const int fk = three_d_ ? (ng + 2 * k) : 0;

        Real sum = 0.0;
        if (three_d_) {
          for (int jj = 0; jj < 2; ++jj) {
            for (int ii = 0; ii < 2; ++ii) {
              const std::uint64_t idx =
                  FaceIndexX3(fi + ii, fj + jj, fk, nout1, nout2);
              sum += b3f[idx];
            }
          }
          sum *= 0.25;
        } else if (multi_d_) {
          sum = 0.25 * (b3f[FaceIndexX3(fi, fj, 0, nout1, nout2)]
                        + b3f[FaceIndexX3(fi + 1, fj, 0, nout1, nout2)]
                        + b3f[FaceIndexX3(fi, fj + 1, 0, nout1, nout2)]
                        + b3f[FaceIndexX3(fi + 1, fj + 1, 0, nout1, nout2)]);
        } else {
          sum = 0.5 * (b3f[FaceIndexX3(fi, 0, 0, nout1, 1)]
                       + b3f[FaceIndexX3(fi + 1, 0, 0, nout1, 1)]);
        }

        const std::uint64_t cidx = FaceIndexX3(i, j, k, nx1_sub, nx2_sub);
        b3c[cidx] = sum;
      }
    }
  }

  const int ncells = nout1 * nout2 * nout3;

  const bool has_energy = (config.nmhd >= 5);
  const bool ideal_eos = (eos_ == "ideal");
  const bool isothermal_eos = (eos_ == "isothermal");

  if (ideal_eos && !has_energy) {
    return false;
  }
  if (!ideal_eos && !isothermal_eos) {
    return false;
  }

  for (int k = 0; k < nx3_sub; ++k) {
    for (int j = 0; j < nx2_sub; ++j) {
      for (int i = 0; i < nx1_sub; ++i) {
        const int fi = ng + 2 * i;
        const int fj = multi_d_ ? (ng + 2 * j) : 0;
        const int fk = three_d_ ? (ng + 2 * k) : 0;

        const std::uint64_t idx0 = CellIndex(fi, fj, fk, nout1, nout2);
        const std::uint64_t idx1 = CellIndex(fi + 1, fj, fk, nout1, nout2);

        Real rho = 0.0;
        Real m1 = 0.0;
        Real m2 = 0.0;
        Real m3 = 0.0;
        Real eng = 0.0;

        auto read_var = [&](int v, std::uint64_t idx) -> Real {
          return mhd[static_cast<std::uint64_t>(v) * ncells + idx];
        };

        if (three_d_) {
          const std::uint64_t idx2 = CellIndex(fi, fj + 1, fk, nout1, nout2);
          const std::uint64_t idx3 = CellIndex(fi + 1, fj + 1, fk, nout1, nout2);
          const std::uint64_t idx4 = CellIndex(fi, fj, fk + 1, nout1, nout2);
          const std::uint64_t idx5 = CellIndex(fi + 1, fj, fk + 1, nout1, nout2);
          const std::uint64_t idx6 = CellIndex(fi, fj + 1, fk + 1, nout1, nout2);
          const std::uint64_t idx7 = CellIndex(fi + 1, fj + 1, fk + 1, nout1, nout2);

          rho = 0.125 * (read_var(0, idx0) + read_var(0, idx1) + read_var(0, idx2)
                         + read_var(0, idx3) + read_var(0, idx4) + read_var(0, idx5)
                         + read_var(0, idx6) + read_var(0, idx7));
          m1 = 0.125 * (read_var(1, idx0) + read_var(1, idx1) + read_var(1, idx2)
                        + read_var(1, idx3) + read_var(1, idx4) + read_var(1, idx5)
                        + read_var(1, idx6) + read_var(1, idx7));
          m2 = 0.125 * (read_var(2, idx0) + read_var(2, idx1) + read_var(2, idx2)
                        + read_var(2, idx3) + read_var(2, idx4) + read_var(2, idx5)
                        + read_var(2, idx6) + read_var(2, idx7));
          m3 = 0.125 * (read_var(3, idx0) + read_var(3, idx1) + read_var(3, idx2)
                        + read_var(3, idx3) + read_var(3, idx4) + read_var(3, idx5)
                        + read_var(3, idx6) + read_var(3, idx7));
          if (has_energy) {
            eng = 0.125 * (read_var(4, idx0) + read_var(4, idx1) + read_var(4, idx2)
                           + read_var(4, idx3) + read_var(4, idx4) + read_var(4, idx5)
                           + read_var(4, idx6) + read_var(4, idx7));
          }
        } else if (multi_d_) {
          const std::uint64_t idx2 = CellIndex(fi, fj + 1, 0, nout1, nout2);
          const std::uint64_t idx3 = CellIndex(fi + 1, fj + 1, 0, nout1, nout2);
          rho = 0.25 * (read_var(0, idx0) + read_var(0, idx1) + read_var(0, idx2)
                        + read_var(0, idx3));
          m1 = 0.25 * (read_var(1, idx0) + read_var(1, idx1) + read_var(1, idx2)
                       + read_var(1, idx3));
          m2 = 0.25 * (read_var(2, idx0) + read_var(2, idx1) + read_var(2, idx2)
                       + read_var(2, idx3));
          m3 = 0.25 * (read_var(3, idx0) + read_var(3, idx1) + read_var(3, idx2)
                       + read_var(3, idx3));
          if (has_energy) {
            eng = 0.25 * (read_var(4, idx0) + read_var(4, idx1) + read_var(4, idx2)
                          + read_var(4, idx3));
          }
        } else {
          rho = 0.5 * (read_var(0, idx0) + read_var(0, idx1));
          m1 = 0.5 * (read_var(1, idx0) + read_var(1, idx1));
          m2 = 0.5 * (read_var(2, idx0) + read_var(2, idx1));
          m3 = 0.5 * (read_var(3, idx0) + read_var(3, idx1));
          if (has_energy) {
            eng = 0.5 * (read_var(4, idx0) + read_var(4, idx1));
          }
        }

        const std::uint64_t b1_left = FaceIndexX1(i, j, k, nx1_sub + 1, nx2_sub);
        const std::uint64_t b1_right = FaceIndexX1(i + 1, j, k, nx1_sub + 1, nx2_sub);
        const std::uint64_t b2_bot = FaceIndexX2(i, j, k, nx1_sub, nx2_sub + 1);
        const std::uint64_t b2_top = FaceIndexX2(i, j + 1, k, nx1_sub, nx2_sub + 1);
        const std::uint64_t b3_bot = FaceIndexX3(i, j, k, nx1_sub, nx2_sub);
        const std::uint64_t b3_top = FaceIndexX3(i, j, k + 1, nx1_sub, nx2_sub);

        const Real bx = 0.5 * (b1c[b1_left] + b1c[b1_right]);
        const Real by = 0.5 * (b2c[b2_bot] + b2c[b2_top]);
        const Real bz = 0.5 * (b3c[b3_bot] + b3c[b3_top]);

        Real v1 = 0.0;
        Real v2 = 0.0;
        Real v3 = 0.0;
        if (rho != 0.0) {
          v1 = m1 / rho;
          v2 = m2 / rho;
          v3 = m3 / rho;
        }

        Real press = 0.0;
        if (ideal_eos) {
          const Real ke = 0.5 * (m1 * m1 + m2 * m2 + m3 * m3) / std::max(rho, 1e-30);
          const Real me = 0.5 * (bx * bx + by * by + bz * bz);
          const Real eint = eng - ke - me;
          press = (gamma_ - 1.0) * eint;
        } else {
          const Real cs2 = iso_sound_speed_ * iso_sound_speed_;
          press = cs2 * rho;
        }

        const std::uint64_t cidx = CellIndex(i, j, k, nx1_sub, nx2_sub);
        (*chunk)[static_cast<std::size_t>(kVarDens) * sub_cells + cidx] =
            static_cast<float>(rho);
        (*chunk)[static_cast<std::size_t>(kVarVel1) * sub_cells + cidx] =
            static_cast<float>(v1);
        (*chunk)[static_cast<std::size_t>(kVarVel2) * sub_cells + cidx] =
            static_cast<float>(v2);
        (*chunk)[static_cast<std::size_t>(kVarVel3) * sub_cells + cidx] =
            static_cast<float>(v3);
        (*chunk)[static_cast<std::size_t>(kVarPress) * sub_cells + cidx] =
            static_cast<float>(press);
        (*chunk)[static_cast<std::size_t>(kVarBcc1) * sub_cells + cidx] =
            static_cast<float>(bx);
        (*chunk)[static_cast<std::size_t>(kVarBcc2) * sub_cells + cidx] =
            static_cast<float>(by);
        (*chunk)[static_cast<std::size_t>(kVarBcc3) * sub_cells + cidx] =
            static_cast<float>(bz);
      }
    }
  }

  return true;
}

bool Downsampler::DownsampleToBinary(const std::string& output_filename) {
  if (!ParseEOS()) {
    return false;
  }
  ComputeCoarseMeshConfig();

  if (!BuildCoarseLogicalLocations()) {
    return false;
  }

  const auto& mesh_indcs = reader_.GetMeshIndcs();
  if ((mesh_indcs.nx1 % 2) != 0 || (multi_d_ && (mesh_indcs.nx2 % 2) != 0)
      || (three_d_ && (mesh_indcs.nx3 % 2) != 0)) {
    if (reader_.GetMyRank() == 0) {
      std::cerr << "ERROR: Global mesh dimensions must be divisible by 2" << std::endl;
    }
    return false;
  }

  const RegionIndcs& mb_indcs = reader_.GetMBIndcs();
  if ((mb_indcs.nx1 % 2) != 0 || (multi_d_ && (mb_indcs.nx2 % 2) != 0)
      || (three_d_ && (mb_indcs.nx3 % 2) != 0)) {
    if (reader_.GetMyRank() == 0) {
      std::cerr << "ERROR: MeshBlock cell counts must be divisible by 2" << std::endl;
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
  std::vector<float> local_data(static_cast<std::size_t>(coarse_nmb_thisrank) * kOutVars
                                * cells, 0.0F);

  std::vector<int> recv_counts(coarse_nmb_thisrank, 0);

  const int fine_gid_start =
      reader_.GetMPIDistribution()->GetStartingMeshBlockID(my_rank);
  const int fine_nmb_thisrank = reader_.GetNMBThisRank();

  const int nx1_sub = nx1 / 2;
  const int nx2_sub = multi_d_ ? (nx2 / 2) : 1;
  const int nx3_sub = three_d_ ? (nx3 / 2) : 1;
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
          static_cast<std::uint64_t>(local_cg) * kOutVars * cells;
      for (int v = 0; v < kOutVars; ++v) {
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

    const std::size_t floats_per_chunk = static_cast<std::size_t>(kOutVars) * sub_cells;
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
          static_cast<std::uint64_t>(local_cg) * kOutVars * cells;
      for (int v = 0; v < kOutVars; ++v) {
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

  const std::vector<std::string> labels = {"dens", "vel1", "vel2", "vel3",
                                           "press", "Bcc1", "Bcc2", "Bcc3"};

  BinaryWriter writer(my_rank, nranks);
  return writer.WriteBinaryFile(output_filename, labels, param_out, reader_.GetTime(),
                               reader_.GetNCycle(), coarse_mesh_size_, coarse_root_level_,
                               nmb_rootx1, nmb_rootx2, nmb_rootx3, coarse_mb_indcs_,
                               coarse_lloc_eachmb_, coarse_dist, local_data);
}
