#include "binary_writer.hpp"

#include <algorithm>
#include <cstdint>
#include <cstring>
#include <iostream>
#include <sstream>

#if MPI_PARALLEL_ENABLED
#include <mpi.h>
#endif

namespace {
constexpr std::uint64_t kMaxMPICount = 2147483648ULL;  // 2^31
}  // namespace

BinaryWriter::BinaryWriter(int my_rank, int nranks) : my_rank_(my_rank), nranks_(nranks) {
}

Real BinaryWriter::LeftEdgeX(int ith, int n, Real xmin, Real xmax) {
  const Real x = static_cast<Real>(ith) / static_cast<Real>(n);
  return (x * xmax - x * xmin) - (0.5 * xmax - 0.5 * xmin) + (0.5 * xmin + 0.5 * xmax);
}

void BinaryWriter::ComputeBlockBounds(const RegionSize& mesh_size, int root_level,
                                      int nmb_rootx1, int nmb_rootx2, int nmb_rootx3,
                                      const LogicalLocation& lloc, RegionSize* mb_size) {
  const std::int32_t lev = lloc.level;

  const std::int32_t nmbx1 = nmb_rootx1 << (lev - root_level);
  if (lloc.lx1 == 0) {
    mb_size->x1min = mesh_size.x1min;
  } else {
    mb_size->x1min = LeftEdgeX(lloc.lx1, nmbx1, mesh_size.x1min, mesh_size.x1max);
  }
  if (lloc.lx1 == nmbx1 - 1) {
    mb_size->x1max = mesh_size.x1max;
  } else {
    mb_size->x1max = LeftEdgeX(lloc.lx1 + 1, nmbx1, mesh_size.x1min, mesh_size.x1max);
  }

  if (nmb_rootx2 <= 1) {
    mb_size->x2min = mesh_size.x2min;
    mb_size->x2max = mesh_size.x2max;
  } else {
    const std::int32_t nmbx2 = nmb_rootx2 << (lev - root_level);
    if (lloc.lx2 == 0) {
      mb_size->x2min = mesh_size.x2min;
    } else {
      mb_size->x2min = LeftEdgeX(lloc.lx2, nmbx2, mesh_size.x2min, mesh_size.x2max);
    }
    if (lloc.lx2 == nmbx2 - 1) {
      mb_size->x2max = mesh_size.x2max;
    } else {
      mb_size->x2max = LeftEdgeX(lloc.lx2 + 1, nmbx2, mesh_size.x2min, mesh_size.x2max);
    }
  }

  if (nmb_rootx3 <= 1) {
    mb_size->x3min = mesh_size.x3min;
    mb_size->x3max = mesh_size.x3max;
  } else {
    const std::int32_t nmbx3 = nmb_rootx3 << (lev - root_level);
    if (lloc.lx3 == 0) {
      mb_size->x3min = mesh_size.x3min;
    } else {
      mb_size->x3min = LeftEdgeX(lloc.lx3, nmbx3, mesh_size.x3min, mesh_size.x3max);
    }
    if (lloc.lx3 == nmbx3 - 1) {
      mb_size->x3max = mesh_size.x3max;
    } else {
      mb_size->x3max = LeftEdgeX(lloc.lx3 + 1, nmbx3, mesh_size.x3min, mesh_size.x3max);
    }
  }
}

bool BinaryWriter::WriteBinaryFile(
    const std::string& filename, const std::vector<std::string>& var_labels,
    const std::string& parameter_string, Real time, int cycle,
    const RegionSize& mesh_size, int root_level, int nmb_rootx1, int nmb_rootx2,
    int nmb_rootx3, const RegionIndcs& mb_indcs,
    const std::vector<LogicalLocation>& lloc_eachmb, const MPIDistribution& mpi_dist,
    const std::vector<float>& local_data) {
  if (!file_.Open(filename.c_str(), IOWrapper::FileMode::write)) {
    if (my_rank_ == 0) {
      std::cerr << "ERROR: Cannot create binary file: " << filename << std::endl;
    }
    return false;
  }

  std::uint64_t header_offset = 0;
  if (my_rank_ == 0) {
    std::stringstream msg;
    msg << "Athena binary output version=1.1\n"
        << "  size of preheader=5\n"
        << "  time=" << time << "\n"
        << "  cycle=" << cycle << "\n"
        << "  size of location=" << sizeof(Real) << "\n"
        << "  size of variable=" << sizeof(float) << "\n"
        << "  number of variables=" << var_labels.size() << "\n"
        << "  variables:  ";
    for (const auto& label : var_labels) {
      msg << label << "  ";
    }
    msg << "\n";

    const std::string preheader = msg.str();
    file_.Write_any_type(preheader.data(), preheader.size(), "byte");

    std::stringstream msg2;
    msg2 << "  header offset=" << parameter_string.size() * sizeof(char) << "\n";
    const std::string header_line = msg2.str();
    file_.Write_any_type(header_line.data(), header_line.size(), "byte");
    file_.Write_any_type(parameter_string.data(), parameter_string.size(), "byte");

    header_offset = file_.GetPosition();
  }

#if MPI_PARALLEL_ENABLED
  MPI_Bcast(&header_offset, sizeof(header_offset), MPI_BYTE, 0, MPI_COMM_WORLD);
#endif

  const int nout_vars = static_cast<int>(var_labels.size());
  const int nx1 = mb_indcs.nx1;
  const int nx2 = (mb_indcs.nx2 > 1) ? mb_indcs.nx2 : 1;
  const int nx3 = (mb_indcs.nx3 > 1) ? mb_indcs.nx3 : 1;
  const std::uint64_t cells =
      static_cast<std::uint64_t>(nx1) * static_cast<std::uint64_t>(nx2)
      * static_cast<std::uint64_t>(nx3);

  const std::uint64_t record_size =
      10ULL * sizeof(std::int32_t) + 6ULL * sizeof(Real)
      + cells * nout_vars * sizeof(float);

  const int ns_mbs = mpi_dist.GetStartingMeshBlockID(my_rank_);
  const int nb_mbs = mpi_dist.GetNumMeshBlocks(my_rank_);

  if (nb_mbs < 0) {
    if (my_rank_ == 0) {
      std::cerr << "ERROR: Invalid local MeshBlock count" << std::endl;
    }
    return false;
  }
  const std::uint64_t expected_floats =
      static_cast<std::uint64_t>(nb_mbs) * nout_vars * cells;
  if (local_data.size() != expected_floats) {
    if (my_rank_ == 0) {
      std::cerr << "ERROR: local_data size mismatch: got " << local_data.size()
                << ", expected " << expected_floats << std::endl;
    }
    return false;
  }

  std::vector<char> buffer(static_cast<std::size_t>(nb_mbs) * record_size);
  for (int m = 0; m < nb_mbs; ++m) {
    const int gid = ns_mbs + m;
    const LogicalLocation& loc = lloc_eachmb.at(gid);
    const std::int32_t phys_level = loc.level - root_level;

    RegionSize mb_size;
    std::memset(&mb_size, 0, sizeof(mb_size));
    ComputeBlockBounds(mesh_size, root_level, nmb_rootx1, nmb_rootx2, nmb_rootx3, loc,
                       &mb_size);

    char* pdata = buffer.data() + static_cast<std::size_t>(m) * record_size;

    const std::int32_t ois = mb_indcs.is;
    const std::int32_t oie = mb_indcs.ie;
    const std::int32_t ojs = mb_indcs.js;
    const std::int32_t oje = mb_indcs.je;
    const std::int32_t oks = mb_indcs.ks;
    const std::int32_t oke = mb_indcs.ke;

    const std::int32_t lx1 = loc.lx1;
    const std::int32_t lx2 = loc.lx2;
    const std::int32_t lx3 = loc.lx3;

    const std::int32_t ints[10] = {ois,  oie,  ojs, oje, oks,
                                   oke,  lx1,  lx2, lx3, phys_level};
    std::memcpy(pdata, ints, sizeof(ints));
    pdata += sizeof(ints);

    const Real coords[6] = {mb_size.x1min, mb_size.x1max, mb_size.x2min,
                            mb_size.x2max, mb_size.x3min, mb_size.x3max};
    std::memcpy(pdata, coords, sizeof(coords));
    pdata += sizeof(coords);

    const std::uint64_t block_float_offset =
        static_cast<std::uint64_t>(m) * static_cast<std::uint64_t>(nout_vars) * cells;
    const float* block_data = local_data.data() + block_float_offset;
    std::memcpy(pdata, block_data,
                static_cast<std::size_t>(nout_vars) * cells * sizeof(float));
  }

  if (record_size * static_cast<std::uint64_t>(nb_mbs) <= kMaxMPICount) {
    const std::uint64_t myoffset =
        header_offset + record_size * static_cast<std::uint64_t>(ns_mbs);
    file_.Write_any_type_at_all(
        buffer.data(), record_size * static_cast<std::uint64_t>(nb_mbs), myoffset,
        "byte");
  } else {
    int noutmbs_min = 0;
    int noutmbs_max = 0;
    mpi_dist.GetRankMinMax(&noutmbs_min, &noutmbs_max);
    for (int m = 0; m < noutmbs_max; ++m) {
      const char* pdata = buffer.data() + static_cast<std::size_t>(m) * record_size;
      const std::uint64_t myoffset =
          header_offset + record_size * static_cast<std::uint64_t>(ns_mbs)
          + record_size * static_cast<std::uint64_t>(m);
      if (m < noutmbs_min) {
        file_.Write_any_type_at_all(pdata, record_size, myoffset, "byte");
      } else if (m < nb_mbs) {
        file_.Write_any_type_at(pdata, record_size, myoffset, "byte");
      }
    }
  }

  file_.Close();
  return true;
}
