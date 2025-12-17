#ifndef BINARY_WRITER_HPP_
#define BINARY_WRITER_HPP_

#include <string>
#include <vector>

#include "common.hpp"
#include "io_wrapper.hpp"
#include "mpi_distribution.hpp"

class BinaryWriter {
 public:
  BinaryWriter(int my_rank, int nranks);
  ~BinaryWriter() = default;

  bool WriteBinaryFile(const std::string& filename,
                       const std::vector<std::string>& var_labels,
                       const std::string& parameter_string, Real time, int cycle,
                       const RegionSize& mesh_size, int root_level,
                       int nmb_rootx1, int nmb_rootx2, int nmb_rootx3,
                       const RegionIndcs& mb_indcs,
                       const std::vector<LogicalLocation>& lloc_eachmb,
                       const MPIDistribution& mpi_dist,
                       const std::vector<float>& local_data);

 private:
  int my_rank_;
  int nranks_;
  IOWrapper file_;

  static Real LeftEdgeX(int ith, int n, Real xmin, Real xmax);
  static void ComputeBlockBounds(const RegionSize& mesh_size, int root_level,
                                 int nmb_rootx1, int nmb_rootx2, int nmb_rootx3,
                                 const LogicalLocation& lloc, RegionSize* mb_size);
};

#endif  // BINARY_WRITER_HPP_

