#ifndef DOWNSAMPLER_HPP_
#define DOWNSAMPLER_HPP_

#include <string>
#include <vector>

#include "common.hpp"
#include "restart_reader.hpp"

class Downsampler {
 public:
  explicit Downsampler(RestartReader& reader, int downsample_factor = 2);
  ~Downsampler() = default;

  bool DownsampleToBinary(const std::string& output_filename);

 private:
  RestartReader& reader_;

  int downsample_factor_ = 2;
  int downsample_levels_ = 1;

  bool three_d_ = false;
  bool multi_d_ = false;
  int nfine_per_coarse_ = 0;

  int coarse_nmb_total_ = 0;
  int coarse_root_level_ = 0;
  RegionSize coarse_mesh_size_;
  RegionIndcs coarse_mesh_indcs_;
  RegionIndcs coarse_mb_indcs_;
  std::vector<LogicalLocation> coarse_lloc_eachmb_;

  std::string eos_;
  Real gamma_ = 5.0 / 3.0;
  Real iso_sound_speed_ = 1.0;

  void ComputeCoarseMeshConfig();
  bool BuildCoarseLogicalLocations();
  bool ParseEOS();

  int ChildOffsetX1(int child_index) const;
  int ChildOffsetX2(int child_index) const;
  int ChildOffsetX3(int child_index) const;

  LogicalLocation CoarseLocation(const LogicalLocation& fine) const;
  static bool SameLocation(const LogicalLocation& a, const LogicalLocation& b);

  static std::uint64_t CellIndex(int i, int j, int k, int nx1, int nx2);
  static std::uint64_t FaceIndexX1(int i, int j, int k, int nx1p1, int nx2);
  static std::uint64_t FaceIndexX2(int i, int j, int k, int nx1, int nx2p1);
  static std::uint64_t FaceIndexX3(int i, int j, int k, int nx1, int nx2);

  bool ComputeDownsampledChunk(int local_fine_mb, std::vector<float>* chunk) const;
};

#endif  // DOWNSAMPLER_HPP_
