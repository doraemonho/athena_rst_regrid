#ifndef RESTART_READER_HPP_
#define RESTART_READER_HPP_
// AthenaK Regridding Tool - RestartReader Class

#include <string>
#include <vector>
#include <memory>

// Standalone includes
#include "io_wrapper.hpp"

// Forward declarations and type definitions
using Real = double;
using IOWrapperSizeT = std::uint64_t;

// Physics variable indices (AthenaK standard)
const int IDN = 0;
const int IM1 = 1;
const int IM2 = 2;
const int IM3 = 3;
const int IEN = 4;

// Utility macros
template<typename T>
inline T SIGN(const T& x) { return (x > 0) - (x < 0); }

struct RegionSize {
  Real x1min, x2min, x3min;
  Real x1max, x2max, x3max;
  Real dx1, dx2, dx3;
};

struct RegionIndcs {
  int ng;
  int nx1, nx2, nx3;
  int is, ie, js, je, ks, ke;
  int cnx1, cnx2, cnx3;
  int cis, cie, cjs, cje, cks, cke;
};

struct LogicalLocation {
  int lx1, lx2, lx3, level;
};

// Must match AthenaK RNG_State structure exactly (from utils/random.hpp)
#define NTAB 32
struct RNG_State {
  int64_t idum;      // 8 bytes
  int64_t idum2;     // 8 bytes  
  int64_t iy;        // 8 bytes
  int64_t iv[NTAB];  // 32 * 8 = 256 bytes
  int iset;          // 4 bytes
  double gset;       // 8 bytes
  // Total: 292 bytes
};

//! \struct RestartData
//! \brief Container for all restart file data
struct RestartData {
  // Header information
  std::string input_params;
  int nmb_total;
  int root_level;
  RegionSize mesh_size;
  RegionIndcs mesh_indcs;
  RegionIndcs mb_indcs;
  Real time;
  Real dt;
  int ncycle;
  
  // MeshBlock structure
  std::vector<LogicalLocation> lloc_eachmb;
  std::vector<float> cost_eachmb;
  
  // Internal state data
  Real z4c_last_output_time;
  std::vector<Real> puncture_positions;  // 3*nco values
  RNG_State turb_rng_state;
  bool has_z4c, has_punctures, has_turb;
  
  // Physics data dimensions
  int nhydro, nmhd, nrad, nforce, nz4c, nadm;
  int nout1, nout2, nout3;  // including ghost zones
  
  // Physics data arrays (stored as 1D for efficiency)
  std::vector<Real> hydro_data;      // [nmb][nhydro][nout3][nout2][nout1]
  std::vector<Real> mhd_data;        // [nmb][nmhd][nout3][nout2][nout1] 
  std::vector<Real> mhd_b1f_data;    // [nmb][nout3][nout2][nout1+1]
  std::vector<Real> mhd_b2f_data;    // [nmb][nout3][nout2+1][nout1]
  std::vector<Real> mhd_b3f_data;    // [nmb][nout3+1][nout2][nout1]
  std::vector<Real> rad_data;        // [nmb][nrad][nout3][nout2][nout1]
  std::vector<Real> force_data;      // [nmb][nforce][nout3][nout2][nout1]
  std::vector<Real> z4c_data;        // [nmb][nz4c][nout3][nout2][nout1]
  std::vector<Real> adm_data;        // [nmb][nadm][nout3][nout2][nout1]
  
  // Helper functions for array indexing
  size_t GetCCIndex(int mb, int var, int k, int j, int i, int nvar) const {
    return mb*(size_t)nvar*nout3*nout2*nout1 + (size_t)var*nout3*nout2*nout1 + 
           (size_t)k*nout2*nout1 + (size_t)j*nout1 + (size_t)i;
  }
  
  size_t GetFCX1Index(int mb, int k, int j, int i) const {
    return mb*(size_t)nout3*nout2*(nout1+1) + (size_t)k*nout2*(nout1+1) + 
           (size_t)j*(nout1+1) + (size_t)i;
  }
  
  size_t GetFCX2Index(int mb, int k, int j, int i) const {
    return mb*(size_t)nout3*(nout2+1)*nout1 + (size_t)k*(nout2+1)*nout1 + 
           (size_t)j*nout1 + (size_t)i;
  }
  
  size_t GetFCX3Index(int mb, int k, int j, int i) const {
    return mb*(size_t)(nout3+1)*nout2*nout1 + (size_t)k*nout2*nout1 + 
           (size_t)j*nout1 + (size_t)i;
  }
};

//! \class RestartReader
//! \brief Handles reading AthenaK restart files
class RestartReader {
 public:
  RestartReader() = default;
  ~RestartReader() = default;
  
  // Main interface functions
  bool ReadRestartFile(const std::string& filename, RestartData& data);
  bool ValidateMHDFile(const RestartData& data) const;
  void PrintFileInfo(const RestartData& data) const;
  
  // Error handling
  const std::string& GetLastError() const { return last_error_; }
  
 private:
  // Internal helper functions
  bool ReadHeader(IOWrapper& file, RestartData& data);
  bool ReadMeshStructure(IOWrapper& file, RestartData& data);
  bool ReadInternalState(IOWrapper& file, RestartData& data);
  bool ReadPhysicsData(IOWrapper& file, RestartData& data);
  
  // Parsing helpers
  bool ParseInputParameters(const std::string& params, RestartData& data);
  void CalculateDataSizes(RestartData& data);
  IOWrapperSizeT CalculateOffset(const RestartData& data, int rank) const;
  
  // Validation helpers
  bool ValidateHeader(const RestartData& data) const;
  bool ValidateMeshStructure(const RestartData& data) const;
  bool ValidateDataConsistency(const RestartData& data) const;
  
  // Error reporting
  void SetError(const std::string& error) const { last_error_ = error; }
  mutable std::string last_error_;
};

#endif // RESTART_READER_HPP_