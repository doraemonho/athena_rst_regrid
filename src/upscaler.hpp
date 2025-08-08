#ifndef UPSCALER_HPP_
#define UPSCALER_HPP_

#include <vector>
#include <memory>
#include "common.hpp"
#include "restart_reader.hpp"
#include "prolongation.hpp"

//----------------------------------------------------------------------------------------
// Upscaler class - transforms restart data from N³ to (2N)³ resolution
//----------------------------------------------------------------------------------------
class Upscaler {
    friend class RestartWriter;  // Allow RestartWriter to access private members
public:
    Upscaler(RestartReader& reader);
    ~Upscaler() = default;
    
    // Main upscaling function
    bool UpscaleRestartFile(const std::string& output_filename);
    
    // Accessor for MPI distribution info (needed by RestartWriter)
    size_t GetFineNMBThisRank() const { return fine_nmb_thisrank_; }
    int GetFineGidsStart() const { return fine_gids_start_; }
    
private:
    RestartReader& reader_;
    
    // Fine mesh configuration
    RegionSize fine_mesh_size_;
    RegionIndcs fine_mesh_indcs_;
    RegionIndcs fine_mb_indcs_;
    size_t fine_nmb_total_;
    std::vector<LogicalLocation> fine_lloc_eachmb_;
    std::vector<float> fine_cost_eachmb_;
    
    // MPI distribution for fine meshblocks
    size_t fine_nmb_thisrank_;     // Number of fine meshblocks on this rank
    int coarse_gids_start_;         // Starting global ID of coarse meshblocks for this rank
    int fine_gids_start_;           // Starting global ID of fine meshblocks for this rank
    int nfine_per_coarse_;          // Number of fine meshblocks per coarse meshblock (2/4/8)
    
    // Helper functions
    void CalculateFineMeshConfiguration();
    void GenerateFineLogicalLocations();
    void SetupFineMPIDistribution();
    int GetGlobalCoarseMBID(int local_mb_id) const;  // Convert local to global coarse MB ID
    int GetLocalFineMBIndex(int local_coarse_mb, int fine_mb_within_coarse) const;  // Get local fine MB index
    
    bool UpscaleMeshBlock(int local_coarse_mb_id, int local_fine_mb_start_id,
                         const Real* coarse_data, Real* fine_data,
                         int nvars, bool is_face_centered = false, int face_dir = 0);
    bool UpscaleFaceCenteredData(int local_coarse_mb_id, int local_fine_mb_start_id,
                                const Real* coarse_x1f, const Real* coarse_x2f, const Real* coarse_x3f,
                                Real* fine_x1f, Real* fine_x2f, Real* fine_x3f);
    
    // Data storage for upscaled data (LOCAL to this rank only)
    std::vector<Real> fine_hydro_data_;
    std::vector<Real> fine_mhd_data_;
    std::vector<Real> fine_rad_data_;
    std::vector<Real> fine_turb_data_;
    std::vector<Real> fine_z4c_data_;
    std::vector<Real> fine_adm_data_;
    
    // Face-centered data for MHD (LOCAL to this rank only)
    std::vector<Real> fine_x1f_data_;
    std::vector<Real> fine_x2f_data_;
    std::vector<Real> fine_x3f_data_;
};

#endif // UPSCALER_HPP_