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
    
private:
    RestartReader& reader_;
    
    // Fine mesh configuration
    RegionSize fine_mesh_size_;
    RegionIndcs fine_mesh_indcs_;
    RegionIndcs fine_mb_indcs_;
    size_t fine_nmb_total_;
    std::vector<LogicalLocation> fine_lloc_eachmb_;
    std::vector<float> fine_cost_eachmb_;
    
    // Helper functions
    void CalculateFineMeshConfiguration();
    void GenerateFineLogicalLocations();
    bool UpscaleMeshBlock(int coarse_mb_id, int fine_mb_start_id,
                         const Real* coarse_data, Real* fine_data,
                         int nvars, bool is_face_centered = false, int face_dir = 0);
    bool UpscaleFaceCenteredData(int coarse_mb_id, int fine_mb_start_id,
                                const Real* coarse_x1f, const Real* coarse_x2f, const Real* coarse_x3f,
                                Real* fine_x1f, Real* fine_x2f, Real* fine_x3f);
    
    // Data storage for upscaled data
    std::vector<Real> fine_hydro_data_;
    std::vector<Real> fine_mhd_data_;
    std::vector<Real> fine_rad_data_;
    std::vector<Real> fine_turb_data_;
    std::vector<Real> fine_z4c_data_;
    std::vector<Real> fine_adm_data_;
    
    // Face-centered data for MHD
    std::vector<Real> fine_x1f_data_;
    std::vector<Real> fine_x2f_data_;
    std::vector<Real> fine_x3f_data_;
};

#endif // UPSCALER_HPP_