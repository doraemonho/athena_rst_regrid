#ifndef MPI_DISTRIBUTION_HPP_
#define MPI_DISTRIBUTION_HPP_

#include <vector>
#include "common.hpp"

//----------------------------------------------------------------------------------------
// MPI rank distribution functionality following AthenaK's load balancing logic
//----------------------------------------------------------------------------------------
class MPIDistribution {
public:
    MPIDistribution(int nranks, int nmb_total);
    
    // Setup distribution using AthenaK's load balancing algorithm
    void SetupDistribution(const std::vector<int>* gids_eachrank = nullptr,
                          const std::vector<int>* nmb_eachrank = nullptr);
    
    // Accessors
    int GetStartingMeshBlockID(int rank) const { return gids_eachrank_[rank]; }
    int GetNumMeshBlocks(int rank) const { return nmb_eachrank_[rank]; }
    void GetRankMinMax(int* noutmbs_min, int* noutmbs_max) const;
    
    // Print distribution info
    void PrintDistribution() const;
    
private:
    int nranks_;
    int nmb_total_;
    std::vector<int> gids_eachrank_;  // Starting MeshBlock ID for each rank
    std::vector<int> nmb_eachrank_;   // Number of MeshBlocks per rank
    
    // AthenaK's load balancing algorithm
    void LoadBalance();
};

#endif // MPI_DISTRIBUTION_HPP_