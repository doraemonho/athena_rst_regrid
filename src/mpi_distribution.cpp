#include "mpi_distribution.hpp"
#include <iostream>
#include <algorithm>

MPIDistribution::MPIDistribution(int nranks, int nmb_total) 
    : nranks_(nranks), nmb_total_(nmb_total) {
    gids_eachrank_.resize(nranks_);
    nmb_eachrank_.resize(nranks_);
}

void MPIDistribution::SetupDistribution(const std::vector<int>* gids_eachrank,
                                       const std::vector<int>* nmb_eachrank) {
    if (gids_eachrank != nullptr && nmb_eachrank != nullptr) {
        // Use provided MPI distribution (from external load balancing)
        gids_eachrank_ = *gids_eachrank;
        nmb_eachrank_ = *nmb_eachrank;
        std::cout << "Using provided MPI distribution" << std::endl;
    } else {
        // Implement AthenaK's load balancing logic
        LoadBalance();
        std::cout << "Using automatic load balancing" << std::endl;
    }
}

void MPIDistribution::LoadBalance() {
    // Following AthenaK's load_balance.cpp:37-76 exactly
    float totalcost = static_cast<float>(nmb_total_);
    float targetcost = totalcost / nranks_;
    
    std::vector<int> rank_eachmb(nmb_total_);
    
    // Create rank list from the end: the master MPI rank should have less load
    int j = nranks_ - 1;
    float mycost = 0.0;
    for (int i = nmb_total_ - 1; i >= 0; i--) {
        // AthenaK's exact error check from load_balance.cpp:52-58
        if (targetcost == 0.0) {
            std::cerr << "### FATAL ERROR: There is at least one process which has no MeshBlock. "
                      << "Decrease the number of processes or use smaller MeshBlocks." << std::endl;
            std::exit(EXIT_FAILURE);
        }
        mycost += 1.0;  // Assume uniform cost per MeshBlock (clist[i] in AthenaK)
        rank_eachmb[i] = j;
        if (mycost >= targetcost && j > 0) {
            j--;
            totalcost -= mycost;
            mycost = 0.0;
            targetcost = totalcost / (j + 1);
        }
    }
    
    // Calculate starting IDs and counts (following AthenaK's logic)
    gids_eachrank_[0] = 0;
    j = 0;
    for (int i = 1; i < nmb_total_; i++) {
        if (rank_eachmb[i] != rank_eachmb[i-1]) {
            nmb_eachrank_[j] = i - gids_eachrank_[j];
            gids_eachrank_[++j] = i;
        }
    }
    nmb_eachrank_[j] = nmb_total_ - gids_eachrank_[j];
}

void MPIDistribution::GetRankMinMax(int* noutmbs_min, int* noutmbs_max) const {
    // Following AthenaK's logic from restart.cpp:121-126
    *noutmbs_max = nmb_eachrank_[0];
    *noutmbs_min = nmb_eachrank_[0];
    for (int i = 0; i < nranks_; ++i) {
        *noutmbs_max = std::max(*noutmbs_max, nmb_eachrank_[i]);
        *noutmbs_min = std::min(*noutmbs_min, nmb_eachrank_[i]);
    }
}

void MPIDistribution::PrintDistribution() const {
    std::cout << "MPI Distribution:" << std::endl;
    for (int i = 0; i < nranks_; i++) {
        std::cout << "  Rank " << i << ": MeshBlocks " << gids_eachrank_[i] 
                  << "-" << (gids_eachrank_[i] + nmb_eachrank_[i] - 1)
                  << " (" << nmb_eachrank_[i] << " total)" << std::endl;
    }
}