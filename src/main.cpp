#include "restart_reader.hpp"
#include "downsampler.hpp"
#include "upscaler.hpp"
#include <iostream>
#include <cstdlib>
#include <string>
#include <cstring>

#if MPI_PARALLEL_ENABLED
#include <mpi.h>
#endif

void PrintUsage(const char* prog_name, int my_rank) {
    if (my_rank == 0) {
        std::cerr << "Usage: " << prog_name << " <restart_file.rst> [options]"
                  << std::endl;
        std::cerr << "Options:" << std::endl;
        std::cerr << "  --upscale <output.rst>  : Upscale from N³ to (2N)³ resolution"
                  << std::endl;
        std::cerr << "  --downsample-bin <output.bin> : Downsample to (N/f)³ and write"
                  << " AthenaK .bin" << std::endl;
        std::cerr << "  --downsample-turb-bin <output.bin> : Downsample turbulence "
                  << "force and write AthenaK .bin" << std::endl;
        std::cerr << "  --downsample-factor <2|4> : Set downsample factor f (default: 2)"
                  << std::endl;
        std::cerr << "Examples:" << std::endl;
        std::cerr << "  Read mode:    " << prog_name
                  << " ../debeg_rst/OutputTest.00000.rst" << std::endl;
        std::cerr << "  Upscale mode: " << prog_name << " input.rst --upscale output.rst"
                  << std::endl;
        std::cerr << "  Downsample:   " << prog_name
                  << " input.rst --downsample-bin output.bin" << std::endl;
        std::cerr << "  Downsample4x: " << prog_name
                  << " input.rst --downsample-bin output.bin --downsample-factor 4"
                  << std::endl;
        std::cerr << "  Down+Turb:    " << prog_name
                  << " input.rst --downsample-bin mhd.bin --downsample-turb-bin turb.bin"
                  << std::endl;
#if MPI_PARALLEL_ENABLED
        std::cerr << "  MPI read:     mpirun -np 4 " << prog_name << " input.rst"
                  << std::endl;
        std::cerr << "  MPI upscale:  mpirun -np 4 " << prog_name
                  << " input.rst --upscale output.rst" << std::endl;
        std::cerr << "  MPI down:     mpirun -np 4 " << prog_name
                  << " input.rst --downsample-bin output.bin" << std::endl;
#endif
    }
}

int main(int argc, char* argv[]) {
    int my_rank = 0;
    int nranks = 1;
    
#if MPI_PARALLEL_ENABLED
    // Initialize MPI
    if (MPI_SUCCESS != MPI_Init(&argc, &argv)) {
        std::cerr << "### FATAL ERROR: MPI Initialization failed." << std::endl;
        return 1;
    }
    
    // Get process rank and total number of ranks
    if (MPI_SUCCESS != MPI_Comm_rank(MPI_COMM_WORLD, &my_rank)) {
        std::cerr << "### FATAL ERROR: MPI_Comm_rank failed." << std::endl;
        MPI_Finalize();
        return 1;
    }
    
    if (MPI_SUCCESS != MPI_Comm_size(MPI_COMM_WORLD, &nranks)) {
        std::cerr << "### FATAL ERROR: MPI_Comm_size failed." << std::endl;
        MPI_Finalize();
        return 1;
    }
#endif
    
    // Parse command line arguments
    if (argc < 2) {
        PrintUsage(argv[0], my_rank);
#if MPI_PARALLEL_ENABLED
        MPI_Finalize();
#endif
        return 1;
    }
    
    std::string filename = argv[1];
    bool upscale_mode = false;
    std::string upscale_output_filename;
    bool downsample_mhd_mode = false;
    bool downsample_turb_mode = false;
    std::string output_mhd_filename;
    std::string output_turb_filename;
    int downsample_factor = 2;
    bool downsample_factor_set = false;
    
    // Check for upscale option
    for (int i = 2; i < argc; i++) {
        if (std::strcmp(argv[i], "--upscale") == 0) {
            if (i + 1 >= argc) {
                if (my_rank == 0) {
                    std::cerr << "Error: --upscale requires an output filename"
                              << std::endl;
                }
                PrintUsage(argv[0], my_rank);
#if MPI_PARALLEL_ENABLED
                MPI_Finalize();
#endif
                return 1;
            }
            upscale_mode = true;
            upscale_output_filename = argv[i + 1];
            i++; // Skip the next argument
        } else if (std::strcmp(argv[i], "--downsample-bin") == 0) {
            if (i + 1 >= argc) {
                if (my_rank == 0) {
                    std::cerr << "Error: --downsample-bin requires an output filename"
                              << std::endl;
                }
                PrintUsage(argv[0], my_rank);
#if MPI_PARALLEL_ENABLED
                MPI_Finalize();
#endif
                return 1;
            }
            downsample_mhd_mode = true;
            output_mhd_filename = argv[i + 1];
            i++; // Skip the next argument
        } else if (std::strcmp(argv[i], "--downsample-turb-bin") == 0) {
            if (i + 1 >= argc) {
                if (my_rank == 0) {
                    std::cerr << "Error: --downsample-turb-bin requires an output "
                              << "filename" << std::endl;
                }
                PrintUsage(argv[0], my_rank);
#if MPI_PARALLEL_ENABLED
                MPI_Finalize();
#endif
                return 1;
            }
            downsample_turb_mode = true;
            output_turb_filename = argv[i + 1];
            i++; // Skip the next argument
        } else if (std::strcmp(argv[i], "--downsample-factor") == 0) {
            if (i + 1 >= argc) {
                if (my_rank == 0) {
                    std::cerr << "Error: --downsample-factor requires an integer value"
                              << std::endl;
                }
                PrintUsage(argv[0], my_rank);
#if MPI_PARALLEL_ENABLED
                MPI_Finalize();
#endif
                return 1;
            }
            downsample_factor = std::atoi(argv[i + 1]);
            downsample_factor_set = true;
            if (downsample_factor != 2 && downsample_factor != 4) {
                if (my_rank == 0) {
                    std::cerr << "Error: --downsample-factor must be 2 or 4"
                              << std::endl;
                }
                PrintUsage(argv[0], my_rank);
#if MPI_PARALLEL_ENABLED
                MPI_Finalize();
#endif
                return 1;
            }
            i++; // Skip the next argument
        }
    }

    if (downsample_factor_set && !downsample_mhd_mode && !downsample_turb_mode) {
        if (my_rank == 0) {
            std::cerr << "Error: --downsample-factor requires a downsample output option"
                      << std::endl;
        }
        PrintUsage(argv[0], my_rank);
#if MPI_PARALLEL_ENABLED
        MPI_Finalize();
#endif
        return 1;
    }
    
    if (my_rank == 0) {
        std::cout << "=== AthenaK Restart Reader with MPI Support ===" << std::endl;
        std::cout << "Running with " << nranks << " MPI rank(s)" << std::endl;
        if (upscale_mode && (downsample_mhd_mode || downsample_turb_mode)) {
            std::cerr << "Error: cannot use --upscale with downsampling options"
                      << std::endl;
#if MPI_PARALLEL_ENABLED
            MPI_Finalize();
#endif
            return 1;
        }
        if (upscale_mode) {
            std::cout << "Mode: UPSCALING from N³ to (2N)³" << std::endl;
            std::cout << "Output file: " << upscale_output_filename << std::endl;
        } else if (downsample_mhd_mode || downsample_turb_mode) {
            std::cout << "Mode: DOWNSAMPLING from N³ to (N/" << downsample_factor
                      << ")³ (binary output)" << std::endl;
            if (downsample_mhd_mode) {
                std::cout << "MHD output file: " << output_mhd_filename << std::endl;
            }
            if (downsample_turb_mode) {
                std::cout << "Turbulence output file: " << output_turb_filename
                          << std::endl;
            }
        } else {
            std::cout << "Mode: READ ONLY" << std::endl;
        }
        std::cout << std::endl;
    }
    
    // Use scope block to ensure reader destructor runs before MPI_Finalize
    {
        RestartReader reader;
        
        // Use automatic load balancing (no external distribution provided)
        if (!reader.ReadRestartFile(filename, my_rank, nranks)) {
            std::cerr << "Failed to read restart file" << std::endl;
            return 1;
        }
        
        if (upscale_mode) {
            // Perform upscaling
            Upscaler upscaler(reader);
            if (!upscaler.UpscaleRestartFile(upscale_output_filename)) {
                std::cerr << "Failed to upscale restart file" << std::endl;
                return 1;
            }
        } else if (downsample_mhd_mode || downsample_turb_mode) {
            Downsampler downsampler(reader, downsample_factor);
            if (downsample_mhd_mode) {
                if (!downsampler.DownsampleToBinary(output_mhd_filename)) {
                    std::cerr << "Failed to downsample restart file (MHD)" << std::endl;
                    return 1;
                }
            }
            if (downsample_turb_mode) {
                if (!downsampler.DownsampleTurbulenceToBinary(output_turb_filename)) {
                    std::cerr << "Failed to downsample restart file (turbulence)"
                              << std::endl;
                    return 1;
                }
            }
        } else {
            // Normal read mode
            reader.DisplayMeshInfo();
            
            if (reader.GetNMBThisRank() > 0) {
                reader.DisplayPhysicsData(0);  // Show data from first local MeshBlock
                reader.ValidateUniformValues();
            } else {
                std::cout << "No MeshBlocks assigned to this rank." << std::endl;
            }
        }
    } // RestartReader destructor called here, before MPI_Finalize
    
#if MPI_PARALLEL_ENABLED
    MPI_Finalize();
#endif
    
    return 0;
}
