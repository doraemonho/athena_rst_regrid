// AthenaK Regridding Tool - Main Program

#include <iostream>
#include <string>
#include <vector>
#include <cstdlib>
#include <getopt.h>

#include "regrid_driver.hpp"

void PrintUsage(const char* program_name);
void PrintVersion();
bool ParseCommandLine(int argc, char* argv[], RegridOptions& options);

//! \fn int main()
//! \brief Main function
int main(int argc, char* argv[]) {
  
  // Parse command line arguments
  RegridOptions options;
  if (!ParseCommandLine(argc, argv, options)) {
    PrintUsage(argv[0]);
    return EXIT_FAILURE;
  }
  
  // Initialize and run regridding
  RegridDriver driver;
  
  try {
    if (!driver.RunRegridding(options)) {
      std::cerr << "ERROR: " << driver.GetLastError() << std::endl;
      return EXIT_FAILURE;
    }
  } catch (const std::exception& e) {
    std::cerr << "EXCEPTION: " << e.what() << std::endl;
    return EXIT_FAILURE;
  } catch (...) {
    std::cerr << "UNKNOWN ERROR: An unexpected error occurred" << std::endl;
    return EXIT_FAILURE;
  }
  
  return EXIT_SUCCESS;
}

//! \fn void PrintUsage()
//! \brief Print usage information
void PrintUsage(const char* program_name) {
  std::cout << "\nAthenaK MHD Restart File Regridding Tool" << std::endl;
  std::cout << "=========================================" << std::endl;
  std::cout << "Usage: " << program_name << " [options] <input_file> <output_file>" << std::endl;
  std::cout << "\nRequired arguments:" << std::endl;
  std::cout << "  input_file          Input AthenaK restart file (.rst)" << std::endl;
  std::cout << "  output_file         Output regridded restart file (.rst)" << std::endl;
  
  std::cout << "\nOptions:" << std::endl;
  std::cout << "  -h, --help          Show this help message" << std::endl;
  std::cout << "  -v, --version       Show version information" << std::endl;
  std::cout << "  -q, --quiet         Suppress progress output" << std::endl;
  std::cout << "  --verbose           Enable detailed output" << std::endl;
  std::cout << "  -b, --benchmark     Show detailed timing information" << std::endl;
  
  std::cout << "\nRegridding options:" << std::endl;
  std::cout << "  -r, --refine=N      Refinement factor (default: 2)" << std::endl;
  
  std::cout << "\nVerification options:" << std::endl;
  std::cout << "  --no-conservation   Skip conservation checks" << std::endl;
  std::cout << "  --no-divergence     Skip ∇·B = 0 verification" << std::endl;
  
  std::cout << "\nOutput options:" << std::endl;
  std::cout << "  --file-number=N     Set output file number" << std::endl;
  std::cout << "  --new-time=T        Set new simulation time" << std::endl;
  
  std::cout << "\nExamples:" << std::endl;
  std::cout << "  " << program_name << " input.rst output.rst" << std::endl;
  std::cout << "  " << program_name << " -r 2 --verbose input.00100.rst output.00100.rst" << std::endl;
  std::cout << "  " << program_name << " --refine=4 --benchmark turb.01000.rst turb_fine.01000.rst" << std::endl;
  
  std::cout << "\nDescription:" << std::endl;
  std::cout << "  This tool regrids AthenaK MHD restart files from N³ to (rN)³ resolution" << std::endl;
  std::cout << "  while preserving the same meshblock size. It uses AthenaK's prolongation" << std::endl;
  std::cout << "  operators to ensure physical accuracy and conservation properties." << std::endl;
  
  std::cout << "\nSupported physics modules:" << std::endl;
  std::cout << "  - MHD (required)" << std::endl;
  std::cout << "  - Hydro" << std::endl;
  std::cout << "  - Radiation" << std::endl;
  std::cout << "  - Turbulence forcing" << std::endl;
  std::cout << "  - Z4c/ADM (general relativity)" << std::endl;
  std::cout << std::endl;
}

//! \fn void PrintVersion()
//! \brief Print version information
void PrintVersion() {
  std::cout << "AthenaK MHD Restart File Regridding Tool" << std::endl;
  std::cout << "Version: " << REGRID_TOOL_VERSION << std::endl;
  std::cout << "Built with AthenaK prolongation operators" << std::endl;
  std::cout << "Supports MHD, hydro, radiation, turbulence, and GR physics modules" << std::endl;
  
#if SINGLE_PRECISION_ENABLED
  std::cout << "Precision: Single (float)" << std::endl;
#else
  std::cout << "Precision: Double (double)" << std::endl;
#endif


  std::cout << "\nFor more information, see develop.md in the source directory." << std::endl;
}

//! \fn bool ParseCommandLine()
//! \brief Parse command line arguments
bool ParseCommandLine(int argc, char* argv[], RegridOptions& options) {
  // Set default values
  options.refinement_factor = 2;
  options.preserve_time = true;
  options.verify_conservation = true;
  options.verify_divergence_b = true;
  options.verbose = false;
  options.show_progress = true;
  options.benchmark = false;
  options.output_file_number = 0;
  options.new_time = 0.0;
  
  // Define long options
  static struct option long_options[] = {
    {"help",            no_argument,       0, 'h'},
    {"version",         no_argument,       0, 'v'},
    {"quiet",           no_argument,       0, 'q'},
    {"verbose",         no_argument,       0, 1001},
    {"benchmark",       no_argument,       0, 'b'},
    {"refine",          required_argument, 0, 'r'},
    {"no-conservation", no_argument,       0, 1003},
    {"no-divergence",   no_argument,       0, 1004},
    {"file-number",     required_argument, 0, 1005},
    {"new-time",        required_argument, 0, 1006},
    {0, 0, 0, 0}
  };
  
  int option_index = 0;
  int c;
  
  while ((c = getopt_long(argc, argv, "hvqbr:", long_options, &option_index)) != -1) {
    switch (c) {
      case 'h':
        PrintUsage(argv[0]);
        exit(EXIT_SUCCESS);
        break;
        
      case 'v':
        PrintVersion();
        exit(EXIT_SUCCESS);
        break;
        
      case 'q':
        options.show_progress = false;
        break;
        
      case 1001: // --verbose
        options.verbose = true;
        break;
        
      case 'b':
        options.benchmark = true;
        break;
        
      case 'r':
        options.refinement_factor = std::atoi(optarg);
        if (options.refinement_factor < 1 || options.refinement_factor > 4) {
          std::cerr << "ERROR: Refinement factor must be between 1 and 4" << std::endl;
          return false;
        }
        break;
        
        
      case 1003: // --no-conservation
        options.verify_conservation = false;
        break;
        
      case 1004: // --no-divergence
        options.verify_divergence_b = false;
        break;
        
      case 1005: // --file-number
        options.output_file_number = std::atoi(optarg);
        break;
        
      case 1006: // --new-time
        options.new_time = std::atof(optarg);
        options.preserve_time = false;
        break;
        
      case '?':
        // getopt_long already printed an error message
        return false;
        
      default:
        std::cerr << "ERROR: Unknown option" << std::endl;
        return false;
    }
  }
  
  // Check for required positional arguments
  if (optind + 2 != argc) {
    std::cerr << "ERROR: Exactly two arguments required: input_file output_file" << std::endl;
    return false;
  }
  
  options.input_file = argv[optind];
  options.output_file = argv[optind + 1];
  
  // Validate file extensions
  if (options.input_file.size() < 4 || 
      options.input_file.substr(options.input_file.size()-4) != ".rst") {
    std::cerr << "WARNING: Input file should have .rst extension" << std::endl;
  }
  
  if (options.output_file.size() < 4 || 
      options.output_file.substr(options.output_file.size()-4) != ".rst") {
    std::cerr << "WARNING: Output file should have .rst extension" << std::endl;
  }
  
  return true;
}