#include "parameter_parser.hpp"
#include <iostream>

PhysicsConfig ParameterParser::ParseParameters(const std::string& param_string) {
    std::cout << "Parsing parameters..." << std::endl;
    
    PhysicsConfig config;
    
    // Look for hydro module
    if (FindParameterBlock(param_string, "hydro")) {
        config.has_hydro = true;
        config.nhydro = 5;  // Default: density, momentum (3), energy
        std::cout << "  Found <hydro> block" << std::endl;
    }
    
    // Look for MHD module - this is the primary physics for most AthenaK simulations
    if (FindParameterBlock(param_string, "mhd")) {
        config.has_mhd = true;
        std::cout << "  Found <mhd> block" << std::endl;
        
        // For MHD simulations, check if hydro is separate or included
        if (!config.has_hydro) {
            // MHD-only: includes hydro variables + magnetic field variables
            config.nmhd = 5;  // ideal MHD : density, momentum(3), energy
            std::cout << "    MHD-only configuration (includes hydro variables)" << std::endl;
        } else {
            // Separate MHD module (less common)
            config.nmhd = 3;  // Just magnetic field components
            std::cout << "    Separate MHD module (magnetic fields only)" << std::endl;
        }
    }
    
    // Look for turbulence driving - multiple possible parameter names
    if (FindParameterName(param_string, "turb_flag") || 
        FindParameterName(param_string, "turb_driving") ||
        FindParameterBlock(param_string, "turb_driving") ||
        FindParameterName(param_string, "driving_type")) {
        config.has_turbulence = true;
        config.nforce = 3;  // 3D force components
        std::cout << "  Found turbulence driving parameters" << std::endl;
    }
    
    // Look for radiation module
    if (FindParameterBlock(param_string, "radiation")) {
        config.has_radiation = true;
        config.nrad = 1;  // Simplified - actual value depends on radiation configuration
        std::cout << "  Found <radiation> block" << std::endl;
    }
    
    // Look for Z4c spacetime evolution
    if (FindParameterBlock(param_string, "z4c")) {
        config.has_z4c = true;
        config.nz4c = 25;  // Standard Z4c variable count
        std::cout << "  Found <z4c> block" << std::endl;
    }
    
    // Look for ADM spacetime evolution  
    if (FindParameterBlock(param_string, "adm")) {
        config.has_adm = true;
        config.nadm = 12;  // Standard ADM variable count
        std::cout << "  Found <adm> block" << std::endl;
    }
    
    // Look for compact objects (black holes, neutron stars)
    if (FindParameterName(param_string, "compact_objects") ||
        FindParameterName(param_string, "num_co")) {
        config.nco = ExtractIntegerParameter(param_string, "num_co");
        if (config.nco > 0) {
            std::cout << "  Found " << config.nco << " compact objects" << std::endl;
        }
    }
    
    // Validate physics configuration
    if (!config.has_hydro && !config.has_mhd) {
        std::cout << "  WARNING: No hydro or MHD physics detected!" << std::endl;
        std::cout << "  Using fallback: MHD + turbulence configuration" << std::endl;
        config.has_mhd = true;
        config.has_turbulence = true;
        config.nmhd = 8;
        config.nforce = 3;
    }
    
    // Final configuration summary
    std::cout << std::endl << "Final physics configuration:" << std::endl;
    std::cout << "  Hydro: " << (config.has_hydro ? "Yes" : "No") << " (nhydro=" << config.nhydro << ")" << std::endl;
    std::cout << "  MHD: " << (config.has_mhd ? "Yes" : "No") << " (nmhd=" << config.nmhd << ")" << std::endl;
    std::cout << "  Turbulence: " << (config.has_turbulence ? "Yes" : "No") << " (nforce=" << config.nforce << ")" << std::endl;
    std::cout << "  Radiation: " << (config.has_radiation ? "Yes" : "No") << " (nrad=" << config.nrad << ")" << std::endl;
    std::cout << "  Z4c: " << (config.has_z4c ? "Yes" : "No") << " (nz4c=" << config.nz4c << ")" << std::endl;
    std::cout << "  ADM: " << (config.has_adm ? "Yes" : "No") << " (nadm=" << config.nadm << ")" << std::endl;
    std::cout << "  Compact Objects: " << config.nco << std::endl;
    
    return config;
}

bool ParameterParser::FindParameterBlock(const std::string& param_string, const std::string& block_name) {
    std::string block_start = "<" + block_name + ">";
    return param_string.find(block_start) != std::string::npos;
}

bool ParameterParser::FindParameterName(const std::string& param_string, const std::string& param_name) {
    return param_string.find(param_name) != std::string::npos;
}

int ParameterParser::ExtractIntegerParameter(const std::string& param_string, const std::string& param_name) {
    size_t pos = param_string.find(param_name);
    if (pos != std::string::npos) {
        size_t eq_pos = param_string.find("=", pos);
        if (eq_pos != std::string::npos) {
            std::string num_str = param_string.substr(eq_pos + 1, 10);
            try {
                return std::stoi(num_str);
            } catch (...) {
                return 0;
            }
        }
    }
    return 0;
}