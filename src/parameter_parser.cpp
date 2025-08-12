#include "parameter_parser.hpp"
#include <iostream>

PhysicsConfig ParameterParser::ParseParameters(const std::string& param_string) {
    std::cout << "Parsing parameters..." << std::endl;
    
    PhysicsConfig config;
    
    // Look for hydro module
    if (FindParameterBlock(param_string, "hydro")) {
        config.has_hydro = true;
        config.nhydro = 5;  // Default: density, momentum (3), energy
    }
    
    // Look for MHD module - this is the primary physics for most AthenaK simulations
    if (FindParameterBlock(param_string, "mhd")) {
        config.has_mhd = true;
        
        // For MHD simulations, check if hydro is separate or included
        if (!config.has_hydro) {
            // MHD-only: includes hydro variables + magnetic field variables
            config.nmhd = 5;  // ideal MHD : density, momentum(3), energy
        } else {
            // Separate MHD module (less common)
            config.nmhd = 3;  // Just magnetic field components
        }
    }
    
    // Look for turbulence driving - multiple possible parameter names
    if (FindParameterName(param_string, "turb_flag") || 
        FindParameterName(param_string, "turb_driving") ||
        FindParameterBlock(param_string, "turb_driving") ||
        FindParameterName(param_string, "driving_type")) {
        config.has_turbulence = true;
        config.nforce = 3;  // 3D force components
    }
    
    // Look for radiation module
    if (FindParameterBlock(param_string, "radiation")) {
        config.has_radiation = true;
        config.nrad = 1;  // Simplified - actual value depends on radiation configuration
    }
    
    // Look for Z4c spacetime evolution
    if (FindParameterBlock(param_string, "z4c")) {
        config.has_z4c = true;
        config.nz4c = 25;  // Standard Z4c variable count
    }
    
    // Look for ADM spacetime evolution  
    if (FindParameterBlock(param_string, "adm")) {
        config.has_adm = true;
        config.nadm = 12;  // Standard ADM variable count
    }
    
    // Look for compact objects (black holes, neutron stars)
    if (FindParameterName(param_string, "compact_objects") ||
        FindParameterName(param_string, "num_co")) {
        config.nco = ExtractIntegerParameter(param_string, "num_co");
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