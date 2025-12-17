#include "parameter_parser.hpp"
#include <iostream>
#include <sstream>
#include <cctype>

PhysicsConfig ParameterParser::ParseParameters(const std::string& param_string) {
    std::cout << "Parsing parameters..." << std::endl;
    
    PhysicsConfig config;
    
    // Look for hydro module
    if (FindParameterBlock(param_string, "hydro")) {
        config.has_hydro = true;
        int nhydro_base = 5;
        auto eos_opt = ExtractStringInBlock(param_string, "hydro", "eos");
        if (eos_opt.has_value() && *eos_opt == "isothermal") {
            nhydro_base = 4;
        }
        int nscalars = 0;
        auto nscalars_opt = ExtractStringInBlock(param_string, "hydro", "nscalars");
        if (nscalars_opt.has_value()) {
            try {
                nscalars = std::stoi(*nscalars_opt);
            } catch (...) {
                nscalars = 0;
            }
        }
        config.nhydro = nhydro_base + nscalars;
    }
    
    // Look for MHD module - this is the primary physics for most AthenaK simulations
    if (FindParameterBlock(param_string, "mhd")) {
        config.has_mhd = true;
        int nmhd_base = 5;
        auto eos_opt = ExtractStringInBlock(param_string, "mhd", "eos");
        if (eos_opt.has_value() && *eos_opt == "isothermal") {
            nmhd_base = 4;
        }
        int nscalars = 0;
        auto nscalars_opt = ExtractStringInBlock(param_string, "mhd", "nscalars");
        if (nscalars_opt.has_value()) {
            try {
                nscalars = std::stoi(*nscalars_opt);
            } catch (...) {
                nscalars = 0;
            }
        }
        config.nmhd = nmhd_base + nscalars;
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

bool ParameterParser::FindParameterBlock(const std::string& param_string,
                                        const std::string& block_name) {
    std::string block_start = "<" + block_name + ">";
    return param_string.find(block_start) != std::string::npos;
}

bool ParameterParser::FindParameterName(const std::string& param_string,
                                       const std::string& param_name) {
    return param_string.find(param_name) != std::string::npos;
}

int ParameterParser::ExtractIntegerParameter(const std::string& param_string,
                                            const std::string& param_name) {
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

std::optional<std::string> ParameterParser::FindBlockText(
    const std::string& param_string, const std::string& block_name) {
  const std::string block_start = "<" + block_name + ">";
  const size_t start_pos = param_string.find(block_start);
  if (start_pos == std::string::npos) {
    return std::nullopt;
  }

  size_t end_pos = param_string.find("\n<", start_pos + block_start.size());
  if (end_pos == std::string::npos) {
    end_pos = param_string.size();
  } else {
    end_pos += 1;  // include newline before next block for line parsing
  }
  return param_string.substr(start_pos, end_pos - start_pos);
}

std::optional<std::string> ParameterParser::ExtractStringInBlock(
    const std::string& param_string, const std::string& block_name,
    const std::string& param_name) {
  auto block_text_opt = FindBlockText(param_string, block_name);
  if (!block_text_opt.has_value()) {
    return std::nullopt;
  }

  std::istringstream block_stream(*block_text_opt);
  std::string line;
  while (std::getline(block_stream, line)) {
    if (line.rfind(param_name, 0) != 0) {
      continue;
    }
    const size_t eq_pos = line.find('=');
    if (eq_pos == std::string::npos) {
      continue;
    }
    size_t value_pos = eq_pos + 1;
    while (value_pos < line.size()
           && std::isspace(static_cast<unsigned char>(line[value_pos])) != 0) {
      ++value_pos;
    }
    size_t value_end = value_pos;
    while (value_end < line.size()
           && std::isspace(static_cast<unsigned char>(line[value_end])) == 0
           && line[value_end] != '#') {
      ++value_end;
    }
    if (value_end <= value_pos) {
      return std::nullopt;
    }
    return line.substr(value_pos, value_end - value_pos);
  }
  return std::nullopt;
}

std::optional<Real> ParameterParser::ExtractRealInBlock(
    const std::string& param_string, const std::string& block_name,
    const std::string& param_name) {
  auto value_opt = ExtractStringInBlock(param_string, block_name, param_name);
  if (!value_opt.has_value()) {
    return std::nullopt;
  }
  try {
    return std::stod(*value_opt);
  } catch (...) {
    return std::nullopt;
  }
}

std::string ParameterParser::ReplaceMeshNx(const std::string& param_string, int nx1,
                                          int nx2, int nx3) {
  std::string updated = param_string;
  const size_t mesh_start = updated.find("<mesh>");
  if (mesh_start == std::string::npos) {
    return updated;
  }

  size_t mesh_end = updated.find("\n<meshblock>", mesh_start);
  if (mesh_end == std::string::npos) {
    mesh_end = updated.find("\n<", mesh_start + 6);
  }
  if (mesh_end == std::string::npos) {
    mesh_end = updated.size();
  }

  const std::string mesh_text = updated.substr(mesh_start, mesh_end - mesh_start);
  std::istringstream mesh_stream(mesh_text);
  std::string rebuilt;
  std::string line;

  auto replace_value = [&](const std::string& name, int value,
                           const std::string& in_line) -> std::string {
    if (in_line.rfind(name, 0) != 0) {
      return in_line;
    }
    const size_t eq_pos = in_line.find('=');
    if (eq_pos == std::string::npos) {
      return in_line;
    }
    size_t value_pos = eq_pos + 1;
    while (value_pos < in_line.size()
           && std::isspace(static_cast<unsigned char>(in_line[value_pos])) != 0) {
      ++value_pos;
    }
    size_t value_end = value_pos;
    while (value_end < in_line.size()
           && std::isspace(static_cast<unsigned char>(in_line[value_end])) == 0
           && in_line[value_end] != '#') {
      ++value_end;
    }
    std::string out_line = in_line;
    out_line.replace(value_pos, value_end - value_pos, std::to_string(value));
    return out_line;
  };

  while (std::getline(mesh_stream, line)) {
    line = replace_value("nx1", nx1, line);
    line = replace_value("nx2", nx2, line);
    line = replace_value("nx3", nx3, line);
    rebuilt.append(line);
    rebuilt.push_back('\n');
  }

  updated.replace(mesh_start, mesh_end - mesh_start, rebuilt);
  return updated;
}
