#ifndef PARAMETER_PARSER_HPP_
#define PARAMETER_PARSER_HPP_

#include <string>
#include "common.hpp"

//----------------------------------------------------------------------------------------
// Parameter parsing functionality
//----------------------------------------------------------------------------------------
class ParameterParser {
public:
    static PhysicsConfig ParseParameters(const std::string& param_string);
    
private:
    static bool FindParameterBlock(const std::string& param_string, const std::string& block_name);
    static bool FindParameterName(const std::string& param_string, const std::string& param_name);
    static int ExtractIntegerParameter(const std::string& param_string, const std::string& param_name);
};

#endif // PARAMETER_PARSER_HPP_