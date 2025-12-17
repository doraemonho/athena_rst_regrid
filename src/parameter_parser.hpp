#ifndef PARAMETER_PARSER_HPP_
#define PARAMETER_PARSER_HPP_

#include <string>
#include <optional>
#include "common.hpp"

//----------------------------------------------------------------------------------------
// Parameter parsing functionality
//----------------------------------------------------------------------------------------
class ParameterParser {
public:
    static PhysicsConfig ParseParameters(const std::string& param_string);

    static std::optional<std::string> ExtractStringInBlock(
        const std::string& param_string, const std::string& block_name,
        const std::string& param_name);

    static std::optional<Real> ExtractRealInBlock(const std::string& param_string,
                                                 const std::string& block_name,
                                                 const std::string& param_name);

    static std::string ReplaceMeshNx(const std::string& param_string, int nx1, int nx2,
                                     int nx3);
    
private:
    static bool FindParameterBlock(const std::string& param_string, const std::string& block_name);
    static bool FindParameterName(const std::string& param_string, const std::string& param_name);
    static int ExtractIntegerParameter(const std::string& param_string, const std::string& param_name);

    static std::optional<std::string> FindBlockText(const std::string& param_string,
                                                   const std::string& block_name);
};

#endif // PARAMETER_PARSER_HPP_
