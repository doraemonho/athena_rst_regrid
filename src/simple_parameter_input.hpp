#ifndef SIMPLE_PARAMETER_INPUT_HPP_
#define SIMPLE_PARAMETER_INPUT_HPP_

#include <string>
#include <list>
#include <cstddef>

//----------------------------------------------------------------------------------------
// Simple parameter input parser for restart file parameter modification
// Based on AthenaK's ParameterInput design but simplified for our needs
//----------------------------------------------------------------------------------------

struct SimpleInputLine {
    std::string param_name;
    std::string param_value;
    std::string param_comment;
    
    SimpleInputLine(const std::string& name, const std::string& value, const std::string& comment = "")
        : param_name(name), param_value(value), param_comment(comment) {}
};

class SimpleInputBlock {
public:
    std::string block_name;
    std::list<SimpleInputLine> lines;
    std::size_t max_len_parname = 0;
    std::size_t max_len_parvalue = 0;
    
    SimpleInputBlock(const std::string& name) : block_name(name) {}
    
    SimpleInputLine* GetLine(const std::string& name);
};

class SimpleParameterInput {
private:
    std::list<SimpleInputBlock> blocks;
    
    SimpleInputBlock* GetBlock(const std::string& name);
    
public:
    void ParseFromString(const std::string& input);
    bool SetInteger(const std::string& block_name, const std::string& param_name, int value);
    std::string DumpToString() const;
};

#endif // SIMPLE_PARAMETER_INPUT_HPP_