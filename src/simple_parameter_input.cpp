#include "simple_parameter_input.hpp"
#include <iostream>
#include <sstream>
#include <algorithm>
#include <iomanip>

SimpleInputLine* SimpleInputBlock::GetLine(const std::string& name) {
    for (auto& line : lines) {
        if (line.param_name == name) {
            return &line;
        }
    }
    return nullptr;
}

SimpleInputBlock* SimpleParameterInput::GetBlock(const std::string& name) {
    for (auto& block : blocks) {
        if (block.block_name == name) {
            return &block;
        }
    }
    return nullptr;
}

void SimpleParameterInput::ParseFromString(const std::string& input) {
    blocks.clear();
    std::istringstream iss(input);
    std::string line;
    SimpleInputBlock* current_block = nullptr;
    
    while (std::getline(iss, line)) {
        // Skip empty lines and comments
        if (line.empty() || line[0] == '#') continue;
        
        // Remove leading/trailing whitespace
        size_t first = line.find_first_not_of(" \t");
        if (first == std::string::npos) continue;
        size_t last = line.find_last_not_of(" \t\r\n");
        line = line.substr(first, last - first + 1);
        
        // Check for block start
        if (line[0] == '<' && line[line.length()-1] == '>') {
            std::string block_name = line.substr(1, line.length()-2);
            
            // Check for end marker
            if (block_name == "par_end") {
                break;
            }
            
            // Skip block end markers
            if (block_name[0] == '/') continue;
            
            // Create new block
            blocks.emplace_back(block_name);
            current_block = &blocks.back();
            continue;
        }
        
        // Parse parameter line
        if (current_block != nullptr) {
            size_t eq_pos = line.find('=');
            if (eq_pos != std::string::npos) {
                std::string name = line.substr(0, eq_pos);
                std::string rest = line.substr(eq_pos + 1);
                
                // Trim whitespace from name
                size_t name_end = name.find_last_not_of(" \t");
                name = name.substr(0, name_end + 1);
                
                // Find value and comment
                std::string value;
                std::string comment;
                
                // Skip leading whitespace in rest
                size_t val_start = rest.find_first_not_of(" \t");
                if (val_start != std::string::npos) {
                    // Find comment start
                    size_t comment_pos = rest.find('#', val_start);
                    if (comment_pos != std::string::npos) {
                        value = rest.substr(val_start, comment_pos - val_start);
                        comment = rest.substr(comment_pos);
                    } else {
                        value = rest.substr(val_start);
                    }
                    
                    // Trim trailing whitespace from value
                    size_t value_end = value.find_last_not_of(" \t");
                    value = value.substr(0, value_end + 1);
                }
                
                // Add to current block
                current_block->lines.emplace_back(name, value, comment);
                
                // Update max lengths for formatting
                current_block->max_len_parname = std::max(current_block->max_len_parname, name.length());
                current_block->max_len_parvalue = std::max(current_block->max_len_parvalue, value.length());
            }
        }
    }
}

bool SimpleParameterInput::SetInteger(const std::string& block_name, const std::string& param_name, int value) {
    SimpleInputBlock* block = GetBlock(block_name);
    if (block == nullptr) return false;
    
    SimpleInputLine* line = block->GetLine(param_name);
    if (line == nullptr) return false;
    
    line->param_value = std::to_string(value);
    block->max_len_parvalue = std::max(block->max_len_parvalue, line->param_value.length());
    return true;
}

std::string SimpleParameterInput::DumpToString() const {
    std::ostringstream oss;
    
    // Write header if needed
    oss << "#------------------------- PAR_DUMP -------------------------" << std::endl;
    
    // Write all blocks
    for (const auto& block : blocks) {
        oss << "<" << block.block_name << ">" << std::endl;
        
        // Write all parameters in the block
        for (const auto& line : block.lines) {
            // Format with padding for alignment
            oss << line.param_name;
            size_t pad_len = block.max_len_parname - line.param_name.length() + 1;
            oss << std::string(pad_len, ' ');
            
            oss << "= " << line.param_value;
            
            // Add padding before comment if exists
            if (!line.param_comment.empty()) {
                pad_len = block.max_len_parvalue - line.param_value.length() + 1;
                oss << std::string(pad_len, ' ');
                oss << line.param_comment;
            }
            
            oss << std::endl;
        }
    }
    
    // Always end with <par_end> and newline
    oss << "<par_end>" << std::endl;
    
    return oss.str();
}