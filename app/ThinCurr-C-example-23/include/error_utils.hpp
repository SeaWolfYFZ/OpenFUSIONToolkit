// error_utils.hpp
// Error handling utilities

#ifndef ERROR_UTILS_HPP
#define ERROR_UTILS_HPP

#include <string>
#include <cstring>
#include <iostream>

inline bool check_error(const char* error_str, const std::string& function_name) {
    if (strlen(error_str) > 0) {
        std::cerr << "ERROR in " << function_name << ": " << error_str << std::endl;
        return true;
    }
    return false;
}

inline void log_info(const std::string& message) {
    std::cout << "< C++ print > " << message << std::endl;
}

inline void log_error(const std::string& message) {
    std::cerr << "< C++ ERROR > " << message << std::endl;
}

#endif // ERROR_UTILS_HPP