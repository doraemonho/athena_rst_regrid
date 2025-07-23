#ifndef IO_WRAPPER_HPP_
#define IO_WRAPPER_HPP_
//========================================================================================
// AthenaK Regridding Tool - Standalone IOWrapper
// Licensed under the 3-clause BSD License (the "LICENSE")
//========================================================================================
//! \file io_wrapper.hpp
//  \brief Standalone IOWrapper class for regridding tool

#include <cstdio>
#include <cstdint>
#include <string>

//----------------------------------------------------------------------------------------
//! Type definitions
using IOWrapperSizeT = std::uint64_t;

//----------------------------------------------------------------------------------------
//! \class IOWrapper
//! \brief Standalone file I/O wrapper class
class IOWrapper {
 public:
  enum class FileMode {
    read,
    write
  };

  IOWrapper() : file_(nullptr), is_open_(false) {}
  ~IOWrapper() { Close(); }

  // File operations
  int Open(const char* filename, FileMode mode);
  void Close();
  bool IsOpen() const { return is_open_; }

  // Read operations
  size_t Read_bytes(void* buf, size_t size, size_t count);
  template<typename T>
  size_t Read_any_type(T* data, size_t count) {
    return Read_bytes(data, sizeof(T), count);
  }

  // Write operations  
  size_t Write_bytes(const void* buf, size_t size, size_t count);
  // AthenaK-compatible interface: count is total bytes, not element count  
  size_t Write_any_type(const void* data, size_t byte_count, const char* type_name) {
    (void)type_name; // Ignore type name in standalone version
    return Write_bytes(data, 1, byte_count);
  }

  // Position operations
  IOWrapperSizeT GetPosition();
  int Seek(IOWrapperSizeT pos);
  int SeekFromEnd(IOWrapperSizeT offset);

  // Error checking
  bool HasError() const;
  std::string GetLastError() const;

 private:
  FILE* file_;
  bool is_open_;
  mutable std::string last_error_;

  void SetError(const std::string& error) const { last_error_ = error; }
};

#endif // IO_WRAPPER_HPP_