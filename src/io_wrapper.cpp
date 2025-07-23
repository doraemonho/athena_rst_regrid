//========================================================================================
// AthenaK Regridding Tool - Standalone IOWrapper Implementation
// Licensed under the 3-clause BSD License (the "LICENSE")
//========================================================================================
//! \file io_wrapper.cpp
//  \brief Implementation of standalone IOWrapper class

#include "io_wrapper.hpp"
#include <iostream>
#include <cstring>
#include <cerrno>

//----------------------------------------------------------------------------------------
//! \fn int IOWrapper::Open()
//! \brief Open a file for reading or writing
int IOWrapper::Open(const char* filename, FileMode mode) {
  Close(); // Close any existing file

  const char* mode_str = (mode == FileMode::read) ? "rb" : "wb";
  file_ = std::fopen(filename, mode_str);
  
  if (file_ == nullptr) {
    SetError("Failed to open file: " + std::string(filename) + 
             " (" + std::string(std::strerror(errno)) + ")");
    is_open_ = false;
    return -1;
  }

  is_open_ = true;
  last_error_.clear();
  return 0;
}

//----------------------------------------------------------------------------------------
//! \fn void IOWrapper::Close()
//! \brief Close the file
void IOWrapper::Close() {
  if (file_ != nullptr) {
    std::fclose(file_);
    file_ = nullptr;
  }
  is_open_ = false;
  last_error_.clear();
}

//----------------------------------------------------------------------------------------
//! \fn size_t IOWrapper::Read_bytes()
//! \brief Read bytes from file
size_t IOWrapper::Read_bytes(void* buf, size_t size, size_t count) {
  if (!is_open_ || file_ == nullptr) {
    SetError("File not open for reading");
    return 0;
  }

  size_t bytes_read = std::fread(buf, size, count, file_);
  
  if (bytes_read != count && !std::feof(file_)) {
    SetError("Read error: expected " + std::to_string(count) + 
             " elements, got " + std::to_string(bytes_read));
  }

  return bytes_read * size; // Return total bytes read
}

//----------------------------------------------------------------------------------------
//! \fn size_t IOWrapper::Write_bytes()  
//! \brief Write bytes to file
size_t IOWrapper::Write_bytes(const void* buf, size_t size, size_t count) {
  if (!is_open_ || file_ == nullptr) {
    SetError("File not open for writing");
    return 0;
  }

  size_t bytes_written = std::fwrite(buf, size, count, file_);
  
  if (bytes_written != count) {
    SetError("Write error: expected " + std::to_string(count) + 
             " elements, wrote " + std::to_string(bytes_written));
  }

  return bytes_written * size; // Return total bytes written
}

//----------------------------------------------------------------------------------------
//! \fn IOWrapperSizeT IOWrapper::GetPosition()
//! \brief Get current file position
IOWrapperSizeT IOWrapper::GetPosition() {
  if (!is_open_ || file_ == nullptr) {
    SetError("File not open");
    return 0;
  }

  long pos = std::ftell(file_);
  if (pos == -1L) {
    SetError("Failed to get file position");
    return 0;
  }

  return static_cast<IOWrapperSizeT>(pos);
}

//----------------------------------------------------------------------------------------
//! \fn int IOWrapper::Seek()
//! \brief Seek to position from beginning of file
int IOWrapper::Seek(IOWrapperSizeT pos) {
  if (!is_open_ || file_ == nullptr) {
    SetError("File not open");
    return -1;
  }

  int result = std::fseek(file_, static_cast<long>(pos), SEEK_SET);
  if (result != 0) {
    SetError("Failed to seek to position " + std::to_string(pos));
  }

  return result;
}

//----------------------------------------------------------------------------------------
//! \fn int IOWrapper::SeekFromEnd()
//! \brief Seek to position from end of file
int IOWrapper::SeekFromEnd(IOWrapperSizeT offset) {
  if (!is_open_ || file_ == nullptr) {
    SetError("File not open");
    return -1;
  }

  int result = std::fseek(file_, -static_cast<long>(offset), SEEK_END);
  if (result != 0) {
    SetError("Failed to seek from end with offset " + std::to_string(offset));
  }

  return result;
}

//----------------------------------------------------------------------------------------
//! \fn bool IOWrapper::HasError()
//! \brief Check if there are any errors
bool IOWrapper::HasError() const {
  return !last_error_.empty();
}

//----------------------------------------------------------------------------------------
//! \fn std::string IOWrapper::GetLastError()
//! \brief Get the last error message
std::string IOWrapper::GetLastError() const {
  return last_error_;
}