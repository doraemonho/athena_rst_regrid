# CMakeLists.txt for AthenaK MHD Regridding Tool
#=========================================================================================
# AthenaK Regridding Tool
# Licensed under the 3-clause BSD License (the "LICENSE")
#=========================================================================================

cmake_minimum_required(VERSION 3.10)
project(athena_regrid VERSION 1.0 LANGUAGES CXX)

# Disable in-source builds
if( "${CMAKE_SOURCE_DIR}" STREQUAL "${CMAKE_BINARY_DIR}" )
  message(FATAL_ERROR "FATAL: In-source builds disabled. Please create a separate build directory.")
endif()

# C++17 standard required
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

# Default build type
if(NOT CMAKE_BUILD_TYPE)
  message(STATUS "Setting build type to 'Release' as none was specified.")
  set(CMAKE_BUILD_TYPE Release)
endif()

#------ Compile options ----------------------------------------------------------------

option(REGRID_SINGLE_PRECISION "Use single precision" OFF)

#------ Precision Configuration --------------------------------------------------------

if(REGRID_SINGLE_PRECISION)
    add_definitions(-DSINGLE_PRECISION_ENABLED=1)
else()
    add_definitions(-DSINGLE_PRECISION_ENABLED=0)
endif()

#------ Include directories ------------------------------------------------------------

include_directories(${CMAKE_CURRENT_SOURCE_DIR}/src)

#------ Source files -------------------------------------------------------------------

set(REGRID_SOURCES
    src/main.cpp
    src/regrid_driver.cpp
    src/restart_reader.cpp
    src/restart_writer.cpp
    src/mesh_regrid.cpp
    src/data_interpolator.cpp
    src/verification.cpp
    src/io_wrapper.cpp
)

#------ Create executables -------------------------------------------------------------

# Main regridding tool
add_executable(athena_regrid ${REGRID_SOURCES})

target_compile_definitions(athena_regrid PRIVATE
    -DREGRID_TOOL_VERSION="1.0"
)

