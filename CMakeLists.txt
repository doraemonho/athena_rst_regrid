cmake_minimum_required(VERSION 3.10)
project(athena_restart_reader)

# Set C++ standard
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Build type configuration
if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE Release CACHE STRING "Build type" FORCE)
endif()

# Compiler flags
set(CMAKE_CXX_FLAGS_DEBUG "-g3 -ggdb -O0 -Wall -Wextra -DDEBUG")
set(CMAKE_CXX_FLAGS_RELEASE "-O3 -DNDEBUG")

# Options
option(ENABLE_MPI "Enable MPI parallelization" ON)

# MPI Configuration
if(ENABLE_MPI)
    find_package(MPI REQUIRED)
    add_compile_definitions(MPI_PARALLEL_ENABLED)
    message(STATUS "MPI support enabled")
    message(STATUS "MPI C++ compiler: ${MPI_CXX_COMPILER}")
else()
    message(STATUS "MPI support disabled")
endif()

# Source files
set(SOURCES
    src/main.cpp
    src/restart_reader.cpp
    src/io_wrapper.cpp
    src/parameter_parser.cpp
    src/mpi_distribution.cpp
    src/physics_reader.cpp
    src/upscaler.cpp
    src/restart_writer.cpp
    src/simple_parameter_input.cpp
    src/downsampler.cpp
    src/binary_writer.cpp
)

# Header files
set(HEADERS
    src/common.hpp
    src/restart_reader.hpp
    src/io_wrapper.hpp
    src/parameter_parser.hpp
    src/mpi_distribution.hpp
    src/physics_reader.hpp
    src/prolongation.hpp
    src/upscaler.hpp
    src/restart_writer.hpp
    src/simple_parameter_input.hpp
    src/downsampler.hpp
    src/binary_writer.hpp
)

# Include directories
include_directories(${CMAKE_CURRENT_SOURCE_DIR}/src)

# Create main executable
add_executable(restart_reader ${SOURCES})

# Create test executable for MPI write (exclude main.cpp from SOURCES)
set(TEST_SOURCES ${SOURCES})
list(REMOVE_ITEM TEST_SOURCES src/main.cpp)

# Create static library (optional, for linking with other projects)
add_library(restart_reader_lib STATIC ${SOURCES})
target_include_directories(restart_reader_lib PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/src)

# MPI linking
if(ENABLE_MPI)
    target_link_libraries(restart_reader PRIVATE MPI::MPI_CXX)
    target_link_libraries(restart_reader_lib PUBLIC MPI::MPI_CXX)
endif()

# Installation
install(TARGETS restart_reader 
        RUNTIME DESTINATION bin)
install(TARGETS restart_reader_lib 
        ARCHIVE DESTINATION lib)
install(FILES ${HEADERS} 
        DESTINATION include/restart_reader)
