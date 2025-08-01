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
set(CMAKE_CXX_FLAGS_DEBUG "-g -O0 -Wall -Wextra")
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
    main.cpp
    restart_reader.cpp
    io_wrapper.cpp
    parameter_parser.cpp
    mpi_distribution.cpp
    physics_reader.cpp
    upscaler.cpp
    restart_writer.cpp
    simple_parameter_input.cpp
)

# Header files
set(HEADERS
    common.hpp
    restart_reader.hpp
    io_wrapper.hpp
    parameter_parser.hpp
    mpi_distribution.hpp
    physics_reader.hpp
    prolongation.hpp
    upscaler.hpp
    restart_writer.hpp
    simple_parameter_input.hpp
)

# Create main executable
add_executable(restart_reader ${SOURCES})

# Create static library (optional, for linking with other projects)
add_library(restart_reader_lib STATIC ${SOURCES})
target_include_directories(restart_reader_lib PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})

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