cmake_minimum_required(VERSION 3.16)
project(linearx CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)
set(CMAKE_INTERPROCEDURAL_OPTIMIZATION ON)  # Enables -flto globally

# Optimization flags
add_compile_options(
    -Winline
    -O3
    -march=native
    -flto
    -ffast-math
    -fno-math-errno
    -funroll-loops
    -fomit-frame-pointer
)

# Optional: disable asserts
add_compile_definitions(NDEBUG)

# Include paths
include_directories(
    /opt/homebrew/include/
    ${PROJECT_SOURCE_DIR}/include
    ${PROJECT_SOURCE_DIR}/external/eigen
)

# Source files
file(GLOB_RECURSE CORE_SRC
    src/alignment/*.cpp
    src/partition/*.cpp
    src/sequence/*.cpp
)

# Copy parameter file
configure_file(${PROJECT_SOURCE_DIR}/data/params/fam_hmm_pars.dat
               ${CMAKE_BINARY_DIR}/fam_hmm_pars.txt COPYONLY)

# Core library
add_library(linearx_core ${CORE_SRC})
target_compile_definitions(linearx_core PRIVATE
    PARAM_FILE_PATH="${CMAKE_BINARY_DIR}/fam_hmm_pars.txt"
)

# Subdirectories
add_subdirectory(external/Catch2)
add_subdirectory(tests)

# Testing
enable_testing()
