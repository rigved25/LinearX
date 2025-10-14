cmake_minimum_required(VERSION 3.16)
project(linearx CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)
set(CMAKE_INTERPROCEDURAL_OPTIMIZATION ON)  # Enables -flto globally

include(CheckIPOSupported)
check_ipo_supported(RESULT result OUTPUT output)
message(STATUS "IPO/LTO support: ${result} (${output})")

# -------- OpenMP Setup --------
find_package(OpenMP)
if(OpenMP_CXX_FOUND)
    message(STATUS "OpenMP found")
    # For Apple Silicon with Homebrew LLVM
    if(APPLE)
        # Use Homebrew LLVM
        set(CMAKE_C_COMPILER "/opt/homebrew/opt/llvm/bin/clang")
        set(CMAKE_CXX_COMPILER "/opt/homebrew/opt/llvm/bin/clang++")

        # Add OpenMP flags manually
        add_compile_options(-Xpreprocessor -fopenmp)
        add_link_options(-L/opt/homebrew/opt/libomp/lib -lomp)

        # Add include path manually
        include_directories(SYSTEM /opt/homebrew/opt/libomp/include)
    else()
        # Standard OpenMP setup for other systems
        set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
        set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}")
    endif()
else()
    message(STATUS "OpenMP not found")
endif()

# Optimization flags
add_compile_options(
    -Winline
    -O3
    -march=native
    -flto
    -ffast-math
    # -fno-rtti
    -fno-math-errno
    -funroll-loops
    -fomit-frame-pointer
    # -Rpass=inline
    # -Rpass-missed=inline
)

# Optional: disable asserts
# add_compile_definitions(NDEBUG)

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
    src/turbofold/*.cpp
)

# Copy parameter file
configure_file(${PROJECT_SOURCE_DIR}/data/params/fam_hmm_pars.dat
               ${CMAKE_BINARY_DIR}/fam_hmm_pars.txt COPYONLY)

# Core library
add_library(linearx_core ${CORE_SRC})
target_compile_definitions(linearx_core PRIVATE
    PARAM_FILE_PATH="${CMAKE_BINARY_DIR}/fam_hmm_pars.txt"
)

# Check GCC version for filesystem library linking
if(CMAKE_CXX_COMPILER_ID STREQUAL "GNU" AND CMAKE_CXX_COMPILER_VERSION VERSION_LESS 9.0)
    target_link_libraries(linearx_core PRIVATE stdc++fs)
endif()

# Subdirectories
add_subdirectory(external/Catch2)
add_subdirectory(tests)

# Testing
enable_testing()

# Executable for linearpartition
add_executable(linearpartition src/partition/main.cpp)
target_link_libraries(linearpartition PRIVATE linearx_core)

# Executable for linearturbofold
add_executable(linearturbofold src/turbofold/main.cpp)
target_link_libraries(linearturbofold PRIVATE linearx_core)
