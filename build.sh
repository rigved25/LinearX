#!/bin/bash
set -e  # Exit immediately on error

# Check if we're in a conda environment
if [[ -z "$CONDA_DEFAULT_ENV" ]]; then
    echo "Warning: No conda environment is active."
    echo "Please activate your conda environment with: conda activate <env_name>"
    exit 1
fi

echo "Using conda environment: $CONDA_DEFAULT_ENV"

# Set boost paths for cmake
export CMAKE_PREFIX_PATH="$CONDA_PREFIX:$CMAKE_PREFIX_PATH"

# Configure build files with boost support
cmake -B build -DCMAKE_BUILD_TYPE=Release \
      -DBoost_ROOT="$CONDA_PREFIX" \
      -DBoost_INCLUDE_DIR="$CONDA_PREFIX/include" \
      -DBoost_LIBRARY_DIR="$CONDA_PREFIX/lib"

# Compile the project
cmake --build build --parallel

# Set runtime library path for execution
echo "Build complete!"
echo "To run the executables, make sure LD_LIBRARY_PATH is set:"
echo "export LD_LIBRARY_PATH=\$CONDA_PREFIX/lib:\$LD_LIBRARY_PATH"
