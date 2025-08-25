#!/bin/bash
set -e  # Exit immediately on error

# Configure build files
cmake -B build -DCMAKE_BUILD_TYPE=Release

# Compile the project
cmake --build build --parallel
