#!/bin/bash
set -e

# Function to check if we're in a conda environment
check_conda_env() {
    if [[ -z "$CONDA_DEFAULT_ENV" ]]; then
        echo "Warning: No conda environment is active."
        echo "Please activate your conda environment or create one with:"
        echo "  conda create -n linearx python=3.8 -c conda-forge"
        echo "  conda activate linearx"
        echo "Then re-run this script."
        exit 1
    fi
    echo "Using conda environment: $CONDA_DEFAULT_ENV"
}

# Check conda environment
check_conda_env

# Install boost if not already installed
echo "Checking for boost..."
if ! conda list | grep -q "^boost "; then
    echo "Installing boost..."
    conda install -c conda-forge boost -y
else
    echo "Boost is already installed."
fi

# Install cmake if not available
if ! command -v cmake &> /dev/null; then
    echo "Installing cmake..."
    conda install -c conda-forge cmake -y
fi

# initialize and update submodules
echo "Initializing submodules..."
git submodule update --init --recursive

# setup RNA Eval
echo "Setting up RNA Eval..."
cd external/rna-eval
chmod +x ./setup_env.sh
./setup_env.sh

# make build script executable
cd ../../
chmod +x ./build.sh

# make python scripts executable
chmod +x ./linearpartition.py
chmod +x ./linearturbofold.py

echo "Installation complete!"
echo "You can now run: ./build.sh"
