#!/bin/bash
set -e

# initialize and update submodules
git submodule update --init --recursive

# setup RNA Eval
cd external/rna-eval
chmod +X ./setup_env.sh
./setup_env.sh

# make build script executable
cd ../../
chmod +x ./build.sh

# make python scripts executable
chmod +x ./linearpartition.py
chmod +x ./linearturbofold.py
