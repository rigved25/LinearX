#!/bin/bash
set -e

# initialize and update submodules
git submodule update --init --recursive

# setup RNA Eval
chmod +X ./external/rna-eval/setup_env.sh
./external/rna-eval/setup_env.sh

# make build script executable
chmod +x ./build.sh
