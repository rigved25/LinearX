#!/bin/bash
set -e

# initialize and update submodules
git submodule update --init --recursive
cd external

# setup RNA Eval
cd rna-eval
chmod +X ./setup_env.sh

# make build script executable
chmod +x ./build.sh
