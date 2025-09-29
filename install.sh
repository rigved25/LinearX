#!/bin/bash
set -e

# initialize and update submodules
git submodule update --init --recursive

# make build script executable
chmod +x ./build.sh