#!/bin/bash

# make sure the script is run from the correct directory
if [ ! -d "build" ]; then
  echo "Please run this script from the root directory of the project."
  exit 1
fi

# run all tests
./build/tests/seq/seq -s --order decl "$@"
