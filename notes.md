cmake -B build        # Configure
cmake --build build   # Build

// to create symlink for tests binary in the root directory
ln -sf build/tests/tests test

// to run test
./test -s -r compact --order decl
(r for compact output format)
