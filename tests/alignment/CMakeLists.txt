# /tests/alignment/CMakeLists.txt

# Gather all test source files
set(ALIGNMENT_TEST_SOURCES
    ${CMAKE_CURRENT_SOURCE_DIR}/main.cpp
)

# Create test binary
add_executable(alignment ${ALIGNMENT_TEST_SOURCES})

# Link Catch2 and core library
target_link_libraries(alignment PRIVATE
    Catch2::Catch2WithMain
    linearx_core
)

# Include headers
target_include_directories(alignment PRIVATE
    ${PROJECT_SOURCE_DIR}/include
)

# Register tests
catch_discover_tests(alignment)
