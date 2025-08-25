# /tests/partition/CMakeLists.txt

# Gather all test source files
set(PARTITION_TEST_SOURCES
    ${CMAKE_CURRENT_SOURCE_DIR}/main.cpp
)

# Create test binary
add_executable(partition ${PARTITION_TEST_SOURCES})

# Link Catch2 and core library
target_link_libraries(partition PRIVATE
    Catch2::Catch2WithMain
    linearx_core
)

# Include headers
target_include_directories(partition PRIVATE
    ${PROJECT_SOURCE_DIR}/include
)

# Register tests
catch_discover_tests(partition)
