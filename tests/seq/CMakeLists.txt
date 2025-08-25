# /tests/seq/CMakeLists.txt

# Gather all test source files
set(SEQ_TEST_SOURCES
    ${CMAKE_CURRENT_SOURCE_DIR}/seq.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/mseq.cpp
)

# Create test binary
add_executable(seq ${SEQ_TEST_SOURCES})

# Link Catch2 and core library
target_link_libraries(seq PRIVATE
    Catch2::Catch2WithMain
    linearx_core
)

# Include headers
target_include_directories(seq PRIVATE
    ${PROJECT_SOURCE_DIR}/include
)

# Register tests
catch_discover_tests(seq)
