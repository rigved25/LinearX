# /tests/CMakeLists.txt

# Enable test support
enable_testing()
include(CTest)

# Strip long paths from __FILE__ in test output
add_compile_options(-ffile-prefix-map=${CMAKE_SOURCE_DIR}=.)

# Add Catch2 test discovery support
include(${PROJECT_SOURCE_DIR}/external/Catch2/extras/Catch.cmake)

# Add individual test modules
add_subdirectory(seq)
add_subdirectory(alignment)
add_subdirectory(partition)
