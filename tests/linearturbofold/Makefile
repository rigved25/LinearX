# Compiler
CXX = ccache g++

# Compiler flags
CXXFLAGS = -Wall -std=c++17 -O3 -MMD -MP -march=native

# Target executable
TARGET = main

# Temporary directory for object files
OBJ_DIR = ./../tmp

# Source files with relative paths
SOURCES = ./../../sequence/seq.cpp ./../../sequence/multi_seq.cpp ./../../energy_model.cpp \
		  ./../../partition/partition.cpp ./../../partition/forward.cpp ./../../partition/backward.cpp \
		  ./../../linear_align/linear_align.cpp ./../../linear_align/backward.cpp ./../../linear_align/phmm.cpp \
		  ./../../turbofold/turbofold.cpp ./../../turbofold/turbo_align.cpp ./../../turbofold/turbo_partition.cpp \
		  ./../../structure/structure.cpp ./main.cpp

# Create object file paths in the temporary directory, preserving directory structure
OBJECTS = $(SOURCES:./../../%.cpp=$(OBJ_DIR)/%.o)

# Default target
all: $(TARGET)

# Create the temporary directory if it doesn't exist
$(OBJ_DIR):
	mkdir -p $(OBJ_DIR)

# Ensure subdirectories exist before compiling
$(OBJ_DIR)/%.o: ./../../%.cpp
	@mkdir -p $(dir $@)  # Ensure the subdirectory exists
	$(CXX) $(CXXFLAGS) -c $< -o $@

# Link object files to create the executable
$(TARGET): $(OBJ_DIR) $(OBJECTS)
	$(CXX) $(CXXFLAGS) -o $(TARGET) $(OBJECTS)

# Clean up object files and executable
clean:
	rm -rf $(OBJ_DIR) $(TARGET)

# Run the tests
run: $(TARGET)
	./$(TARGET)