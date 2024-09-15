# Compiler
CXX = ccache g++

# Compiler flags
CXXFLAGS = -Wall -std=c++17 -O3 -MMD -MP -march=native

# Target executable
TARGET = main

# Temporary directory for object files
OBJ_DIR = ./../tmp

# Source files with relative paths
SOURCES = ./../../sequence/seq.cpp ./../../energy_model.cpp ./../../partition/partition.cpp ./../../partition/forward.cpp ./../../partition/backward.cpp ./main.cpp

# Create object file paths in the temporary directory
OBJECTS = $(patsubst %.cpp, $(OBJ_DIR)/%.o, $(notdir $(SOURCES)))

# Default target
all: $(TARGET)

# Create the temporary directory if it doesn't exist
$(OBJ_DIR):
	mkdir -p $(OBJ_DIR)

# Link object files to create the executable
$(TARGET): $(OBJ_DIR) $(OBJECTS)
	$(CXX) $(CXXFLAGS) -o $(TARGET) $(OBJECTS)

# Compile .cpp files to object files in the temporary directory
$(OBJ_DIR)/%.o: $(SOURCES)
	$(CXX) $(CXXFLAGS) -c $(filter %$*.cpp,$(SOURCES)) -o $@

# Clean up object files and executable
clean:
	rm -rf $(OBJ_DIR) $(TARGET)

# Run the tests
run: $(TARGET)
	./$(TARGET)