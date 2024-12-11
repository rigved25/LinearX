# Compiler
CXX = ccache g++

# Compiler flags
CXXFLAGS = -Wall -std=c++17 -O3 -MMD -MP -march=native

# Targets
TARGETS = main1 main2

# Temporary directory for object files
OBJ_DIR = ./../tmp

# Source files shared by both targets
COMMON_SOURCES = ./../../sequence/seq.cpp ./../../sequence/multi_seq.cpp ./../../energy_model.cpp \
                 ./../../partition/partition.cpp ./../../partition/forward.cpp ./../../partition/backward.cpp \
				 ./../../structure/structure.cpp
				 

# Per-target sources
MAIN1_SOURCES = ./main1.cpp
MAIN2_SOURCES = ./main2.cpp

# Object files for each target
MAIN1_OBJECTS = $(COMMON_SOURCES:./../../%.cpp=$(OBJ_DIR)/%.o) $(MAIN1_SOURCES:./%.cpp=$(OBJ_DIR)/%.o)
MAIN2_OBJECTS = $(COMMON_SOURCES:./../../%.cpp=$(OBJ_DIR)/%.o) $(MAIN2_SOURCES:./%.cpp=$(OBJ_DIR)/%.o)

# Default target (build both executables)
all: $(TARGETS)

# Create the temporary directory if it doesn't exist
$(OBJ_DIR):
	mkdir -p $(OBJ_DIR)

# Build main1
main1: $(OBJ_DIR) $(MAIN1_OBJECTS)
	$(CXX) $(CXXFLAGS) -o $@ $(MAIN1_OBJECTS)

# Build main2
main2: $(OBJ_DIR) $(MAIN2_OBJECTS)
	$(CXX) $(CXXFLAGS) -o $@ $(MAIN2_OBJECTS)

# Compile .cpp files to object files in the temporary directory
$(OBJ_DIR)/%.o: ./../../%.cpp
	@mkdir -p $(dir $@)  # Ensure the subdirectory exists
	$(CXX) $(CXXFLAGS) -c $< -o $@

$(OBJ_DIR)/%.o: ./%.cpp
	@mkdir -p $(dir $@)  # Ensure the subdirectory exists
	$(CXX) $(CXXFLAGS) -c $< -o $@

# Clean up object files and executables
clean:
	rm -rf $(OBJ_DIR) $(TARGETS)

# Run main1
run_main1: main1
	./main1

# Run main2
run_main2: main2
	./main2
	