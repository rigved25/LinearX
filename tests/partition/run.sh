#!/bin/bash

# example run: ./run.sh seqs tmp_no_lzy 0 3>&1

# Path to the compiled C++ executable
CPP_EXECUTABLE="./main1"  # Replace with the actual path to your C++ executable

# Ensure the correct number of arguments are passed
if [ "$#" -ne 3 ]; then
    echo "Usage: $0 <input_dir> <output_dir> <lazy_outside>"
    exit 1
fi

# Parse arguments
INPUT_DIR=$1
OUTPUT_DIR=$2
LAZY_OUTSIDE=$3

# Check if the C++ program exists and is executable
if [ ! -x "$CPP_EXECUTABLE" ]; then
    echo "Error: C++ program $CPP_EXECUTABLE does not exist or is not executable."
    exit 1
fi

# Check if input directory exists
if [ ! -d "$INPUT_DIR" ]; then
    echo "Error: Input directory $INPUT_DIR does not exist."
    exit 1
fi

# Create the output directory if it does not exist
mkdir -p "$OUTPUT_DIR/stdout"
mkdir -p "$OUTPUT_DIR/stderr"

# Process each file in the input directory
for INPUT_FILE in "$INPUT_DIR"/*.fasta*; do
    if [ -f "$INPUT_FILE" ]; then
        echo "Processing $INPUT_FILE"

        # Extract the base name of the file
        BASE_NAME=$(basename "$INPUT_FILE")
        
        # Define the output file path
        OUTPUT_FILE_STDOUT="$OUTPUT_DIR/stdout/$BASE_NAME"
        OUTPUT_FILE_STDERR="$OUTPUT_DIR/stderr/$BASE_NAME"

        # Execute the C++ program and store the output
        "$CPP_EXECUTABLE" "$INPUT_FILE" "$LAZY_OUTSIDE" > "$OUTPUT_FILE_STDOUT" 2> "$OUTPUT_FILE_STDERR"
        if [ $? -eq 0 ]; then
            echo "Processed $INPUT_FILE -> $OUTPUT_FILE_STDOUT"
            echo "Processed $INPUT_FILE -> $OUTPUT_FILE_STDERR"
        else
            echo "Error processing $INPUT_FILE"
        fi
    fi
done
