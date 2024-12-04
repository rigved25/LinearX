#!/bin/bash

# Example run: ./run.sh ./test_seqs/ ./tmp 0 1

CPP_EXECUTABLE="./main"

# Ensure the correct number of arguments
if [ "$#" -ne 4 ]; then
    echo "Usage: $0 <input_dir> <output_dir> <energy_params> <lazy_outside>"
    exit 1
fi

INPUT_DIR=$1
OUTPUT_DIR=$2
ENERGY_PARAMS=$3
LAZY_OUTSIDE=$4

# Check if the C++ executable exists and is executable
if [ ! -x "$CPP_EXECUTABLE" ]; then
    echo "Error: C++ program $CPP_EXECUTABLE does not exist or is not executable."
    exit 1
fi

# Check if the input directory exists
if [ ! -d "$INPUT_DIR" ]; then
    echo "Error: Input directory $INPUT_DIR does not exist."
    exit 1
fi

# Create necessary directories
mkdir -p "$OUTPUT_DIR/stdout"
mkdir -p "$OUTPUT_DIR/stderr"
mkdir -p "$OUTPUT_DIR/logs"

# Process each file in the input directory
for INPUT_FILE in "$INPUT_DIR"/*.fasta*; do
    if [ -f "$INPUT_FILE" ]; then
        echo "Processing $INPUT_FILE"

        BASE_NAME=$(basename "$INPUT_FILE")
        OUTPUT_FILE_STDOUT="$OUTPUT_DIR/stdout/$BASE_NAME"
        OUTPUT_FILE_STDERR="$OUTPUT_DIR/stderr/$BASE_NAME"
        LOG_FILE="$OUTPUT_DIR/logs/$BASE_NAME"

        # Temporary file to capture time output
        TIME_LOG=$(mktemp)

        # Initialize exit status
        EXIT_STATUS=0

        # Run the executable with `time` and redirect output
        /usr/bin/time -l -p -o "$TIME_LOG" "$CPP_EXECUTABLE" "$INPUT_FILE" "$ENERGY_PARAMS" "3" "$LAZY_OUTSIDE" > "$OUTPUT_FILE_STDOUT" 2> "$OUTPUT_FILE_STDERR" || EXIT_STATUS=$?

        # Check the exit status
        if [ "$EXIT_STATUS" -ne 0 ]; then
            echo "Error processing $INPUT_FILE"
            {
                echo "Input File: $INPUT_FILE"
                echo "Error: Program execution failed."
                echo "Exit Status: $EXIT_STATUS"
            } > "$LOG_FILE"
            rm -f "$TIME_LOG"
            continue
        fi

        # Extract elapsed time from time command output
        ELAPSED_TIME=$(grep real "$TIME_LOG" | awk '{print $2}')

        # Clean up temporary file
        rm -f "$TIME_LOG"

        # Write the log file
        {
            echo "Input File: $INPUT_FILE"
            echo "Output File (stdout): $OUTPUT_FILE_STDOUT"
            echo "Output File (stderr): $OUTPUT_FILE_STDERR"
            echo "Elapsed Time (s): $ELAPSED_TIME"
            echo "Exit Status: $EXIT_STATUS"
        } > "$LOG_FILE"

        echo "Processed $INPUT_FILE -> $OUTPUT_FILE_STDOUT"
        echo "Processed $INPUT_FILE -> $OUTPUT_FILE_STDERR"
        echo "Logged results to $LOG_FILE"
    fi
done
