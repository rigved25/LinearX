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

# Detect the OS
OS_TYPE=$(uname)

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

        # Run the executable with `time` and capture its PID
        /usr/bin/time -p -o "$TIME_LOG" "$CPP_EXECUTABLE" "$INPUT_FILE" "$ENERGY_PARAMS" "3" "$LAZY_OUTSIDE" > "$OUTPUT_FILE_STDOUT" 2> "$OUTPUT_FILE_STDERR" &
        PID=$!

        PEAK_RSS=0
        PEAK_VSIZE=0
        PEAK_SWAP=0
        PEAK_SHARED=0

        while ps -p $PID > /dev/null 2>&1; do
            if [[ "$OS_TYPE" == "Darwin" ]]; then
                # macOS: Use ps to get memory stats
                MEMORY_INFO=$(ps -o rss,vsz -p $PID | tail -n 1)
                CURRENT_RSS=$(echo "$MEMORY_INFO" | awk '{print $1}')
                CURRENT_VSIZE=$(echo "$MEMORY_INFO" | awk '{print $2}')

                # macOS Swap and Shared Memory
                SWAP_USED=$(memory_pressure | awk '/Swap Used:/ {print $3}' | tr -d 'K')
                SHARED_MEMORY=$(vm_stat | awk '/Pages shared/ {print $3 * 4096 / 1024}')

            elif [[ "$OS_TYPE" == "Linux" ]]; then
                # Linux: Use /proc/<PID>/status
                STATUS_FILE="/proc/$PID/status"
                if [ -f "$STATUS_FILE" ]; then
                    CURRENT_RSS=$(awk '/VmRSS/ {print $2}' $STATUS_FILE)
                    CURRENT_VSIZE=$(awk '/VmSize/ {print $2}' $STATUS_FILE)
                    CURRENT_SWAP=$(awk '/VmSwap/ {print $2}' $STATUS_FILE)
                    SMAPS_FILE="/proc/$PID/smaps"
                    if [ -f "$SMAPS_FILE" ]; then
                        CURRENT_SHARED=$(awk '/Shared_Clean/ {shared += $2} END {print shared}' $SMAPS_FILE)
                    fi
                fi
            else
                echo "Unsupported OS: $OS_TYPE"
                exit 1
            fi

            # Update peak values
            PEAK_RSS=$((CURRENT_RSS > PEAK_RSS ? CURRENT_RSS : PEAK_RSS))
            PEAK_VSIZE=$((CURRENT_VSIZE > PEAK_VSIZE ? CURRENT_VSIZE : PEAK_VSIZE))
            PEAK_SWAP=$((CURRENT_SWAP > PEAK_SWAP ? CURRENT_SWAP : PEAK_SWAP))
            PEAK_SHARED=$((CURRENT_SHARED > PEAK_SHARED ? CURRENT_SHARED : PEAK_SHARED))

            sleep 0.1
        done

        wait $PID
        EXIT_STATUS=$?

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
            echo "Peak RSS (kB): $PEAK_RSS"
            echo "Peak Virtual Memory (VmSize, kB): $PEAK_VSIZE"
            echo "Peak Swap Usage (kB): $PEAK_SWAP"
            echo "Peak Shared Memory (kB): $PEAK_SHARED"
            echo "Exit Status: $EXIT_STATUS"
        } > "$LOG_FILE"

        if [ $EXIT_STATUS -eq 0 ]; then
            echo "Processed $INPUT_FILE -> $OUTPUT_FILE_STDOUT"
            echo "Processed $INPUT_FILE -> $OUTPUT_FILE_STDERR"
            echo "Logged results to $LOG_FILE"
        else
            echo "Error processing $INPUT_FILE"
        fi
    fi
done
