import os
import re
import matplotlib.pyplot as plt
from collections import defaultdict
import difflib
import argparse

# Extract time from log file
def extract_time_from_log(log_file):
    with open(log_file, "r") as f:
        content = f.readlines()

    # Join all lines into a single string for regex matching
    content_str = "".join(content)

    # Extract time in both possible formats
    real_time_match = re.search(r"(\d+\.\d+)\s+real", content_str)
    elapsed_time_match = re.search(r"Elapsed Time \(s\):\s*(\d+\.\d+)", content_str)

    # Use the first match found
    if real_time_match:
        return float(real_time_match.group(1))
    elif elapsed_time_match:
        return float(elapsed_time_match.group(1))
    else:
        return 0  # Default to 0 if no match is found

# Process log files in a directory
def process_log_dir(log_dir):
    time_data = defaultdict(lambda: {"real": 0, "count": 0, "values": []})
    for log_file in os.listdir(log_dir):
        if log_file.endswith(".log") or log_file.endswith(".fasta"):
            family = log_file.split(".")[0]  # Extract family name
            real_time = extract_time_from_log(os.path.join(log_dir, log_file))
            time_data[family]["real"] += real_time
            time_data[family]["count"] += 1
            time_data[family]["values"].append(real_time)
    return time_data

# Function to find the closest match for a family from the log files
def get_closest_family_name(correct_family, log_families):
    closest_match = difflib.get_close_matches(correct_family, log_families, n=1)
    return closest_match[0] if closest_match else None

def main(log_dir1, log_dir2):
    # Process each directory
    time_data1 = process_log_dir(log_dir1)
    time_data2 = process_log_dir(log_dir2)

    # Calculate average time
    avg_time1 = {
        family: data["real"] / data["count"] for family, data in time_data1.items()
    }
    avg_time2 = {
        family: data["real"] / data["count"] for family, data in time_data2.items()
    }

    # List of correct family names and sequence lengths
    correct_families = [
        "tRNA",
        "5S",
        "RNaseP",
        "tmRNA",
        "Group 1",
        "telomerase",
        "16S",
    ]
    seq_lengths = [77.1, 116.2, 360.0, 367.4, 428.5, 444.9, 1419.2]

    # Ensure log family names match the correct ones using closest match
    log_families1 = list(avg_time1.keys())
    log_families2 = list(avg_time2.keys())
    avg_time_values1 = []
    avg_time_values2 = []
    time_scatter_values1 = []
    time_scatter_values2 = []

    for family in correct_families:
        closest_family1 = get_closest_family_name(family, log_families1)
        closest_family2 = get_closest_family_name(family, log_families2)

        if closest_family1:
            avg_time_values1.append(avg_time1[closest_family1])
            time_scatter_values1.append(time_data1[closest_family1]["values"])
        else:
            avg_time_values1.append(0)
            time_scatter_values1.append([])

        if closest_family2:
            avg_time_values2.append(avg_time2[closest_family2])
            time_scatter_values2.append(time_data2[closest_family2]["values"])
        else:
            avg_time_values2.append(0)
            time_scatter_values2.append([])

    # Plotting the data
    fig, ax1 = plt.subplots()

    # Set the size of the plot
    fig.set_size_inches(10, 6)

    # Plot average time for algorithm 1 with a blue line
    ax1.set_xlabel("Sequence Length")
    ax1.set_xticks(seq_lengths)
    ax1.set_xticklabels([f"{length}" for length in seq_lengths])
    ax1.set_ylabel("Avg Time (s)")
    ax1.plot(
        seq_lengths,
        avg_time_values1,
        color="tab:blue",
        label="LTF I (Vienna, Normal Outside)",
    )

    # Scatter actual time data points for algorithm 1
    for i, (length, times) in enumerate(zip(seq_lengths, time_scatter_values1)):
        ax1.scatter([length] * len(times), times, color="tab:blue", alpha=0.4)

    # Plot average time for algorithm 2 with a red line
    ax1.plot(
        seq_lengths,
        avg_time_values2,
        color="tab:red",
        label="LTF II (Vienna, Lazy Outside)",
    )

    # Scatter actual time data points for algorithm 2
    for i, (length, times) in enumerate(zip(seq_lengths, time_scatter_values2)):
        ax1.scatter([length] * len(times), times, color="tab:red", alpha=0.4)

    # Add legend
    ax1.legend()

    # Add title and layout adjustments
    plt.title("Average Time Comparison per Family")
    plt.tight_layout()
    plt.grid()

    # Show the plot
    plt.show()


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Compare average time data from two log directories.")
    parser.add_argument("log_dir1", type=str, help="Path to the first log directory.")
    parser.add_argument("log_dir2", type=str, help="Path to the second log directory.")

    args = parser.parse_args()
    main(args.log_dir1, args.log_dir2)
    