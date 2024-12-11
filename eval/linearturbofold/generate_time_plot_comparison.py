import os
import re
import matplotlib.pyplot as plt
from collections import defaultdict
import difflib

# Directories for log files
log_dir1 = "./output/log_ltf1_rnastraln"
log_dir2 = "./output/logs_ltf2_rnastraln_lazy_ltf1_params"

# log_dir = "./output/logs_ltf2_rnastraln_lazy_ltf1_params"
# log_dir = "./output/log_ltf1_rnastraln"

# Initialize dictionaries to store total time and count for averaging
time_data1 = defaultdict(lambda: {"real": 0, "count": 0, "values": []})
time_data2 = defaultdict(lambda: {"real": 0, "count": 0, "values": []})


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


# Process each log file for algorithm 1
for log_file in os.listdir(log_dir1):
    if log_file.endswith(".log") or log_file.endswith(".fasta"):
        family = log_file.split(".")[0]  # Extract family name
        real_time = extract_time_from_log(os.path.join(log_dir1, log_file))
        time_data1[family]["real"] += real_time
        time_data1[family]["count"] += 1
        time_data1[family]["values"].append(real_time)

# Process each log file for algorithm 2
for log_file in os.listdir(log_dir2):
    if log_file.endswith(".log") or log_file.endswith(".fasta"):
        family = log_file.split(".")[0]  # Extract family name
        real_time = extract_time_from_log(os.path.join(log_dir2, log_file))
        time_data2[family]["real"] += real_time
        time_data2[family]["count"] += 1
        time_data2[family]["values"].append(real_time)

# Calculate average time
avg_time1 = {
    family: data["real"] / data["count"] for family, data in time_data1.items()
}
avg_time2 = {
    family: data["real"] / data["count"] for family, data in time_data2.items()
}

# List of correct family names sorted by sequence length, excluding 23S
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


# Function to find the closest match for a family from the log files
def get_closest_family_name(correct_family, log_families):
    closest_match = difflib.get_close_matches(correct_family, log_families, n=1)
    return closest_match[0] if closest_match else None


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

# Create multi-line labels to include avg seq length in brackets
family_labels = [
    f"{family}\n({length})" for family, length in zip(correct_families, seq_lengths)
]

# Plotting the data
fig, ax1 = plt.subplots()

# Set the size of the plot
fig.set_size_inches(10, 6)

# Plot average time for algorithm 1 with a blue line
ax1.set_xlabel("Family")
ax1.set_xticks(range(len(correct_families)))
ax1.set_xticklabels(family_labels)
ax1.set_ylabel("Avg Time (s)")
ax1.plot(
    range(len(correct_families)),
    avg_time_values1,
    color="tab:blue",
    label="LTF I (Vienna, Normal Outside)",
)

# Scatter actual time data points for algorithm 1
for i, times in enumerate(time_scatter_values1):
    ax1.scatter([i] * len(times), times, color="tab:blue", alpha=0.4)

# Plot average time for algorithm 2 with a red line
ax1.plot(
    range(len(correct_families)),
    avg_time_values2,
    color="tab:red",
    label="LTF II (Vienna, Lazy Outside)",
)

# Scatter actual time data points for algorithm 2
for i, times in enumerate(time_scatter_values2):
    ax1.scatter([i] * len(times), times, color="tab:red", alpha=0.4)

# Add legend
ax1.legend()

# Add title and layout adjustments
plt.title("Average Time Comparison per Family")
plt.tight_layout()
plt.grid()

# Show the plot
plt.show()
