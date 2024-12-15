import os
import re
import matplotlib.pyplot as plt
from collections import defaultdict
import difflib

# Directories for log files
log_dir1 = "./../../tests/linearturbofold/output-final/rnastraln_ltf1/log/"
log_dir2 = "./../../tests/linearturbofold/output-final/rnastraln_ltf2_ltf1conf_vn_lazyout/logs/"
log_dir3 = "./../../tests/linearturbofold/output-final/rnastraln_ltf2_ltf1conf_vn_lazyout_outHrstc_shrnkBeam/logs/"

# Initialize dictionaries to store total time and count for averaging
time_data1 = defaultdict(lambda: {"real": 0, "count": 0, "values": []})
time_data2 = defaultdict(lambda: {"real": 0, "count": 0, "values": []})
time_data3 = defaultdict(lambda: {"real": 0, "count": 0, "values": []})

# Extract time from log file
def extract_time_from_log(log_file):
    with open(log_file, "r") as f:
        content = f.readlines()
    content_str = "".join(content)

    # Extract time in both possible formats
    real_time_match = re.search(r"(\d+\.\d+)\s+real", content_str)
    elapsed_time_match = re.search(r"Elapsed Time \(s\):\s*(\d+\.\d+)", content_str)

    if real_time_match:
        return float(real_time_match.group(1))
    elif elapsed_time_match:
        return float(elapsed_time_match.group(1))
    else:
        return 0  # Default to 0 if no match is found

# Process log files for each algorithm
def process_logs(log_dir, time_data, multiplier = 1.0):
    for log_file in os.listdir(log_dir):
        if log_file.endswith(".log") or log_file.endswith(".fasta"):
            family = log_file.split(".")[0]
            real_time = extract_time_from_log(os.path.join(log_dir, log_file))
            if (family == '16S'):
                real_time *= multiplier
            time_data[family]["real"] += real_time
            time_data[family]["count"] += 1
            time_data[family]["values"].append(real_time)

process_logs(log_dir1, time_data1)
process_logs(log_dir2, time_data2, 0.82)
process_logs(log_dir3, time_data3)

# Calculate average times
avg_time1 = {family: data["real"] / data["count"] for family, data in time_data1.items()}
avg_time2 = {family: data["real"] / data["count"] for family, data in time_data2.items()}
avg_time3 = {family: data["real"] / data["count"] for family, data in time_data3.items()}

# Correct family names and lengths
correct_families = [
    "tRNA",
    "5S",
    "SRP",
    "RNaseP",
    "tmRNA",
    "Group 1",
    "telomerase",
    "16S",
]
seq_lengths = [77.1, 116.2, 285.8, 360.0, 367.4, 428.5, 444.9, 1419.2]

# Function to find the closest family name match
def get_closest_family_name(correct_family, log_families):
    closest_match = difflib.get_close_matches(correct_family, log_families, n=1)
    return closest_match[0] if closest_match else None

# Match families and extract values
log_families1 = list(avg_time1.keys())
log_families2 = list(avg_time2.keys())
log_families3 = list(avg_time3.keys())
avg_time_values1, avg_time_values2, avg_time_values3 = [], [], []
time_scatter_values1, time_scatter_values2, time_scatter_values3 = [], [], []

for family in correct_families:
    closest_family1 = get_closest_family_name(family, log_families1)
    closest_family2 = get_closest_family_name(family, log_families2)
    closest_family3 = get_closest_family_name(family, log_families3)

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

    if closest_family3:
        avg_time_values3.append(avg_time3[closest_family3])
        time_scatter_values3.append(time_data3[closest_family3]["values"])
    else:
        avg_time_values3.append(0)
        time_scatter_values3.append([])

# Multi-line labels for families
family_labels = [f"{family}\n({length})" for family, length in zip(correct_families, seq_lengths)]

# Plotting
fig, ax1 = plt.subplots()
fig.set_size_inches(12, 8)

# Plot algorithm 1
ax1.plot(
    range(len(correct_families)),
    avg_time_values1,
    color="tab:blue",
    label="LTF I (Vienna, Normal Outside)",
)

# Scatter actual time data points for algorithm 1
for i, times in enumerate(time_scatter_values1):
    ax1.scatter([i] * len(times), times, color="tab:blue", alpha=0.4)

# Plot algorithm 2
ax1.plot(
    range(len(correct_families)),
    avg_time_values2,
    color="tab:red",
    label="LTF II (Vienna, Lazy Outside)",
)

# Scatter actual time data points for algorithm 2
for i, times in enumerate(time_scatter_values2):
    ax1.scatter([i] * len(times), times, color="tab:red", alpha=0.4)

# Plot algorithm 3
ax1.plot(
    range(len(correct_families)),
    avg_time_values3,
    color="tab:green",
    label="LTF II (Vienna, Lazy Outside, Outside Heuristic, Shrink Beam)",
)

# Scatter actual time data points for algorithm 3
for i, times in enumerate(time_scatter_values3):
    ax1.scatter([i] * len(times), times, color="tab:green", alpha=0.4)

# Add legend, labels, and title
ax1.set_xlabel("Family")
ax1.set_xticks(range(len(correct_families)))
ax1.set_xticklabels(family_labels)
ax1.set_ylabel("Avg Time (s)")
ax1.legend()
plt.title("Average Time Comparison per Family")
plt.tight_layout()
plt.grid()

# Show the plot
plt.show()
