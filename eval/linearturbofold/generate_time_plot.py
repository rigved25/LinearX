import os
import re
import matplotlib.pyplot as plt
from collections import defaultdict
import difflib

# Path to your log files
# log_dir = "./output/logs_ltf2_rnastraln_lazy_ltf1_params"
log_dir = "./output/log_ltf1_rnastraln"

# Initialize dictionaries to store total time and count for averaging
time_data = defaultdict(lambda: {"real": 0, "count": 0, "values": []})


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


# Process each log file
for log_file in os.listdir(log_dir):
    if log_file.endswith(".log") or log_file.endswith(".fasta"):
        family = log_file.split(".")[0]  # Extract family name

        # Extract time data from log
        real_time = extract_time_from_log(os.path.join(log_dir, log_file))

        # Update time data for the family
        time_data[family]["real"] += real_time
        time_data[family]["count"] += 1
        time_data[family]["values"].append(real_time)

# Calculate average time
avg_time = {family: data["real"] / data["count"] for family, data in time_data.items()}

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
log_families = list(avg_time.keys())  # These are the family names from the logs
avg_time_values = []
time_scatter_values = []

for family in correct_families:
    closest_family = get_closest_family_name(family, log_families)
    if closest_family:
        avg_time_values.append(avg_time[closest_family])
        time_scatter_values.append(
            time_data[closest_family]["values"]
        )  # All real time values
    else:
        print(f"Warning: No match found for family '{family}'")

# Create multi-line labels to include avg seq length in brackets
family_labels = [
    f"{family}\n({length})" for family, length in zip(correct_families, seq_lengths)
]

# Plotting the data
fig, ax1 = plt.subplots()

# Set the size of the plot
fig.set_size_inches(10, 6)

# Plot average time with a line
ax1.set_xlabel("Family")
ax1.set_xticks(range(len(correct_families)))
ax1.set_xticklabels(family_labels)
ax1.set_ylabel("Avg Time (s)", color="tab:blue")
ax1.plot(
    range(len(correct_families)), avg_time_values, color="tab:blue", label="Avg Time"
)

# Scatter actual time data points with lower opacity
for i, times in enumerate(time_scatter_values):
    ax1.scatter([i] * len(times), times, color="tab:blue", alpha=0.4)

ax1.tick_params(axis="y", labelcolor="tab:blue")

# Add title and layout adjustments
plt.title("Average Time per Family")
plt.tight_layout()

ax1.grid(True)

# Show the plot
plt.show()
