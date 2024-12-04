import os
import re
import matplotlib.pyplot as plt
from collections import defaultdict
import difflib

# Path to your log files
log_dir = './output/log_ltf2_rnastraln'

# Initialize dictionaries to store total time, memory and count for averaging
time_data = defaultdict(lambda: {'real': 0, 'count': 0, 'values': []})
memory_data = defaultdict(lambda: {'memory': 0, 'count': 0, 'values': []})

# Extract time and maximum resident set size (RSS) from log file
def extract_data_from_log(log_file):
    with open(log_file, 'r') as f:
        content = f.readlines()
    
    # Extract time (real time)
    real_time_match = re.search(r'(\d+\.\d+)\s+real', ''.join(content))
    real_time = float(real_time_match.group(1)) if real_time_match else 0
    
    # Extract memory usage (maximum resident set size)
    rss_match = re.search(r'(\d+)\s+maximum resident set size', ''.join(content))
    max_rss = int(rss_match.group(1)) if rss_match else 0
    
    return real_time, max_rss

# Process each log file
for log_file in os.listdir(log_dir):
    if log_file.endswith('.log'):
        family = log_file.split('.')[0]  # Extract family name
        
        # Extract time and memory (RSS) data from log
        real_time, max_rss = extract_data_from_log(os.path.join(log_dir, log_file))
        
        # Update time and memory data for the family
        time_data[family]['real'] += real_time
        time_data[family]['count'] += 1
        time_data[family]['values'].append(real_time)
        
        memory_data[family]['memory'] += max_rss
        memory_data[family]['count'] += 1
        memory_data[family]['values'].append(max_rss)

# Convert memory to GB and calculate average time and memory
avg_time = {family: data['real'] / data['count'] for family, data in time_data.items()}
avg_memory = {family: (data['memory'] / data['count']) / (1024 ** 3) for family, data in memory_data.items()}  # Convert bytes to GB

# List of correct family names sorted by sequence length, excluding 23S
correct_families = ['tRNA', '5S', 'SRP', 'RNaseP', 'tmRNA', 'Group 1', 'telomerase', '16S']
seq_lengths = [77.1, 116.2, 285.8, 360.0, 367.4, 428.5, 444.9, 1419.2]

# Function to find the closest match for a family from the log files
def get_closest_family_name(correct_family, log_families):
    closest_match = difflib.get_close_matches(correct_family, log_families, n=1)
    return closest_match[0] if closest_match else None

# Ensure log family names match the correct ones using closest match
log_families = list(avg_time.keys())  # These are the family names from the logs
avg_time_values = []
avg_memory_values = []
time_scatter_values = []
memory_scatter_values = []

for family in correct_families:
    closest_family = get_closest_family_name(family, log_families)
    if closest_family:
        avg_time_values.append(avg_time[closest_family])
        avg_memory_values.append(avg_memory[closest_family])
        time_scatter_values.append(time_data[closest_family]['values'])  # All real time values
        memory_scatter_values.append([x / (1024 ** 3) for x in memory_data[closest_family]['values']])  # Convert memory to GB
    else:
        print(f"Warning: No match found for family '{family}'")

# Create multi-line labels to include avg seq length in brackets
family_labels = [f'{family}\n({length})' for family, length in zip(correct_families, seq_lengths)]

# Plotting the data
fig, ax1 = plt.subplots()

# set the size of the plot
fig.set_size_inches(10, 6)

# Plot average time with a line
ax1.set_xlabel('Family')
ax1.set_xticks(range(len(correct_families)))
ax1.set_xticklabels(family_labels)
ax1.set_ylabel('Avg Time (s)', color='tab:blue')
ax1.plot(range(len(correct_families)), avg_time_values, color='tab:blue', label='Avg Time')

# Scatter actual time data points with lower opacity
for i, times in enumerate(time_scatter_values):
    ax1.scatter([i] * len(times), times, color='tab:blue', alpha=0.4)

ax1.tick_params(axis='y', labelcolor='tab:blue')

# Create another y-axis for memory
ax2 = ax1.twinx()
ax2.set_ylabel('Avg Memory (GB)', color='tab:red')
ax2.plot(range(len(correct_families)), avg_memory_values, color='tab:red', label='Avg Memory')

# Scatter actual memory (RSS) data points with lower opacity
for i, memories in enumerate(memory_scatter_values):
    ax2.scatter([i] * len(memories), memories, color='tab:red', alpha=0.2)

ax2.tick_params(axis='y', labelcolor='tab:red')

# Add title and layout adjustments
plt.title('Average Time and Maximum Resident Set Size (Memory) per Family')
plt.tight_layout()

ax1.grid(True)

# Show the plot
plt.show()
