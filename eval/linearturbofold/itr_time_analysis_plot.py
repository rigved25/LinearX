import argparse
import re
import matplotlib.pyplot as plt
from matplotlib.ticker import MultipleLocator

def parse_file(file_path):
    """Parse the input file and extract ALIGNMENT and FOLDING runtimes."""
    mode = None
    data = {}

    with open(file_path, 'r') as file:
        for line in file:
            # Extract the MODE
            if line.startswith("MODE:"):
                mode = line.strip().split("MODE:")[1].strip()
                data[mode] = {"FOLDING": [], "ALIGNMENT": []}

            # Extract ALIGNMENT times
            match_align = re.match(r"\[ALIGNMENT\] Total Time taken for iteration (\d+): (\d+)ms", line)
            if match_align:
                time = int(match_align[2])
                data[mode]["ALIGNMENT"].append(time)

            # Extract FOLDING times
            match_fold = re.match(r"\[FOLDING\] Total Time taken for iteration (\d+): (\d+)ms", line)
            if match_fold:
                time = int(match_fold[2])
                data[mode]["FOLDING"].append(time)

    return data

def percentage_reduction(value1, value2):
    """Calculate percentage reduction from value1 to value2."""
    return ((value1 - value2) / value1) * 100 if value1 != 0 else 0

def plot_data(data):
    """Plot FOLDING and ALIGNMENT curves with arrows and percentage reduction labels."""
    x_labels = ['f0', 'a1', 'f1', 'a2', 'f2', 'a3', 'f3']  # Alternating x-axis labels
    colors = ['blue', 'green']  # Colors for the two curves

    fig, ax = plt.subplots(figsize=(8, 6))  # Adjust the figure size
    all_points = []  # To store points for arrows

    for idx, (mode, values) in enumerate(data.items()):
        combined_runtime = []
        markers = []

        # Prepare combined runtimes and markers
        for i in range(len(x_labels)):
            if i % 2 == 0:  # Even indices are FOLDING
                combined_runtime.append(values["FOLDING"][i // 2])
                markers.append('o')  # Circle marker
            else:           # Odd indices are ALIGNMENT
                combined_runtime.append(values["ALIGNMENT"][i // 2])
                markers.append('s')  # Square marker

        # Store all points for arrow drawing
        all_points.append(combined_runtime)

        # Plot the combined data with proper color, markers, and line
        for j, (x, y) in enumerate(zip(x_labels, combined_runtime)):
            if markers[j] == 'o':  # Filled circle for FOLDING
                ax.plot(x, y, 'o', color=colors[idx], markersize=8)
            else:  # Hollow square for ALIGNMENT
                ax.plot(x, y, 's', color=colors[idx], markerfacecolor='none', markersize=8)

        # Connect the points with a dotted line
        ax.plot(x_labels, combined_runtime, linestyle=':', linewidth=2, color=colors[idx], label=mode)

    # Flatten all y-values to determine dynamic tick intervals
    all_y_values = [y for sublist in all_points for y in sublist]
    y_min, y_max = min(all_y_values), max(all_y_values)

    # Determine tick intervals dynamically
    y_range = y_max - y_min    
    major_tick_interval = 6 ** (len(str(int(y_range))) - 1)  # Dynamic major tick interval
    minor_tick_interval = major_tick_interval / 2            # Minor ticks are half the major interval

    # Apply dynamic ticks
    ax.yaxis.set_major_locator(MultipleLocator(major_tick_interval))
    ax.yaxis.set_minor_locator(MultipleLocator(minor_tick_interval))

    # Add grid lines: solid for major, dotted for minor
    ax.grid(which='major', linestyle='-', linewidth=0.5, color='gray')
    ax.grid(which='minor', linestyle=':', linewidth=0.5, color='gray')

    # Draw arrows and percentage reduction labels
    for i, x in enumerate(x_labels):
        y1 = all_points[0][i]  # First mode
        y2 = all_points[1][i]  # Second mode
        reduction = percentage_reduction(y1, y2)
        if (reduction > 5):
                ax.annotate(
                "", xy=(x, y2), xytext=(x, y1),
                arrowprops=dict(arrowstyle="->", color='orange', linewidth=1)
                )
                ax.text(
                x, (y1 + y2) / 2, f"{reduction:.1f}%", color='red',
                fontsize=10, ha='center', va='center', fontweight='bold'
                )

    # Move the legend to the top
    ax.legend(loc='upper center', bbox_to_anchor=(0.5, 1.15), ncol=2, fontsize=10)

    # Plot formatting
    ax.set_xlabel("Iterations", fontsize=12)
    ax.set_ylabel("Runtime (ms)", fontsize=12)
    ax.set_title("Iteration-wise Alignment and Folding Runtimes", fontsize=14)
    plt.tight_layout()

    # Show the plot
    plt.show()

def main():
    # Argument parser for input file
    parser = argparse.ArgumentParser(description="Plot combined ALIGNMENT and FOLDING runtimes.")
    parser.add_argument('input_file', type=str, help='Path to the input text file.')
    args = parser.parse_args()

    # Parse the input file and extract data
    data = parse_file(args.input_file)
    print(data)
    
    # Plot the extracted data
    plot_data(data)

if __name__ == "__main__":
    main()
