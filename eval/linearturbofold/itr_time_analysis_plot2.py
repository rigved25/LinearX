import argparse
import re
import matplotlib.pyplot as plt
import numpy as np

def parse_file(file_path):
    """Parse the input file and extract INSIDE and OUTSIDE runtimes."""
    mode = None
    data = {}

    with open(file_path, 'r') as file:
        for line in file:
            # Extract the MODE
            if line.startswith("MODE:"):
                mode = line.strip().split("MODE:")[1].strip()
                data[mode] = {"FOLDING_INSIDE": [], "FOLDING_OUTSIDE": [],
                              "ALIGNMENT_INSIDE": [], "ALIGNMENT_OUTSIDE": []}

            # Extract ALIGNMENT inside and outside times
            match_align_inside = re.match(r"\[ALIGNMENT\] Total inside time for iteration \d+: (\d+)ms", line)
            match_align_outside = re.match(r"\[ALIGNMENT\] Total outside time for iteration \d+: (\d+)ms", line)
            if match_align_inside:
                time = int(match_align_inside[1])
                data[mode]["ALIGNMENT_INSIDE"].append(time)
            if match_align_outside:
                time = int(match_align_outside[1])
                data[mode]["ALIGNMENT_OUTSIDE"].append(time)

            # Extract FOLDING inside and outside times
            match_fold_inside = re.match(r"\[FOLDING\] Total inside time for iteration \d+: (\d+)ms", line)
            match_fold_outside = re.match(r"\[FOLDING\] Total outside time for iteration \d+: (\d+)ms", line)
            if match_fold_inside:
                time = int(match_fold_inside[1])
                data[mode]["FOLDING_INSIDE"].append(time)
            if match_fold_outside:
                time = int(match_fold_outside[1])
                data[mode]["FOLDING_OUTSIDE"].append(time)

    return data

def plot_comparison(data):
    """Plot comparison of INSIDE and OUTSIDE times using bar plots."""
    x_labels = ['f0', 'a1', 'f1', 'a2', 'f2', 'a3', 'f3']  # Iteration labels
    width = 0.2  # Width of each bar
    x = np.arange(len(x_labels))  # X positions for the groups

    fig, ax = plt.subplots(figsize=(8, 6))  # Adjust the figure size

    # Extract data for the two modes
    modes = list(data.keys())
    mode1 = modes[0]
    mode2 = modes[1]

    print( data[mode1]["FOLDING_INSIDE"] + data[mode1]["ALIGNMENT_INSIDE"])
    print( data[mode2]["FOLDING_INSIDE"] + data[mode2]["ALIGNMENT_INSIDE"])

    def join_data_alternate(lst1, lst2):
        """Join two lists alternately."""
        return [val for pair in zip(lst1, lst2) for val in pair] + [lst1[-1]]

    # Plot bars for MODE 1
    ax.bar(x - 1.5 * width, join_data_alternate(data[mode1]["FOLDING_INSIDE"], data[mode1]["ALIGNMENT_INSIDE"]),
           width, label=f"{mode1} - Inside", color='blue', alpha=1)
    ax.bar(x - 0.5 * width, join_data_alternate(data[mode1]["FOLDING_OUTSIDE"], data[mode1]["ALIGNMENT_OUTSIDE"]),
           width, label=f"{mode1} - Outside", color='blue', alpha=0.5)

    # Plot bars for MODE 2
    ax.bar(x + 0.5 * width, join_data_alternate(data[mode2]["FOLDING_INSIDE"], data[mode2]["ALIGNMENT_INSIDE"]),
           width, label=f"{mode2} - Inside", color='green', alpha=1)
    ax.bar(x + 1.5 * width, join_data_alternate(data[mode2]["FOLDING_OUTSIDE"], data[mode2]["ALIGNMENT_OUTSIDE"]),
           width, label=f"{mode2} - Outside", color='green', alpha=0.5)

    # Format the plot
    ax.set_xlabel("Iterations", fontsize=12)
    ax.set_ylabel("Runtime (ms)", fontsize=12)
    ax.set_title("Iteration Wise Alignment and Folding Runtimes", fontsize=14)
    ax.set_xticks(x)
    ax.set_xticklabels(x_labels, fontsize=10)

    
    ax.legend(loc='upper center', bbox_to_anchor=(0.5, 1.2), ncol=2, fontsize=10)

    plt.tight_layout()
    plt.grid(axis='y', linestyle=':', alpha=0.5)  # Add horizontal grid lines
    plt.show()

def main():
    # Argument parser for input file
    parser = argparse.ArgumentParser(description="Plot comparison of INSIDE and OUTSIDE times for two modes.")
    parser.add_argument('input_file', type=str, help='Path to the input text file.')
    args = parser.parse_args()

    # Parse the input file and extract data
    data = parse_file(args.input_file)
    
    # Plot the extracted data
    plot_comparison(data)

if __name__ == "__main__":
    main()
