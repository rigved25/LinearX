import matplotlib.pyplot as plt
import numpy as np
from matplotlib.ticker import MaxNLocator, MultipleLocator

# Function to read probabilities from file
def read_probs(filename):
    probs = {}
    with open(filename, 'r') as f:
        for line in f:
            # Extracting the indices and the probability value
            parts = line.strip().split(' = ')
            indices = tuple(map(int, parts[0][2:-1].split(', ')))
            prob = float(parts[1])
            probs[indices] = prob
    return probs

# Read data from the two files
data_p1 = read_probs('p1.txt')
data_p2 = read_probs('p2.txt')

# Extract keys that are present in both
common_keys = set(data_p1.keys()).intersection(data_p2.keys())

# Create lists for X and Y values
x_vals = [data_p1[key] for key in common_keys]
y_vals = [data_p2[key] for key in common_keys]

# Plot the data with a line representing perfect agreement (y=x) and red dots
plt.figure(figsize=(8, 8))
plt.scatter(x_vals, y_vals, s=5, color='red')  # Red dots for probabilities
plt.plot([0, 1], [0, 1], color='blue', linestyle='--', label='y = x')  # Diagonal line for perfect agreement

# Set axis labels and title
plt.xlabel('Regular Outside')
plt.ylabel('Lazy Outside')
plt.title('Probabilities Deviation Analysis (Red: probs)')

# Enable the grid with both major and minor lines
plt.grid(True, which='both', linestyle='--', linewidth=0.5)

# Set custom tick locators to ensure visibility when zooming
ax = plt.gca()

# Ensure ticks for major and minor axes remain visible
ax.xaxis.set_major_locator(MaxNLocator(10))  # Adjust maximum number of major ticks to 10
ax.yaxis.set_major_locator(MaxNLocator(10))
ax.xaxis.set_minor_locator(MultipleLocator(0.01))  # Minor ticks every 0.01
ax.yaxis.set_minor_locator(MultipleLocator(0.01))

# Customize the grid for minor ticks
plt.grid(which='minor', linestyle=':', linewidth=0.5)  # Minor grid lines with dotted style

# Ensure ticks remain visible by setting tick parameters
ax.tick_params(which='both', direction='in', top=True, right=True, length=6)

plt.legend()
plt.show()
