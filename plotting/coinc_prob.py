import numpy as np
import matplotlib.pyplot as plt
import matplotlib.transforms as transforms


def plot_prob_heatmap(ax, data, seq1_len, seq2_len, xlabel='Sequence 2 Indices', ylabel='Sequence 1 Indices', rotation_angle=0):
    # Create an empty matrix for the heatmap
    heatmap = np.zeros((seq1_len, seq2_len))

    # Load the data from the list of tuples
    for (idx1, idx2, prob) in data:
        heatmap[idx1, idx2] = prob

    # Apply rotation to the entire axis using a transformation
    trans = transforms.Affine2D().rotate_deg(rotation_angle) + ax.transData
    cax = ax.imshow(heatmap, cmap='Purples', interpolation='nearest', aspect='auto', origin='lower', transform=trans)

    # Move y-axis labels and ticks to the right
    ax.yaxis.tick_right()
    ax.yaxis.set_label_position('right')

    # Set x and y labels
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)

    # Add a colorbar
    cbar = plt.colorbar(cax, ax=ax, label='Probability', pad=0.12)
    cbar_pos = cbar.ax.get_position()
    cbar.ax.set_position([cbar_pos.x0, cbar_pos.y0, cbar_pos.width, cbar_pos.height])
