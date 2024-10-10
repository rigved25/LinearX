import matplotlib.pyplot as plt
from matplotlib.patches import Arc
import matplotlib.transforms as transforms
import numpy as np


def plot_rna_bpp_arcs(ax, data, seq_len, title="", xlabel="", rotation_angle=0):
    """
    This function creates a linear plot with arcs representing RNA base pairs.
    The arc intensity is given by the base pairing probability (BPP).
    
    Args:
    data (list of tuples): Each tuple contains (i, j, prob) where i and j are
                           the base pair indices, and prob is the probability.
    seq_len (int): The length of the RNA sequence.
    rotation_angle (float): The angle (in degrees) to rotate the plot content.
    
    Returns:
    fig: The figure object with the RNA base pair plot.
    """
    
    # Function to rotate points (i, j) around the origin by a given angle
    def rotate_point(x, y, angle_deg):
        angle_rad = np.deg2rad(angle_deg)
        x_new = x * np.cos(angle_rad) - y * np.sin(angle_rad)
        y_new = x * np.sin(angle_rad) + y * np.cos(angle_rad)
        return x_new, y_new
    
    # Rotate the linear RNA sequence as a black line
    x1, y1 = rotate_point(0, 0, rotation_angle)
    x2, y2 = rotate_point(seq_len - 1, 0, rotation_angle)
    ax.plot([x1, x2], [y1, y2], color='black', lw=1)  # Straight black line for the RNA sequence

    # Set limits to accommodate the rotation
    ax.set_xlim(-seq_len - 1, seq_len - 1)
    ax.set_ylim(-seq_len - 1, seq_len - 1)
    ax.get_yaxis().set_visible(False)
    
    # Draw arcs for each base pair, rotating each arc's center point
    for (i, j, prob) in data:
        if i < j:
            # Calculate the arc's height based on the distance between the base pairs
            arc_radius = (j - i) / 2
            arc_center = (i + j) / 2
            
            # Set the arc color intensity based on the probability
            color_intensity = plt.cm.Purples(prob)  # You can change the colormap
            
            # Rotate the center of the arc
            arc_center_x, arc_center_y = rotate_point(arc_center, 0, rotation_angle)
            
            # Add an arc between the base pairs with the rotated center
            arc = Arc([arc_center_x, arc_center_y], width=(j - i), height=arc_radius,
                      angle=rotation_angle, theta1=0, theta2=180, color=color_intensity, lw=2, alpha=0.9, zorder= prob)
            arc.set_clip_on(False) 
            ax.add_patch(arc)

    # Remove the spines (box) around the plot
    for spine in ax.spines.values():
        spine.set_visible(False)

    # Hide x-axis ticks
    ax.set_xticks([])

    # Label the x-axis
    ax.set_xlabel(xlabel)

    # Set the title
    ax.set_title(title)
