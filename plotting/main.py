import argparse
import matplotlib.pyplot as plt
import subprocess
import os

from utils import read_probs_data
from coinc_prob import plot_prob_heatmap
from bpp import plot_rna_bpp_arcs
from PIL import Image

plt.rcParams.update({'font.size': 16})  # Example: Set font size to 14

def main():
    # Set up argument parser
    parser = argparse.ArgumentParser(description='Generate the BPPs and coincident probabilities visualization.')
    parser.add_argument('prob_dir_path', type=str, help='Path to the probabilities directory')
    parser.add_argument('itr', type=str, help='Iteration number')
    parser.add_argument('seq1_len', type=int, help='Length of sequence 1')
    parser.add_argument('seq2_len', type=int, help='Length of sequence 2')
    parser.add_argument('-o', '--output', type=str, help='Output file directory', required=False)

    # Parse the arguments
    args = parser.parse_args()

    # ./vb_info/1_aln_0_1.bpp.txt ./vb_info/1_pf_seq1.bpp.txt ./vb_info/1_pf_seq2.bpp.txt

    coinc_prob_path = f"{args.prob_dir_path}/{args.itr}_aln_0_1.bpp.txt"
    seq1_bpp_path = f"{args.prob_dir_path}/{args.itr}_pf_seq1.bpp.txt"
    seq2_bpp_path = f"{args.prob_dir_path}/{args.itr}_pf_seq2.bpp.txt"


    # Read the data from the files
    coinc_prob_data = read_probs_data(coinc_prob_path)
    seq_1_bpp_data = read_probs_data(seq1_bpp_path)
    seq_2_bpp_data = read_probs_data(seq2_bpp_path)

    fig = plt.figure(figsize=(12, 10))

    # Manually adding axes with custom positions and sizes (move and resize)
    ax1 = fig.add_axes([-0.381, -0.528, 1.4, 1.214])  # [left, bottom, width, height]
    ax2 = fig.add_axes([-0.197, -0.020, 1.036, 1.4])  # Smaller size for second plot
    ax3 = fig.add_axes([0.327, 0.086, 0.7, 0.6])  # Different position and size

    ax1.set_box_aspect(1)  # Forces the axes to maintain equal aspect ratio, avoiding potential clipping
    ax2.set_box_aspect(1)  # Forces the axes to maintain equal aspect ratio, avoiding potential clipping
    ax1.set_zorder(2)
    ax2.set_zorder(1)
    ax3.set_zorder(3)

    plot_rna_bpp_arcs(ax1, seq_1_bpp_data, args.seq1_len, rotation_angle=90)
    plot_rna_bpp_arcs(ax2, seq_2_bpp_data, args.seq2_len, rotation_angle=0)
    plot_prob_heatmap(ax3, coinc_prob_data, args.seq1_len, args.seq2_len) 

    # Save and crop the PNG
    if args.output:
        if not os.path.exists(args.output):
            os.makedirs(args.output)

        dpi = 500

        # Save the plot as PNG with tight layout and transparent background
        png_output_path = os.path.join(args.output, f"itr_{args.itr}_ltf_plot.png")
        plt.tight_layout()
        plt.savefig(png_output_path, bbox_inches='tight', pad_inches=0, transparent=True, dpi=dpi)

        # Open the saved PNG image with Pillow
        with Image.open(png_output_path) as img:
            # Convert the image to RGBA (if not already in that format)
            img = img.convert("RGBA")
            
            # Get the bounding box of non-transparent pixels
            bbox = img.getbbox()

            # Crop the image using the bounding box
            cropped_img = img.crop(bbox)
            
            # Save the cropped image
            cropped_png_output_path = os.path.join(args.output, f"itr_{args.itr}_ltf_plot_cropped.png")
            cropped_img.save(cropped_png_output_path, dpi=(dpi, dpi))

        print(f"Cropped PNG saved to {cropped_png_output_path}")
    
    plt.show()

if __name__ == "__main__":
    main()