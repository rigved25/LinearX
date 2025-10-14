#!/usr/bin/env python3
"""
Script to visualize coincidence probabilities from .bpp.txt files
Creates matrix-style plots showing probabilistic alignment constraints across iterations
- Groups files by type (_aln_prob_, _aln_, _mct_aln_prob_, etc.)
- Creates separate visualizations for each file type and sequence pair
- Shows present positions as shaded and missing positions as white/unshaded
- Generates plots for ALL sequence pairs automatically

Usage:
    python3 visualize_coincidence_probs.py [folder_path] [--output-dir OUTPUT_DIR] [--threshold 1e-6]
    
Examples:
    python3 visualize_coincidence_probs.py ./vb_info/
    python3 visualize_coincidence_probs.py ./vb_info/ --output-dir ./plots/
    python3 visualize_coincidence_probs.py ./vb_info/ --threshold 1e-4 -o ./results/
"""

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import os
import glob
import argparse
import re
from pathlib import Path

def parse_bpp_file(filepath):
    """Parse a .bpp.txt file and return DataFrame with alignment probabilities"""
    data = []
    
    if not os.path.exists(filepath):
        print(f"File not found: {filepath}")
        return pd.DataFrame()
    
    with open(filepath, 'r') as f:
        lines = f.readlines()
    
    for line in lines:
        line = line.strip()
        if not line:
            continue
            
        # Parse format: i=0, j=2, probs=3.167837e-04
        try:
            parts = line.split(', ')
            if len(parts) >= 3:
                i_part = parts[0].split('=')[1]
                j_part = parts[1].split('=')[1] 
                prob_part = parts[2].split('=')[1]
                
                i = int(i_part)
                j = int(j_part)
                prob = float(prob_part)
                
                data.append({
                    'i': i,
                    'j': j, 
                    'prob': prob
                })
        except (ValueError, IndexError) as e:
            # Skip malformed lines
            continue
    
    return pd.DataFrame(data)

def extract_file_info(filename):
    """Extract iteration number, file type, and sequence pair from filename"""
    # Try different patterns:
    # 1_aln_prob_0_10.bpp.txt -> iter=1, type=aln_prob, seq1=0, seq2=10
    # 0_msa_pct_alnprobs_0_1.bpp.txt -> iter=0, type=msa_pct_alnprobs, seq1=0, seq2=1
    # 2_mct_aln_prob_1_5.bpp.txt -> iter=2, type=mct_aln_prob, seq1=1, seq2=5
    
    # Pattern to extract: iteration_type_seq1_seq2.bpp.txt
    pattern = r'(\d+)_(.+?)_(\d+)_(\d+)\.bpp\.txt'
    match = re.match(pattern, filename)
    
    if match:
        iteration = int(match.group(1))
        file_type = match.group(2)
        seq1 = int(match.group(3))
        seq2 = int(match.group(4))
        return iteration, file_type, seq1, seq2
    
    return None, None, None, None

def group_files_by_type_and_sequence_pair(folder_path):
    """Group .bpp.txt files by file type and sequence pairs"""
    
    # Find all .bpp.txt files
    pattern = os.path.join(folder_path, '*.bpp.txt')
    files = glob.glob(pattern)
    
    if not files:
        print(f"No .bpp.txt files found in {folder_path}")
        return {}
    
    # Group by file type, then by sequence pairs
    grouped_data = {}
    
    for filepath in files:
        filename = os.path.basename(filepath)
        iteration, file_type, seq1, seq2 = extract_file_info(filename)
        
        if iteration is not None and file_type is not None and seq1 is not None and seq2 is not None:
            # Create consistent pair key (always smaller number first)
            pair_key = f"{min(seq1, seq2)}_{max(seq1, seq2)}"
            
            # Group by file type first
            if file_type not in grouped_data:
                grouped_data[file_type] = {}
            
            if pair_key not in grouped_data[file_type]:
                grouped_data[file_type][pair_key] = []
            
            grouped_data[file_type][pair_key].append({
                'filepath': filepath,
                'filename': filename,
                'iteration': iteration,
                'file_type': file_type,
                'seq1': seq1,
                'seq2': seq2
            })
        else:
            print(f"Could not parse filename: {filename}")
    
    # Sort files within each pair by iteration
    for file_type in grouped_data:
        for pair_key in grouped_data[file_type]:
            grouped_data[file_type][pair_key].sort(key=lambda x: x['iteration'])
    
    return grouped_data

def create_matrix_visualization_for_pair(files_info, file_type, pair_key, threshold, output_dir):
    """Create matrix-style visualization for a specific sequence pair and file type"""
    
    n_iterations = len(files_info)
    if n_iterations == 0:
        return
    
    # Determine subplot layout - arrange in 2 columns
    if n_iterations == 1:
        fig, axes = plt.subplots(1, 1, figsize=(10, 8))
        axes = [axes]
    elif n_iterations == 2:
        fig, axes = plt.subplots(1, 2, figsize=(20, 8))
    else:
        # For more than 2 iterations, use 2-column grid layout
        rows = int(np.ceil(n_iterations / 2))
        cols = 2
        fig, axes = plt.subplots(rows, cols, figsize=(20, 8*rows))
        axes = axes.flatten() if rows > 1 else axes
    
    seq1, seq2 = pair_key.split('_')
    
    # Find global max dimensions across all iterations for consistent matrix size
    global_max_i, global_max_j = 0, 0
    all_data = []
    
    for file_info in files_info:
        df = parse_bpp_file(file_info['filepath'])
        if not df.empty:
            global_max_i = max(global_max_i, df['i'].max())
            global_max_j = max(global_max_j, df['j'].max())
        all_data.append(df)
    
    # Create matrix plots
    for idx, (file_info, df) in enumerate(zip(files_info, all_data)):
        ax = axes[idx]
        
        if df.empty:
            ax.text(0.5, 0.5, 'No data', ha='center', va='center', transform=ax.transAxes)
            ax.set_title(f"Iteration {file_info['iteration']}")
            continue
        
        # Create matrix with white background (missing positions)
        matrix = np.full((global_max_i + 1, global_max_j + 1), np.nan)
        
        # Fill in the available positions
        for _, row in df.iterrows():
            if row['prob'] > threshold:
                matrix[int(row['i']), int(row['j'])] = row['prob']
        
        # Create colormap with white for NaN values (missing positions)
        cmap = plt.cm.viridis.copy()
        cmap.set_bad(color='white')
        
        # Create the plot
        im = ax.imshow(matrix, cmap=cmap, aspect='equal', origin='lower', 
                      vmin=threshold, interpolation='nearest')
        
        # Formatting
        ax.set_xlabel(f'Sequence {seq2} Position (j)')
        ax.set_ylabel(f'Sequence {seq1} Position (i)')
        ax.set_title(f"Iteration {file_info['iteration']}\n{len(df[df['prob'] > threshold])} positions > {threshold}")
        
        # Add diagonal reference line if sequences have similar lengths
        if global_max_i > 0 and global_max_j > 0:
            diagonal_length = min(global_max_i, global_max_j)
            ax.plot([0, diagonal_length], [0, diagonal_length], 'r--', alpha=0.7, linewidth=1)
        
        # Add colorbar for each subplot
        cbar = plt.colorbar(im, ax=ax, shrink=0.8)
        cbar.set_label('α(i,j)·β(i,j)/Q(x) > θ', rotation=270, labelpad=15)
    
    # Hide unused subplots
    if n_iterations < len(axes):
        for idx in range(n_iterations, len(axes)):
            axes[idx].set_visible(False)
    
    plt.suptitle(f'Probabilistic Alignment Constraint: {file_type.upper()}\n'
                 f'Sequences {seq1} ↔ {seq2}, Threshold θ = {threshold}', 
                 fontsize=14, y=0.98)
    plt.tight_layout()
    
    # Save the plot
    output_file = os.path.join(output_dir, f'{file_type}_alignment_matrix_{pair_key}.png')
    plt.savefig(output_file, dpi=300, bbox_inches='tight')
    
    print(f"Matrix visualization saved: {output_file}")
    plt.close()  # Close to save memory

def create_all_visualizations(grouped_data, threshold, output_dir):
    """Create matrix visualizations for all file types and sequence pairs"""
    
    total_plots = 0
    for file_type in grouped_data:
        for pair_key in grouped_data[file_type]:
            total_plots += 1
    
    print(f"Creating {total_plots} visualizations...")
    
    plot_count = 0
    for file_type in grouped_data:
        print(f"\nProcessing file type: {file_type}")
        
        for pair_key in grouped_data[file_type]:
            plot_count += 1
            print(f"  [{plot_count}/{total_plots}] Creating plot for sequence pair {pair_key}")
            
            files_info = grouped_data[file_type][pair_key]
            create_matrix_visualization_for_pair(files_info, file_type, pair_key, threshold, output_dir)
    
    print(f"\nCompleted all {total_plots} visualizations!")

def print_summary(grouped_data):
    """Print summary of found files"""
    print("=== FOUND FILES SUMMARY ===")
    
    for file_type, sequence_pairs in grouped_data.items():
        print(f"\nFile Type: {file_type}")
        for pair, files_info in sequence_pairs.items():
            print(f"  Sequence Pair {pair}:")
            for file_info in files_info:
                print(f"    Iteration {file_info['iteration']}: {file_info['filename']}")

def main():
    parser = argparse.ArgumentParser(description='Visualize alignment probabilities as matrix plots for all file types and sequence pairs')
    parser.add_argument('folder_path', nargs='?', default='./vb_info/',
                       help='Path to folder containing .bpp.txt files')
    parser.add_argument('--output-dir', '-o', default=None,
                       help='Output directory for saving plots (default: same as input folder)')
    parser.add_argument('--threshold', type=float, default=1e-100,
                       help='Probability threshold for visualization (default: 1e-6)')
    
    args = parser.parse_args()
    
    folder_path = os.path.abspath(args.folder_path)
    output_dir = os.path.abspath(args.output_dir) if args.output_dir else folder_path
    threshold = args.threshold
    
    print("=== ALIGNMENT PROBABILITY MATRIX VISUALIZER ===")
    print(f"Input folder: {folder_path}")
    print(f"Output folder: {output_dir}")
    print(f"Threshold: {threshold}")
    print()
    
    if not os.path.exists(folder_path):
        print(f"Error: Directory '{folder_path}' does not exist!")
        return
    
    # Group files by file type and sequence pairs
    grouped_data = group_files_by_type_and_sequence_pair(folder_path)
    
    if not grouped_data:
        print("No valid .bpp.txt files found!")
        return
    
    # Create output directory if it doesn't exist
    os.makedirs(output_dir, exist_ok=True)
    
    # Print summary
    print_summary(grouped_data)
    print()
    
    # Create all visualizations
    create_all_visualizations(grouped_data, threshold, output_dir)
    
    print("\nAnalysis complete! All matrix visualizations have been saved.")

if __name__ == "__main__":
    main()