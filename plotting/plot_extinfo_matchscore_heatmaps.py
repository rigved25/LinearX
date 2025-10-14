#!/usr/bin/env python3
"""
Generate heatmaps for extrinsic info and match scores across iterations
"""

import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import argparse
from pathlib import Path

def load_data(filepath):
    """Load CSV file and convert to numeric types"""
    if not Path(filepath).exists():
        print(f"Warning: {filepath} not foute nd")
        return None
    
    df = pd.read_csv(filepath)
    # Convert all columns to numeric, dropping rows that can't convert (headers)
    for col in df.columns:
        df[col] = pd.to_numeric(df[col], errors='coerce')
    df = df.dropna()
    return df

def plot_extrinsic_info_heatmaps(ext_df, output_dir="./plots"):
    """Plot extrinsic info heatmaps: one PNG per sequence, all iterations as subplots"""
    Path(output_dir).mkdir(parents=True, exist_ok=True)
    
    sequences = sorted(ext_df['seq_id'].unique())
    iterations = sorted(ext_df['iteration'].unique())
    
    print(f"\n📊 Generating Extrinsic Info Heatmaps:")
    print(f"   Sequences: {len(sequences)}, Iterations: {len(iterations)}")
    
    for seq_id in sequences:
        seq_data = ext_df[ext_df['seq_id'] == seq_id]
        
        # Create subplot grid: one column per iteration
        n_iters = len(iterations)
        fig, axes = plt.subplots(1, n_iters, figsize=(6*n_iters, 5))
        if n_iters == 1:
            axes = [axes]
        
        for idx, iteration in enumerate(iterations):
            iter_data = seq_data[seq_data['iteration'] == iteration]
            
            if len(iter_data) == 0:
                axes[idx].text(0.5, 0.5, 'No data', ha='center', va='center')
                axes[idx].set_title(f'Iter {int(iteration)}')
                continue
            
            # Create matrix
            max_pos = int(max(iter_data['i'].max(), iter_data['j'].max()) + 1)
            matrix = np.zeros((max_pos, max_pos))
            
            for _, row in iter_data.iterrows():
                i, j = int(row['i']), int(row['j'])
                matrix[i, j] = row['ext_info']
            
            # Plot heatmap
            im = axes[idx].imshow(matrix, cmap='hot', aspect='auto', origin='lower')
            axes[idx].set_title(f'Iteration {int(iteration)}')
            axes[idx].set_xlabel('Position j')
            if idx == 0:
                axes[idx].set_ylabel('Position i')
            plt.colorbar(im, ax=axes[idx], label='Extrinsic Info')
        
        fig.suptitle(f'Extrinsic Info: Sequence {int(seq_id)} - e_{int(seq_id)}(i,j)^t', fontsize=14, y=1.02)
        plt.tight_layout()
        
        output_file = Path(output_dir) / f'extrinsic_info_seq{int(seq_id)}_all_iters.png'
        plt.savefig(output_file, dpi=150, bbox_inches='tight')
        plt.close()
        print(f"   ✓ Saved: {output_file}")

def plot_pairwise_match_score_heatmaps(pairwise_df, output_dir="./plots"):
    """Plot pairwise match score heatmaps: one PNG per sequence pair, all iterations as subplots"""
    Path(output_dir).mkdir(parents=True, exist_ok=True)
    
    # Get unique sequence pairs
    pairs = pairwise_df[['seq1_id', 'seq2_id']].drop_duplicates().values
    iterations = sorted(pairwise_df['iteration'].unique())
    
    print(f"\n📊 Generating Pairwise Match Score Heatmaps:")
    print(f"   Sequence pairs: {len(pairs)}, Iterations: {len(iterations)}")
    
    for seq1_id, seq2_id in pairs:
        pair_data = pairwise_df[(pairwise_df['seq1_id'] == seq1_id) & 
                                (pairwise_df['seq2_id'] == seq2_id)]
        
        # Create subplot grid: one column per iteration
        n_iters = len(iterations)
        fig, axes = plt.subplots(1, n_iters, figsize=(6*n_iters, 5))
        if n_iters == 1:
            axes = [axes]
        
        for idx, iteration in enumerate(iterations):
            iter_data = pair_data[pair_data['iteration'] == iteration]
            
            if len(iter_data) == 0:
                axes[idx].text(0.5, 0.5, 'No data', ha='center', va='center')
                axes[idx].set_title(f'Iter {int(iteration)}')
                continue
            
            # Create matrix
            max_pos1 = int(iter_data['pos1'].max() + 1)
            max_pos2 = int(iter_data['pos2'].max() + 1)
            matrix = np.zeros((max_pos1, max_pos2))
            
            for _, row in iter_data.iterrows():
                i, j = int(row['pos1']), int(row['pos2'])
                matrix[i, j] = row['match_score']
            
            # Plot heatmap
            im = axes[idx].imshow(matrix, cmap='viridis', aspect='auto', origin='lower')
            axes[idx].set_title(f'Iteration {int(iteration)}')
            axes[idx].set_xlabel(f'Seq {int(seq2_id)} pos')
            if idx == 0:
                axes[idx].set_ylabel(f'Seq {int(seq1_id)} pos')
            plt.colorbar(im, ax=axes[idx], label='Match Score')
        
        fig.suptitle(f'Match Scores: Seq {int(seq1_id)} × Seq {int(seq2_id)} - m_{int(seq1_id)},{int(seq2_id)}(i,j)^t', 
                     fontsize=14, y=1.02)
        plt.tight_layout()
        
        output_file = Path(output_dir) / f'match_scores_seq{int(seq1_id)}_seq{int(seq2_id)}_all_iters.png'
        plt.savefig(output_file, dpi=150, bbox_inches='tight')
        plt.close()
        print(f"   ✓ Saved: {output_file}")

def plot_match_scores_single_seq_heatmaps(match_df, output_dir="./plots"):
    """Plot match scores (upstream/downstream) as heatmaps for each sequence"""
    Path(output_dir).mkdir(parents=True, exist_ok=True)
    
    sequences = sorted(match_df['seq_id'].unique())
    iterations = sorted(match_df['iteration'].unique())
    
    print(f"\n📊 Generating Single-Sequence Match Score Heatmaps:")
    print(f"   Sequences: {len(sequences)}, Iterations: {len(iterations)}")
    
    for seq_id in sequences:
        seq_data = match_df[match_df['seq_id'] == seq_id]
        
        # Create subplot grid: 2 rows (upstream/downstream) × iterations columns
        n_iters = len(iterations)
        fig, axes = plt.subplots(2, n_iters, figsize=(5*n_iters, 8))
        if n_iters == 1:
            axes = axes.reshape(2, 1)
        
        for idx, iteration in enumerate(iterations):
            iter_data = seq_data[seq_data['iteration'] == iteration].sort_values('pos')
            
            if len(iter_data) == 0:
                continue
            
            max_pos = int(iter_data['pos'].max() + 1)
            
            # Upstream heatmap (as 1D array shown as 2D)
            upstream = np.zeros((20, max_pos))  # 20 rows for visual height
            for _, row in iter_data.iterrows():
                pos = int(row['pos'])
                upstream[:, pos] = row['upstream']
            
            im1 = axes[0, idx].imshow(upstream, cmap='Reds', aspect='auto', vmin=0, vmax=1)
            axes[0, idx].set_title(f'Iter {int(iteration)}')
            axes[0, idx].set_xlabel('Position')
            if idx == 0:
                axes[0, idx].set_ylabel('Upstream')
            axes[0, idx].set_yticks([])
            if idx == n_iters - 1:
                plt.colorbar(im1, ax=axes[0, idx])
            
            # Downstream heatmap
            downstream = np.zeros((20, max_pos))
            for _, row in iter_data.iterrows():
                pos = int(row['pos'])
                downstream[:, pos] = row['downstream']
            
            im2 = axes[1, idx].imshow(downstream, cmap='Blues', aspect='auto', vmin=0, vmax=1)
            axes[1, idx].set_xlabel('Position')
            if idx == 0:
                axes[1, idx].set_ylabel('Downstream')
            axes[1, idx].set_yticks([])
            if idx == n_iters - 1:
                plt.colorbar(im2, ax=axes[1, idx])
        
        fig.suptitle(f'Match Scores: Sequence {int(seq_id)} (Upstream/Downstream)', fontsize=14, y=0.995)
        plt.tight_layout()
        
        output_file = Path(output_dir) / f'match_scores_single_seq{int(seq_id)}_all_iters.png'
        plt.savefig(output_file, dpi=150, bbox_inches='tight')
        plt.close()
        print(f"   ✓ Saved: {output_file}")

def main():
    parser = argparse.ArgumentParser(description='Generate heatmaps for LinearTurboFold iteration data')
    parser.add_argument('--ext-info', default='./vb_info/extrinsic_info.csv', 
                       help='Path to extrinsic info CSV')
    parser.add_argument('--match-scores', default='./vb_info/match_scores.csv', 
                       help='Path to match scores CSV (single seq)')
    parser.add_argument('--pairwise-match', default='./vb_info/pairwise_match_scores.csv', 
                       help='Path to pairwise match scores CSV')
    parser.add_argument('--output-dir', default='./plots/heatmaps', 
                       help='Output directory for plots')
    parser.add_argument('--skip-pairwise', action='store_true',
                       help='Skip pairwise match score plots (can be slow for large data)')
    
    args = parser.parse_args()
    
    print("=" * 70)
    print("LinearTurboFold Heatmap Generator")
    print("=" * 70)
    
    # Load data
    print("\n[1] Loading data...")
    ext_df = load_data(args.ext_info)
    match_df = load_data(args.match_scores)
    pairwise_df = None if args.skip_pairwise else load_data(args.pairwise_match)
    
    if ext_df is not None:
        print(f"    Loaded {len(ext_df)} extrinsic info entries")
    if match_df is not None:
        print(f"    Loaded {len(match_df)} match score entries")
    if pairwise_df is not None:
        print(f"    Loaded {len(pairwise_df)} pairwise match score entries")
    
    # Generate plots
    print("\n[2] Generating heatmaps...")
    
    if ext_df is not None:
        plot_extrinsic_info_heatmaps(ext_df, args.output_dir)
    
    if match_df is not None:
        plot_match_scores_single_seq_heatmaps(match_df, args.output_dir)
    
    if pairwise_df is not None:
        plot_pairwise_match_score_heatmaps(pairwise_df, args.output_dir)
    
    print(f"\n[3] All plots saved to: {args.output_dir}")
    print("=" * 70)

if __name__ == '__main__':
    main()

