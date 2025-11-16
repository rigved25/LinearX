#!/usr/bin/env python3

import os
import re
import csv
from pathlib import Path
from collections import defaultdict
import statistics
import sys

def extract_family_name(dirname):
    """Extract family name from directory like '5S.RNAStrAlign.k_30_1'"""
    # Pattern: FAMILY.RNAStrAlign.k_30_NUM
    match = re.match(r'^([^.]+)\.', dirname)
    if match:
        return match.group(1)
    return None

def parse_alignment_log(filepath):
    """Parse alignment itr_1.txt and extract timing metrics"""
    metrics = {}
    
    try:
        with open(filepath, 'r') as f:
            lines = f.readlines()
        
        # Parse lines 5-8 (index 4-7, 0-based)
        content = ''.join(lines)
        
        # Extract Best Execution Time
        best_match = re.search(r'Best Execution Time \(ms\):\s+(\d+)\s+ms', content)
        if best_match:
            metrics['best_time_ms'] = int(best_match.group(1))
        
        # Extract Inside Execution Time
        inside_match = re.search(r'Inside Execution Time \(ms\):\s+(\d+)\s+ms', content)
        if inside_match:
            metrics['inside_time_ms'] = int(inside_match.group(1))
        
        # Extract Outside Execution Time
        outside_match = re.search(r'Outside Execution Time \(ms\):\s+(\d+)\s+ms', content)
        if outside_match:
            metrics['outside_time_ms'] = int(outside_match.group(1))
        
        # Extract Coincidence Probs - average all occurrences
        coincidence_times = re.findall(r'Coincidence Probs Execution Time:\s+([\d.]+)\s+ms', content)
        if coincidence_times:
            metrics['coincidence_time_ms'] = statistics.mean([float(x) for x in coincidence_times])
        
        # Extract number of sequence pairs
        seq_identities = re.findall(r'Sequence Identity:\s+([\d.]+)', content)
        if seq_identities:
            metrics['num_pairs'] = len(seq_identities)
    
    except Exception as e:
        print(f"Error parsing alignment log {filepath}: {e}", file=sys.stderr)
    
    return metrics

def parse_partition_log(filepath):
    """Parse partition itr_0.txt and extract timing metrics"""
    metrics = {}
    
    try:
        with open(filepath, 'r') as f:
            content = f.read()
        
        # Extract Inside Execution Time
        inside_match = re.search(r'Inside Execution Time \(ms\):\s+(\d+)\s+ms', content)
        if inside_match:
            metrics['inside_time_ms'] = int(inside_match.group(1))
        
        # Extract Outside Execution Time
        outside_match = re.search(r'Outside Execution Time \(ms\):\s+(\d+)\s+ms', content)
        if outside_match:
            metrics['outside_time_ms'] = int(outside_match.group(1))
        
        # Extract BPP Execution Time
        bpp_match = re.search(r'BPP Execution Time \(ms\):\s+(\d+)\s+ms', content)
        if bpp_match:
            metrics['bpp_time_ms'] = int(bpp_match.group(1))
    
    except Exception as e:
        print(f"Error parsing partition log {filepath}: {e}", file=sys.stderr)
    
    return metrics

def process_run_directory(base_dir):
    """Process all families in a run directory"""
    
    if not os.path.exists(base_dir):
        print(f"Error: Directory not found: {base_dir}")
        return None
    
    # Group datasets by family
    family_data = defaultdict(lambda: {'alignment': [], 'partition': []})
    
    # List all dataset directories
    dataset_dirs = sorted([d for d in os.listdir(base_dir) 
                          if os.path.isdir(os.path.join(base_dir, d))])
    
    print(f"Processing {len(dataset_dirs)} datasets in {base_dir}\n")
    
    for dataset in dataset_dirs:
        family = extract_family_name(dataset)
        if not family:
            continue
        
        dataset_path = os.path.join(base_dir, dataset)
        logs_path = os.path.join(dataset_path, 'logs')
        
        if not os.path.exists(logs_path):
            continue
        
        # Parse alignment log (itr_1.txt)
        align_file = os.path.join(logs_path, 'alignment', 'itr_1.txt')
        if os.path.exists(align_file):
            metrics = parse_alignment_log(align_file)
            if metrics:
                family_data[family]['alignment'].append(metrics)
        
        # Parse partition log (itr_0.txt)
        part_file = os.path.join(logs_path, 'partition', 'itr_0.txt')
        if os.path.exists(part_file):
            metrics = parse_partition_log(part_file)
            if metrics:
                family_data[family]['partition'].append(metrics)
    
    return family_data

def calculate_family_averages(family_data):
    """Calculate average metrics per family"""
    
    output_data = []
    
    for family in sorted(family_data.keys()):
        row = {'family': family}
        
        # Alignment metrics
        align_metrics = family_data[family]['alignment']
        if align_metrics:
            row['align_num_runs'] = len(align_metrics)
            row['align_best_time_avg_ms'] = statistics.mean([m.get('best_time_ms', 0) for m in align_metrics if 'best_time_ms' in m])
            row['align_inside_time_avg_ms'] = statistics.mean([m.get('inside_time_ms', 0) for m in align_metrics if 'inside_time_ms' in m])
            row['align_outside_time_avg_ms'] = statistics.mean([m.get('outside_time_ms', 0) for m in align_metrics if 'outside_time_ms' in m])
            row['align_coincidence_time_avg_ms'] = statistics.mean([m.get('coincidence_time_ms', 0) for m in align_metrics if 'coincidence_time_ms' in m])
            row['align_num_pairs_avg'] = statistics.mean([m.get('num_pairs', 0) for m in align_metrics if 'num_pairs' in m])
        
        # Partition metrics
        part_metrics = family_data[family]['partition']
        if part_metrics:
            row['part_num_runs'] = len(part_metrics)
            row['part_inside_time_avg_ms'] = statistics.mean([m.get('inside_time_ms', 0) for m in part_metrics if 'inside_time_ms' in m])
            row['part_outside_time_avg_ms'] = statistics.mean([m.get('outside_time_ms', 0) for m in part_metrics if 'outside_time_ms' in m])
            row['part_bpp_time_avg_ms'] = statistics.mean([m.get('bpp_time_ms', 0) for m in part_metrics if 'bpp_time_ms' in m])
        
        output_data.append(row)
    
    return output_data

def write_output(data, output_file):
    """Write per-family data to CSV"""
    
    if not data:
        print("No data to write!")
        return False
    
    # Get all column names
    fieldnames = ['family']
    for row in data:
        for key in row.keys():
            if key not in fieldnames:
                fieldnames.append(key)
    
    with open(output_file, 'w', newline='') as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(data)
    
    return True

def print_summary(data):
    """Print formatted summary"""
    if not data:
        return
    
    print("\n" + "="*100)
    print("FAMILY-WISE TIMING SUMMARY")
    print("="*100 + "\n")
    
    for row in data:
        print(f"Family: {row['family']}")
        print(f"  Alignment (n={row.get('align_num_runs', 0)})")
        if 'align_best_time_avg_ms' in row:
            print(f"    Best Time: {row['align_best_time_avg_ms']:.2f} ms")
        if 'align_inside_time_avg_ms' in row:
            print(f"    Inside Time: {row['align_inside_time_avg_ms']:.2f} ms")
        if 'align_outside_time_avg_ms' in row:
            print(f"    Outside Time: {row['align_outside_time_avg_ms']:.2f} ms")
        if 'align_coincidence_time_avg_ms' in row:
            print(f"    Coincidence Prob Time: {row['align_coincidence_time_avg_ms']:.2f} ms")
        if 'align_num_pairs_avg' in row:
            print(f"    Avg Num Pairs: {row['align_num_pairs_avg']:.1f}")
        
        print(f"  Partition (n={row.get('part_num_runs', 0)})")
        if 'part_inside_time_avg_ms' in row:
            print(f"    Inside Time: {row['part_inside_time_avg_ms']:.2f} ms")
        if 'part_outside_time_avg_ms' in row:
            print(f"    Outside Time: {row['part_outside_time_avg_ms']:.2f} ms")
        if 'part_bpp_time_avg_ms' in row:
            print(f"    BPP Time: {row['part_bpp_time_avg_ms']:.2f} ms")
        print()

if __name__ == '__main__':
    if len(sys.argv) < 2:
        print("Usage: python3 extract_family_timing.py <base_directory>")
        print("Example: python3 extract_family_timing.py /path/to/ltf2_lz_rs_AStar")
        sys.exit(1)
    
    base_dir = sys.argv[1]
    run_name = os.path.basename(base_dir.rstrip('/'))
    output_file = f'{run_name}_family_timing.csv'
    
    print(f"Extracting family-wise timing data from: {base_dir}\n")
    
    # Process directory
    family_data = process_run_directory(base_dir)
    
    if not family_data:
        print("No data extracted!")
        sys.exit(1)
    
    # Calculate averages
    output_data = calculate_family_averages(family_data)
    
    # Write to CSV
    if write_output(output_data, output_file):
        print(f"\n✓ Family timing data written to: {output_file}")
        print_summary(output_data)
    else:
        print("Failed to write output!")
        sys.exit(1)

