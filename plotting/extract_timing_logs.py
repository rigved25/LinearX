#!/usr/bin/env python3

import os
import re
import csv
from pathlib import Path
from collections import defaultdict
import statistics

# Directories to process
base_dirs = [
    '/nfs/hpc/share/naukarkr/LiangHuang/LTF3/LinearX/outputs/ltf2_lz_rs_AStar/',
    '/nfs/hpc/share/naukarkr/LiangHuang/LTF3/LinearX/outputs/ltf2_lz_rs/'
]

# Data structure to store metrics per iteration
alignment_data = defaultdict(list)  # itr_num -> list of metrics
partition_data = defaultdict(list)  # itr_num -> list of metrics

def parse_alignment_logs(filepath):
    """Parse alignment log files (itr_1.txt, etc.) and extract timing data"""
    metrics = {}
    
    try:
        with open(filepath, 'r') as f:
            content = f.read()
            
        # Extract Best Execution Time (ms)
        best_match = re.search(r'Best Execution Time \(ms\):\s+(\d+)\s+ms', content)
        if best_match:
            metrics['best_time'] = int(best_match.group(1))
        
        # Extract Inside Execution Time (ms)
        inside_match = re.search(r'Inside Execution Time \(ms\):\s+(\d+)\s+ms', content)
        if inside_match:
            metrics['inside_time'] = int(inside_match.group(1))
        
        # Extract Outside Execution Time (ms) - just the first occurrence
        outside_match = re.search(r'Outside Execution Time \(ms\):\s+(\d+)\s+ms', content)
        if outside_match:
            metrics['outside_time'] = int(outside_match.group(1))
        
        # Extract Coincidence Probs Execution Time (ms) - average across all pairs
        coincidence_times = re.findall(r'Coincidence Probs Execution Time:\s+([\d.]+)\s+ms', content)
        if coincidence_times:
            metrics['coincidence_time'] = statistics.mean([float(x) for x in coincidence_times])
        
        # Extract Sequence Identities - average
        seq_identities = re.findall(r'Sequence Identity:\s+([\d.]+)', content)
        if seq_identities:
            metrics['avg_seq_identity'] = statistics.mean([float(x) for x in seq_identities])
            metrics['num_pairs'] = len(seq_identities)
    
    except Exception as e:
        print(f"Error parsing {filepath}: {e}")
    
    return metrics

def parse_partition_logs(filepath):
    """Parse partition log files (itr_0.txt, etc.) and extract timing data"""
    metrics = {}
    
    try:
        with open(filepath, 'r') as f:
            content = f.read()
            
        # Extract Inside Execution Time (ms)
        inside_match = re.search(r'Inside Execution Time \(ms\):\s+(\d+)\s+ms', content)
        if inside_match:
            metrics['inside_time'] = int(inside_match.group(1))
        
        # Extract Outside Execution Time (ms)
        outside_match = re.search(r'Outside Execution Time \(ms\):\s+(\d+)\s+ms', content)
        if outside_match:
            metrics['outside_time'] = int(outside_match.group(1))
        
        # Extract BPP Execution Time (ms)
        bpp_match = re.search(r'BPP Execution Time \(ms\):\s+(\d+)\s+ms', content)
        if bpp_match:
            metrics['bpp_time'] = int(bpp_match.group(1))
    
    except Exception as e:
        print(f"Error parsing {filepath}: {e}")
    
    return metrics

def process_directories():
    """Process all directories and extract timing data"""
    
    for base_dir in base_dirs:
        if not os.path.exists(base_dir):
            print(f"Directory not found: {base_dir}")
            continue
        
        # List all dataset directories
        dataset_dirs = [d for d in os.listdir(base_dir) 
                       if os.path.isdir(os.path.join(base_dir, d))]
        
        print(f"\nProcessing {len(dataset_dirs)} datasets in {base_dir}")
        
        for dataset in dataset_dirs:
            dataset_path = os.path.join(base_dir, dataset)
            logs_path = os.path.join(dataset_path, 'logs')
            
            if not os.path.exists(logs_path):
                continue
            
            # Process partition logs
            partition_dir = os.path.join(logs_path, 'partition')
            if os.path.exists(partition_dir):
                for log_file in sorted(os.listdir(partition_dir)):
                    if log_file.startswith('itr_') and log_file.endswith('.txt'):
                        itr_num = int(log_file.split('_')[1].split('.')[0])
                        filepath = os.path.join(partition_dir, log_file)
                        metrics = parse_partition_logs(filepath)
                        if metrics:
                            partition_data[itr_num].append(metrics)
            
            # Process alignment logs
            alignment_dir = os.path.join(logs_path, 'alignment')
            if os.path.exists(alignment_dir):
                for log_file in sorted(os.listdir(alignment_dir)):
                    if log_file.startswith('itr_') and log_file.endswith('.txt'):
                        itr_num = int(log_file.split('_')[1].split('.')[0])
                        filepath = os.path.join(alignment_dir, log_file)
                        metrics = parse_alignment_logs(filepath)
                        if metrics:
                            alignment_data[itr_num].append(metrics)

def calculate_averages():
    """Calculate averages for each iteration"""
    
    output_data = []
    
    # Get all iteration numbers
    all_itrs = set(alignment_data.keys()) | set(partition_data.keys())
    
    for itr in sorted(all_itrs):
        row = {'iteration': itr}
        
        # Alignment metrics
        if itr in alignment_data:
            align_metrics = alignment_data[itr]
            if align_metrics:
                row['align_num_datasets'] = len(align_metrics)
                row['align_best_time_avg'] = statistics.mean([m.get('best_time', 0) for m in align_metrics if 'best_time' in m])
                row['align_inside_time_avg'] = statistics.mean([m.get('inside_time', 0) for m in align_metrics if 'inside_time' in m])
                row['align_outside_time_avg'] = statistics.mean([m.get('outside_time', 0) for m in align_metrics if 'outside_time' in m])
                row['align_coincidence_time_avg'] = statistics.mean([m.get('coincidence_time', 0) for m in align_metrics if 'coincidence_time' in m])
                row['align_seq_identity_avg'] = statistics.mean([m.get('avg_seq_identity', 0) for m in align_metrics if 'avg_seq_identity' in m])
        
        # Partition metrics
        if itr in partition_data:
            part_metrics = partition_data[itr]
            if part_metrics:
                row['part_num_datasets'] = len(part_metrics)
                row['part_inside_time_avg'] = statistics.mean([m.get('inside_time', 0) for m in part_metrics if 'inside_time' in m])
                row['part_outside_time_avg'] = statistics.mean([m.get('outside_time', 0) for m in part_metrics if 'outside_time' in m])
                row['part_bpp_time_avg'] = statistics.mean([m.get('bpp_time', 0) for m in part_metrics if 'bpp_time' in m])
        
        output_data.append(row)
    
    return output_data

def write_output(data):
    """Write aggregated data to CSV file"""
    
    output_file = '/nfs/hpc/share/naukarkr/LiangHuang/LTF3/LinearX/timing_analysis.csv'
    
    if not data:
        print("No data to write!")
        return
    
    # Get all column names from data
    fieldnames = ['iteration']
    for row in data:
        for key in row.keys():
            if key not in fieldnames:
                fieldnames.append(key)
    
    with open(output_file, 'w', newline='') as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(data)
    
    print(f"\nTiming analysis written to: {output_file}")
    print(f"Total iterations processed: {len(data)}")
    
    # Print summary
    print("\nSummary of extracted data:")
    for row in data[:5]:
        print(row)
    if len(data) > 5:
        print(f"... and {len(data) - 5} more iterations")

if __name__ == '__main__':
    print("Starting timing data extraction...")
    process_directories()
    print(f"\nAlignment iterations found: {sorted(alignment_data.keys())}")
    print(f"Partition iterations found: {sorted(partition_data.keys())}")
    
    data = calculate_averages()
    write_output(data)
    print("\nDone!")

