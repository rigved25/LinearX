#!/usr/bin/env python3
import os
import re
from collections import OrderedDict
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np


BASE_DIR = "/nfs/hpc/share/naukarkr/LiangHuang/LTF2/LinearX/tests/linearturbofold"

# Config name -> (subdirectory, display_name)
CONFIG_DIRS = OrderedDict([
    ("noflag", ("outputbasic_noflag", "No Flag")),
    ("res", ("outputbasic_res", "Posterior Pruning")),
    ("lazy", ("outputbasic_lazy", "Lazy")),
    ("lazy_res", ("outputbasic_lazy_res", "Lazy + Posterior Pruning")),
])

LOG_SUBDIR = "logs"

# Regex to capture the phase and milliseconds
RE_LINE = re.compile(r"\[Multi Sequence Alignment\]\s+Total Time taken for\s+(.*?):\s+(\d+)ms")

# Normalize phase labels to consistent keys
PHASE_MAP = {
    "max expected accuracy calculation": "MEA",
    "multi consistency transformation": "PCT",
    "compute tree": "ComputeTree",
    "process tree": "ProcessTree",
    "iterative refinement": "Iterative",
}

PLOT_PHASES = ["MEA", "PCT", "Tree", "Iterative"]

# Colors for each phase
PHASE_COLORS = {
    "MEA": "#2E8B57",        # Sea Green
    "PCT": "#FF6347",        # Tomato
    "Tree": "#4682B4",       # Steel Blue
    "Iterative": "#DAA520",  # Goldenrod
}


def parse_log_file(filepath: str) -> dict:
    """Parse a log file and extract phase timings."""
    phases = {}
    with open(filepath, "r", errors="ignore") as f:
        for line in f:
            m = RE_LINE.search(line)
            if not m:
                continue
            phase_raw = m.group(1).strip().lower()
            ms = int(m.group(2))
            key = PHASE_MAP.get(phase_raw)
            if key:
                phases[key] = ms
    
    # Combine compute + process tree if available
    if "ComputeTree" in phases or "ProcessTree" in phases:
        phases["Tree"] = phases.get("ComputeTree", 0) + phases.get("ProcessTree", 0)
    
    return phases


def collect_config_data(config_dirname: str) -> dict:
    """Collect timing data for a specific config."""
    logs_dir = os.path.join(BASE_DIR, config_dirname, LOG_SUBDIR)
    results = {}
    
    if not os.path.isdir(logs_dir):
        return results
    
    for name in os.listdir(logs_dir):
        # Skip RNAStrAlign files
        if "RNAStrAlign" in name:
            continue
        
        path = os.path.join(logs_dir, name)
        if not os.path.isfile(path):
            continue
        
        phases = parse_log_file(path)
        if phases:
            key = os.path.splitext(name)[0]
            results[key] = phases
    
    return results


def aggregate_config_stats(config_data: dict) -> dict:
    """Aggregate statistics for a config across all datasets."""
    if not config_data:
        return {phase: 0 for phase in PLOT_PHASES}
    
    # Collect all values for each phase
    phase_values = {phase: [] for phase in PLOT_PHASES}
    
    for dataset_name, phases in config_data.items():
        for phase in PLOT_PHASES:
            if phase in phases:
                phase_values[phase].append(phases[phase])
    
    # Calculate mean for each phase
    stats = {}
    for phase in PLOT_PHASES:
        if phase_values[phase]:
            stats[phase] = np.mean(phase_values[phase])
        else:
            stats[phase] = 0
    
    return stats


def create_phase_timing_plot(all_stats: dict, out_dir: str):
    """Create a clean bar chart showing phase-wise timing for each configuration."""
    
    configs = list(CONFIG_DIRS.keys())
    config_labels = [CONFIG_DIRS[c][1] for c in configs]
    
    # Set up the plot
    plt.style.use("default")
    fig, ax = plt.subplots(figsize=(12, 8))
    
    # Bar parameters
    bar_width = 0.6
    x_pos = np.arange(len(configs))
    
    # Create stacked bars
    bottoms = np.zeros(len(configs))
    
    for phase in PLOT_PHASES:
        values = [all_stats[config][phase] for config in configs]
        
        bars = ax.bar(x_pos, values, bar_width, bottom=bottoms, 
                     color=PHASE_COLORS[phase], label=phase,
                     edgecolor='white', linewidth=1.5)
        
        # Add value labels on bars if they're significant
        for i, (bar, value) in enumerate(zip(bars, values)):
            if value > 50:  # Only label if > 50ms
                height = bottoms[i] + value/2
                ax.text(bar.get_x() + bar.get_width()/2, height, 
                       f'{int(value)}', ha='center', va='center', 
                       fontweight='bold', fontsize=9, color='white')
        
        bottoms += values
    
    # Add total time labels on top of each bar
    for i, total in enumerate(bottoms):
        if total > 0:
            ax.text(i, total + max(20, total * 0.02), f'{int(total)}ms', 
                   ha='center', va='bottom', fontweight='bold', fontsize=11)
    
    # Customize the plot
    ax.set_xlabel('Configuration', fontsize=12, fontweight='bold')
    ax.set_ylabel('Average Time (ms)', fontsize=12, fontweight='bold')
    ax.set_title('MSA Phase-wise Timing Comparison\n(Average across all datasets)', 
                 fontsize=14, fontweight='bold', pad=20)
    
    # Set x-axis
    ax.set_xticks(x_pos)
    ax.set_xticklabels(config_labels, fontsize=11)
    
    # Add grid for better readability
    ax.grid(True, axis='y', alpha=0.3, linestyle='--')
    ax.set_axisbelow(True)
    
    # Legend
    ax.legend(bbox_to_anchor=(1.05, 1), loc='upper left', fontsize=10, framealpha=0.9)
    
    # Adjust layout
    plt.tight_layout()
    
    # Save the plot
    out_path = os.path.join(out_dir, "msa_phase_timing_comparison.png")
    plt.savefig(out_path, dpi=300, bbox_inches='tight')
    plt.close()
    
    print(f"[INFO] Created phase timing comparison: {out_path}")
    
    # Print summary statistics
    print("\n=== Phase Timing Summary (Average across datasets) ===")
    print(f"{'Config':<25} {'MEA':<8} {'PCT':<8} {'Tree':<8} {'Iterative':<10} {'Total':<8}")
    print("-" * 75)
    
    for config in configs:
        display_name = CONFIG_DIRS[config][1]
        stats = all_stats[config]
        total = sum(stats.values())
        print(f"{display_name:<25} {stats['MEA']:<8.0f} {stats['PCT']:<8.0f} "
              f"{stats['Tree']:<8.0f} {stats['Iterative']:<10.0f} {total:<8.0f}")


def ensure_plot_dir():
    """Ensure output directory exists."""
    out_dir = os.path.join(BASE_DIR, "timing_plots")
    os.makedirs(out_dir, exist_ok=True)
    return out_dir


def main():
    """Main function to generate the phase timing comparison."""
    print("Collecting MSA phase timing data...")
    
    out_dir = ensure_plot_dir()
    
    # Collect data for all configs
    all_config_data = {}
    all_stats = {}
    
    for config_name, (dirname, display_name) in CONFIG_DIRS.items():
        print(f"Processing {display_name}...")
        config_data = collect_config_data(dirname)
        all_config_data[config_name] = config_data
        all_stats[config_name] = aggregate_config_stats(config_data)
        
        if config_data:
            print(f"  Found {len(config_data)} datasets")
        else:
            print(f"  No data found in {dirname}")
    
    # Create the comparison plot
    create_phase_timing_plot(all_stats, out_dir)


if __name__ == "__main__":
    main()