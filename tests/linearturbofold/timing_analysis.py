#!/usr/bin/env python3
import os
import glob
import re
import numpy as np
import matplotlib.pyplot as plt

# Order for plotting
ITERATION_ORDER = ['F0', 'A1', 'F1', 'A2', 'F2', 'A3', 'F3', 'A4', 'MSA4']
ITERATION_LABELS = ['F0', 'A1', 'F1', 'A2', 'F2', 'A3', 'F3', 'A4', 'MSA4']

# Directories for each configuration; each should contain a logs/ subdirectory
CONFIGS = {
    'noflag': {
        'path': '/nfs/hpc/share/naukarkr/LiangHuang/LTF2/LinearX/tests/linearturbofold/outputbasic_noflag',
        'label': 'noflags', 'color': 'blue', 'linestyle': '--', 'marker': 'o'
    },
    'pos': {
        'path': '/nfs/hpc/share/naukarkr/LiangHuang/LTF2/LinearX/tests/linearturbofold/outputbasic_lazy_pos',
        'label': 'lazy + UBH', 'color': 'green', 'linestyle': '-.', 'marker': 's'
    },
    'res': {
        'path': '/nfs/hpc/share/naukarkr/LiangHuang/LTF2/LinearX/tests/linearturbofold/outputbasic_res',
        'label': 'posterior pruning', 'color': 'red', 'linestyle': ':', 'marker': '^'
    },
    'lazy': {
        'path': '/nfs/hpc/share/naukarkr/LiangHuang/LTF2/LinearX/tests/linearturbofold/outputbasic_lazy',
        'label': 'lazy', 'color': 'purple', 'linestyle': '-', 'marker': 'd'
    },
    'lazy_res': {
        'path': '/nfs/hpc/share/naukarkr/LiangHuang/LTF2/LinearX/tests/linearturbofold/outputbasic_lazy_res',
        'label': 'lazy + posterior pruning', 'color': 'orange', 'linestyle': '-', 'marker': 'x'
    },
}

OUTPUT_PNG = '/nfs/hpc/share/naukarkr/LiangHuang/LTF2/LinearX/tests/linearturbofold/timing_comparison.png'

def parse_single_log(log_path):
    segment_times = {}
    with open(log_path, 'r', encoding='utf-8') as f:
        for line in f:
            line = line.strip()

            m = re.search(r'\[FOLDING\].*iteration\s+(\d+).*?(\d+)ms', line, re.IGNORECASE)
            if m:
                iteration = int(m.group(1))
                segment_times[f'F{iteration}'] = int(m.group(2))
                continue

            m = re.search(r'\[ALIGNMENT\].*iteration\s+(\d+).*?(\d+)ms', line, re.IGNORECASE)
            if m:
                iteration = int(m.group(1))
                segment_times[f'A{iteration}'] = int(m.group(2))
                continue

            m = re.search(r'\[Multi Sequence Alignment\].*iteration\s+(\d+).*?(\d+)ms', line, re.IGNORECASE)
            if m:
                iteration = int(m.group(1))
                segment_times[f'MSA{iteration}'] = int(m.group(2))
                continue
    return segment_times

def aggregate_logs_for_config(config_dir):
    log_files = sorted(glob.glob(os.path.join(config_dir, 'logs', '*.fasta')))
    if not log_files:
        return {}

    bucket = {seg: [] for seg in ITERATION_ORDER}
    for lf in log_files:
        segs = parse_single_log(lf)
        for seg, val in segs.items():
            if seg in bucket:
                bucket[seg].append(val)

    aggregated = {}
    for seg, values in bucket.items():
        if values:
            aggregated[seg] = int(round(float(np.mean(values))))
    return aggregated

def collect_data():
    data = {}
    for key, cfg in CONFIGS.items():
        agg = aggregate_logs_for_config(cfg['path'])
        if agg:
            data[key] = agg
    return data

def create_timing_chart():
    data = collect_data()

    if not data:
        raise RuntimeError('No timing data found. Check CONFIGS paths and logs/*.fasta presence.')

    plt.figure(figsize=(16, 8))

    for key, cfg in CONFIGS.items():
        if key not in data:
            continue
        times = [data[key].get(seg, np.nan) / 1000.0 for seg in ITERATION_ORDER]  # Convert ms to seconds
        plt.plot(
            range(len(ITERATION_ORDER)), times,
            color=cfg['color'], linestyle=cfg['linestyle'],
            marker=cfg['marker'], label=cfg['label'],
            linewidth=2, markersize=8
        )

    plt.xlabel('Iterations', fontsize=12, fontweight='bold')
    plt.ylabel('Runtime (seconds)', fontsize=12, fontweight='bold')
    plt.title('LinearTurboFold: Segment-wise Timing Comparison)',
              fontsize=14, fontweight='bold', pad=20)
    plt.xticks(range(len(ITERATION_ORDER)), ITERATION_LABELS)
    plt.grid(True, alpha=0.3, linestyle='-')
    # Place legend outside the axes on the right
    plt.legend(loc='upper left', bbox_to_anchor=(1.02, 1.0),
               frameon=True, fancybox=True, shadow=True, borderaxespad=0)
    # Leave room on the right for the legend
    plt.tight_layout(rect=[0, 0, 0.78, 1])

    plt.savefig(OUTPUT_PNG, dpi=300, bbox_inches='tight')
    print(f"Chart saved as PNG: {OUTPUT_PNG}")

    # Summary
    print("\n=== TIMING ANALYSIS SUMMARY (mean across logs) ===")
    for key, cfg in CONFIGS.items():
        if key not in data:
            continue
        total_ms = sum(v for v in data[key].values() if isinstance(v, (int, float)))
        print(f"{cfg['label']}: {total_ms/1000.0:.2f}s")

    if 'noflag' in data:
        print("\n=== SEGMENT-WISE VS noflags ===")
        for seg in ITERATION_ORDER:
            base = data['noflag'].get(seg)
            for cmp_key in ('pos', 'res', 'lazy', 'lazy_res'):
                if cmp_key in data and base and seg in data[cmp_key]:
                    t = data[cmp_key][seg]
                    speedup = float(base) / float(t) if t else np.nan
                    print(f"{seg} - {CONFIGS[cmp_key]['label']}: {base}ms → {t}ms (x{speedup:.2f})")

    return data

if __name__ == '__main__':
    create_timing_chart()
    print("\nChart saved as 'timing_comparison.png'")
