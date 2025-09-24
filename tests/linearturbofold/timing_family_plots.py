#!/usr/bin/env python3
import os
import glob
import re
import numpy as np
import matplotlib.pyplot as plt

# Order for plotting
ITERATION_ORDER = ['F0', 'A1', 'F1', 'A2', 'F2', 'A3', 'F3', 'A4', 'MSA4']
ITERATION_LABELS = ['F0', 'A1', 'F1', 'A2', 'F2', 'A3', 'F3', 'A4', 'MSA4']

# Final run directories for each configuration; each contains a logs/ subdirectory
CONFIGS = {
    'noflag':   {'path': '/nfs/hpc/share/naukarkr/LiangHuang/LTF2/LinearX/tests/linearturbofold/outputfinal_rnastraln_noflag',   'label': 'noflag',                  'color': 'blue',   'linestyle': '--', 'marker': 'o'},
    'res':      {'path': '/nfs/hpc/share/naukarkr/LiangHuang/LTF2/LinearX/tests/linearturbofold/outputfinal_rnastraln_res',       'label': 'posterior pruning',       'color': 'red',    'linestyle': ':',  'marker': '^'},
    'lazy':     {'path': '/nfs/hpc/share/naukarkr/LiangHuang/LTF2/LinearX/tests/linearturbofold/outputfinal_rnastraln_lazy',      'label': 'lazy',                     'color': 'purple', 'linestyle': '-',  'marker': 'd'},
    'lazy_res': {'path': '/nfs/hpc/share/naukarkr/LiangHuang/LTF2/LinearX/tests/linearturbofold/outputfinal_rnastraln_lazy_res',  'label': 'lazy + posterior pruning', 'color': 'orange', 'linestyle': '-',  'marker': 'x'},
}

OUTPUT_DIR = '/nfs/hpc/share/naukarkr/LiangHuang/LTF2/LinearX/tests/linearturbofold/plots'


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


def aggregate_logs_for_config_by_family(config_dir):
    log_files = sorted(glob.glob(os.path.join(config_dir, 'logs', '*.fasta')))
    if not log_files:
        return {}
    family_to_bucket = {}
    for lf in log_files:
        family = os.path.basename(lf).split('.')[0]
        segs = parse_single_log(lf)
        if family not in family_to_bucket:
            family_to_bucket[family] = {seg: [] for seg in ITERATION_ORDER}
        for seg, val in segs.items():
            if seg in family_to_bucket[family]:
                family_to_bucket[family][seg].append(val)
    family_to_means = {}
    for family, bucket in family_to_bucket.items():
        family_to_means[family] = {}
        for seg, values in bucket.items():
            if values:
                family_to_means[family][seg] = int(round(float(np.mean(values))))
    return family_to_means


def collect_family_data():
    # family -> config_key -> seg -> mean_ms
    out = {}
    for key, cfg in CONFIGS.items():
        per_family = aggregate_logs_for_config_by_family(cfg['path'])
        for family, segs in per_family.items():
            out.setdefault(family, {})[key] = segs
    return out


def create_timing_charts_per_family():
    data = collect_family_data()
    if not data:
        raise RuntimeError('No timing data found. Check CONFIGS paths and logs/*.fasta presence.')
    os.makedirs(OUTPUT_DIR, exist_ok=True)

    for family, per_cfg in data.items():
        plt.figure(figsize=(16, 8))
        for key in ('noflag', 'res', 'lazy', 'lazy_res'):
            if key not in CONFIGS or key not in per_cfg:
                continue
            cfg = CONFIGS[key]
            times = [per_cfg[key].get(seg, np.nan) / 1000.0 for seg in ITERATION_ORDER]
            plt.plot(
                range(len(ITERATION_ORDER)), times,
                color=cfg['color'], linestyle=cfg['linestyle'],
                marker=cfg['marker'], label=cfg['label'],
                linewidth=2, markersize=8
            )
        plt.xlabel('Iterations', fontsize=12, fontweight='bold')
        plt.ylabel('Runtime (seconds)', fontsize=12, fontweight='bold')
        plt.title(f'LinearTurboFold: {family} Segment-wise Timing', fontsize=14, fontweight='bold', pad=20)
        plt.xticks(range(len(ITERATION_ORDER)), ITERATION_LABELS)
        plt.grid(True, alpha=0.3, linestyle='-')
        plt.legend(loc='upper left', bbox_to_anchor=(1.02, 1.0), frameon=True, fancybox=True, shadow=True, borderaxespad=0)
        plt.tight_layout(rect=[0, 0, 0.78, 1])

        out_png = os.path.join(OUTPUT_DIR, f'timing_{family}.png')
        plt.savefig(out_png, dpi=300, bbox_inches='tight')
        plt.close()
        print(f"Saved: {out_png}")


if __name__ == '__main__':
    create_timing_charts_per_family()


