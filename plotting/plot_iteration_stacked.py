#!/usr/bin/env python3
"""Generate iteration-wise stacked bar charts for alignment and partition logs.

Given an output experiment directory (e.g. .../outputs/<config>/<dataset>),
this script parses the per-iteration log files under `logs/alignment/` and
`logs/partition/`, sums the relevant timing metrics across all sequence pairs
or sequences, and produces stacked bar charts showing the contribution of each
metric per iteration.
"""

import argparse
import re
from pathlib import Path
from typing import Dict, Iterable, Mapping, Tuple

import matplotlib.pyplot as plt
import numpy as np


MetricTotals = Dict[str, float]
IterationTotals = Dict[int, MetricTotals]


ALIGNMENT_METRICS: Mapping[str, Tuple[str, re.Pattern[str]]] = {
    "best": (
        "Best",
        re.compile(r"Best Execution Time \(ms\):\s*([0-9.]+)\s*ms", re.IGNORECASE),
    ),
    "inside": (
        "Inside",
        re.compile(r"Inside Execution Time \(ms\):\s*([0-9.]+)\s*ms", re.IGNORECASE),
    ),
    "outside": (
        "Outside",
        re.compile(r"Outside Execution Time \(ms\):\s*([0-9.]+)\s*ms", re.IGNORECASE),
    ),
    "coincidence": (
        "Coincidence",
        re.compile(
            r"Coincidence Probs Execution Time:\s*([0-9.]+)\s*ms", re.IGNORECASE
        ),
    ),
}

PARTITION_METRICS: Mapping[str, Tuple[str, re.Pattern[str]]] = {
    "inside": (
        "Inside",
        re.compile(r"Inside Execution Time \(ms\):\s*([0-9.]+)\s*ms", re.IGNORECASE),
    ),
    "outside": (
        "Outside",
        re.compile(r"Outside Execution Time \(ms\):\s*([0-9.]+)\s*ms", re.IGNORECASE),
    ),
    "bpp": (
        "BPP",
        re.compile(r"BPP Execution Time \(ms\):\s*([0-9.]+)\s*ms", re.IGNORECASE),
    ),
}

ALIGNMENT_COLORS = {
    "best": "#4E79A7",
    "inside": "#F28E2C",
    "outside": "#E15759",
    "coincidence": "#76B7B2",
}

PARTITION_COLORS = {
    "inside": "#59A14F",
    "outside": "#EDC948",
    "bpp": "#AF7AA1",
}


def parse_iteration_file(
    file_path: Path, metrics: Mapping[str, Tuple[str, re.Pattern[str]]]
) -> MetricTotals:
    totals: MetricTotals = {key: 0.0 for key in metrics}
    with file_path.open("r", encoding="utf-8") as handle:
        for raw_line in handle:
            line = raw_line.strip()
            if not line:
                continue
            for key, (_, pattern) in metrics.items():
                match = pattern.search(line)
                if match:
                    totals[key] += float(match.group(1))
                    break
    return totals


def collect_iteration_totals(
    directory: Path, metrics: Mapping[str, Tuple[str, re.Pattern[str]]]
) -> IterationTotals:
    if not directory.exists():
        raise FileNotFoundError(f"Missing logs directory: {directory}")

    totals: IterationTotals = {}
    iteration_files = sorted(
        directory.glob("itr_*.txt"),
        key=lambda path: int(re.search(r"itr_(\d+)\.txt", path.name).group(1)),
    )

    for file_path in iteration_files:
        match = re.search(r"itr_(\d+)\.txt", file_path.name)
        if not match:
            continue
        iteration = int(match.group(1))
        totals[iteration] = parse_iteration_file(file_path, metrics)

    if not totals:
        raise RuntimeError(f"No iteration logs found in {directory}")

    return totals


def plot_stacked_bars(
    ax: plt.Axes,
    iterations: Iterable[int],
    totals: IterationTotals,
    metric_keys: Iterable[str],
    metric_labels: Mapping[str, str],
    colors: Mapping[str, str],
    title: str,
) -> None:
    xs = np.arange(len(list(iterations)))
    bottoms = np.zeros(len(xs), dtype=float)
    iterations = list(iterations)

    for metric in metric_keys:
        values = np.array([totals[it].get(metric, 0.0) for it in iterations], dtype=float)
        ax.bar(
            xs,
            values,
            bottom=bottoms,
            color=colors.get(metric, None),
            edgecolor="black",
            linewidth=0.3,
            label=metric_labels[metric],
        )
        bottoms += values

    ax.set_xticks(xs, [str(it) for it in iterations])
    ax.set_xlabel("Iteration")
    ax.set_ylabel("Time (ms)")
    ax.set_title(title)
    ax.grid(axis="y", linestyle="--", alpha=0.4)
    ax.legend()


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="Create iteration-wise stacked bar charts from LinearX logs."
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        required=True,
        help="Path to the experiment output directory (contains logs/).",
    )
    parser.add_argument(
        "--save-dir",
        type=Path,
        default=Path("/nfs/hpc/share/naukarkr/LiangHuang/LTF3/LinearX/Graphs"),
        help="Directory where the plot image will be saved.",
    )
    parser.add_argument(
        "--filename",
        type=str,
        default=None,
        help="Optional filename for the output plot. Defaults to "
        "'<dataset>_iteration_stacked.png'.",
    )
    return parser


def main() -> None:
    parser = build_parser()
    args = parser.parse_args()

    output_dir: Path = args.output_dir.resolve()
    save_dir: Path = args.save_dir.resolve()

    logs_dir = output_dir / "logs"
    alignment_dir = logs_dir / "alignment"
    partition_dir = logs_dir / "partition"

    alignment_totals = collect_iteration_totals(alignment_dir, ALIGNMENT_METRICS)
    partition_totals = collect_iteration_totals(partition_dir, PARTITION_METRICS)

    alignment_iters = sorted(alignment_totals)
    partition_iters = sorted(partition_totals)

    dataset_label = output_dir.name
    save_dir.mkdir(parents=True, exist_ok=True)

    filename = (
        args.filename
        if args.filename
        else f"{dataset_label}_iteration_stacked.png"
    )
    output_path = save_dir / filename

    plt.style.use("seaborn-v0_8-colorblind")
    fig, axes = plt.subplots(1, 2, figsize=(14, 6), sharey=True)

    plot_stacked_bars(
        axes[0],
        alignment_iters,
        alignment_totals,
        ["best", "inside", "outside", "coincidence"],
        {key: label for key, (label, _) in ALIGNMENT_METRICS.items()},
        ALIGNMENT_COLORS,
        f"Alignment Metrics ({dataset_label})",
    )

    plot_stacked_bars(
        axes[1],
        partition_iters,
        partition_totals,
        ["inside", "outside", "bpp"],
        {key: label for key, (label, _) in PARTITION_METRICS.items()},
        PARTITION_COLORS,
        f"Partition Metrics ({dataset_label})",
    )

    fig.suptitle("Iteration-wise Timing Breakdown", fontsize=16, fontweight="bold")
    fig.tight_layout(rect=[0, 0, 1, 0.95])
    fig.savefig(output_path, dpi=300)

    print(f"Saved stacked bar chart to: {output_path}")


if __name__ == "__main__":
    main()

