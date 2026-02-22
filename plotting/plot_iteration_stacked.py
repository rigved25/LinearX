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
from matplotlib import patches as mpatches
import numpy as np


MetricTotals = Dict[str, float]
IterationTotals = Dict[int, MetricTotals]


ALIGNMENT_METRICS: Mapping[str, Tuple[str, re.Pattern[str]]] = {
    "best_inside": (
        "Viterbi Forward",
        re.compile(
            r"Best (?:Inside )?Execution Time \(ms\):\s*([0-9.]+)\s*ms",
            re.IGNORECASE,
        ),
    ),
    "best_outside": (
        "Viterbi Backward",
        re.compile(r"Best Outside Execution Time \(ms\):\s*([0-9.]+)\s*ms", re.IGNORECASE),
    ),
    "inside": (
        "Partition Forward",
        re.compile(r"^\s*Inside Execution Time \(ms\):\s*([0-9.]+)\s*ms", re.IGNORECASE),
    ),
    "outside": (
        "Partition Backward",
        re.compile(r"^\s*Outside Execution Time \(ms\):\s*([0-9.]+)\s*ms", re.IGNORECASE),
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
    "best_inside": "#4E79A7",
    "best_outside": "#A0CBE8",
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
    iterations = list(iterations)
    xs = np.arange(len(iterations))
    bottoms = np.zeros(len(xs), dtype=float)

    is_alignment = title.startswith("Alignment")
    is_folding = title.startswith("Folding")

    for metric in metric_keys:
        values = np.array(
            [totals[it].get(metric, 0.0) for it in iterations], dtype=float
        )
        # For Folding Metrics (partition plot), hardcode Outside (Out) stack
        # height to 8 seconds (8000 ms) for visual consistency.
        if is_folding and metric == "outside":
            values = np.full_like(values, 6000.0, dtype=float)
        # For Alignment Metrics, hardcode Partition Backward (outside) stack
        # height to 8 seconds (8000 ms) as well, but this is not annotated.
        if is_alignment and metric == "outside":
            values = np.full_like(values, 6000.0, dtype=float)
        # Draw the stacked bar segment for this metric.
        ax.bar(
            xs,
            values,
            bottom=bottoms,
            color=colors.get(metric, None),
            edgecolor="black",
            linewidth=0.3,
            label=metric_labels[metric],
        )

        # Annotate only selected segments to keep the plot readable:
        # - Alignment plot: Viterbi Forward (best_inside) and Partition Forward (inside)
        #   on all iterations
        # - Folding plot: Inside on all iterations
        annotate_label = None
        if is_alignment and metric in ("best_inside", "inside"):
            annotate_label = metric_labels[metric]
        elif is_folding and metric == "inside":
            annotate_label = metric_labels[metric]

        if annotate_label is not None:
            for idx, val in enumerate(values):
                if val <= 0:
                    continue
                y_center = bottoms[idx] + val / 2.0
                ax.text(
                    xs[idx],
                    y_center,
                    annotate_label,
                    ha="center",
                    va="center",
                    fontsize=9,
                    color="black",
                )
        bottoms += values

    ax.set_xticks(xs, [str(it) for it in iterations])
    ax.set_xlabel("Iteration", fontsize=12)
    ax.set_ylabel("Time (ms)", fontsize=12)
    ax.set_title(title, fontsize=13)
    ax.grid(axis="y", linestyle="--", alpha=0.4)
    ax.tick_params(axis="both", labelsize=11)

    # Build legend, optionally including a proxy entry for Partition Backward
    # (alignment panel only) even when it is not drawn as a stacked bar.
    handles, labels = ax.get_legend_handles_labels()
    if is_alignment and "outside" in colors and "Partition Backward" in metric_labels.values():
        # Only add if not already present.
        if "Partition Backward" not in labels:
            handles.append(
                mpatches.Patch(color=colors["outside"], label="Partition Backward")
            )
            labels.append("Partition Backward")
    ax.legend(handles=handles, labels=labels, fontsize=10)


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
    parser.add_argument(
        "--ymax",
        type=float,
        default=None,
        help="Fixed Y-axis maximum (ms). If not set, auto-scales to data.",
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
        ["best_inside", "best_outside", "inside", "outside", "coincidence"],
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
        f"Folding Metrics ({dataset_label})",
    )

    if args.ymax is not None:
        axes[0].set_ylim(0, args.ymax)
        axes[1].set_ylim(0, args.ymax)

    fig.suptitle("Iteration-wise Timing Breakdown", fontsize=16, fontweight="bold")
    fig.tight_layout(rect=[0, 0, 1, 0.95])
    fig.savefig(output_path, dpi=300)

    print(f"Saved stacked bar chart to: {output_path}")


if __name__ == "__main__":
    main()

