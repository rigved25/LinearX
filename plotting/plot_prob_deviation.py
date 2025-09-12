import argparse
import math
import matplotlib.pyplot as plt
import os
import re

def read_probs(filepath):
    probs = {}
    with open(filepath) as f:
        for line in f:
            line = line.strip()
            if not line:
                continue

            # Format 1: labeled format like i=1, j=360, probs=9.9906e-01
            if "i=" in line and "j=" in line and "probs=" in line:
                match = re.match(r"i=(\d+),\s*j=(\d+),\s*probs=([\deE\.\+-]+)", line)
                if match:
                    i, j, p = int(match[1]), int(match[2]), float(match[3])
                    probs[(i, j)] = p
                continue

            # Format 2: whitespace-separated i j prob
            parts = line.split()
            if len(parts) == 3:
                try:
                    i, j, p = int(parts[0]), int(parts[1]), float(parts[2])
                    probs[(i, j)] = p
                except ValueError:
                    continue  # skip bad lines

    return probs


def main():
    parser = argparse.ArgumentParser(
        description="Compare two BPP files and plot XY scatter with stats"
    )
    parser.add_argument("file1", help="Path to first BPP file")
    parser.add_argument("file2", help="Path to second BPP file")
    parser.add_argument("--xlabel", help="X-axis label")
    parser.add_argument("--ylabel", help="Y-axis label")
    parser.add_argument("--title", help="Plot title")
    parser.add_argument(
        "--out", help="Save plot as image (e.g., 'plot.png') instead of displaying"
    )
    parser.add_argument(
        "--colored", action="store_true", help="Color scatter by deviation"
    )
    args = parser.parse_args()

    p1 = read_probs(args.file1)
    p2 = read_probs(args.file2)
    print(
        f"[Info] Read {len(p1)} entries from {args.file1}"
        f" and {len(p2)} entries from {args.file2}"
    )

    all_keys = set(p1.keys()).union(p2.keys())
    x_vals, y_vals, diffs, squared_diffs = [], [], [], []

    for key in all_keys:
        v1 = p1.get(key, 0.0)
        v2 = p2.get(key, 0.0)
        x_vals.append(v1)
        y_vals.append(v2)
        diff = abs(v1 - v2)
        diffs.append(diff)
        squared_diffs.append((v1 - v2) ** 2)

    avg_diff = sum(diffs) / len(diffs)
    rmsd = math.sqrt(sum(squared_diffs) / len(squared_diffs))

    # === Labels ===
    x_file = os.path.basename(args.file1).rsplit(".", 1)[0]
    y_file = os.path.basename(args.file2).rsplit(".", 1)[0]
    xlabel = args.xlabel or f"Probs in file {x_file}"
    ylabel = args.ylabel or f"Probs in file {y_file}"
    title = args.title or f"{xlabel} vs {ylabel}"

    # === Plot ===
    fig, ax = plt.subplots(figsize=(8, 8))
    if args.colored:
        sc = ax.scatter(
            x_vals, y_vals, c=diffs, cmap="coolwarm", s=12, alpha=0.6, edgecolors="none"
        )
        plt.colorbar(sc, label="|Δprobability|")
    else:
        ax.scatter(x_vals, y_vals, s=10, alpha=0.4)

    ax.plot([0, 1], [0, 1], "r--", linewidth=1, label="Perfect Match")

    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    ax.set_xticks([i * 0.1 for i in range(11)])
    ax.set_yticks([i * 0.1 for i in range(11)])
    ax.grid(True, linestyle="--", alpha=0.3)

    ax.set_xlabel(xlabel, fontsize=12, fontweight="bold")
    ax.set_ylabel(ylabel, fontsize=12, fontweight="bold")
    ax.set_title(title, fontsize=14, fontweight="bold")
    ax.legend(loc="lower center")

    # === Stats box 1: unique entry counts ===
    only_in_1 = len(set(p1.keys()) - set(p2.keys()))
    only_in_2 = len(set(p2.keys()) - set(p1.keys()))
    total = len(all_keys)

    stats_text = (
        f"Total entries: {total}\n"
        f"Unique to Y: {only_in_1}\n"
        f"Unique to X: {only_in_2}"
    )
    ax.text(
        0.02,
        0.98,
        stats_text,
        transform=ax.transAxes,
        verticalalignment="top",
        horizontalalignment="left",
        fontsize=10,
        bbox=dict(
            facecolor="white", edgecolor="black", boxstyle="round,pad=0.4", alpha=0.85
        ),
    )

    # === Stats box 2: deviation numbers ===
    dev_text = f"[Avg Δ] {avg_diff:.5f}\n" f"[RMSD]  {rmsd:.5f}"
    ax.text(
        0.98,
        0.02,
        dev_text,
        transform=ax.transAxes,
        verticalalignment="bottom",
        horizontalalignment="right",
        fontsize=10,
        bbox=dict(
            facecolor="white", edgecolor="black", boxstyle="round,pad=0.4", alpha=0.85
        ),
    )

    plt.tight_layout()

    # === Console Output ===
    print("\n================ Probability Deviation Report ================\n")
    print(f"Average Probability Deviation : {avg_diff:.6f}")
    print(f"Root Mean Square Deviation    : {rmsd:.6f}")
    print("\nInterpretation Guide:")
    print("-------------------------------------------------------------")
    print("Value     | Avg Dev / RMSD     | Interpretation")
    print("----------|---------------------|----------------------------------------")
    print("~0.00     | Perfect match        | Identical BPPs ✅")
    print("<  0.02   | Very similar         | Only slight variations")
    print("0.02-0.05 | Moderate             | Noticeable differences")
    print(">  0.10   | Large                | Very different pairing probabilities ❗")
    print("-------------------------------------------------------------\n")

    if args.out:
        plt.savefig(args.out, bbox_inches="tight", dpi=300)
        print(f"[Saved] Plot saved to {args.out}")
    else:
        plt.show()


if __name__ == "__main__":
    main()
