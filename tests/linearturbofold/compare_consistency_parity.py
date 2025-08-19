#!/usr/bin/env python3
"""
Compare posterior/consistency matrices written as sparse text files across two folders
by making per-file parity (x vs y with y=x identity line) plots and summary statistics.

Example filenames:
  1_aln_0_1.bpp.txt

Example lines (sparse):
  i=0, j=0, probs=9.7953e-01

USAGE:
  python compare_bpp_parity.py --dir1 /path/to/folderA --dir2 /path/to/folderB --out ./parity_out
  # optional flags:
  #   --missing-policy zero|skip  (default: zero -> treat missing entries as 0)
  #   --aggregate               (also produce an aggregate parity plot across all files)
  #   --max-aggregate-points N  (subsample points for the aggregate plot if huge; default 300000)

Notes:
- We match files by EXACT same basename in both directories.
- We parse all (i,j,prob) triples; values not present in a file are treated per --missing-policy.
- Metrics reported: Pearson r, Spearman rho, MAE, RMSE, count.

python compare_consistency_parity.py \
  --dir1 vb_info_new \
  --dir2 vb_info_old \
  --out ./parity_out \
  --missing-policy zero \
  --aggregate

"""

import argparse
import re
from pathlib import Path
from typing import Dict, Tuple, List
import math

import numpy as np
import matplotlib.pyplot as plt
from scipy import stats

TRIPLE_RE = re.compile(r"\s*i\s*=\s*(\d+)\s*,\s*j\s*=\s*(\d+)\s*,\s*probs\s*=\s*([eE0-9\.\+\-]+)\s*$")

def parse_bpp_file(fp: Path) -> Dict[Tuple[int,int], float]:
    d: Dict[Tuple[int,int], float] = {}
    with fp.open("r", encoding="utf-8", errors="ignore") as f:
        for line in f:
            m = TRIPLE_RE.match(line)
            if not m:
                continue
            i = int(m.group(1)); j = int(m.group(2)); p = float(m.group(3))
            d[(i,j)] = p
    return d

def build_xy(a: Dict[Tuple[int,int],float],
             b: Dict[Tuple[int,int],float],
             missing_policy: str):
    if missing_policy == "zero":
        keys = set(a.keys()) | set(b.keys())
        xs = np.fromiter((a.get(k, 0.0) for k in keys), dtype=float)
        ys = np.fromiter((b.get(k, 0.0) for k in keys), dtype=float)
    elif missing_policy == "skip":
        keys = set(a.keys()) & set(b.keys())
        xs = np.fromiter((a[k] for k in keys), dtype=float)
        ys = np.fromiter((b[k] for k in keys), dtype=float)
    else:
        raise ValueError("missing_policy must be 'zero' or 'skip'")
    return xs, ys

def metrics(xs: np.ndarray, ys: np.ndarray):
    if xs.size == 0:
        return {"pearson_r": np.nan, "spearman_rho": np.nan, "mae": np.nan, "rmse": np.nan, "n": 0}
    # Safe stats
    try:
        pearson_r = stats.pearsonr(xs, ys).statistic
    except Exception:
        pearson_r = np.nan
    try:
        spearman_rho = stats.spearmanr(xs, ys).correlation
    except Exception:
        spearman_rho = np.nan
    mae = float(np.mean(np.abs(xs - ys)))
    rmse = float(math.sqrt(np.mean((xs - ys)**2)))
    return {"pearson_r": pearson_r, "spearman_rho": spearman_rho, "mae": mae, "rmse": rmse, "n": int(xs.size)}

def parity_plot(xs: np.ndarray, ys: np.ndarray, title: str, out_path: Path, annotate: dict | None = None):
    plt.figure(figsize=(7,7))
    plt.scatter(xs, ys, s=6, alpha=0.5, edgecolors='none')
    # identity line
    minv = float(min(xs.min() if xs.size else 0.0, ys.min() if ys.size else 0.0, 0.0))
    maxv = float(max(xs.max() if xs.size else 1.0, ys.max() if ys.size else 1.0, 1.0))
    plt.plot([minv, maxv], [minv, maxv])
    plt.xlabel("LTF 2")
    plt.ylabel("LTF 1")
    plt.title(title)
    if annotate:
        txt = (f"n={annotate.get('n')}\n"
               f"Pearson r={annotate.get('pearson_r'):.4f}\n"
               f"Spearman ρ={annotate.get('spearman_rho'):.4f}\n"
               f"MAE={annotate.get('mae'):.4e}\n"
               f"RMSE={annotate.get('rmse'):.4e}")
        # place annotation inside axes
        plt.gca().text(0.05, 0.95, txt, transform=plt.gca().transAxes, va='top', ha='left')
    plt.tight_layout()
    plt.savefig(out_path, dpi=150)
    plt.close()

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--dir1", required=True, type=Path, help="Folder A with .bpp.txt files")
    ap.add_argument("--dir2", required=True, type=Path, help="Folder B with .bpp.txt files")
    ap.add_argument("--out", required=True, type=Path, help="Output folder for plots and CSV")
    ap.add_argument("--missing-policy", choices=["zero","skip"], default="zero",
                    help="How to handle (i,j) present in one file but not the other")
    ap.add_argument("--aggregate", action="store_true", help="Also generate an aggregate parity plot across all files")
    ap.add_argument("--max-aggregate-points", type=int, default=300_000,
                    help="Subsample the aggregate scatter if more points than this")
    args = ap.parse_args()

    args.out.mkdir(parents=True, exist_ok=True)

    files1 = {p.name: p for p in args.dir1.glob("*.bpp.txt")}
    files2 = {p.name: p for p in args.dir2.glob("*.bpp.txt")}
    common = sorted(set(files1.keys()) & set(files2.keys()))
    if not common:
        print("No matching *.bpp.txt filenames between the two folders.", flush=True)
        return

    import csv
    summary_csv = args.out / "summary_metrics.csv"
    with summary_csv.open("w", newline="") as fcsv:
        w = csv.writer(fcsv)
        w.writerow(["filename", "n", "pearson_r", "spearman_rho", "mae", "rmse"])

        agg_xs: List[np.ndarray] = []
        agg_ys: List[np.ndarray] = []

        for name in common:
            a = parse_bpp_file(files1[name])
            b = parse_bpp_file(files2[name])
            xs, ys = build_xy(a, b, args.missing_policy)
            m = metrics(xs, ys)
            # per-file plot
            out_png = args.out / f"{name}.parity.png"
            parity_plot(xs, ys, f"{name} (parity plot)", out_png, annotate=m)
            # write metrics
            w.writerow([name, m["n"], f"{m['pearson_r']:.6f}", f"{m['spearman_rho']:.6f}",
                        f"{m['mae']:.6e}", f"{m['rmse']:.6e}"])
            # accumulate
            if args.aggregate:
                agg_xs.append(xs); agg_ys.append(ys)

        if args.aggregate and agg_xs:
            X = np.concatenate(agg_xs)
            Y = np.concatenate(agg_ys)
            # subsample if too big
            if X.size > args.max_aggregate_points:
                idx = np.random.default_rng(0).choice(X.size, size=args.max_aggregate_points, replace=False)
                X = X[idx]; Y = Y[idx]
            m = metrics(X, Y)
            parity_plot(X, Y, f"Aggregate parity plot ({len(common)} files)", args.out / "aggregate_parity.png", annotate=m)

    print(f"Wrote parity plots to: {args.out}")
    print(f"Summary CSV: {summary_csv}")

if __name__ == "__main__":
    main()