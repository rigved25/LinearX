#!/usr/bin/env python3

import argparse
import subprocess
import os
import sys


def run_linearturbofold(msa_file, out_dir, args):
    cmd = [
        "./build/linearturbofold",
        msa_file,
        out_dir,
        str(args.energy_model),
        str(args.num_iterations),
        str(int(args.use_lazy_outside)),
        str(int(args.restrict_search)),
        str(int(args.astar_viterbi)),
        str(int(args.max_marginal)),
        str(int(args.verbose)),
        str(int(args.save_logs)),
        str(int(args.save_probs)),
        str(args.alignment_threshold),
        str(args.folding_threshold),
        str(args.max_marginal_pruning_threshold),
    ]

    print("[Executing] " + " ".join(cmd))

    # Control OpenMP parallelism. If the user has explicitly set OMP_NUM_THREADS
    # in their environment, respect it. Otherwise:
    #   - with --parallel/-parallel: use all available cores
    #   - without --parallel: force single-threaded execution.
    env = os.environ.copy()
    if "OMP_NUM_THREADS" not in env:
        if args.parallel:
            try:
                env["OMP_NUM_THREADS"] = str(os.cpu_count() or 1)
            except Exception:
                env["OMP_NUM_THREADS"] = "1"
        else:
            env["OMP_NUM_THREADS"] = "1"

    try:
        subprocess.run(cmd, check=True, env=env)
    except subprocess.CalledProcessError as e:
        print(f"[Error] LinearTurboFold failed with code {e.returncode}")
        sys.exit(e.returncode)


def main():
    parser = argparse.ArgumentParser(
        description="Run the LinearTurboFold executable with command-line arguments."
    )

    parser.add_argument("msa_path", help="Path to input MSA FASTA file or directory of FASTA files")
    parser.add_argument(
        "out_dir",
        nargs="?",
        default="",
        help="Output directory (if a directory is passed, each file gets its own subdirectory)",
    )

    parser.add_argument("--energy_model", type=int, choices=[0, 1], default=0,
                        help="Energy model (0 = Vienna, 1 = BL*), default: 0")
    parser.add_argument("--num-iterations", "-it", type=int, default=3,
                        help="Number of TurboFold iterations (default: 3)")
    parser.add_argument("--use-lazy-outside", "-lz", action="store_true", default=False,
                        help="Use lazy outside computation (default: False)")
    parser.add_argument("--restrict-search", "-rs", action="store_true", default=False,
                        help="Restrict search space (default: False)")
    parser.add_argument("--astar-viterbi", "-av", action="store_true", default=False,
                        help="Use A* Viterbi search in iterations 2+ (default: False)")
    parser.add_argument("--max-marginal", "-mm", action="store_true", default=False,
                        help="Enable max-marginal: compute and save best-mode beta, prune partition alpha (default: False)")
    parser.add_argument("--verbose", "-v", action="store_true", default=False,
                        help="Enable verbose output (default: False)")
    parser.add_argument("--save-logs", "-sl", action="store_true", default=False,
                        help="Save execution logs (default: False)")
    parser.add_argument("--save-probs", "-sp", action="store_true", default=False,
                        help="Save BPP and coincidence probabilities (default: False)")

    parser.add_argument("--parallel", "-parallel", action="store_true", default=False,
                        help="Enable OpenMP-based parallel execution inside LinearTurboFold (default: False)")

    parser.add_argument("--alignment-threshold", "-at", type=float, default=9.91152,
                        help="Alignment pruning threshold (default: 9.91152, same as DEVIATION_THRESHOLD)")
    parser.add_argument("--folding-threshold", "-ft", type=float, default=19.82304,
                        help="Folding pruning threshold (default: 19.82304, 2 * DEVIATION_THRESHOLD)")
    parser.add_argument("--max-marginal-pruning-threshold", "-mmpt", type=float, default=9.91152,
                        help="Max-marginal / BEST outside pruning threshold for grid search (default: 9.91152)")

    args = parser.parse_args()

    if os.path.isfile(args.msa_path):
        # Single FASTA file -> run once
        run_linearturbofold(args.msa_path, args.out_dir, args)

    elif os.path.isdir(args.msa_path):
        # Directory of FASTA files -> run for each
        fasta_files = [f for f in os.listdir(args.msa_path)
                       if f.lower().endswith((".fa", ".fasta"))]

        if not fasta_files:
            print(f"[Warning] No FASTA files found in {args.msa_path}")
            return

        for f in fasta_files:
            msa_file = os.path.join(args.msa_path, f)
            if len(args.out_dir):
                sub_out_dir = os.path.join(args.out_dir, os.path.splitext(f)[0])
            else:
                sub_out_dir = ""
            run_linearturbofold(msa_file, sub_out_dir, args)

    else:
        print(f"[Error] {args.msa_path} is not a valid file or directory")
        sys.exit(1)


if __name__ == "__main__":
    main()
    