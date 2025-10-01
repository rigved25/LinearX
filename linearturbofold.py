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
        str(int(args.use_prev_itr_beta)),
        str(int(args.restrict_search)),
        str(int(args.verbose)),
        str(int(args.save_logs)),
        str(int(args.save_probs)),
    ]

    print("[Executing] " + " ".join(cmd))

    try:
        subprocess.run(cmd, check=True)
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
    parser.add_argument("--use-prev-itr-beta", "-pb", action="store_true", default=False,
                        help="Use beta from previous iteration (default: False)")
    parser.add_argument("--restrict-search", "-rs", action="store_true", default=False,
                        help="Restrict search space (default: False)")
    parser.add_argument("--verbose", "-v", action="store_true", default=False,
                        help="Enable verbose output (default: False)")
    parser.add_argument("--save-logs", "-sl", action="store_true", default=False,
                        help="Save execution logs (default: False)")
    parser.add_argument("--save-probs", "-sp", action="store_true", default=False,
                        help="Save BPP and coincidence probabilities (default: False)")

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
    