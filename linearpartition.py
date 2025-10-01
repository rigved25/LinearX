#!/usr/bin/env python3

import argparse
import subprocess
import os
import sys


def run_linearpartition(input_seq_or_file, args):
    # Always forward 9 arguments expected by C++ binary
    cmd = [
        "./build/linearpartition",
        input_seq_or_file,
        str(args.energy_model),
        str(int(args.use_lazy_outside)),
        str(int(args.mfe_mode)),
        str(int(args.verbose)),
        str(int(args.threshknot)),
        str(int(args.bpp)),
        args.bpp_path if args.bpp_path else "none",
    ]

    print("[Executing] " + " ".join(cmd))

    try:
        subprocess.run(cmd, check=True)
    except subprocess.CalledProcessError as e:
        print(f"[Error] LinearPartition failed with code {e.returncode}")
        sys.exit(e.returncode)


def main():
    parser = argparse.ArgumentParser(
        description="Run the LinearPartition executable with command-line arguments."
    )

    parser.add_argument(
        "input",
        nargs="?",
        help="RNA sequence string or path to FASTA file. "
             "If omitted, input will be read from stdin."
    )
    parser.add_argument("--energy_model", "-em", type=int, choices=[0, 1], default=0,
                        help="Energy model (0 = Vienna, 1 = BL*), default: 0")
    parser.add_argument("--use_lazy_outside", "-lz", action="store_true", default=False,
                        help="Use lazy outside computation (default: False)")
    parser.add_argument("--mfe_mode", "-mfe", action="store_true", default=False,
                        help="Use MFE mode (default: False, i.e. partition function mode)")
    parser.add_argument("--verbose", "-v", action="store_true", default=False,
                        help="Enable verbose output (default: False)")

    # New options
    parser.add_argument("--threshknot", "-tk", action="store_true", default=False,
                        help="Compute ThreshKnot structure (default: False)")
    parser.add_argument("--bpp", "-bpp", action="store_true", default=False,
                        help="Compute base pairing probabilities (default: False)")
    parser.add_argument("--bpp-path", default="",
                        help="Directory to save BPP matrix (used only if --bpp is set)")

    args = parser.parse_args()

    # If no positional input given, read from stdin
    if not args.input:
        args.input = sys.stdin.read().strip()

    if not args.input:
        print("[Error] No input sequence or file provided")
        sys.exit(1)

    run_linearpartition(args.input, args)


if __name__ == "__main__":
    main()
