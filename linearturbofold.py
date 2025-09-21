#!/usr/bin/env python3

import argparse
import subprocess


def main():
    parser = argparse.ArgumentParser(
        description="Run the LinearTurboFold executable with command-line arguments."
    )

    parser.add_argument("msa_file", help="Path to input MSA FASTA file")
    parser.add_argument(
        "out_dir",
        nargs="?",
        default="",
        help="Output directory",
    )

    parser.add_argument(
        "--energy_model",
        type=int,
        choices=[0, 1],
        default=0,
        help="Energy model (0 = Vienna, 1 = BL*), default: 0",
    )
    parser.add_argument(
        "--num_iterations",
        "-it",
        type=int,
        default=3,
        help="Number of TurboFold iterations (default: 3)",
    )

    parser.add_argument(
        "--use_lazy_outside",
        "-lz",
        action="store_true",
        default=False,
        help="Use lazy outside computation (default: False)",
    )
    parser.add_argument(
        "--use_prev_itr_beta",
        "-pb",
        action="store_true",
        default=False,
        help="Use beta from previous iteration (default: False)",
    )
    parser.add_argument(
        "--restrict_search",
        "-rs",
        action="store_true",
        default=False,
        help="Restrict search space (default: False)",
    )
    parser.add_argument(
        "--verbose",
        "-v",
        action="store_true",
        default=False,
        help="Enable verbose output (default: False)",
    )
    parser.add_argument(
        "--save_logs",
        "-sl",
        action="store_true",
        default=False,
        help="Save execution logs (default: False)",
    )
    parser.add_argument(
        "--save_probs",
        "-sp",
        action="store_true",
        default=False,
        help="Save BPP and coincidence probabilities (default: False)",
    )
    args = parser.parse_args()

    cmd = [
        "./build/linearturbofold",
        args.msa_file,
        args.out_dir,
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
        exit(e.returncode)


if __name__ == "__main__":
    main()
