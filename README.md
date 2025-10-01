# LinearX

LinearX is a high-performance C++ framework for RNA secondary structure prediction and related algorithms.  
It implements beam search–based approximations for classical RNA dynamic programming algorithms, providing scalability to long sequences.  

Currently included programs:  
- **LinearPartition** → Fast approximation of the RNA partition function and minimum free energy (MFE) structure.  
- **LinearTurboFold** → Fast iterative alignment and folding of multiple RNA sequences (TurboFold-style).  

---

## Features

### LinearPartition
- Approximate computation of RNA partition function and ensemble free energy  
- MFE mode for minimum free energy structure prediction  
- Optional ThreshKnot postprocessing to extract consensus structures from base-pair probabilities  
- Support for Vienna and BL* energy models  
- Scales efficiently to long RNA sequences using beam search  

### LinearTurboFold
- Iterative folding and alignment across multiple RNA sequences  
- Lazy outside computation and restricted search options  
- Produces base-pair probability matrices and alignment-aware predictions  

---

## Prerequisites

- **Python**: 3.8 or newer
- **Built executables**:
  - `./build/linearpartition`
  - `./build/linearturbofold`
- **POSIX shell tools** (optional for examples):
  - `bash`, `echo`, `cat`

> If you have not built the project yet, see *Build & Install* below.

---

## Build & Install (one-time)

1) Install dependencies:
```bash
chmod +x install.sh
./install.sh
```

2) Configure and build:
```bash
chmod +x build.sh
./build.sh
```
This will produce the binaries under `./build/`.

3) (Optional) Make Python wrappers executable:
```bash
chmod +x scripts/linearpartition.py
chmod +x scripts/linearturbofold.py
```

---

## Directory Layout (relevant pieces)

```
linearx/
├─ build/                      # compiled executables
│  ├─ linearpartition
│  └─ linearturbofold
├── include/linearx/       # Header files for core modules
│   ├── energy/            # Energy model definitions
│   ├── partition/         # LinearPartition algorithm
│   ├── sequence/          # Sequence and ms utilities
│   └── turbofold/         # LinearTurboFold algorithm
├── src/                   # Source files for executables
│   ├── partition/         # LinearPartition main
│   └── turbofold/         # LinearTurboFold main
├── scripts/               # Python wrappers and utilities
├── install.sh             # Script to install dependencies
├── build.sh               # Script to configure and build
└── README.md
```

---

## Running with Python Wrappers

The Python wrappers simplify execution by validating arguments and calling the correct binary under `./build/`. You can pass a **raw sequence**, a **FASTA file path**, or **stdin** (piped input).


### 1) LinearPartition (Python)

**Script**: `./scripts/linearpartition.py`

**Purpose**:
- Partition function (ensemble free energy)
- MFE (minimum free energy) structure
- Optional ThreshKnot postprocessing
- Optional base-pair probability (BPP) computation and saving

#### Command-line synopsis
```bash
./scripts/linearpartition.py [INPUT] [OPTIONS]
```

- `INPUT`: (optional) Either a raw RNA sequence (e.g. `ACGUACGU`) or a FASTA file path.
  - If omitted, the script reads from **stdin**.
- `OPTIONS`:
  - `--energy_model {0,1}`: Energy model (default: `0`)
    - `0`: Vienna
    - `1`: BL*
  - `--use_lazy_outside`, `-lz`: Use lazy outside computation (default: off)
  - `--mfe_mode`, `-mfe`: Enable MFE mode (default: partition function mode)
  - `--verbose`, `-v`: Verbose output (default: off)
  - `--tk`: Also compute ThreshKnot structure (default: off)
  - `--bpp`: Compute base-pair probabilities (default: off)
  - `--bpp-path PATH`: Directory where BPP files will be saved (used only if `--bpp` is set)

> Internally the Python wrapper converts flags to the fixed-arity C++ interface:
> ```
> ./build/linearpartition <sequence_or_fasta> <energy_model> <use_lazy_outside> <mfe_mode> <verbose> <tk> <bpp> <bpp_path>
> ```
> When `--bpp` is not set or `--bpp-path` is not provided, it passes `bpp=0` and `bpp_path="none"`.

#### Examples

**Raw sequence (MFE, ThreshKnot, BPP saved)**
```bash
./scripts/linearpartition.py "ACGUACGUAC" --mfe --tk --bpp --bpp-path ./outputs
```

**From FASTA file (partition function, verbose)**
```bash
./scripts/linearpartition.py data/example.fa --verbose
```

**From stdin (pipe a raw sequence)**
```bash
echo ACGUACGU | ./scripts/linearpartition.py --mfe
```

**From stdin (pipe a tiny FASTA)**
```bash
echo -e ">seq1\nACGUACGU" | ./scripts/linearpartition.py --bpp --bpp-path ./bpp_out
```

#### Expected outputs

- **Partition function mode**: prints free energy of the ensemble to stdout.
- **MFE mode**: prints MFE structure (dot-bracket) to stdout.
- `--tk`: prints an additional “ThreshKnot Structure: …” line.
- `--bpp --bpp-path DIR`: writes BPP matrices to `DIR` (filenames are determined by the C++ code), and prints a confirmation line.
- With `--verbose`, the binary will also print timing and pruning stats to stderr.

#### Common issues & tips

- **“No input sequence or file provided”**: Make sure you pass `INPUT` or pipe data via stdin.
- **BPP not written**: Provide `--bpp` **and** `--bpp-path` (a directory). The wrapper will pass `"none"` to the C++ binary if `--bpp-path` is missing.
- **Performance**: For large sequences, use a release build (default in `build.sh`) and consider running without `--verbose` for cleaner logs.


### 2) LinearTurboFold (Python)

**Script**: `./scripts/linearturbofold.py`

**Purpose**:
- TurboFold-style iterative folding/alignment across multiple sequences
- Supports both single FASTA files and directories of FASTA files

#### Command-line synopsis
```bash
./scripts/linearturbofold.py <ms_path> [out_dir] [OPTIONS]
```

- `ms_path`: Path to an input FASTA file **or** a directory of FASTA files.
- `out_dir`: (optional) Output directory. If `ms_path` is a directory, each file will get its own subdirectory.
- `OPTIONS`:
  - `--energy_model {0,1}` (default: `0`): 0=Vienna, 1=BL*
  - `--num_iterations, -it N` (default: `3`): Number of TurboFold iterations
  - `--use_lazy_outside, -lz`: Use lazy outside computation
  - `--use_prev_itr_beta, -pb`: Use beta from previous iteration
  - `--restrict_search, -rs`: Restrict search space
  - `--verbose, -v`: Verbose output
  - `--save_logs, -sl`: Save execution logs
  - `--save_probs, -sp`: Save BPP and coincidence probabilities

#### Examples

**Single FASTA input**
```bash
./scripts/linearturbofold.py data/ms.fa ./outputs --num_iterations 3 --lz --verbose --save_probs
```

**Directory of FASTA files**
```bash
./scripts/linearturbofold.py data/ms_dir ./outputs --it 4 -v
```

The script will enumerate `*.fa`/`*.fasta` files and invoke the binary for each one, creating subdirectories under `./outputs/`.

#### Expected outputs

- Iteration logs (optionally saved with `--save_logs`)
- BPP matrices and coincidence probabilities (if `--save_probs` is set)
- Console logs indicating progress and execution time

---

## FAQ

**Q: Do I have to set executable bits for the Python scripts?**  
A: Not strictly—running `python3 scripts/linearpartition.py …` works without `chmod +x`. But setting the executable bit allows `./scripts/linearpartition.py …`.

**Q: Can I pass FASTA content directly as the first argument to LinearPartition?**  
A: If the argument begins with `>`, the C++ program treats it as FASTA content. Alternatively, pipe FASTA via stdin to the Python wrapper.

**Q: What energy model should I pick?**  
A: Start with `--energy_model 0` (Vienna). Use `1` (BL*) if you need those parameters specifically.

**Q: Where are BPP files written?**  
A: To the directory provided by `--bpp-path`. The naming is managed by the C++ implementation.

**Q: What is ThreshKnot?**  
A: A postprocessing step that turns base-pair probabilities into a single structure prediction. Enable it with `--tk` after running inside/outside.

---

## Troubleshooting

- Ensure `./build/linearpartition` and/or `./build/linearturbofold` exist and are executable.
- Use absolute paths for inputs if relative paths are failing.
- Run with `--verbose` to get more diagnostics (stderr) from the binaries.
- Check your Python version with `python3 --version`.
- If the wrapper prints a nonzero exit code, it propagates the binary’s error—check the stderr logs printed above the error.

---

## Versioning

- Python wrappers are thin shims and rely on the exact argument order and behavior of the compiled binaries. If you change the C++ CLI, update the corresponding wrapper to match.

---

## License

This project is released under the MIT License. See `LICENSE` for details.
