#!/usr/bin/env bash
set -euo pipefail

# Build and run external/LinearTurboFold on a directory of multi-seq FASTA files.
# Writes each input FASTA's results into its own output subdirectory.

# ./run_external_linearturbofold.sh input_tRNA the_ltf1_testrun --pf --bpp


INPUT_DIR="${1:-./input_paper}"
OUT_BASE="${2:-./outputs/external_linearturbofold}"
shift $(( $#>0 ? 1 : 0 )) || true
shift $(( $#>0 ? 1 : 0 )) || true
EXTRA_ARGS=( "$@" ) # forwarded to external/LinearTurboFold/linearturbofold (e.g. -it 3 -b1 100 -b2 100 --pf --bpp)

ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
LTF_DIR="${ROOT_DIR}/external/LinearTurboFold"
LTF_WRAPPER="${LTF_DIR}/linearturbofold"
LTF_BIN="${LTF_DIR}/bin/linearturbofold"

abspath() {
  python - "$1" <<'PY'
import os, sys
print(os.path.abspath(sys.argv[1]))
PY
}

INPUT_DIR_ABS="$(abspath "${INPUT_DIR}")"
OUT_BASE_ABS="$(abspath "${OUT_BASE}")"

if [[ ! -d "${INPUT_DIR_ABS}" ]]; then
  echo "Error: input directory not found: ${INPUT_DIR_ABS}" >&2
  exit 1
fi
mkdir -p "${OUT_BASE_ABS}"

if [[ ! -x "${LTF_WRAPPER}" ]]; then
  echo "Error: wrapper not found or not executable: ${LTF_WRAPPER}" >&2
  exit 1
fi

# Build if needed
if [[ ! -x "${LTF_BIN}" ]]; then
  echo "Building LinearTurboFold in ${LTF_DIR} ..." >&2
  ( cd "${LTF_DIR}" && make -j 8 )
fi

shopt -s nullglob
fasta_files=( "${INPUT_DIR_ABS}"/*.fa "${INPUT_DIR_ABS}"/*.fasta "${INPUT_DIR_ABS}"/*.fna )
if (( ${#fasta_files[@]} == 0 )); then
  echo "Error: no FASTA files found in ${INPUT_DIR_ABS} (expected *.fa/*.fasta/*.fna)" >&2
  exit 1
fi

# Overwrite safety: set FORCE=1 to overwrite existing non-empty output dirs.
FORCE="${FORCE:-0}"

for fasta in "${fasta_files[@]}"; do
  base="$(basename "$fasta")"
  stem="${base%.*}"
  out_dir="${OUT_BASE_ABS}/${stem}"

  if [[ -d "${out_dir}" ]] && [[ "${FORCE}" != "1" ]]; then
    # If directory exists and is non-empty, skip to avoid overwriting.
    if [[ -n "$(ls -A "${out_dir}" 2>/dev/null || true)" ]]; then
      echo "Skipping ${base}: output exists (set FORCE=1 to overwrite): ${out_dir}" >&2
      continue
    fi
  fi
  mkdir -p "${out_dir}"

  echo ""
  echo "========== external/LinearTurboFold: ${base} =========="
  echo "Output: ${out_dir}"

  # Run and capture logs.
  (
    cd "${LTF_DIR}"
    "${LTF_WRAPPER}" -i "${fasta}" -o "${out_dir}" "${EXTRA_ARGS[@]}" \
      > "${out_dir}/ltf.stdout.log" 2> "${out_dir}/ltf.stderr.log"
  )
done

echo ""
echo "Done. Results in: ${OUT_BASE_ABS}"

