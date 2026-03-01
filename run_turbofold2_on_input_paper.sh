#!/usr/bin/env bash

# run ex - ./run_turbofold2_on_input_paper.sh ./input_paper ./outputs/turbofold2

set -euo pipefail

INPUT_DIR="${1:-./input_paper}"
OUT_BASE="${2:-./outputs/turbofold2}"

ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
RNASTRUCTURE_DIR="${ROOT_DIR}/deletelater/RNAstructure"
TURBOFOLD_EXE="${RNASTRUCTURE_DIR}/exe/TurboFold"
DATAPATH_DIR="${RNASTRUCTURE_DIR}/data_tables"

if [[ ! -d "${INPUT_DIR}" ]]; then
  echo "Error: input directory not found: ${INPUT_DIR}" >&2
  exit 1
fi

mkdir -p "${OUT_BASE}"

if [[ ! -x "${TURBOFOLD_EXE}" ]]; then
  echo "TurboFold binary not found; compiling in ${RNASTRUCTURE_DIR} ..." >&2
  ( cd "${RNASTRUCTURE_DIR}" && make TurboFold -j 8 )
fi

if [[ ! -d "${DATAPATH_DIR}" ]]; then
  echo "Error: DATAPATH directory not found: ${DATAPATH_DIR}" >&2
  exit 1
fi
export DATAPATH="${DATAPATH_DIR}"

abspath() {
  local p="$1"
  python - "$p" <<'PY'
import os, sys
print(os.path.abspath(sys.argv[1]))
PY
}

count_fasta_records() {
  local fasta="$1"
  python - "$fasta" <<'PY'
import sys
path = sys.argv[1]
n = 0
with open(path, "r", encoding="utf-8", errors="replace") as f:
    for line in f:
        if line.startswith(">"):
            n += 1
print(n)
PY
}

make_group_list() {
  local prefix="$1" suffix="$2" n="$3"
  local i name
  printf "{"
  for ((i=1; i<=n; i++)); do
    name="${prefix}$(printf "%03d" "$i")${suffix}"
    printf "%s;" "$name"
  done
  printf "}"
}

shopt -s nullglob
fasta_files=( "${INPUT_DIR}"/*.fa "${INPUT_DIR}"/*.fasta "${INPUT_DIR}"/*.fna )
if (( ${#fasta_files[@]} == 0 )); then
  echo "Error: no FASTA files found in ${INPUT_DIR} (expected *.fa/*.fasta/*.fna)" >&2
  exit 1
fi

for fasta in "${fasta_files[@]}"; do
  fasta_abs="$(abspath "$fasta")"
  base="$(basename "$fasta")"
  stem="${base%.*}"
  out_dir="${OUT_BASE}/${stem}"
  out_dir_abs="$(abspath "$out_dir")"
  mkdir -p "${out_dir}"

  # Skip if this run already completed (resume-friendly)
  if [[ -f "${out_dir}/Output.aln" ]]; then
    echo "Skipping ${base}: already complete (Output.aln exists)" >&2
    continue
  fi

  nseq="$(count_fasta_records "$fasta_abs")"
  if [[ "${nseq}" -lt 2 ]]; then
    echo "Skipping ${base}: TurboFold requires >=2 sequences (found ${nseq})" >&2
    continue
  fi

  save0="$(make_group_list "iter0_" ".pfs" "${nseq}")"
  saveN="$(make_group_list "final_" ".pfs" "${nseq}")"
  outct="$(make_group_list "seq_" ".ct" "${nseq}")"

  conf="${out_dir}/turbofold.conf"
  cat > "${conf}" <<EOF
Mode = MEA
InFasta = ${fasta_abs}
OutCT = ${outct}
StartingSaveFiles = ${save0}
SaveFiles = ${saveN}
OutAln = Output.aln
AlnFormat = Fasta
ColumnNumber = 60
Gamma = 0.3
Iterations = 3
EOF

  echo ""
  echo "========== TurboFold: ${base} (${nseq} seqs) =========="
  echo "Output: ${out_dir_abs}"

  (
    cd "${out_dir}"
    "${TURBOFOLD_EXE}" "turbofold.conf" > turbofold.stdout.log 2> turbofold.stderr.log
  )
done

echo ""
echo "Done. Results in: ${OUT_BASE}"

