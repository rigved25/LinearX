#!/usr/bin/env bash
#
#!/usr/bin/env bash
#
# Run all LTF2 optimization combinations on a directory of multi-seq FASTA files,
# split across two tmux sessions.
#
# Usage:
#   On machine 1 (tmux session ltf2):
#     ./run_paper_optimization_tests.sh ltf2 [INPUT_DIR [OUT_BASE]]
#   On machine 2 (tmux session ltf3):
#     ./run_paper_optimization_tests.sh ltf3 [INPUT_DIR [OUT_BASE]]
#
# Examples:
#   # Original paper inputs (backwards compatible):
#   ./run_paper_optimization_tests.sh ltf2
#   ./run_paper_optimization_tests.sh ltf3
#
#   # Scale-over-k experiments:
#   ./run_paper_optimization_tests.sh ltf2 ./input_scale_over_k ./outputs/scale_over_k
#   ./run_paper_optimization_tests.sh ltf3 ./input_scale_over_k ./outputs/scale_over_k
#
#   # Scale-over-n experiments:
#   ./run_paper_optimization_tests.sh ltf2 ./input_scale_over_n ./outputs/scale_over_n
#   ./run_paper_optimization_tests.sh ltf3 ./input_scale_over_n ./outputs/scale_over_n
#
# Workload is split so slower runs (fewer optimizations) are balanced between
# the two sessions. Common flags: -sl -sp for all runs.
#

set -e

# All 16 optimization combinations (name, flags).
# Fewer flags = slower; distribution below balances slow/medium/fast across ltf2 and ltf3.
declare -a NAMES
declare -a FLAGS

# 0 opts (slowest)
NAMES+=( "ltf2_noflags" );    FLAGS+=( "" )
# 1 opt
NAMES+=( "ltf2_lazy" );        FLAGS+=( "-lz" )
NAMES+=( "ltf2_postprun" );   FLAGS+=( "-rs" )
NAMES+=( "ltf2_maxprun" );    FLAGS+=( "-mm" )
NAMES+=( "ltf2_astarvit" );   FLAGS+=( "-av" )
# 2 opts
NAMES+=( "ltf2_lazy_postprun" );      FLAGS+=( "-lz -rs" )
NAMES+=( "ltf2_lazy_maxprun" );       FLAGS+=( "-lz -mm" )
NAMES+=( "ltf2_lazy_astarvit" );      FLAGS+=( "-lz -av" )
NAMES+=( "ltf2_postprun_maxprun" );   FLAGS+=( "-rs -mm" )
NAMES+=( "ltf2_postprun_astarvit" );  FLAGS+=( "-rs -av" )
NAMES+=( "ltf2_maxprun_astarvit" );   FLAGS+=( "-mm -av" )
# 3 opts
NAMES+=( "ltf2_lazy_postprun_maxprun" );   FLAGS+=( "-lz -rs -mm" )
NAMES+=( "ltf2_lazy_postprun_astarvit" );  FLAGS+=( "-lz -rs -av" )
NAMES+=( "ltf2_lazy_maxprun_astarvit" );   FLAGS+=( "-lz -mm -av" )
NAMES+=( "ltf2_postprun_maxprun_astarvit" ); FLAGS+=( "-rs -mm -av" )
# 4 opts (fastest)
NAMES+=( "ltf2_all" );        FLAGS+=( "-lz -rs -mm -av" )

Indices for each session (balanced: mix of slow/medium/fast).
# ltf2: 8 runs — 3 slow, 3 medium, 2 fast
LTF2_IDX=( 0 1 3 5 9 10 11 13 )
# ltf3: 8 runs — 2 slow, 3 medium, 3 fast
LTF3_IDX=( 2 4 6 7 8 12 14 15 )

run_one() {
  local name=$1
  shift
  local out_dir="${OUT_BASE}/${name}"
  echo ""
  echo "========== ${name} =========="
  ./linearturbofold.py "${INPUT_DIR}" "${out_dir}" "$@" ${COMMON}
}

main() {
  local session=${1:-}
  local input_dir=${2:-./input_paper}
  local out_base=${3:-./outputs/paper}

  INPUT_DIR="${input_dir}"
  OUT_BASE="${out_base}"
  COMMON="-sl -sp"
  if [[ "$session" != "ltf2" && "$session" != "ltf3" ]]; then
    echo "Usage: $0 <ltf2|ltf3> [INPUT_DIR [OUT_BASE]]"
    echo "  ltf2     - run first half of tests (tmux session ltf2 on machine 1)"
    echo "  ltf3     - run second half of tests (tmux session ltf3 on machine 2)"
    echo "  INPUT_DIR (optional, default ./input_paper)"
    echo "  OUT_BASE  (optional, default ./outputs/paper)"
    exit 1
  fi

  if [[ ! -d "$INPUT_DIR" ]]; then
    echo "Error: input directory not found: $INPUT_DIR"
    exit 1
  fi

  local idx
  if [[ "$session" == "ltf2" ]]; then
    for idx in "${LTF2_IDX[@]}"; do
      run_one "${NAMES[$idx]}" ${FLAGS[$idx]}
    done
  else
    for idx in "${LTF3_IDX[@]}"; do
      run_one "${NAMES[$idx]}" ${FLAGS[$idx]}
    done
  fi

  echo ""
  echo "Done: session $session"
}

main "$@"
