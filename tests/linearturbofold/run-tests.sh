#!/bin/bash

# Execute the commands one by one

# 1. LTF 2: Vienna + Lazy Outside
./run.sh ./../../eval/rnastralign/data/v1/no_aln/ output-final/rnastraln_ltf2_ltf1conf_vn_lazyout 0 1 0 0

# 2. LTF 1
./../../../../final-laf/other_programs/LinearTurboFold/run_samples.sh ./../../eval/rnastralign/data/v1/no_aln/ ./output-final/rnastraln_ltf1

# 3. LTF 2: Vienna + Lazy Outside + Outside Heuristic + Shrink Beam 
./run.sh ./../../eval/rnastralign/data/v1/no_aln/ output-final/rnastraln_ltf2_ltf1conf_vn_lazyout_outHrstc_shrnkBeam 0 1 1 1

# 4. LTF 2: Vienna + Lazy Outside + Outside Heuristic
# ./run.sh ./../../eval/rnastralign/data/v1/no_aln/ output-final/rnastraln_ltf2_ltf1conf_vn_lazyout_outHrstc 0 1 1 0

# 5. LTF 2: Vienna + Lazy Outside + Shrink Beam
# ./run.sh ./../../eval/rnastralign/data/v1/no_aln/ output-final/rnastraln_ltf2_ltf1conf_vn_lazyout_shrnkBeam 0 1 0 1