#include "energy_model.hpp"
#include <iostream>

void EnergyModel::init_tetra_hex_tri(std::vector<std::string> &seq_MSA_no_gap,
                                     std::vector<std::vector<int> > &if_tetraloops_MSA,
                                     std::vector<std::vector<int> > &if_hexaloops_MSA,
                                     std::vector<std::vector<int> > &if_triloops_MSA) {
  for (int s = 0; s < if_tetraloops_MSA.size(); s++) {
    std::string seq = seq_MSA_no_gap[s];

    int seq_length = seq.size();

    // TetraLoops
    if_tetraloops_MSA[s].resize(seq_length - 5 < 0 ? 0 : seq_length - 5, -1);
    for (int i = 0; i < seq_length - 5; ++i) {
      if (!(seq[i] == 'C' && seq[i + 5] == 'G'))
        continue;
      const char *ts;
      if ((ts = strstr(Tetraloops, seq.substr(i, 6).c_str())))
        if_tetraloops_MSA[s][i] = (ts - Tetraloops) / 7;
    }

    // Triloops
    if_triloops_MSA[s].resize(seq_length - 4 < 0 ? 0 : seq_length - 4, -1);
    for (int i = 0; i < seq_length - 4; ++i) {
      if (!((seq[i] == 'C' && seq[i + 4] == 'G') || (seq[i] == 'G' && seq[i + 4] == 'C')))
        continue;
      const char *ts;
      if ((ts = strstr(Triloops, seq.substr(i, 5).c_str())))
        if_triloops_MSA[s][i] = (ts - Triloops) / 6;
    }

    // Hexaloops
    if_hexaloops_MSA[s].resize(seq_length - 7 < 0 ? 0 : seq_length - 7, -1);
    for (int i = 0; i < seq_length - 7; ++i) {
      if (!(seq[i] == 'A' && seq[i + 7] == 'U'))
        continue;
      const char *ts;
      if ((ts = strstr(Hexaloops, seq.substr(i, 8).c_str())))
        if_hexaloops_MSA[s][i] = (ts - Hexaloops) / 9;
    }
  }
  return;
}