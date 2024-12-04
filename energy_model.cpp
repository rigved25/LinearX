#include "energy_model.hpp"

#include <iostream>

void EnergyModel::init_tetra_hex_tri(std::string &seq, int seq_length, std::vector<int> &if_tetraloops,
                                     std::vector<int> &if_hexaloops, std::vector<int> &if_triloops) {
    // TetraLoops
    if_tetraloops.resize(seq_length - 5 < 0 ? 0 : seq_length - 5, -1);
    for (int i = 0; i < seq_length - 5; ++i) {
        if (!(seq[i] == 'C' && seq[i + 5] == 'G')) continue;
        const char *ts;
        if ((ts = strstr(Tetraloops, seq.substr(i, 6).c_str()))) if_tetraloops[i] = (ts - Tetraloops) / 7;
    }

    // Triloops
    if_triloops.resize(seq_length - 4 < 0 ? 0 : seq_length - 4, -1);
    for (int i = 0; i < seq_length - 4; ++i) {
        if (!((seq[i] == 'C' && seq[i + 4] == 'G') || (seq[i] == 'G' && seq[i + 4] == 'C'))) continue;
        const char *ts;
        if ((ts = strstr(Triloops, seq.substr(i, 5).c_str()))) if_triloops[i] = (ts - Triloops) / 6;
    }

    // Hexaloops
    if_hexaloops.resize(seq_length - 7 < 0 ? 0 : seq_length - 7, -1);
    for (int i = 0; i < seq_length - 7; ++i) {
        if (!(seq[i] == 'A' && seq[i + 7] == 'U')) continue;
        const char *ts;
        if ((ts = strstr(Hexaloops, seq.substr(i, 8).c_str()))) if_hexaloops[i] = (ts - Hexaloops) / 9;
    }
    return;
}
