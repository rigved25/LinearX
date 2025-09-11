// src/partition/forward.cpp
#include <linearx/partition/linear_partition.hpp>

using namespace linearx::utils;
using namespace linearx::math;
using namespace linearx::utils::io;
using namespace linearx::constants::energy;
using namespace linearx::constants::math;
using namespace std;

// unsigned long LinearPartition::beam_prune(const int j, const StateType type, std::unordered_map<int, State> &beamstep) {
//     if (run_beam_size_ == 0 || beamstep.size() <= run_beam_size_) {
//         return 0;
//     }
//     unsigned long num_pruned = 0;
//     beam_scores.clear();
//     beam_scores.reserve(beamstep.size());
//     for (auto &item : beamstep) {
//         const int i = item.first;
//         const State &cand = item.second;
//         const int k = i - 1;
//         const value_type newalpha = ((k >= 0 ? bestC[k].alpha : 0) + cand.alpha);
//         beam_scores.emplace_back(newalpha, i);
//     }
//     const value_type threshold =
//         linearx::utils::quickselect(beam_scores, 0, beam_scores.size() - 1, beam_scores.size() - run_beam_size_);
//     for (auto &p : beam_scores) {
//         if (p.first < threshold) {
//             beamstep.erase(p.second);
//             num_pruned++;
//         }
//     }
//     return num_pruned;
// }
// void LinearPartition::update_score(const StateType type, const int i, const int j, State &state, int new_score,
//                                    value_type prev_score) {
//     if (run_mode_ == BEST) {
//         state.alpha = std::max(state.alpha, prev_score + new_score);
//     } else {
//         state.alpha = LOG_SUM(state.alpha, prev_score + (new_score * linearx::constants::energy::INV_KT));
//     }
// }

PartitionInsideLog LinearPartition::compute_inside(const Mode mode, const unsigned beam_size,
                                                   const bool verbose_output) {
    const auto start_time = chrono::high_resolution_clock::now();
    run_mode_ = mode;
    run_beam_size_ = beam_size;
    if (verbose_output) {
        fprintf(stderr, "[LinearPartition] Running Inside Algorithm (ID: %d | Length: %zu | Name: %s):\n", seq.id,
                seq.length(), seq.name.c_str());
    }
    reset_beams();
    unsigned long nodes_pruned = 0;
    for (int j = 0; j < seq_length; j++) {
        if (verbose_output) {
            showProgressBar(j, seq_length - 1);
        }
        // beam of H
        nodes_pruned += beamstep_H(j);
        if (j == 0) {
            continue;
        }
        // beam of Multi
        nodes_pruned += beamstep_Multi(j);
        // beam of P
        nodes_pruned += beamstep_P(j);
        // beam of M2
        nodes_pruned += beamstep_M2(j);
        // beam of M
        nodes_pruned += beamstep_M(j);
        // beam of C
        beamstep_C(j);  // beam of C
    }
    // update/print time stats
    const auto end_time = chrono::high_resolution_clock::now();
    const double execution_time = chrono::duration_cast<chrono::milliseconds>(end_time - start_time).count();
    _last_inside_exec_time = execution_time;  // store the last inside execution time
    const value_type score = mode == Mode::BEST ? bestC[seq_length - 1].alpha / -100 : get_ensemble_energy();
    if (verbose_output) {
        fprintf(stderr, "  - Execution Time: %.3f ms\n", execution_time);
        fprintf(stderr, "  - Nodes Pruned: %lu\n", nodes_pruned);
        if (mode == Mode::BEST) {
            fprintf(stderr, "  - MFE (Minimum Free Energy): %.5f kcal/mol\n", score);
        } else {
            fprintf(stderr, "  - Free Energy of the Ensemble: %.5f kcal/mol\n", score);
        }
    }
    return PartitionInsideLog{score, execution_time, beam_size, nodes_pruned};
}

unsigned long LinearPartition::beamstep_H(const int j) {
    unsigned long nodes_pruned = beam_prune(j, StateType::H, bestH[j]);
    int jnext = next_pair[seq.enc[j]][j];
    while (!allow_sharp_turn && jnext < seq_length && (jnext - j) < 4) {
        jnext = next_pair[seq.enc[j]][jnext];
    }
    if (jnext < seq_length) {
        int tetra_hex_tri = -1;
        if (jnext - j - 1 == 4)  // 6:tetra
            tetra_hex_tri = if_tetraloops[j];
        else if (jnext - j - 1 == 6)  // 8:hexa
            tetra_hex_tri = if_hexaloops[j];
        else if (jnext - j - 1 == 3)  // 5:tri
            tetra_hex_tri = if_triloops[j];
        const int score = -energy_model.score_hairpin(j, jnext, seq.enc[j], seq.enc[j + 1], seq.enc[jnext - 1],
                                                      seq.enc[jnext], tetra_hex_tri);
        update_score(StateType::H, j, jnext, bestH[jnext][j], score);
    }
    // for every state h in H[j]
    //   1. extend H(i, j) to H(i, jnext)
    //   2. generate P(i, j)
    for (auto &item : bestH[j]) {
        const int i = item.first;
        const State &state = item.second;
        // 1. extend H(i, j) to H(i, jnext)
        jnext = next_pair[seq.enc[i]][j];
        if (jnext < seq_length) {
            int tetra_hex_tri = -1;
            if (jnext - i - 1 == 4)  // 6:tetra
                tetra_hex_tri = if_tetraloops[i];
            else if (jnext - i - 1 == 6)  // 8:hexa
                tetra_hex_tri = if_hexaloops[i];
            else if (jnext - i - 1 == 3)  // 5:tri
                tetra_hex_tri = if_triloops[i];
            const int score = -energy_model.score_hairpin(i, jnext, seq.enc[i], seq.enc[i + 1], seq.enc[jnext - 1],
                                                          seq.enc[jnext], tetra_hex_tri);
            update_score(StateType::H, i, jnext, bestH[jnext][i], score);
        }
        // 2. generate P(i, j)
        update_score(StateType::P, i, j, bestP[j][i], 0, state.alpha);
    }
    return nodes_pruned;
}

unsigned long LinearPartition::beamstep_Multi(const int j) {
    unsigned long nodes_pruned = beam_prune(j, StateType::Multi, bestMulti[j]);
    for (auto &item : bestMulti[j]) {
        const int i = item.first;
        const State &state = item.second;
        // 1. extend Multi(i, j) to Multi(i, jnext)
        const int jnext = next_pair[seq.enc[i]][j];
        if (jnext < seq_length) {
            const int new_score = -energy_model.score_multi_unpaired(j, jnext);
            update_score(StateType::Multi, i, jnext, bestMulti[jnext][i], new_score, state.alpha);
        }
        // 2. generate P(i, j)
        const int new_score =
            -energy_model.score_multi(i, j, seq.enc[i], seq.enc[i + 1], seq.enc[j - 1], seq.enc[j], seq_length);
        update_score(StateType::P, i, j, bestP[j][i], new_score, state.alpha);
    }
    return nodes_pruned;
}

unsigned long LinearPartition::beamstep_P(const int j) {
    unsigned long nodes_pruned = beam_prune(j, StateType::P, bestP[j]);
    // for every state in P[j]
    //   1. generate new P (helix/bulge)
    //   2. M = P
    //   3. M2 = M + P
    //   4. C = C + P
    const int nucj1 = (j + 1 < seq_length ? seq.enc[j + 1] : -1);
    for (auto &item : bestP[j]) {
        const int i = item.first;
        const State &state = item.second;
        // 1. generate new P (helix/bulge)
        int q = seq_length;
        for (int p = i - 1; p >= max(0, i - MAXLOOPSIZE); --p) {
            q = next_pair[seq.enc[p]][j];
            while (q < seq_length && ((i - p) + (q - j) - 2) <= MAXLOOPSIZE) {
                // current shape is: p...i (pair) j...q
                const int new_score =
                    -energy_model.score_single_loop(p, q, i, j, seq.enc[p], seq.enc[p + 1], seq.enc[q - 1], seq.enc[q],
                                                    seq.enc[i - 1], seq.enc[i], seq.enc[j], seq.enc[j + 1]);
                update_score(StateType::P, p, q, bestP[q][p], new_score, state.alpha);
                q = next_pair[seq.enc[p]][q];
            }
        }
        // 2. M = P
        if (i > 0 && j < seq_length - 1) {
            const int new_score =
                -energy_model.score_M1(i, j, j, seq.enc[i - 1], seq.enc[i], seq.enc[j], nucj1, seq_length);
            update_score(StateType::M, i, j, bestM[j][i], new_score, state.alpha);
        }
        // 3. M2 = M + P
        const int h = i - 1;
        if (h > 0 && !bestM[h].empty()) {
            const int new_score =
                -energy_model.score_M1(i, j, j, seq.enc[i - 1], seq.enc[i], seq.enc[j], nucj1, seq_length);
            for (auto &m_item : bestM[h]) {
                const int g = m_item.first;
                const State &m_state = m_item.second;
                update_score(StateType::M2, g, j, bestM2[j][g], new_score, m_state.alpha + state.alpha);
            }
        }
        // 4. C = C + P
        if (h >= 0) {
            const State &prefix_C = bestC[h];
            const int new_score =
                -energy_model.score_external_paired(i, j, seq.enc[h], seq.enc[h + 1], seq.enc[j], nucj1, seq_length);
            update_score(StateType::C, 0, j, bestC[j], new_score, prefix_C.alpha + state.alpha);
        } else {
            const int new_score =
                -energy_model.score_external_paired(0, j, -1, seq.enc[0], seq.enc[j], nucj1, seq_length);
            update_score(StateType::C, 0, j, bestC[j], new_score, state.alpha);
        }
    }
    return nodes_pruned;
}

unsigned long LinearPartition::beamstep_M2(const int j) {
    unsigned long nodes_pruned = beam_prune(j, StateType::M2, bestM2[j]);
    for (auto &item : bestM2[j]) {
        const int i = item.first;
        const State &state = item.second;
        // 1. multi-loop = M2
        int q = seq_length;
        for (int p = i - 1; p >= max(0, i - MAXLOOPSIZE); --p) {
            q = next_pair[seq.enc[p]][j];
            if (q < seq_length) {
                const int new_score =
                    -(energy_model.score_multi_unpaired(p, i - 1) + energy_model.score_multi_unpaired(j, q - 1));
                update_score(StateType::Multi, p, q, bestMulti[q][p], new_score, state.alpha);
            }
        }
        // 2. M = M2
        update_score(StateType::M, i, j, bestM[j][i], 0, state.alpha);
    }
    return nodes_pruned;
}

unsigned long LinearPartition::beamstep_M(const int j) {
    unsigned long nodes_pruned = beam_prune(j, StateType::M, bestM[j]);
    for (auto &item : bestM[j]) {
        const int i = item.first;
        const State &state = item.second;
        if (j < seq_length - 1) {
            const int new_score = -energy_model.score_multi_unpaired(j, j + 1);
            update_score(StateType::M, i, j + 1, bestM[j + 1][i], new_score, state.alpha);
        }
    }
    return nodes_pruned;
}

void LinearPartition::beamstep_C(const int j) {
    if (j < seq_length - 1) {
        const int new_score = -energy_model.score_external_unpaired(j + 1, j + 1);
        update_score(StateType::C, 0, j + 1, bestC[j + 1], new_score, bestC[j].alpha);
    }
}
