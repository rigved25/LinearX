// src/partition/forward.cpp
#include <linearx/partition/linear_partition.hpp>

using namespace linearx::utils;
using namespace linearx::math;
using namespace linearx::utils::io;
using namespace linearx::constants::energy;
using namespace linearx::constants::math;
using namespace std;

unsigned LinearPartition::beam_prune(unordered_map<int, State> &beamstep, const int beam_size) {
    if (beam_size == 0 || beamstep.size() <= beam_size) {
        return 0;
    }
    unsigned num_pruned = 0;
    beam_scores.clear();
    beam_scores.reserve(beamstep.size());
    for (auto &item : beamstep) {
        const int i = item.first;
        const State &cand = item.second;
        const int k = i - 1;
        const value_type newalpha = ((k >= 0 ? bestC[k].alpha : 0) + cand.alpha);
        beam_scores.emplace_back(newalpha, i);
    }
    const value_type threshold = quickselect(beam_scores, 0, beam_scores.size() - 1, beam_scores.size() - beam_size);
    for (auto &p : beam_scores) {
        if (p.first < threshold) {
            beamstep.erase(p.second);
            num_pruned++;
        }
    }
    return num_pruned;
}

inline __attribute__((always_inline)) bool LinearPartition::check_state(const StateType type, const int i,
                                                                        const int j) {
    return true;
}

inline __attribute__((always_inline)) void LinearPartition::update_score(const Mode mode, State &state,
                                                                         const int new_score,
                                                                         const value_type prev_score) {
    if (mode == BEST) {
        state.alpha = max(state.alpha, prev_score + new_score);
    } else {
        state.alpha = LOG_SUM(state.alpha, prev_score + (((value_type)new_score) * INV_KT));
    }
}

PartitionInsideLog LinearPartition::compute_inside(const Mode mode, const unsigned beam_size,
                                                   const bool verbose_output) {
    const auto start_time = chrono::high_resolution_clock::now();
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
        nodes_pruned += beam_prune(bestH[j], beam_size);
        beamstep_H(j, mode);
        if (j == 0) {
            continue;
        }
        // beam of Multi
        nodes_pruned += beam_prune(bestMulti[j], beam_size);
        beamstep_Multi(j, mode);
        // beam of P
        nodes_pruned += beam_prune(bestP[j], beam_size);
        beamstep_P(j, mode);
        // beam of M2
        nodes_pruned += beam_prune(bestM2[j], beam_size);
        beamstep_M2(j, mode);
        // beam of M
        nodes_pruned += beam_prune(bestM[j], beam_size);
        beamstep_M(j, mode);
        // beam of C
        beamstep_C(j, mode);  // beam of C
    }
    // update/print time stats
    const auto end_time = chrono::high_resolution_clock::now();
    const double execution_time = chrono::duration_cast<chrono::milliseconds>(end_time - start_time).count();
    _last_inside_exec_time = execution_time;  // store the last inside execution time
    if (verbose_output) {
        fprintf(stderr, "  - Execution Time: %.3f ms\n", execution_time);
        fprintf(stderr, "  - Nodes Pruned: %lu\n", nodes_pruned);
        if (mode == Mode::BEST) {
            fprintf(stderr, "  - MFE (Minimum Free Energy): %.5f kcal/mol\n", bestC[seq_length - 1].alpha / -100.0);
        } else {
            fprintf(stderr, "  - Free Energy of the Ensemble: %.5f kcal/mol\n", get_ensemble_energy());
        }
        fprintf(stderr, "\n");
    }
    return PartitionInsideLog{bestC[seq_length - 1].alpha, execution_time, beam_size, nodes_pruned};
}

void LinearPartition::beamstep_H(const int j, const Mode mode) {
    int jnext = next_pair[seq.enc[j]][j];
    while (!allow_sharp_turn && jnext < seq_length && (jnext - j) < 4) {
        jnext = next_pair[seq.enc[j]][jnext];
    }

    if (jnext < seq_length && check_state(StateType::H, j, jnext)) {
        int tetra_hex_tri = -1;
        if (jnext - j - 1 == 4)  // 6:tetra
            tetra_hex_tri = if_tetraloops[j];
        else if (jnext - j - 1 == 6)  // 8:hexa
            tetra_hex_tri = if_hexaloops[j];
        else if (jnext - j - 1 == 3)  // 5:tri
            tetra_hex_tri = if_triloops[j];

        const int score = -energy_model.score_hairpin(j, jnext, seq.enc[j], seq.enc[j + 1], seq.enc[jnext - 1],
                                                      seq.enc[jnext], tetra_hex_tri);
        update_score(mode, bestH[jnext][j], score);
    }

    // for every state h in H[j]
    //   1. extend H(i, j) to H(i, jnext)
    //   2. generate P(i, j)
    for (auto &item : bestH[j]) {
        const int i = item.first;
        const State &state = item.second;

        // 1. extend H(i, j) to H(i, jnext)
        jnext = next_pair[seq.enc[i]][j];
        if (jnext < seq_length && check_state(StateType::H, i, jnext)) {
            int tetra_hex_tri = -1;
            if (jnext - i - 1 == 4)  // 6:tetra
                tetra_hex_tri = if_tetraloops[i];
            else if (jnext - i - 1 == 6)  // 8:hexa
                tetra_hex_tri = if_hexaloops[i];
            else if (jnext - i - 1 == 3)  // 5:tri
                tetra_hex_tri = if_triloops[i];

            const int score = -energy_model.score_hairpin(i, jnext, seq.enc[i], seq.enc[i + 1], seq.enc[jnext - 1],
                                                          seq.enc[jnext], tetra_hex_tri);
            update_score(mode, bestH[jnext][i], score);
        }

        // 2. generate P(i, j)
        if (check_state(StateType::P, i, j)) {
            update_score(mode, bestP[j][i], 0, state.alpha);
        }
    }
}

void LinearPartition::beamstep_Multi(const int j, const Mode mode) {
    for (auto &item : bestMulti[j]) {
        const int i = item.first;
        const State &state = item.second;

        // 1. extend Multi(i, j) to Multi(i, jnext)
        const int jnext = next_pair[seq.enc[i]][j];
        if (jnext < seq_length && check_state(StateType::Multi, i, jnext)) {
            const int new_score = -energy_model.score_multi_unpaired(j, jnext);
            update_score(mode, bestMulti[jnext][i], new_score, state.alpha);
        }

        // 2. generate P(i, j)
        if (check_state(StateType::P, i, j)) {
            const int new_score =
                -energy_model.score_multi(i, j, seq.enc[i], seq.enc[i + 1], seq.enc[j - 1], seq.enc[j], seq_length);
            update_score(mode, bestP[j][i], new_score, state.alpha);
        }
    }
}

void LinearPartition::beamstep_P(const int j, const Mode mode) {
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
                if (check_state(StateType::P, p, q)) {
                    const int new_score = -energy_model.score_single_loop(p, q, i, j, seq.enc[p], seq.enc[p + 1],
                                                                          seq.enc[q - 1], seq.enc[q], seq.enc[i - 1],
                                                                          seq.enc[i], seq.enc[j], seq.enc[j + 1]);
                    update_score(mode, bestP[q][p], new_score, state.alpha);
                }
                q = next_pair[seq.enc[p]][q];
            }
        }

        // 2. M = P
        if (i > 0 && j < seq_length - 1 && check_state(StateType::M, i, j)) {
            const int new_score =
                -energy_model.score_M1(i, j, j, seq.enc[i - 1], seq.enc[i], seq.enc[j], nucj1, seq_length);
            update_score(mode, bestM[j][i], new_score, state.alpha);
        }

        // 3. M2 = M + P
        const int h = i - 1;
        if (h > 0 && !bestM[h].empty()) {
            const int new_score =
                -energy_model.score_M1(i, j, j, seq.enc[i - 1], seq.enc[i], seq.enc[j], nucj1, seq_length);
            for (auto &m_item : bestM[h]) {
                const int g = m_item.first;
                if (check_state(StateType::M2, g, j)) {
                    const State &m_state = m_item.second;
                    update_score(mode, bestM2[j][g], new_score, m_state.alpha + state.alpha);
                }
            }
        }

        // 4. C = C + P
        if (h >= 0) {
            const State &prefix_C = bestC[h];
            const int new_score =
                -energy_model.score_external_paired(i, j, seq.enc[h], seq.enc[h + 1], seq.enc[j], nucj1, seq_length);
            update_score(mode, bestC[j], new_score, prefix_C.alpha + state.alpha);
        } else {
            const int new_score =
                -energy_model.score_external_paired(0, j, -1, seq.enc[0], seq.enc[j], nucj1, seq_length);
            update_score(mode, bestC[j], new_score, state.alpha);
        }
    }
}

void LinearPartition::beamstep_M2(const int j, const Mode mode) {
    for (auto &item : bestM2[j]) {
        const int i = item.first;
        const State &state = item.second;

        // 1. multi-loop = M2
        int q = seq_length;
        for (int p = i - 1; p >= max(0, i - MAXLOOPSIZE); --p) {
            q = next_pair[seq.enc[p]][j];
            if (q < seq_length && check_state(StateType::Multi, p, q)) {
                const int new_score =
                    -(energy_model.score_multi_unpaired(p, i - 1) + energy_model.score_multi_unpaired(j, q - 1));
                update_score(mode, bestMulti[q][p], new_score, state.alpha);
            }
        }

        // 2. M = M2
        if (check_state(StateType::M, i, j)) {
            update_score(mode, bestM[j][i], 0, state.alpha);
        }
    }
}

void LinearPartition::beamstep_M(const int j, const Mode mode) {
    for (auto &item : bestM[j]) {
        const int i = item.first;
        const State &state = item.second;

        if (j < seq_length - 1 && check_state(StateType::M, i, j + 1)) {
            const int new_score = -energy_model.score_multi_unpaired(j, j + 1);
            update_score(mode, bestM[j + 1][i], new_score, state.alpha);
        }
    }
}

void LinearPartition::beamstep_C(const int j, const Mode mode) {
    if (j < seq_length - 1) {
        const int new_score = -energy_model.score_external_unpaired(j + 1, j + 1);
        update_score(mode, bestC[j + 1], new_score, bestC[j].alpha);
    }
}
