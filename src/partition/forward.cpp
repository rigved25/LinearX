// src/partition/forward.cpp
#include <linearx/partition/config.hpp>
#include <linearx/partition/linear_partition.hpp>

using namespace linearx::utils;
using namespace linearx::math;
using namespace linearx::utils::io;
using namespace linearx::constants::energy;
using namespace linearx::constants::math;
using namespace std;

template <typename T>
PartitionLog LinearPartitionInterface<T>::compute_inside(const Mode mode, const unsigned beam_size,
                                                         const bool verbose_output) {
    const auto start_time = chrono::high_resolution_clock::now();
    if (verbose_output) {
        fprintf(stderr, "[LinearPartition] Running Inside Algorithm (ID: %d | Length: %zu | Name: %s):\n", seq.id,
                seq.length(), seq.name.c_str());
    }
    unsigned long nodes_pruned = 0;
    for (int j = 0; j < seq_length; j++) {
        if (verbose_output) {
            showProgressBar(j, seq_length - 1);
        }
        // beam of H
        nodes_pruned += beam_prune(bestH[j], beam_size, StateType::H, j);
        beamstep_H<PARTITION_INSIDE>(j);
        if (j == 0) {
            continue;
        }
        // beam of Multi
        nodes_pruned += beam_prune(bestMulti[j], beam_size, StateType::Multi, j);
        beamstep_Multi<PARTITION_INSIDE>(j);
        // beam of P
        nodes_pruned += beam_prune(bestP[j], beam_size, StateType::P, j);
        beamstep_P<PARTITION_INSIDE>(j);
        // beam of M2
        nodes_pruned += beam_prune(bestM2[j], beam_size, StateType::M2, j);
        beamstep_M2<PARTITION_INSIDE>(j);
        // beam of M
        nodes_pruned += beam_prune(bestM[j], beam_size, StateType::M, j);
        beamstep_M<PARTITION_INSIDE>(j);
        // beam of C
        beamstep_C<PARTITION_INSIDE>(j);
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
    beam_scores.clear();
    return PartitionLog{get_ensemble_energy(),
                        bestC[seq_length - 1].alpha,
                        LOG_ZERO,
                        execution_time,
                        -LOG_ZERO,
                        float(beam_size),
                        "Inside",
                        0,
                        nodes_pruned,
                        0,
                        0};
}

template <typename T>
template <Mode mode>
unsigned long LinearPartitionInterface<T>::beamstep_H(const int j) {
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

        if constexpr (mode != Mode::PARTITION_OUTSIDE) {
            const int new_score = -energy_model.score_hairpin(j, jnext, seq.enc[j], seq.enc[j + 1], seq.enc[jnext - 1],
                                                              seq.enc[jnext], tetra_hex_tri);
            if constexpr (mode == Mode::BEST) {
                bestH[jnext][j].alpha = std::max(bestH[jnext][j].alpha, static_cast<value_type>(new_score));
            } else {
                State* next_state = check_state(StateType::H, j, jnext);
                if (next_state) {
                    HEdge(new_score, nullptr, nullptr).update_state_alpha(*next_state);
                }
            }
        }
    }
    // for every state h in H[j]
    //   1. extend H(i, j) to H(i, jnext)
    //   2. generate P(i, j)
    for (auto& item : bestH[j]) {
        const int i = item.first;
        State& state = item.second;
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

            // update_score(StateType::H, i, jnext, score);
            if constexpr (mode == Mode::PARTITION_OUTSIDE) {
                continue;  // do nothing
            } else {
                const int new_score = -energy_model.score_hairpin(i, jnext, seq.enc[i], seq.enc[i + 1],
                                                                  seq.enc[jnext - 1], seq.enc[jnext], tetra_hex_tri);
                if constexpr (mode == Mode::BEST) {
                    bestH[jnext][i].alpha = std::max(bestH[jnext][i].alpha, static_cast<value_type>(new_score));
                } else {
                    State* next_state = check_state(StateType::H, i, jnext);
                    if (next_state) {
                        HEdge(new_score, nullptr, nullptr).update_state_alpha(*next_state);
                    }
                }
            }
        }
        // 2. H(i, j) -> P(i, j)
        if constexpr (mode == Mode::PARTITION_OUTSIDE) {
            continue;  // do nothing
        } else {
            const int new_score = 0;
            if constexpr (mode == Mode::BEST) {
                bestP[j][i].alpha = std::max(bestP[j][i].alpha, state.alpha + new_score);
            } else {
                State* next_state = check_state(StateType::P, i, j);
                if (next_state) {
                    HEdge(new_score, &state, nullptr).update_state_alpha(*next_state);
                }
            }
        }
    }
    return 0;
}

template <typename T>
template <Mode mode>
unsigned long LinearPartitionInterface<T>::beamstep_Multi(const int j) {
    for (auto& item : bestMulti[j]) {
        const int i = item.first;
        State& state = item.second;
        // 1. Multi(i, j) -> Multi(i, jnext)
        const int jnext = next_pair[seq.enc[i]][j];
        if (jnext < seq_length) {
            if constexpr (mode == Mode::PARTITION_OUTSIDE) {
                auto it = bestMulti[jnext].find(i);
                if (it != bestMulti[jnext].end()) {
                    const int new_score = -energy_model.score_multi_unpaired(j, jnext);
                    update_state_beta(HEdge(new_score, &state, nullptr), &it->second, StateType::Multi, i, jnext);
                }
            } else {
                const int new_score = -energy_model.score_multi_unpaired(j, jnext);
                if constexpr (mode == Mode::BEST) {
                    bestMulti[jnext][i].alpha = std::max(bestMulti[jnext][i].alpha, state.alpha + new_score);
                } else {
                    State* next_state = check_state(StateType::Multi, i, jnext);
                    if (next_state) {
                        HEdge(new_score, &state, nullptr).update_state_alpha(*next_state);
                    }
                }
            }
        }

        // 2. generate P(i, j)
        if constexpr (mode == Mode::PARTITION_OUTSIDE) {
            auto it = bestP[j].find(i);
            if (it != bestP[j].end()) {
                const int new_score =
                    -energy_model.score_multi(i, j, seq.enc[i], seq.enc[i + 1], seq.enc[j - 1], seq.enc[j], seq_length);
                update_state_beta(HEdge(new_score, &state, nullptr), &it->second, StateType::P, i, j);
            }
        } else {
            const int new_score =
                -energy_model.score_multi(i, j, seq.enc[i], seq.enc[i + 1], seq.enc[j - 1], seq.enc[j], seq_length);
            if constexpr (mode == Mode::BEST) {
                bestP[j][i].alpha = std::max(bestP[j][i].alpha, state.alpha + new_score);
            } else {
                State* next_state = check_state(StateType::P, i, j);
                if (next_state) {
                    HEdge(new_score, &state, nullptr).update_state_alpha(*next_state);
                }
            }
        }
    }
    return 0;
}

template <typename T>
template <Mode mode>
unsigned long LinearPartitionInterface<T>::beamstep_P(const int j) {
    const int nucj1 = (j + 1 < seq_length ? seq.enc[j + 1] : -1);
    for (auto& item : bestP[j]) {
        const int i = item.first;
        State& state = item.second;
        // 1. P(i, j) -> P(p, q) [scan left & jump right, generate new P (helix/bulge)]
        int q = seq_length;
        for (int p = i - 1; p >= max(0, i - MAXLOOPSIZE); --p) {
            q = next_pair[seq.enc[p]][j];
            while (q < seq_length && ((i - p) + (q - j) - 2) <= MAXLOOPSIZE) {
                // current shape is: p...i (pair) j...q
                if constexpr (mode == Mode::PARTITION_OUTSIDE) {
                    auto it = bestP[q].find(p);
                    if (it != bestP[q].end()) {
                        const int new_score =
                            -energy_model.score_single_loop(p, q, i, j, seq.enc[p], seq.enc[p + 1], seq.enc[q - 1],
                                                            seq.enc[q], seq.enc[i - 1], seq.enc[i], seq.enc[j], nucj1);
                        update_state_beta(HEdge(new_score, &state, nullptr), &it->second, StateType::P, p, q);
                    }
                } else {
                    const int new_score =
                        -energy_model.score_single_loop(p, q, i, j, seq.enc[p], seq.enc[p + 1], seq.enc[q - 1],
                                                        seq.enc[q], seq.enc[i - 1], seq.enc[i], seq.enc[j], nucj1);
                    if constexpr (mode == Mode::BEST) {
                        bestP[q][p].alpha = std::max(bestP[q][p].alpha, state.alpha + new_score);
                    } else {
                        State* next_state = check_state(StateType::P, p, q);
                        if (next_state) {
                            HEdge(new_score, &state, nullptr).update_state_alpha(*next_state);
                        }
                    }
                }
                q = next_pair[seq.enc[p]][q];
            }
        }
        // 2. P -> M
        if (i > 0 && j < seq_length - 1) {
            if constexpr (mode == Mode::PARTITION_OUTSIDE) {
                auto it = bestM[j].find(i);
                if (it != bestM[j].end()) {
                    const int new_score =
                        -energy_model.score_M1(i, j, j, seq.enc[i - 1], seq.enc[i], seq.enc[j], nucj1, seq_length);
                    update_state_beta(HEdge(new_score, &state, nullptr), &it->second, StateType::M, i, j);
                }
            } else {
                const int new_score =
                    -energy_model.score_M1(i, j, j, seq.enc[i - 1], seq.enc[i], seq.enc[j], nucj1, seq_length);
                if constexpr (mode == Mode::BEST) {
                    bestM[j][i].alpha = std::max(bestM[j][i].alpha, state.alpha + new_score);
                } else {
                    State* next_state = check_state(StateType::M, i, j);
                    if (next_state) {
                        HEdge(new_score, &state, nullptr).update_state_alpha(*next_state);
                    }
                }
            }
        }
        // 3. M + P -> M2
        const int h = i - 1;
        if (h > 0 && !bestM[h].empty()) {
            const int new_score =
                -energy_model.score_M1(i, j, j, seq.enc[i - 1], seq.enc[i], seq.enc[j], nucj1, seq_length);
            for (auto& m_item : bestM[h]) {
                const int g = m_item.first;
                State& m_state = m_item.second;
                if constexpr (mode == Mode::PARTITION_OUTSIDE) {
                    auto it = bestM2[j].find(g);
                    if (it != bestM2[j].end()) {
                        update_state_beta(HEdge(new_score, &m_state, &state), &it->second, StateType::M2, g, j);
                    }
                } else {
                    if constexpr (mode == Mode::BEST) {
                        bestM2[j][g].alpha = std::max(bestM2[j][g].alpha, m_state.alpha + state.alpha + new_score);
                    } else {
                        State* next_state = check_state(StateType::M2, g, j);
                        if (next_state) {
                            HEdge(new_score, &m_state, &state).update_state_alpha(*next_state);
                        }
                    }
                }
            }
        }

        // // 4. C + P -> C
        if constexpr (mode == Mode::PARTITION_OUTSIDE) {
            const int new_score = -energy_model.score_external_paired(i, j, i > 0 ? seq.enc[h] : -1, seq.enc[i],
                                                                      seq.enc[j], nucj1, seq_length);
            update_state_beta(HEdge(new_score, &bestC[h], &state), &bestC[j], StateType::C, 0, j);
        } else {
            const int new_score = -energy_model.score_external_paired(i, j, i > 0 ? seq.enc[h] : -1, seq.enc[i],
                                                                      seq.enc[j], nucj1, seq_length);
            if constexpr (mode == Mode::BEST) {
                bestC[j].alpha = std::max(bestC[j].alpha, bestC[h].alpha + state.alpha + new_score);
            } else {
                HEdge(new_score, &bestC[h], &state).update_state_alpha(bestC[j]);
            }
        }
    }
    return 0;
}

template <typename T>
template <Mode mode>
unsigned long LinearPartitionInterface<T>::beamstep_M2(const int j) {
    for (auto& item : bestM2[j]) {
        const int i = item.first;
        State& state = item.second;
        // 1. M2 -> Multi
        int q = seq_length;
        for (int p = i - 1; p >= max(0, i - MAXLOOPSIZE); --p) {
            q = next_pair[seq.enc[p]][j];
            if (q < seq_length) {
                if constexpr (mode == Mode::PARTITION_OUTSIDE) {
                    auto it = bestMulti[q].find(p);
                    if (it != bestMulti[q].end()) {
                        const int new_score = -(energy_model.score_multi_unpaired(p, i - 1) +
                                                energy_model.score_multi_unpaired(j, q - 1));
                        update_state_beta(HEdge(new_score, &state, nullptr), &it->second, StateType::Multi, p, q);
                    }
                } else {
                    const int new_score =
                        -(energy_model.score_multi_unpaired(p, i - 1) + energy_model.score_multi_unpaired(j, q - 1));
                    if constexpr (mode == Mode::BEST) {
                        bestMulti[q][p].alpha = std::max(bestMulti[q][p].alpha, state.alpha + new_score);
                    } else {
                        State* next_state = check_state(StateType::Multi, p, q);
                        if (next_state) {
                            HEdge(new_score, &state, nullptr).update_state_alpha(*next_state);
                        }
                    }
                }
            }
        }
        // 2. M2 -> M
        if constexpr (mode == Mode::PARTITION_OUTSIDE) {
            auto it = bestM[j].find(i);
            if (it != bestM[j].end()) {
                const int new_score = 0;
                update_state_beta(HEdge(new_score, &state, nullptr), &it->second, StateType::M, i, j);
            }
        } else {
            const int new_score = 0;
            if constexpr (mode == Mode::BEST) {
                bestM[j][i].alpha = std::max(bestM[j][i].alpha, state.alpha + new_score);
            } else {
                State* next_state = check_state(StateType::M, i, j);
                if (next_state) {
                    HEdge(new_score, &state, nullptr).update_state_alpha(*next_state);
                }
            }
        }
    }
    return 0;
}

template <typename T>
template <Mode mode>
unsigned long LinearPartitionInterface<T>::beamstep_M(const int j) {
    // M -> M + U
    for (auto& item : bestM[j]) {
        const int i = item.first;
        State& state = item.second;
        if (j < seq_length - 1) {
            if constexpr (mode == Mode::PARTITION_OUTSIDE) {
                auto it = bestM[j + 1].find(i);
                if (it != bestM[j + 1].end()) {
                    const int new_score = -energy_model.score_multi_unpaired(j, j + 1);
                    update_state_beta(HEdge(new_score, &state, nullptr), &it->second, StateType::M, i, j + 1);
                }
            } else {
                const int new_score = -energy_model.score_multi_unpaired(j, j + 1);
                if constexpr (mode == Mode::BEST) {
                    bestM[j + 1][i].alpha = std::max(bestM[j + 1][i].alpha, state.alpha + new_score);
                } else {
                    State* next_state = check_state(StateType::M, i, j + 1);
                    if (next_state) {
                        HEdge(new_score, &state, nullptr).update_state_alpha(*next_state);
                    }
                }
            }
        }
    }
    return 0;
}

template <typename T>
template <Mode mode>
void LinearPartitionInterface<T>::beamstep_C(const int j) {
    // C -> C + U
    if (j < seq_length - 1) {
        const int new_score = -energy_model.score_external_unpaired(j + 1, j + 1);
        if constexpr (mode == Mode::PARTITION_OUTSIDE) {
            update_state_beta(HEdge(new_score, &bestC[j], nullptr), &bestC[j + 1], StateType::C, 0, j + 1);
        } else {
            if constexpr (mode == Mode::BEST) {
                bestC[j + 1].alpha = std::max(bestC[j + 1].alpha, bestC[j].alpha + new_score);
            } else {
                HEdge(new_score, &bestC[j], nullptr).update_state_alpha(bestC[j + 1]);
            }
        }
    }
}

// Instantiate templates for LinearPartitionInterface with desired types
#define INSTANTIATE_FUNCS_WITH_MODE(FUNC, TYPE)                                               \
    template unsigned long LinearPartitionInterface<TYPE>::FUNC<Mode::BEST>(int);             \
    template unsigned long LinearPartitionInterface<TYPE>::FUNC<Mode::PARTITION_INSIDE>(int); \
    template unsigned long LinearPartitionInterface<TYPE>::FUNC<Mode::PARTITION_OUTSIDE>(int);

#define DECLARE_FUNCS(TYPE)                                                                 \
    template class LinearPartitionInterface<TYPE>;                                          \
    template void LinearPartitionInterface<TYPE>::beamstep_C<Mode::BEST>(int);              \
    template void LinearPartitionInterface<TYPE>::beamstep_C<Mode::PARTITION_INSIDE>(int);  \
    template void LinearPartitionInterface<TYPE>::beamstep_C<Mode::PARTITION_OUTSIDE>(int); \
    INSTANTIATE_FUNCS_WITH_MODE(beamstep_H, TYPE)                                           \
    INSTANTIATE_FUNCS_WITH_MODE(beamstep_Multi, TYPE)                                       \
    INSTANTIATE_FUNCS_WITH_MODE(beamstep_P, TYPE)                                           \
    INSTANTIATE_FUNCS_WITH_MODE(beamstep_M2, TYPE)                                          \
    INSTANTIATE_FUNCS_WITH_MODE(beamstep_M, TYPE)

#define X(TYPE) DECLARE_FUNCS(TYPE)
LP_TEMPLATE_TYPES
#undef X
#undef DECLARE_FUNCS
#undef INSTANTIATE_FUNCS_WITH_MODE
