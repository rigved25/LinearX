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
void LinearPartitionInterface<T>::compute_inside(const Mode mode, const unsigned beam_size, const bool verbose_output) {
    const auto start_time = chrono::high_resolution_clock::now();
    if (verbose_output) {
        fprintf(stderr, "[LinearPartition] Running Inside Algorithm (ID: %d | Length: %zu | Name: %s):\n", seq.id,
                seq.length(), seq.name.c_str());
    }
    unsigned long states_pruned = 0;
    for (int j = 0; j < seq_length; j++) {
        if (verbose_output) {
            showProgressBar(j, seq_length - 1);
        }
        // beam of H
        states_pruned += beam_prune(bestH[j], beam_size, StateType::H, j);
        beamstep_H<PARTITION_INSIDE>(j);
        if (j == 0) {
            continue;
        }
        // beam of Multi
        states_pruned += beam_prune(bestMulti[j], beam_size, StateType::Multi, j);
        beamstep_Multi<PARTITION_INSIDE>(j);
        // beam of P
        states_pruned += beam_prune(bestP[j], beam_size, StateType::P, j);
        beamstep_P<PARTITION_INSIDE>(j);
        // beam of M2
        states_pruned += beam_prune(bestM2[j], beam_size, StateType::M2, j);
        beamstep_M2<PARTITION_INSIDE>(j);
        // beam of M
        states_pruned += beam_prune(bestM[j], beam_size, StateType::M, j);
        beamstep_M<PARTITION_INSIDE>(j);
        // beam of C
        beamstep_C<PARTITION_INSIDE>(j);
    }
    // update/print time stats
    const auto end_time = chrono::high_resolution_clock::now();
    const double execution_time = chrono::duration_cast<chrono::milliseconds>(end_time - start_time).count();
    const value_type score = mode == Mode::BEST ? bestC[seq_length - 1].alpha / -100 : get_ensemble_energy();
    if (verbose_output) {
        fprintf(stderr, "  - Execution Time: %.3f ms\n", execution_time);
        fprintf(stderr, "  - Nodes Pruned: %lu\n", states_pruned);
        if (mode == Mode::BEST) {
            fprintf(stderr, "  - MFE (Minimum Free Energy): %.5f kcal/mol\n", score);
        } else {
            fprintf(stderr, "  - Free Energy of the Ensemble: %.5f kcal/mol\n", score);
        }
    }
    beam_scores.clear();

    // update logs
    log.free_energy_of_ensemble = get_ensemble_energy();
    log.total_energy.first = bestC[seq_length - 1].alpha;
    log.exec_time.first = execution_time;
    log.states_pruned.first = states_pruned;
    log.effective_beam_size.first = beam_size;
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

        const auto get_score = [this, j, jnext, tetra_hex_tri]() constexpr -> int {
            return -energy_model.score_hairpin(j, jnext, seq.enc[j], seq.enc[j + 1], seq.enc[jnext - 1], seq.enc[jnext],
                                               tetra_hex_tri);
        };
        update_state<mode>(get_score, nullptr, nullptr, StateType::H, j, jnext);
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

            const auto get_score = [this, i, jnext, tetra_hex_tri]() constexpr -> int {
                return -energy_model.score_hairpin(i, jnext, seq.enc[i], seq.enc[i + 1], seq.enc[jnext - 1],
                                                   seq.enc[jnext], tetra_hex_tri);
            };
            update_state<mode>(get_score, nullptr, nullptr, StateType::H, i, jnext);
        }
        // 2. H(i, j) -> P(i, j)
        const auto get_score = []() constexpr -> int { return 0; };
        update_state<mode>(get_score, &state, nullptr, StateType::P, i, j);
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
            const auto get_score = [this, j, jnext]() constexpr -> int {
                return -energy_model.score_multi_unpaired(j, jnext);
            };
            update_state<mode>(get_score, &state, nullptr, StateType::Multi, i, jnext);
        }

        // 2. generate P(i, j)
        const auto get_score = [this, i, j]() constexpr -> int {
            return -energy_model.score_multi(i, j, seq.enc[i], seq.enc[i + 1], seq.enc[j - 1], seq.enc[j], seq_length);
        };
        update_state<mode>(get_score, &state, nullptr, StateType::P, i, j);
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
                // lambda function
                const auto get_score = [this, p, q, i, j, nucj1]() constexpr -> int {
                    return -energy_model.score_single_loop(p, q, i, j, seq.enc[p], seq.enc[p + 1], seq.enc[q - 1],
                                                           seq.enc[q], seq.enc[i - 1], seq.enc[i], seq.enc[j], nucj1);
                };
                // current shape is: p...i (pair) j...q
                update_state<mode>(get_score, &state, nullptr, StateType::P, p, q);
                q = next_pair[seq.enc[p]][q];
            }
        }
        // 2. P -> M
        if (i > 0 && j < seq_length - 1) {
            const auto get_score = [this, i, j, nucj1]() constexpr -> int {
                return -energy_model.score_M1(i, j, j, seq.enc[i - 1], seq.enc[i], seq.enc[j], nucj1, seq_length);
            };
            update_state<mode>(get_score, &state, nullptr, StateType::M, i, j);
        }
        // 3. M + P -> M2
        const int h = i - 1;
        if (h > 0 && !bestM[h].empty()) {
            const auto get_score = [this, i, j, nucj1]() constexpr -> int {
                return -energy_model.score_M1(i, j, j, seq.enc[i - 1], seq.enc[i], seq.enc[j], nucj1, seq_length);
            };
            for (auto& m_item : bestM[h]) {
                const int g = m_item.first;
                State& m_state = m_item.second;
                update_state<mode>(get_score, &m_state, &state, StateType::M2, g, j);
            }
        }

        // 4. C + P -> C
        const auto get_score = [this, h, i, j, nucj1]() constexpr -> int {
            return -energy_model.score_external_paired(i, j, i > 0 ? seq.enc[h] : -1, seq.enc[i], seq.enc[j], nucj1,
                                                       seq_length);
        };
        update_state<mode>(get_score, &bestC[h], &state, StateType::C, 0, j);
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
                const auto get_score = [this, p, i, j, q]() constexpr -> int {
                    return -(energy_model.score_multi_unpaired(p, i - 1) + energy_model.score_multi_unpaired(j, q - 1));
                };
                update_state<mode>(get_score, &state, nullptr, StateType::Multi, p, q);
            }
        }
        // 2. M2 -> M
        const auto get_score = []() constexpr -> int { return 0; };
        update_state<mode>(get_score, &state, nullptr, StateType::M, i, j);
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
            const auto get_score = [this, j]() constexpr -> int {
                return -energy_model.score_multi_unpaired(j, j + 1);
            };
            update_state<mode>(get_score, &state, nullptr, StateType::M, i, j + 1);
        }
    }
    return 0;
}

template <typename T>
template <Mode mode>
void LinearPartitionInterface<T>::beamstep_C(const int j) {
    // C -> C + U
    if (j < seq_length - 1) {
        const auto get_score = [this, j]() constexpr -> int {
            return -energy_model.score_external_unpaired(j + 1, j + 1);
        };
        update_state<mode>(get_score, &bestC[j], nullptr, StateType::C, 0, j + 1);
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
