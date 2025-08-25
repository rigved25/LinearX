// src/partition/backward.cpp
#include <linearx/partition/linear_partition.hpp>

using namespace linearx::math;
using namespace linearx::utils::io;
using namespace linearx::constants::energy;

using namespace std;

inline __attribute__((always_inline)) void LinearPartition::update_best_trace(const HEdge &new_hedge,
                                                                              const TraceInfo &new_trace) {
    const double best_value = best_hedge.weight + (best_hedge.left ? best_hedge.left->alpha : 0) +
                              (best_hedge.right ? best_hedge.right->alpha : 0);
    const double new_value = new_hedge.weight + (new_hedge.left ? new_hedge.left->alpha : 0) +
                             (new_hedge.right ? new_hedge.right->alpha : 0);
    if (new_value >= best_value) {
        best_hedge = new_hedge;
        best_trace = new_trace;
    }
}

void LinearPartition::mfe_backtrack(const int i, const int j, const StateType type, Structure &structure) {
    if (i >= j) return;

    switch (type) {
        case H:
            return;
        case Multi:
            get_incoming_hedges_Multi<Mode::BEST>(i, j);
            break;
        case P:
            structure.addPair(i, j);
            get_incoming_hedges_P<Mode::BEST>(i, j);
            break;
        case M2:
            get_incoming_hedges_M2<Mode::BEST>(i, j);
            break;
        case M:
            get_incoming_hedges_M<Mode::BEST>(i, j);
            break;
        case C:
            get_incoming_hedges_C<Mode::BEST>(j);
            break;
    }

    mfe_backtrack(best_trace.i, best_trace.t, best_trace.type_left, structure);
    mfe_backtrack(best_trace.t + 1, best_trace.j, best_trace.type_right, structure);
}

PartitionOutsideLog LinearPartition::compute_outside(double deviation_threshold, const bool verbose_output) {
    const double global_threshold = bestC[seq_length - 1].alpha - deviation_threshold;
    unsigned long total_nodes = 0, nodes_visited = 0;
    unsigned long edges_saved = 0, edges_pruned = 0;
    incoming_hedges.reserve(4 * bestP[seq_length - 1].size());
    auto process_beam = [&](const int j, unordered_map<int, State> &beam, const StateType type) {
        for (auto &item : beam) {
            const int i = item.first;
            State &state = item.second;
            // if (state.beta > -deviation_threshold) {    // Major Bug Here (fixed!)
            if (state.beta > LOG_OF_ZERO && state.alpha + state.beta > global_threshold) {
                const double edge_threshold = global_threshold - state.beta;
                const auto local_edges_info = backward_update(i, j, state, type, edge_threshold);
                edges_saved += local_edges_info.first;
                edges_pruned += local_edges_info.second;
                nodes_visited += 1;
            }
            total_nodes += 1;
        }
    };
    if (verbose_output) {
        fprintf(stderr, "[LinearPartition] Running Outside Algorithm (ID: %d | Length: %zu | Name: %s):\n", seq.id,
                seq.length(), seq.name.c_str());
    }
    const auto start_time = chrono::high_resolution_clock::now();
    for (int j = seq_length - 1; j >= 0; --j) {
        // reverse topological order: C->M->M2->P->Multi
        if (verbose_output) {
            showProgressBar(seq_length - 1 - j, seq_length - 1);
        }
        // if (bestC[j].beta > -deviation_threshold) {    // Major Bug Here (fixed!)
        if (bestC[j].alpha + bestC[j].beta > global_threshold) {
            const double edge_threshold = global_threshold - bestC[j].beta;
            backward_update(0, j, bestC[j], StateType::C, edge_threshold);
        }
        process_beam(j, bestM[j], StateType::M);
        process_beam(j, bestM2[j], StateType::M2);
        process_beam(j, bestP[j], StateType::P);
        process_beam(j, bestMulti[j], StateType::Multi);
    }
    const auto end_time = chrono::high_resolution_clock::now();
    const double execution_time = chrono::duration_cast<chrono::milliseconds>(end_time - start_time).count();
    const double effective_beam_size = double(nodes_visited) / (4 * seq_length);  // 4 beams (M, M2, P, Multi) per j
    if (verbose_output) {
        fprintf(stderr, "  - Execution Time: %.2f ms (%.2f%% of inside time)\n", execution_time,
                100.0 * execution_time / max(_last_inside_exec_time, 1.0));
        fprintf(stderr, "  - Visited Edges: %lu (saved) + %lu (pruned)\n", edges_saved, edges_pruned);
        fprintf(stderr, "  - Visited Nodes (%.2f%%): %lu (visited) / %lu (total)\n",
                100.0 * nodes_visited / total_nodes, nodes_visited, total_nodes);
        fprintf(stderr, "  - Effective Beam Size: %.2f\n", effective_beam_size);
        fprintf(stderr, "  - Alpha(C(n)): %.5f | Beta(C(0)): %.5f\n", bestC[seq_length - 1].alpha, bestC[-1].beta);
        fprintf(stderr, "\n");
    }
    return PartitionOutsideLog{bestC[-1].beta,      execution_time, deviation_threshold,
                               effective_beam_size, nodes_visited,  total_nodes - nodes_visited,
                               edges_saved,         edges_pruned};
}

inline __attribute__((always_inline)) void LinearPartition::get_incoming_edges_state(const int i, const int j,
                                                                                     const StateType type) {
    switch (type) {
        case H:
            break;
        case Multi:
            get_incoming_hedges_Multi<Mode::PARTITION>(i, j);
            break;
        case P:
            get_incoming_hedges_P<Mode::PARTITION>(i, j);
            break;
        case M2:
            get_incoming_hedges_M2<Mode::PARTITION>(i, j);
            break;
        case M:
            get_incoming_hedges_M<Mode::PARTITION>(i, j);
            break;
        case C:
            get_incoming_hedges_C<Mode::PARTITION>(j);
            break;
    }
    for (auto &hedge : incoming_hedges) {
        hedge.weight *= INV_KT;
    }
}

pair<unsigned long, unsigned long> LinearPartition::backward_update(const int i, const int j, State &state,
                                                                    const StateType type, const double edge_threshold) {
    get_incoming_edges_state(i, j, type);
    if (incoming_hedges.empty()) {
        return make_pair(0, 0);
    }

    saved_hedges.clear();
    saved_hedges.reserve(incoming_hedges.size());
    HEdge *best_hedge = nullptr;

    double best_inside = xlog(0.0);
    double saved_inside = xlog(0.0);

    unsigned long num_local_edges_pruned = 0;
    unsigned long num_local_edges_saved = 0;

    for (auto &hedge : incoming_hedges) {
        double edge_inside = hedge.weight + hedge.left->alpha + (hedge.right ? hedge.right->alpha : 0);
        if (edge_inside > edge_threshold) {  // keep the edge
            Fast_LogPlusEquals(saved_inside, edge_inside);
            saved_hedges.push_back(&hedge);
        } else {  // prune the edge
            num_local_edges_pruned++;
            if (saved_hedges.empty() && edge_inside > best_inside) {
                best_inside = edge_inside;
                best_hedge = &hedge;
            }
        }
    }

    double delta;  // scaling factor to compensate for edge pruning
    if (!saved_hedges.empty()) {
        delta = state.alpha - saved_inside;
    } else {
        delta = state.alpha - best_inside;
        saved_hedges.push_back(best_hedge);
        num_local_edges_pruned -= 1;  // one more edge recovered
    }

    for (auto &hedge : saved_hedges) {
        State *left = hedge->left, *right = hedge->right;
        if (!right) {
            Fast_LogPlusEquals(left->beta, state.beta + hedge->weight + delta);
        } else {
            Fast_LogPlusEquals(left->beta, state.beta + right->alpha + hedge->weight + delta);
            Fast_LogPlusEquals(right->beta, state.beta + left->alpha + hedge->weight + delta);
        }
    }

    num_local_edges_saved += saved_hedges.size();
    return make_pair(num_local_edges_saved, num_local_edges_pruned);
}

template <Mode mode>
void LinearPartition::get_incoming_hedges_C(const int j) {
    if constexpr (mode == Mode::BEST) {
        best_hedge.reset();
        best_trace.reset();
    } else {
        incoming_hedges.clear();
    }

    // C = C + U
    int new_score = 0;
    if constexpr (mode == Mode::BEST) {
        const HEdge new_hedge(new_score, &bestC[j - 1], nullptr);
        const TraceInfo new_trace(0, j - 1, j - 1, StateType::C, StateType::C);
        update_best_trace(new_hedge, new_trace);
    } else {
        incoming_hedges.emplace_back(new_score, &bestC[j - 1], nullptr);
    }

    // C = C + P
    for (auto &item : bestP[j]) {
        const int i = item.first;
        State &state = item.second;
        new_score = -energy_model.score_external_paired(i, j, (i > 0 ? seq.enc[i - 1] : -1), seq.enc[i], seq.enc[j],
                                                        (j + 1 < seq_length ? seq.enc[j + 1] : -1), seq_length);
        if constexpr (mode == Mode::BEST) {
            const HEdge new_hedge(new_score, &bestC[i - 1], &state);
            const TraceInfo new_trace(0, j, i - 1, StateType::C, StateType::P);
            update_best_trace(new_hedge, new_trace);
        } else {
            incoming_hedges.emplace_back(new_score, &bestC[i - 1], &state);
        }
    }
}

template <Mode mode>
void LinearPartition::get_incoming_hedges_P(const int i, const int j) {
    if constexpr (mode == Mode::BEST) {
        best_hedge.reset();
        best_trace.reset();
    } else {
        incoming_hedges.clear();
    }

    // P = H
    int new_score = 0;
    auto &mp1 = bestH[j];
    auto itr = mp1.find(i);
    if (itr != mp1.end()) {
        if constexpr (mode == Mode::BEST) {
            const HEdge new_hedge(new_score, &itr->second, nullptr);
            const TraceInfo new_trace(i, j, j, StateType::H, StateType::H);
            update_best_trace(new_hedge, new_trace);
        } else {
            incoming_hedges.emplace_back(new_score, &itr->second, nullptr);
        }
    }

    // P = P (scan left & jump right)
    for (int p = i + 1; (p < j - 1) && (p - i <= MAXLOOPSIZE); ++p) {
        int q = prev_pair[seq.enc[p]][j];
        while ((q != -1) && (p < q) && ((p - i) + (j - q) - 2 <= MAXLOOPSIZE)) {
            auto &mp2 = bestP[q];
            itr = mp2.find(p);
            if (itr != mp2.end()) {
                // current shape is: i...p (pair) q...j
                new_score =
                    -energy_model.score_single_loop(i, j, p, q, seq.enc[i], seq.enc[i + 1], seq.enc[j - 1], seq.enc[j],
                                                    seq.enc[p - 1], seq.enc[p], seq.enc[q], seq.enc[q + 1]);
                if constexpr (mode == Mode::BEST) {
                    const HEdge new_hedge(new_score, &itr->second, nullptr);
                    const TraceInfo new_trace(p, q, q, StateType::P, StateType::P);
                    update_best_trace(new_hedge, new_trace);
                } else {
                    incoming_hedges.emplace_back(new_score, &itr->second, nullptr);
                }
            }
            q = prev_pair[seq.enc[p]][q];
        }
    }

    // P = Multi
    auto &mp3 = bestMulti[j];
    itr = mp3.find(i);
    if (itr != mp3.end()) {
        new_score = -energy_model.score_multi(i, j, seq.enc[i], seq.enc[i + 1], seq.enc[j - 1], seq.enc[j], seq_length);
        if constexpr (mode == Mode::BEST) {
            const HEdge new_hedge(new_score, &itr->second, nullptr);
            const TraceInfo new_trace(i, j, j, StateType::Multi, StateType::Multi);
            update_best_trace(new_hedge, new_trace);
        } else {
            incoming_hedges.emplace_back(new_score, &itr->second, nullptr);
        }
    }
}

template <Mode mode>
void LinearPartition::get_incoming_hedges_M(const int i, const int j) {
    if constexpr (mode == Mode::BEST) {
        best_hedge.reset();
        best_trace.reset();
    } else {
        incoming_hedges.clear();
    }

    // M = M + U
    if (j > 0) {
        auto &mp1 = bestM[j - 1];
        auto itr = mp1.find(i);
        if (itr != mp1.end()) {
            int new_score = -energy_model.score_multi_unpaired(j - 1, j);
            if constexpr (mode == Mode::BEST) {
                const HEdge new_hedge(new_score, &itr->second, nullptr);
                const TraceInfo new_trace(i, j - 1, j - 1, StateType::M, StateType::M);
                update_best_trace(new_hedge, new_trace);
            } else {
                incoming_hedges.emplace_back(new_score, &itr->second, nullptr);
            }
        }
    }

    // M = P
    auto &mp2 = bestP[j];
    auto itr = mp2.find(i);
    if (itr != mp2.end()) {
        int new_score = -energy_model.score_M1(i, j, j, seq.enc[i - 1], seq.enc[i], seq.enc[j],
                                               (j + 1 < seq_length ? seq.enc[j + 1] : -1), seq_length);
        if constexpr (mode == Mode::BEST) {
            const HEdge new_hedge(new_score, &itr->second, nullptr);
            const TraceInfo new_trace(i, j, j, StateType::P, StateType::P);
            update_best_trace(new_hedge, new_trace);
        } else {
            incoming_hedges.emplace_back(new_score, &itr->second, nullptr);
        }
    }

    // M = M2
    auto &mp3 = bestM2[j];
    itr = mp3.find(i);
    if (itr != mp3.end()) {
        int new_score = 0;
        if constexpr (mode == Mode::BEST) {
            const HEdge new_hedge(new_score, &itr->second, nullptr);
            const TraceInfo new_trace(i, j, j, StateType::M2, StateType::M2);
            update_best_trace(new_hedge, new_trace);
        } else {
            incoming_hedges.emplace_back(new_score, &itr->second, nullptr);
        }
    }
}

template <Mode mode>
void LinearPartition::get_incoming_hedges_M2(const int i, const int j) {
    if constexpr (mode == Mode::BEST) {
        best_hedge.reset();
        best_trace.reset();
    } else {
        incoming_hedges.clear();
    }

    // [TODO] sort P?

    // M2 = M + P
    for (auto &item : bestP[j]) {
        const int t = item.first;
        State &state = item.second;
        if (t > i) {
            auto &mp = bestM[t - 1];
            auto itr = mp.find(i);
            if (itr != mp.end()) {
                int new_score = -energy_model.score_M1(t, j, j, seq.enc[t - 1], seq.enc[t], seq.enc[j],
                                                       (j + 1 < seq_length ? seq.enc[j + 1] : -1), seq_length);
                if constexpr (mode == Mode::BEST) {
                    const HEdge new_hedge(new_score, &itr->second, &state);
                    const TraceInfo new_trace(i, t - 1, t - 1, StateType::M, StateType::P);
                    update_best_trace(new_hedge, new_trace);
                } else {
                    incoming_hedges.emplace_back(new_score, &itr->second, &state);
                }
            }
        }
    }
}

template <Mode mode>
void LinearPartition::get_incoming_hedges_Multi(const int i, const int j) {
    if constexpr (mode == Mode::BEST) {
        best_hedge.reset();
        best_trace.reset();
    } else {
        incoming_hedges.clear();
    }

    int jprev = prev_pair[seq.enc[i]][j];
    if (jprev == -1) {
        return;
    }

    // Multi = Multi (jump right)
    auto &mp1 = bestMulti[jprev];
    auto itr = mp1.find(i);
    if (itr != mp1.end()) {
        int new_score = -energy_model.score_multi_unpaired(jprev, j);
        if constexpr (mode == Mode::BEST) {
            const HEdge new_hedge(new_score, &itr->second, nullptr);
            const TraceInfo new_trace(i, jprev, jprev, StateType::Multi, StateType::Multi);
            update_best_trace(new_hedge, new_trace);
        } else {
            incoming_hedges.emplace_back(new_score, &itr->second, nullptr);
        }
    }

    // Multi = M2 (scan left & jump right)
    for (int q = j - 1; q >= jprev; --q) {
        auto &mp2 = bestM2[q];
        for (int p = i + 1; (p < q) && (p - i <= MAXLOOPSIZE); ++p) {
            auto itr = mp2.find(p);
            if (itr != mp2.end()) {
                // the current shape is i..p M2 q..j
                int new_score =
                    -(energy_model.score_multi_unpaired(i, p - 1) + energy_model.score_multi_unpaired(q, j - 1));
                if constexpr (mode == Mode::BEST) {
                    const HEdge new_hedge(new_score, &itr->second, nullptr);
                    const TraceInfo new_trace(p, q, q, StateType::M2, StateType::M2);
                    update_best_trace(new_hedge, new_trace);
                } else {
                    incoming_hedges.emplace_back(new_score, &itr->second, nullptr);
                }
            }
        }
    }
}

template void LinearPartition::get_incoming_hedges_C<Mode::BEST>(int);
template void LinearPartition::get_incoming_hedges_C<Mode::PARTITION>(int);

template void LinearPartition::get_incoming_hedges_M<Mode::BEST>(int, int);
template void LinearPartition::get_incoming_hedges_M<Mode::PARTITION>(int, int);

template void LinearPartition::get_incoming_hedges_M2<Mode::BEST>(int, int);
template void LinearPartition::get_incoming_hedges_M2<Mode::PARTITION>(int, int);

template void LinearPartition::get_incoming_hedges_P<Mode::BEST>(int, int);
template void LinearPartition::get_incoming_hedges_P<Mode::PARTITION>(int, int);

template void LinearPartition::get_incoming_hedges_Multi<Mode::BEST>(int, int);
template void LinearPartition::get_incoming_hedges_Multi<Mode::PARTITION>(int, int);
