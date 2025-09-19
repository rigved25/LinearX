// src/alignment/backward.cpp
#include <iomanip>
#include <linearx/alignment/config.hpp>
#include <linearx/alignment/linear_align.hpp>

using namespace std;
using namespace linearx::utils;
using namespace linearx::math;

template <typename T>
void LinearAlignmentInterface<T>::update_best_trace(const AlnEdge& new_edge, const HStateType& new_trace) {
    const value_type new_score = new_edge.weight + (new_edge.prev ? new_edge.prev->alpha : 0);
    const value_type best_score = best_edge.weight + (best_edge.prev ? best_edge.prev->alpha : 0);
    if (new_score >= best_score) {
        best_edge = new_edge;
        best_trace = new_trace;
    }
}

template <typename T>
MultiSeq LinearAlignmentInterface<T>::get_alignment() {
    int i = seq1.length();
    int j = seq2.length();
    get_incoming_edges<Mode::BEST>(i + 1, j + 1, HStateType::ALN);
    HStateType h = best_trace;

    Sequence s1("", seq1.name, seq1.id);
    Sequence s2("", seq2.name, seq2.id);
    while (i > 0 || j > 0) {
        get_incoming_edges<Mode::BEST>(i, j, h);
        HStateType h_prev = best_trace;
        switch (h) {
            case HStateType::ALN:
                i -= 1;
                j -= 1;
                s1.add_nuc(seq1[i]);
                s2.add_nuc(seq2[j]);
                break;

            case HStateType::INS1:
                i -= 1;
                s1.add_nuc(seq1[i]);
                s2.add_nuc('-');
                break;

            case HStateType::INS2:
                j -= 1;
                s1.add_nuc('-');
                s2.add_nuc(seq2[j]);
                break;
        }
        h = h_prev;
    }

    s1.reverse();
    s2.reverse();

    MultiSeq alignment;
    alignment.add_sequence(s1);
    alignment.add_sequence(s2);
    return alignment;
}

template <typename T>
AlignmentLog LinearAlignmentInterface<T>::compute_outside(bool use_lazy_outside, value_type deviation_threshold,
                                                          bool verbose_output) {
    if (!use_lazy_outside) {
        return run_normal_outside(verbose_output);
    }

    const value_type global_threshold =
        bestALN[seq_len_sum + 2][{seq1.length() + 1, seq2.length() + 1}].alpha - deviation_threshold;
    unsigned long total_states = 0, nodes_visited = 0;
    unsigned long edges_saved = 0, edges_pruned = 0;
    auto process_beam = [&](unordered_map<pair<int, int>, HState, PairHash>& beam, const HStateType type) {
        for (auto& item : beam) {
            const int i = item.first.first;
            const int j = item.first.second;
            HState& state = item.second;
            if (state.alpha + state.beta > global_threshold) {
                value_type edge_threshold = global_threshold - state.beta;
                pair<int, int> local_edges_info = backward_update(i, j, state, type, edge_threshold);
                edges_saved += local_edges_info.first;
                edges_pruned += local_edges_info.second;
                nodes_visited += 1;
            }
            total_states += 1;
        }
    };

    if (verbose_output) {
        cerr << "[LinearAlignment] Running Outside Algorithm:" << endl;
    }
    auto start_time = chrono::high_resolution_clock::now();
    process_beam(bestALN[seq_len_sum + 2], HStateType::ALN);
    for (int s = seq_len_sum; s > 0; --s) {
        if (verbose_output) {
            linearx::utils::io::showProgressBar(seq_len_sum - s, seq_len_sum - 1);
        }
        // reverse topological order: ALN->INS2->INS1
        process_beam(bestALN[s], HStateType::ALN);
        process_beam(bestINS2[s], HStateType::INS2);
        process_beam(bestINS1[s], HStateType::INS1);
    }
    auto end_time = chrono::high_resolution_clock::now();
    const value_type execution_time = chrono::duration_cast<chrono::milliseconds>(end_time - start_time).count();
    const float effective_beam_size = float(nodes_visited) / (3 * seq_len_sum);
    if (verbose_output) {
        fprintf(stderr, "  - Execution Time: %.2f ms (%.2f%% of inside time)\n", execution_time,
                100.0 * execution_time / max(_last_inside_exec_time, 1.0));
        fprintf(stderr, "  - Visited Edges: %lu (saved) + %lu (pruned)\n", edges_saved, edges_pruned);
        fprintf(stderr, "  - Visited Nodes (%.2f%%): %lu (visited) / %lu (total)\n",
                100.0 * nodes_visited / total_states, nodes_visited, total_states);
        fprintf(stderr, "  - Effective Beam Size: %.2f\n", effective_beam_size);
        fprintf(stderr, "  - Alpha(ALN(n)): %.5f | Beta(ALN(0)): %.5f\n",
                bestALN[seq_len_sum + 2][{seq1.length() + 1, seq2.length() + 1}].alpha, bestALN[0][{0, 0}].beta);
    }
    // clean up
    incoming_edges.clear();
    incoming_edges.shrink_to_fit();
    saved_edges.clear();
    saved_edges.shrink_to_fit();
    return AlignmentLog{-1,
                        bestALN[seq_len_sum + 2][{seq1.length() + 1, seq2.length() + 1}].alpha,
                        bestALN[0][{0, 0}].beta,
                        _last_inside_exec_time,
                        execution_time,
                        effective_beam_size,
                        "Lazy Outside",
                        nodes_visited,
                        total_states - nodes_visited,
                        edges_saved,
                        edges_pruned};
}

template <typename T>
pair<unsigned long, unsigned long> LinearAlignmentInterface<T>::backward_update(const int i, const int j,
                                                                                const HState& state,
                                                                                const HStateType type,
                                                                                const value_type edge_threshold) {
    if ((i == 0 || j == 0) && type == HStateType::ALN) {
        return make_pair(0, 0);
    }

    get_incoming_edges<Mode::PARTITION_OUTSIDE>(i, j, type);
    if (incoming_edges.empty()) {
        return make_pair(0, 0);
    }

    saved_edges.clear();
    saved_edges.reserve(incoming_edges.size());
    AlnEdge* best_edge = nullptr;

    value_type best_inside = LOG_ZERO;
    value_type saved_inside = LOG_ZERO;

    unsigned long edges_pruned = 0;
    unsigned long edges_saved = 0;

    for (auto& edge : incoming_edges) {
        value_type edge_inside = LOG_MUL(edge.prev->alpha, edge.weight);  // LOG_MUL(a, b) -> a + b
        if (edge_inside > edge_threshold) {                               // keep the edge
            saved_inside = LOG_SUM(saved_inside, edge_inside);
            saved_edges.push_back(&edge);
        } else {  // prune the edge
            edges_pruned++;
            if (saved_edges.empty() && edge_inside > best_inside) {
                best_inside = edge_inside;
                best_edge = &edge;
            }
        }
    }

    value_type delta;  // scaling factor to compensate for edge pruning
    if (!saved_edges.empty()) {
        delta = LOG_DIV(state.alpha, saved_inside);  // LOG_DIV(a, b) -> a - b
    } else {
        delta = LOG_DIV(state.alpha, best_inside);  // state.alpha - best_inside
        saved_edges.push_back(best_edge);
        edges_pruned -= 1;  // one more edge recovered
    }

    for (auto& edge : saved_edges) {
        edge->prev->beta = LOG_SUM(edge->prev->beta, state.beta + edge->weight + delta);
    }

    edges_saved += saved_edges.size();
    return make_pair(edges_saved, edges_pruned);
}

template <typename T>
template <Mode mode>
void LinearAlignmentInterface<T>::get_incoming_edges(const int i, const int j, const HStateType type) {
    if constexpr (mode == Mode::BEST) {
        best_edge.reset();
    } else {
        incoming_edges.clear();
    }

    int p = i, q = j;
    switch (type) {
        case HStateType::ALN:
            p = i - 1;
            q = j - 1;
            break;
        case HStateType::INS1:
            p = i - 1;
            break;
        case HStateType::INS2:
            q = j - 1;
            break;
    }

    for (const HStateType h_prev : hstate_types) {
        vector<unordered_map<pair<int, int>, HState, PairHash>>& beam = get_beams(h_prev);
        auto& mp = beam[p + q];
        auto it = mp.find({p, q});
        if (it != mp.end()) {
            value_type edge_weight = get_trans_emit_prob(i, j, type, h_prev);
            if constexpr (mode == Mode::BEST) {
                if (type == HStateType::ALN) {
                    edge_weight = LOG_MUL(edge_weight, get_match_score(i - 1, j - 1));
                }
                AlnEdge new_edge(&(it->second), edge_weight);
                update_best_trace(new_edge, h_prev);
            } else {
                incoming_edges.emplace_back(&(it->second), edge_weight);
            }
        }
    }
}

template <typename T>
AlignmentLog LinearAlignmentInterface<T>::run_normal_outside(bool verbose_output) {
    const auto start_time = chrono::high_resolution_clock::now();
    if (verbose_output) {
        cerr << "[LinearAlignment] Running Outside Algorithm:" << endl;
    }
    unsigned long nodes_visited = 0;
    unsigned long edges_visited = 0;
    for (int s = seq_len_sum; s >= 0; --s) {
        if (verbose_output) {
            linearx::utils::io::showProgressBar(seq_len_sum - s, seq_len_sum);
        }
        for (const HStateType h : hstate_types) {
            vector<unordered_map<pair<int, int>, HState, PairHash>>& beam = get_beams(h);
            for (const auto& item : beam[s]) {
                const int i = item.first.first;
                const int j = item.first.second;
                HState& state = beam[s][{i, j}];
                nodes_visited += 1;
                // INS1
                if (i < seq1.length() && j <= seq2.length()) {
                    auto it = bestINS1[s + 1].find({i + 1, j});
                    if ((it != bestINS1[s + 1].end())) {
                        edges_visited += 1;
                        const value_type new_score = get_trans_emit_prob(i + 1, j, HStateType::INS1, h);
                        AlnEdge(&state, new_score).update_state_beta(it->second);
                    }
                }
                // INS2
                if (i <= seq1.length() && j < seq2.length()) {
                    auto it = bestINS2[s + 1].find({i, j + 1});
                    if ((it != bestINS2[s + 1].end())) {
                        edges_visited += 1;
                        const value_type new_score = get_trans_emit_prob(i, j + 1, HStateType::INS2, h);
                        AlnEdge(&state, new_score).update_state_beta(it->second);
                    }
                }
                // ALN
                const bool end_check = (i == seq1.length() && j == seq2.length());
                if ((i < seq1.length() && j < seq2.length()) || end_check) {
                    auto it = bestALN[s + 2].find({i + 1, j + 1});
                    if (it != bestALN[s + 2].end()) {
                        edges_visited += 1;
                        value_type new_score = get_trans_emit_prob(i + 1, j + 1, HStateType::ALN, h);
                        new_score = LOG_MUL(new_score, get_match_score(i, j));
                        AlnEdge(&state, new_score).update_state_beta(it->second);
                    }
                }
            }
        }
    }
    // update/print time stats
    auto end_time = chrono::high_resolution_clock::now();
    const value_type execution_time = chrono::duration_cast<chrono::milliseconds>(end_time - start_time).count();
    const float effective_beam_size = float(nodes_visited) / (3 * seq_len_sum);
    if (verbose_output) {
        fprintf(stderr, "  - Execution Time: %.2f ms (%.2f%% of inside time)\n", execution_time,
                100.0 * execution_time / max(_last_inside_exec_time, 1.0));
        fprintf(stderr, "  - Visited Edges: %lu (saved) + %lu (pruned)\n", edges_visited, (unsigned long)0);
        fprintf(stderr, "  - Visited Nodes (%.2f%%): %lu (visited) / %lu (total)\n",
                100.0 * nodes_visited / nodes_visited, nodes_visited, nodes_visited);
        fprintf(stderr, "  - Effective Beam Size: %.2f\n", effective_beam_size);
        fprintf(stderr, "  - Alpha(ALN(n)): %.5f | Beta(ALN(0)): %.5f\n",
                bestALN[seq_len_sum + 2][{seq1.length() + 1, seq2.length() + 1}].alpha, bestALN[0][{0, 0}].beta);
    }
    return AlignmentLog{-1,
                        bestALN[seq_len_sum + 2][{seq1.length() + 1, seq2.length() + 1}].alpha,
                        bestALN[0][{0, 0}].beta,
                        _last_inside_exec_time,
                        execution_time,
                        0,
                        "Regular Outside",
                        nodes_visited,
                        0,
                        edges_visited,
                        0};
}

// instantiate templates for LinearAlignmentInterface with desired types
#define DECLARE_FUNCS(TYPE)                                                                             \
    template class LinearAlignmentInterface<TYPE>;                                                      \
    template void LinearAlignmentInterface<TYPE>::get_incoming_edges<Mode::BEST>(int, int, HStateType); \
    template void LinearAlignmentInterface<TYPE>::get_incoming_edges<Mode::PARTITION_OUTSIDE>(int, int, HStateType);

#define X(TYPE) DECLARE_FUNCS(TYPE)
LA_TEMPLATE_TYPES
#undef X
#undef DECLARE_FUNCS
