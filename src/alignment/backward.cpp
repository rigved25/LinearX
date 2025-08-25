// src/alignment/backward.cpp
#include <linearx/alignment/linear_align.hpp>

using namespace std;
using namespace linearx::utils;

inline __attribute__((always_inline)) void LinearAlignment::update_best_trace(const AlnEdge &new_edge,
                                                                              const HStateType &new_trace) {
    const double new_score = new_edge.weight + (new_edge.prev ? new_edge.prev->alpha : 0);
    const double best_score = best_edge.weight + (best_edge.prev ? best_edge.prev->alpha : 0);
    if (new_score >= best_score) {
        best_edge = new_edge;
        best_trace = new_trace;
    }
}

MultiSeq LinearAlignment::get_alignment() {
    int i = seq1.length();
    int j = seq2.length();
    const bool use_match_score = (pm1 != nullptr && pm2 != nullptr);

    get_incoming_edges<Mode::BEST>(i + 1, j + 1, HStateType::ALN, use_match_score);
    HStateType h = best_trace;

    string aln1 = "";
    string aln2 = "";

    while (i > 0 || j > 0) {
        get_incoming_edges<Mode::BEST>(i, j, h, use_match_score);
        HStateType h_prev = best_trace;
        switch (h) {
            case HStateType::ALN:
                i -= 1;
                j -= 1;
                aln1 += to_string(seq1.enc[i]);
                aln2 += to_string(seq2.enc[j]);
                break;

            case HStateType::INS1:
                i -= 1;
                aln1 += to_string(seq1.enc[i]);
                aln2 += "-";
                break;

            case HStateType::INS2:
                j -= 1;
                aln1 += "-";
                aln2 += to_string(seq2.enc[j]);
                break;
        }
        h = h_prev;
    }

    reverse(aln1.begin(), aln1.end());
    reverse(aln2.begin(), aln2.end());

    MultiSeq alignment;
    alignment.add_sequence(Sequence(aln1, seq1.name, seq1.id));
    alignment.add_sequence(Sequence(aln2, seq2.name, seq2.id));
    return alignment;
}

inline __attribute__((always_inline)) void LinearAlignment::update_state_beta(HState &state, const double new_score) {
    state.beta = LOG_SUM(state.beta, new_score);
}

void LinearAlignment::compute_outside(bool use_lazy_outside, double deviation_threshold, bool verbose_output) {
    if (!use_lazy_outside) {
        run_normal_outside(verbose_output);
        return;
    }

    const double global_threshold =
        bestALN[seq_len_sum + 2][{seq1.length() + 1, seq2.length() + 1}].alpha + deviation_threshold;

    unsigned long total_states = 0, states_visited = 0;
    unsigned long edges_saved = 0, edges_pruned = 0;

    auto process_beam = [&](const int s, unordered_map<pair<int, int>, HState, PairHash> &beam, const HStateType type) {
        for (auto &item : beam) {
            const int i = item.first.first;
            const int j = item.first.second;
            HState &state = item.second;
            if (state.alpha + state.beta > global_threshold) {
                double edge_threshold = global_threshold - state.beta;
                pair<int, int> local_edges_info = backward_update(i, j, state, type, edge_threshold);
                edges_saved += local_edges_info.first;
                edges_pruned += local_edges_info.second;
                states_visited += 1;
            }
            total_states += 1;
        }
    };

    if (verbose_output) {
        cerr << "[LinearAlignment] Running Outside Algorithm:" << endl;
    }
    auto start_time = chrono::high_resolution_clock::now();
    process_beam(seq_len_sum + 2, bestALN[seq_len_sum + 2], HStateType::ALN);
    for (int s = seq_len_sum; s > 0; --s) {
        if (verbose_output) {
            linearx::utils::io::showProgressBar(seq_len_sum - s, seq_len_sum - 1);
        }

        // reverse topological order: ALN->INS2->INS1
        process_beam(s, bestALN[s], HStateType::ALN);
        process_beam(s, bestINS2[s], HStateType::INS2);
        process_beam(s, bestINS1[s], HStateType::INS1);
    }
    auto end_time = chrono::high_resolution_clock::now();
    const double execution_time =
        chrono::duration_cast<chrono::microseconds>(end_time - start_time).count() / 1000.0;  // in milliseconds

    if (verbose_output) {
        fprintf(stderr, "  - Execution Time: %.2f ms (%.2f%% of inside time)\n", execution_time,
                100.0 * execution_time / max(_last_inside_exec_time, 1.0));
        fprintf(stderr, "  - Visited Edges: %lu (saved) + %lu (pruned)\n", edges_saved, edges_pruned);
        fprintf(stderr, "  - Visited Nodes (%.2f%%): %lu (visited) / %lu (total)\n\n",
                100.0 * states_visited / total_states, states_visited, total_states);
    }
}

pair<unsigned long, unsigned long> LinearAlignment::backward_update(const int i, const int j, const HState &state,
                                                                    const HStateType type,
                                                                    const double edge_threshold) {
    if ((i == 0 || j == 0) && type == HStateType::ALN) {
        return make_pair(0, 0);
    }

    get_incoming_edges<Mode::PARTITION>(i, j, type, false);
    if (incoming_edges.empty()) {
        return make_pair(0, 0);
    }

    saved_edges.clear();
    saved_edges.reserve(incoming_edges.size());
    AlnEdge *best_edge = nullptr;

    double best_inside = LOG(0.0);
    double saved_inside = LOG(0.0);

    unsigned long edges_pruned = 0;
    unsigned long edges_saved = 0;

    for (auto &edge : incoming_edges) {
        double edge_inside = LOG_MUL(edge.prev->alpha, edge.weight);  // LOG_MUL(a, b) -> a + b
        if (edge_inside > edge_threshold) {                           // keep the edge
            saved_inside = LOG_SUM(saved_inside, edge_inside);        // Fast_LogPlusEquals(saved_inside, edge_inside);
            saved_edges.push_back(&edge);
        } else {  // prune the edge
            edges_pruned++;
            if (saved_edges.empty() && edge_inside > best_inside) {
                best_inside = edge_inside;
                best_edge = &edge;
            }
        }
    }

    double delta;  // scaling factor to compensate for edge pruning
    if (!saved_edges.empty()) {
        delta = LOG_DIV(state.alpha, saved_inside);  // LOG_DIV(a, b) -> a - b
    } else {
        delta = LOG_DIV(state.alpha, best_inside);  // state.alpha - best_inside
        saved_edges.push_back(best_edge);
        edges_pruned -= 1;  // one more edge recovered
    }

    for (auto &edge : saved_edges) {
        edge->prev->beta = LOG_SUM(edge->prev->beta, state.beta + edge->weight + delta);
    }

    edges_saved += saved_edges.size();
    return make_pair(edges_saved, edges_pruned);
}

template <Mode mode>
void LinearAlignment::get_incoming_edges(const int i, const int j, const HStateType type, bool use_match_score) {
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
        vector<unordered_map<pair<int, int>, HState, PairHash>> &beam = get_beam(h_prev);
        auto &mp = beam[p + q];
        auto it = mp.find({p, q});
        if (it != mp.end()) {
            double edge_weight = get_trans_emit_prob(i, j, type, h_prev);
            if (use_match_score && type == HStateType::ALN) {
                double match_score = get_match_score(i - 1, j - 1);
                edge_weight = LOG_MUL(edge_weight, match_score);
            }
            // new_edge.set(&(it->second), edge_weight);
            // new_trace = h_prev;
            // edges_trace_update_helper(incoming_edges, new_edge, best_edge, new_trace, best_trace);
            if constexpr (mode == Mode::BEST) {
                AlnEdge new_edge(&(it->second), edge_weight);
                HStateType new_trace = h_prev;
                update_best_trace(new_edge, new_trace);
            } else {
                incoming_edges.emplace_back(&(it->second), edge_weight);
            }
        }
    }
}

template void LinearAlignment::get_incoming_edges<Mode::BEST>(int, int, HStateType, bool);
template void LinearAlignment::get_incoming_edges<Mode::PARTITION>(int, int, HStateType, bool);

void LinearAlignment::run_normal_outside(bool verbose_output) {
    const auto start_time = chrono::high_resolution_clock::now();
    if (verbose_output) {
        cerr << "[LinearAlignment] Running Outside Algorithm:" << endl;
    }
    bool use_match_score = (pm1 != nullptr && pm2 != nullptr);
    int states_visited = 0;
    for (int s = seq_len_sum; s >= 0; --s) {
        if (verbose_output) {
            linearx::utils::io::showProgressBar(seq_len_sum - s, seq_len_sum);
        }
        for (const HStateType h : hstate_types) {
            vector<unordered_map<pair<int, int>, HState, PairHash>> &beam = get_beam(h);
            for (const auto &item : beam[s]) {
                int i = item.first.first;
                int j = item.first.second;
                HState &state = beam[s][{i, j}];
                states_visited += 1;
                // INS1
                if (i < seq1.length() && j <= seq2.length()) {
                    auto it = bestINS1[s + 1].find({i + 1, j});
                    if ((it != bestINS1[s + 1].end())) {
                        double prob = get_trans_emit_prob(i + 1, j, HStateType::INS1, h);
                        double score = LOG_MUL(prob, it->second.beta);
                        update_state_beta(state, score);
                    }
                }
                // INS2
                if (i <= seq1.length() && j < seq2.length()) {
                    auto it = bestINS2[s + 1].find({i, j + 1});
                    if ((it != bestINS2[s + 1].end())) {
                        double prob = get_trans_emit_prob(i, j + 1, HStateType::INS2, h);
                        double score = LOG_MUL(prob, it->second.beta);
                        update_state_beta(state, score);
                    }
                }
                // ALN
                const bool end_check = (i == seq1.length() && j == seq2.length());
                if ((i < seq1.length() && j < seq2.length()) || end_check) {
                    auto it = bestALN[s + 2].find({i + 1, j + 1});
                    if (it != bestALN[s + 2].end()) {
                        double prob = get_trans_emit_prob(i + 1, j + 1, HStateType::ALN, h);
                        double score = LOG_MUL(prob, it->second.beta);
                        if (use_match_score) {
                            double match_score = get_match_score(i, j);
                            score = LOG_MUL(score, match_score);
                        }
                        update_state_beta(state, score);
                    }
                }
            }
        }
    }
    // update/print time stats
    const auto end_time = chrono::high_resolution_clock::now();
    const double exec_time = chrono::duration_cast<chrono::milliseconds>(end_time - start_time).count();
    if (verbose_output) {
        fprintf(stderr, "  - Execution Time: %.2f ms (%.2f%% of inside time)\n", exec_time,
                100.0 * exec_time / max(_last_inside_exec_time, 1.0));
        fprintf(stderr, "  - Visited Nodes: %d\n\n", states_visited);
    }
}
