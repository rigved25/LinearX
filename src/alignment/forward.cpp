// src/alignment/forward.cpp
#include <linearx/alignment/linear_align.hpp>

using namespace std;
using namespace linearx::utils;

template <Mode mode>
AlignmentLog LinearAlignment::compute_inside(const unsigned beam_size, bool verbose_output) {
    run_beam_size_ = beam_size;
    const auto start_time = chrono::high_resolution_clock::now();
    if (verbose_output) {
        cerr << "[LinearAlignment] Running Inside Algorithm:" << endl;
    }
    unsigned long nodes_pruned = 0;
    for (int s = 0; s <= seq_len_sum; ++s) {
        if (verbose_output) {
            linearx::utils::io::showProgressBar(s, seq_len_sum);
        }
        for (const HStateType h : hstate_types) {
            vector<unordered_map<pair<int, int>, HState, PairHash>> &beam = get_beams(h);
            nodes_pruned += beam_prune(beam[s]);
            for (auto &item : beam[s]) {
                const int i = item.first.first;
                const int j = item.first.second;
                HState &state = item.second;

                // INS1
                if (i < seq1.length() && j <= seq2.length()) {
                    if constexpr (mode == Mode::BEST) {
                        const value_type new_score = get_trans_emit_prob(i + 1, j, HStateType::INS1, h);
                        HState &next_state = bestINS1[s + 1][{i + 1, j}];
                        next_state.alpha = max(next_state.alpha, LOG_MUL(state.alpha, new_score));
                    } else {
                        HState *next_state = check_state(HStateType::INS1, i + 1, j);
                        if (next_state) {
                            const value_type new_score = get_trans_emit_prob(i + 1, j, HStateType::INS1, h);
                            AlnEdge(&state, new_score).update_state_alpha(*next_state);
                        }
                    }
                }

                // INS2
                if (i <= seq1.length() && j < seq2.length()) {
                    if constexpr (mode == Mode::BEST) {
                        const value_type new_score = get_trans_emit_prob(i, j + 1, HStateType::INS2, h);
                        HState &next_state = bestINS2[s + 1][{i, j + 1}];
                        next_state.alpha = max(next_state.alpha, LOG_MUL(state.alpha, new_score));
                    } else {
                        HState *next_state = check_state(HStateType::INS2, i, j + 1);
                        if (next_state) {
                            const value_type new_score = get_trans_emit_prob(i, j + 1, HStateType::INS2, h);
                            AlnEdge(&state, new_score).update_state_alpha(*next_state);
                        }
                    }
                }

                // ALN
                const bool end_check = (i == seq1.length() && j == seq2.length());
                if ((i < seq1.length() && j < seq2.length() || end_check)) {
                    if constexpr (mode == Mode::BEST) {
                        value_type new_score = get_trans_emit_prob(i + 1, j + 1, HStateType::ALN, h);
                        new_score = LOG_MUL(new_score, get_match_score(i, j));
                        HState &next_state = bestALN[s + 2][{i + 1, j + 1}];
                        next_state.alpha = max(next_state.alpha, LOG_MUL(state.alpha, new_score));
                    } else {
                        HState *next_state = check_state(HStateType::ALN, i + 1, j + 1);
                        if (next_state) {
                            value_type new_score = get_trans_emit_prob(i + 1, j + 1, HStateType::ALN, h);
                            new_score = LOG_MUL(new_score, get_match_score(i, j));
                            AlnEdge(&state, new_score).update_state_alpha(*next_state);
                        }
                    }
                }
            }
        }
    }
    // update/print time stats
    const auto end_time = chrono::high_resolution_clock::now();
    const value_type execution_time = chrono::duration_cast<chrono::milliseconds>(end_time - start_time).count();
    _last_inside_exec_time = execution_time;
    if (verbose_output) {
        fprintf(stderr, "  - Execution Time: %.3f ms\n", execution_time);
        fprintf(stderr, "  - Nodes Pruned: %lu\n", nodes_pruned);
        fprintf(stderr, "  - Score: %.4f\n", bestALN[seq_len_sum + 2][{seq1.length() + 1, seq2.length() + 1}].alpha);
    }
    return AlignmentLog{-1.0,
                        bestALN[seq_len_sum + 2][{seq1.length() + 1, seq2.length() + 1}].alpha,
                        linearx::math::LOG_ZERO,
                        execution_time,
                        -linearx::math::LOG_ZERO,
                        float(beam_size),
                        "Inside",
                        0,
                        nodes_pruned,
                        0,
                        0};
}

// explicit template instantiation
template AlignmentLog LinearAlignment::compute_inside<Mode::BEST>(const unsigned beam_size, bool verbose_output);
template AlignmentLog LinearAlignment::compute_inside<Mode::PARTITION_INSIDE>(const unsigned beam_size,
                                                                              bool verbose_output);
