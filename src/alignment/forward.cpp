// src/alignment/forward.cpp
#include <linearx/alignment/linear_align.hpp>

using namespace std;
using namespace linearx::utils;
using namespace linearx::constants::math;

unsigned LinearAlignment::beam_prune(unordered_map<pair<int, int>, HState, PairHash> &beamstep, const int beam_size) {
    if (beam_size == 0 || beamstep.size() <= beam_size) {
        return 0;
    }
    unsigned num_pruned = 0;
    beam_scores.clear();
    beam_scores.reserve(beamstep.size());
    for (auto &item : beamstep) {
        const pair<int, int> aij = item.first;
        HState &cand = item.second;
        beam_scores.emplace_back(cand.alpha, aij);
    }
    double threshold = quickselect(beam_scores, 0, beam_scores.size() - 1, beam_scores.size() - beam_size);
    for (auto &p : beam_scores) {
        if (p.first < threshold) {
            beamstep.erase(p.second);
            num_pruned++;
        }
    }
    return num_pruned;
}

inline __attribute__((always_inline)) void LinearAlignment::update_state_alpha(const Mode mode, HState &state,
                                                                               const double new_score) {
    if (mode == Mode::BEST) {
        if (new_score > state.alpha) {
            state.alpha = new_score;
        }
    } else {
        state.alpha = LOG_SUM(state.alpha, new_score);
    }
}

void LinearAlignment::compute_inside(const Mode mode, const unsigned beam_size, bool verbose_output) {
    const auto start_time = chrono::high_resolution_clock::now();
    if (verbose_output) {
        cerr << "[LinearAlignment] Running Inside Algorithm:" << endl;
    }
    unsigned long nodes_pruned = 0;
    const bool use_match_score = (pm1 != nullptr && pm2 != nullptr);
    for (int s = 0; s <= seq_len_sum; ++s) {
        if (verbose_output) {
            linearx::utils::io::showProgressBar(s, seq_len_sum);
        }
        for (const HStateType h : hstate_types) {
            vector<unordered_map<pair<int, int>, HState, PairHash>> &beam = get_beam(h);
            auto &current_beam = beam[s];
            nodes_pruned += beam_prune(current_beam, beam_size);
            for (const auto &item : current_beam) {
                const int i = item.first.first;
                const int j = item.first.second;
                const HState &state = item.second;

                // INS1
                if (i < seq1.length() && j <= seq2.length() && check_state(i + 1, j, HStateType::INS1)) {
                    const double prob = get_trans_emit_prob(i + 1, j, HStateType::INS1, h);
                    const double new_score = LOG_MUL(state.alpha, prob);
                    update_state_alpha(mode, bestINS1[s + 1][{i + 1, j}], new_score);
                }

                // INS2
                if (i <= seq1.length() && j < seq2.length() && check_state(i, j + 1, HStateType::INS2)) {
                    const double prob = get_trans_emit_prob(i, j + 1, HStateType::INS2, h);
                    const double new_score = LOG_MUL(state.alpha, prob);
                    update_state_alpha(mode, bestINS2[s + 1][{i, j + 1}], new_score);
                }

                // ALN
                const bool end_check = (i == seq1.length() && j == seq2.length());
                if ((i < seq1.length() && j < seq2.length() && check_state(i + 1, j + 1, HStateType::ALN)) ||
                    end_check) {
                    const double prob = get_trans_emit_prob(i + 1, j + 1, HStateType::ALN, h);
                    double new_score = LOG_MUL(state.alpha, prob);
                    if (use_match_score) {
                        const double match_score = get_match_score(i, j);
                        new_score = LOG_MUL(new_score, match_score);
                    }
                    update_state_alpha(mode, bestALN[s + 2][{i + 1, j + 1}], new_score);
                }
            }
        }
    }
    // update/print time stats
    const auto end_time = chrono::high_resolution_clock::now();
    const double execution_time = chrono::duration_cast<chrono::milliseconds>(end_time - start_time).count();
    _last_inside_exec_time = execution_time;
    if (verbose_output) {
        fprintf(stderr, "  - Execution Time: %.3f ms\n", execution_time);
        fprintf(stderr, "  - Nodes Pruned: %lu\n", nodes_pruned);
        fprintf(stderr, "  - Score: %.4f\n", bestALN[seq_len_sum + 2][{seq1.length() + 1, seq2.length() + 1}].alpha);
        fprintf(stderr, "\n");
    }
}
