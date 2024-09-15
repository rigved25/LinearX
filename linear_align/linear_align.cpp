#include "./linear_align.hpp"

void LinearAlign::update_state_alpha(HState &state, const double new_score, const HStateType h, const HStateType pre,
                                     const bool best_only = false) {
    if (best_only) {
        if (new_score > state.alpha) {
            state.alpha = new_score;
            state.pre = pre;
        }
    } else {
        state.alpha = xlog_sum(state.alpha, new_score);
    }
}

void LinearAlign::update_state_beta(HState &state, const double new_score) {
    state.beta = xlog_sum(state.beta, new_score);
}

void LinearAlign::update_coinc_prob(const std::pair<int, int> &pos, const double prob) {
    coinc_prob[pos] = xlog_sum(coinc_prob[pos], prob);
}

double LinearAlign::get_trans_emit_prob(const int i, const int j, const HStateType h, const HStateType h_prev) {
    int emit_idx;

    if (i == seq1->size() + 1 && j == seq2->size() + 1) {
        emit_idx = 26;
    } else {
        // seq encoding is (0:N, 1:A, 2:C, 3:G, 4:U, 5:T)
        // lookup table encoding is (0:A, 1:C, 2:G, 3:U, 3:T, 4:.)
        // gap (-) is encoded as 4 in the lookup table
        int nuci = 4, nucj = 4;
        if (h == HStateType::INS1) {
            nuci = seq1->at(i - 1);
        } else if (h == HStateType::INS2) {
            nucj = seq2->at(j - 1);
        } else {
            nuci = seq1->at(i - 1);
            nucj = seq2->at(j - 1);
        }
        emit_idx = nuci * 5 + nucj;
    }

    int prev_h = static_cast<int>(h_prev);
    int curr_h = static_cast<int>(h);

    return xlog_mul((*trans_probs_)[prev_h][curr_h], (*emit_probs_)[emit_idx][curr_h]);
}

float LinearAlign::get_match_score(const int i, const int j) {
    if (i >= seq1->size() || j >= seq2->size()) {
        // std::cout << i << " " << seq1->size() << " " << j << " " << seq2->size() << std::endl;
        return 0;
    }

    float t1 = sqrt(pm1->upstrm.at(i) * pm2->upstrm.at(j));
    float t2 = sqrt(pm1->dwnstrm.at(i) * pm2->dwnstrm.at(j));
    float t3 = sqrt((1 - pm1->upstrm.at(i) - pm1->dwnstrm.at(i)) * (1 - pm2->upstrm.at(j) - pm2->dwnstrm.at(j)));
    // std::cout << "i: " << i << " j: " << j << " t1: " << t1 << " t2: " << t2 << " t3: " << t3 << " "
    //           << pm1->upstrm.at(i) << " " << pm1->dwnstrm.at(i) << " " << pm2->upstrm.at(j) << " " <<
    //           pm2->dwnstrm.at(j)
    //           << std::endl;
    return ((t1 + t2) * alpha1) + (t3 * alpha2) + (alpha3);
}

void LinearAlign::set_prob_accm(ProbAccm &prob_accm1, ProbAccm &prob_accm2) {
    pm1 = &prob_accm1;
    pm2 = &prob_accm2;
}

void LinearAlign::compute_inside() {
    auto start_time = std::chrono::high_resolution_clock::now();
    if (verbose_output) {
        std::cout << "\n[LinearAlign] Running Inside Algorithm:" << std::endl;
    }
    bool use_match_score = (pm1 != nullptr && pm2 != nullptr);
    for (int s = 0; s <= seq_len_sum; ++s) {
        if (verbose_output) {
            Utility::showProgressBar(s, seq_len_sum);
        }
        for (const HStateType h : hstate_types) {
            std::unordered_map<std::pair<int, int>, HState, PairHash> *beam = get_beam(h);
            for (const auto &item : beam[s]) {
                int i = item.first.first;
                int j = item.first.second;
                HState &state = beam[s][{i, j}];

                // INS1
                if (i < seq1->size() && j <= seq2->size()) {
                    double prob = get_trans_emit_prob(i + 1, j, HStateType::INS1, h);
                    double new_score = xlog_mul(state.alpha, prob);
                    update_state_alpha(bestINS1[s + 1][{i + 1, j}], new_score, HStateType::INS1, h);
                }

                // INS2
                if (i <= seq1->size() && j < seq2->size()) {
                    double prob = get_trans_emit_prob(i, j + 1, HStateType::INS2, h);
                    double new_score = xlog_mul(state.alpha, prob);
                    update_state_alpha(bestINS2[s + 1][{i, j + 1}], new_score, HStateType::INS2, h);
                }

                // ALN
                const bool end_check = (i == seq1->size() && j == seq2->size());
                if ((i < seq1->size() && j < seq2->size()) || end_check) {
                    double prob = get_trans_emit_prob(i + 1, j + 1, HStateType::ALN, h);
                    double new_score = xlog_mul(state.alpha, prob);
                    if (use_match_score) {
                        float match_score = get_match_score(i, j);
                        new_score = xlog_mul(new_score, match_score);
                    }
                    update_state_alpha(bestALN[s + 2][{i + 1, j + 1}], new_score, HStateType::ALN, h);
                }
            }
        }
    }
    // update/print time stats
    auto end_time = std::chrono::high_resolution_clock::now();
    inside_execution_time = std::chrono::duration_cast<std::chrono::milliseconds>(end_time - start_time).count();
    if (verbose_output) {
        std::cout << "  - Execution Time: " << inside_execution_time << " ms" << std::endl;
    }
}

// legacy methods below
// ------------------------------------------------------------------------------------------------------------------------
void LinearAlign::run_backward_phase(bool compute_coinc_probs = true) {
    float p_xy = bestALN[seq_len_sum + 2][{seq1->size() + 1, seq2->size() + 1}].alpha;
    if (compute_coinc_probs) {
        coinc_prob[{0, 0}] = p_xy;
    }

    auto start_time = std::chrono::high_resolution_clock::now();
    if (verbose_output) {
        std::cout << "\n[LinearAlign] Running Outside Algorithm:" << std::endl;
    }
    bool use_match_score = (pm1 != nullptr && pm2 != nullptr);
    for (int s = seq_len_sum; s >= 0; --s) {
        if (verbose_output) {
            Utility::showProgressBar(seq_len_sum - s, seq_len_sum);
        }
        for (const HStateType h : hstate_types) {
            std::unordered_map<std::pair<int, int>, HState, PairHash> *beam = get_beam(h);
            for (const auto &item : beam[s]) {
                int i = item.first.first;
                int j = item.first.second;

                HState &state = beam[s][{i, j}];

                const bool end_check = (i == seq1->size() && j == seq2->size());
                // ALN
                if ((i < seq1->size() && j < seq2->size()) || end_check) {
                    double prob = get_trans_emit_prob(i + 1, j + 1, HStateType::ALN, h);
                    double score = xlog_mul(prob, bestALN[s + 2][{i + 1, j + 1}].beta);
                    if (use_match_score) {
                        float match_score = get_match_score(i, j);
                        score = xlog_mul(score, match_score);
                    }
                    update_state_beta(state, score);
                }

                // INS1
                if (i < seq1->size() && j <= seq2->size()) {
                    double prob = get_trans_emit_prob(i + 1, j, HStateType::INS1, h);
                    double score = xlog_mul(prob, bestINS1[s + 1][{i + 1, j}].beta);
                    update_state_beta(state, score);
                }

                // INS2
                if (i <= seq1->size() && j < seq2->size()) {
                    double prob = get_trans_emit_prob(i, j + 1, HStateType::INS2, h);
                    double score = xlog_mul(prob, bestINS2[s + 1][{i, j + 1}].beta);
                    update_state_beta(state, score);
                }

                // coincidence probability
                if (compute_coinc_probs) {
                    double prob = xlog_div(xlog_mul(state.alpha, state.beta), p_xy);
                    update_coinc_prob({i, j}, prob);
                }
            }
        }
    }
    // update/print time stats
    auto end_time = std::chrono::high_resolution_clock::now();
    outside_execution_time = std::chrono::duration_cast<std::chrono::milliseconds>(end_time - start_time).count();
    if (verbose_output) {
        printf("  - Execution Time: %.2f ms (%.2f%% of inside time)\n", outside_execution_time,
               100.0 * outside_execution_time / std::max(inside_execution_time, 1.0f));
    }
}

void LinearAlign::traceback() {
    HStateType h = bestALN[seq_len_sum + 2][{seq1->size() + 1, seq2->size() + 1}].pre;

    int i = seq1->size();
    int j = seq2->size();

    std::string aln1 = "";
    std::string aln2 = "";

    while (i != 0 || j != 0) {
        switch (h) {
        case ALN:
            // std::cout << "ALN" << " " << i << " " << j << std::endl;
            h = bestALN[i + j][{i, j}].pre;
            i -= 1;
            j -= 1;
            aln1 += std::to_string(seq1->at(i));
            aln2 += std::to_string(seq2->at(j));
            break;

        case INS1:
            // std::cout << "INS1" << " " << i << " " << j << std::endl;
            h = bestINS1[i + j][{i, j}].pre;
            i -= 1;
            aln1 += std::to_string(seq1->at(i));
            aln2 += "-";
            break;

        case INS2:
            // std::cout << "INS2" << " " << i << " " << j << std::endl;
            h = bestINS2[i + j][{i, j}].pre;
            j -= 1;
            aln1 += "-";
            aln2 += std::to_string(seq2->at(j));
            break;
        }
    }

    std::reverse(aln1.begin(), aln1.end());
    std::reverse(aln2.begin(), aln2.end());

    std::cout << "\nAlignment: " << std::endl;
    std::cout << aln1 << std::endl;
    std::cout << aln2 << std::endl;
}