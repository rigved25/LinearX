#include "./linear_align.hpp"

void AlignBeam::reset(bool force) {
    if (!has_data && !force) {
        return;
    }

    bestALN.assign(seq_len_sum + 3, std::unordered_map<std::pair<int, int>, HState, PairHash>());
    bestINS1.assign(seq_len_sum + 1, std::unordered_map<std::pair<int, int>, HState, PairHash>());
    bestINS2.assign(seq_len_sum + 1, std::unordered_map<std::pair<int, int>, HState, PairHash>());

    has_data = false;
}

void AlignBeam::save(std::unordered_map<std::pair<int, int>, HState, PairHash> *bestALN,
                     std::unordered_map<std::pair<int, int>, HState, PairHash> *bestINS1,
                     std::unordered_map<std::pair<int, int>, HState, PairHash> *bestINS2) {
    for (int s = 0; s < seq_len_sum + 3; ++s) {
        this->bestALN[s] = bestALN[s];
    }
    for (int s = 0; s < seq_len_sum + 1; ++s) {
        this->bestINS1[s] = bestINS1[s];
        this->bestINS2[s] = bestINS2[s];
    }
    has_data = true;
}

void AlignBeam::save(LinearAlign &la) {
    save(la.bestALN, la.bestINS1, la.bestINS2);
    has_data = true;
}

void LinearAlign::reset_beams() {
    delete[] bestALN;
    delete[] bestINS1;
    delete[] bestINS2;

    bestALN = new std::unordered_map<std::pair<int, int>, HState, PairHash>[seq_len_sum + 3];
    bestINS1 = new std::unordered_map<std::pair<int, int>, HState, PairHash>[seq_len_sum + 1];
    bestINS2 = new std::unordered_map<std::pair<int, int>, HState, PairHash>[seq_len_sum + 1];

    bestALN[0][{0, 0}].alpha = LOG(1.0);
    bestALN[seq_len_sum + 2][{seq1->size() + 1, seq2->size() + 1}].beta = LOG(1.0);
}

double LinearAlign::beam_prune(std::unordered_map<std::pair<int, int>, HState, PairHash> &beamstep, HStateType h,
                               int beam_size) {
    std::vector<std::pair<double, std::pair<int, int>>> scores;
    scores.reserve(beamstep.size());
    for (auto &item : beamstep) {
        std::pair<int, int> aij = item.first;
        HState &cand = item.second;
        scores.push_back(std::make_pair(cand.alpha, aij));
    }
    if (scores.size() <= beam_size) return NEG_INF;
    double threshold = Utility::quickselect(scores, 0, scores.size() - 1, scores.size() - beam_size);
    for (auto &p : scores) {
        if (p.first < threshold) beamstep.erase(p.second);
    }
    return threshold;
}

void LinearAlign::update_state_alpha(HState &state, const double new_score, const HStateType h, const HStateType pre,
                                     const bool best_only = false) {
    if (best_only) {
        if (new_score > state.alpha) {
            state.alpha = new_score;
        }
    } else {
        state.alpha = LOG_SUM(state.alpha, new_score);
    }
}

void LinearAlign::update_state_beta(HState &state, const double new_score) {
    state.beta = LOG_SUM(state.beta, new_score);
}

double LinearAlign::get_trans_emit_prob(const int i, const int j, const HStateType h, const HStateType h_prev) {
    int emit_idx;

    if (i == seq1->size() + 1 && j == seq2->size() + 1) {
        emit_idx = 26;
    } else {
        // seq encoding is: N->0, A->1, C->2, G->3, U->4
        // lookup table encoding is: A->0, C->1, G->2, U->3, .->4
        // gap (-) is encoded as 4 in the lookup table
        int nuci = 4, nucj = 4;
        if (h == HStateType::INS1) {
            nuci = seq1->at(i - 1) - 1;
        } else if (h == HStateType::INS2) {
            nucj = seq2->at(j - 1) - 1;
        } else {
            nuci = seq1->at(i - 1) - 1;
            nucj = seq2->at(j - 1) - 1;
        }
        emit_idx = nuci * 5 + nucj;
    }

    int prev_h = static_cast<int>(h_prev);
    int curr_h = static_cast<int>(h);

    const double tp_val = phmm->get_trans_prob(prev_h, curr_h);
    const double ep_val = phmm->get_emit_prob(emit_idx, curr_h);
    const double score = LOG_MUL(tp_val, ep_val);
    return score;
}

double LinearAlign::get_match_score(const int i, const int j) {
    if (i >= seq1->size() || j >= seq2->size()) {
        return 0;
    }

    const double t1 = sqrt(pm1->upstrm.at(i) * pm2->upstrm.at(j));
    const double t2 = sqrt(pm1->dwnstrm.at(i) * pm2->dwnstrm.at(j));
    const double t3 = sqrt(std::max(1 - pm1->upstrm.at(i) - pm1->dwnstrm.at(i), 0.0) *
                           std::max(1 - pm2->upstrm.at(j) - pm2->dwnstrm.at(j), 0.0));

    const double output = ((t1 + t2) * alpha1) + (t3 * alpha2) + (alpha3);
    return LOG(output);
}

void LinearAlign::set_prob_accm(ProbAccm &prob_accm1, ProbAccm &prob_accm2) {
    pm1 = &prob_accm1;
    pm2 = &prob_accm2;
}

void LinearAlign::compute_inside(bool best_only, int beam_size, bool verbose_output) {
    auto start_time = std::chrono::high_resolution_clock::now();
    if (verbose_output) {
        std::cerr << "[LinearAlign] Running Inside Algorithm:" << std::endl;
    }
    bool use_match_score = (pm1 != nullptr && pm2 != nullptr);
    for (int s = 0; s <= seq_len_sum; ++s) {
        if (verbose_output) {
            Utility::showProgressBar(s, seq_len_sum);
        }
        for (const HStateType h : hstate_types) {
            std::unordered_map<std::pair<int, int>, HState, PairHash> *beam = get_beam(h);
            auto &current_beam = beam[s];
            if (beam_size > 0 && current_beam.size() > beam_size) {
                beam_prune(current_beam, h, beam_size);
            }
            for (const auto &item : current_beam) {
                int i = item.first.first;
                int j = item.first.second;
                HState &state = current_beam[{i, j}];

                // INS1
                if (i < seq1->size() && j <= seq2->size()) {
                    double prob = get_trans_emit_prob(i + 1, j, HStateType::INS1, h);
                    double new_score = LOG_MUL(state.alpha, prob);
                    update_state_alpha(bestINS1[s + 1][{i + 1, j}], new_score, HStateType::INS1, h, best_only);
                }

                // INS2
                if (i <= seq1->size() && j < seq2->size()) {
                    double prob = get_trans_emit_prob(i, j + 1, HStateType::INS2, h);
                    double new_score = LOG_MUL(state.alpha, prob);
                    update_state_alpha(bestINS2[s + 1][{i, j + 1}], new_score, HStateType::INS2, h, best_only);
                }

                // ALN
                const bool end_check = (i == seq1->size() && j == seq2->size());
                if ((i < seq1->size() && j < seq2->size()) || end_check) {
                    double prob = get_trans_emit_prob(i + 1, j + 1, HStateType::ALN, h);
                    double new_score = LOG_MUL(state.alpha, prob);
                    if (use_match_score) {
                        double match_score = get_match_score(i, j);
                        new_score = LOG_MUL(new_score, match_score);
                    }
                    update_state_alpha(bestALN[s + 2][{i + 1, j + 1}], new_score, HStateType::ALN, h, best_only);
                }
            }
        }
    }
    // update/print time stats
    auto end_time = std::chrono::high_resolution_clock::now();
    inside_execution_time = std::chrono::duration_cast<std::chrono::milliseconds>(end_time - start_time).count();
    if (verbose_output) {
        std::cerr << "  - Execution Time: " << inside_execution_time << " ms\n" << std::endl;
    }
}

void LinearAlign::dump_coinc_probs(const std::string &filepath, const float threshold) const {
    if (coinc_prob == nullptr) {
        throw std::runtime_error(
            "[LinearAlign Error] Coincidence probabilities not computed yet! You must run "
            "compute_coincidence_probabilities() first.");
    }

    // open the file for writing
    std::ofstream file(filepath);
    if (!file) {
        std::cerr
            << "[Hint] The directory for the output file may not exist. Please create it before running the method."
            << std::endl;
        throw std::runtime_error("[LinearAlign Error] Unable to open the file " + filepath +
                                 " for writing coincidence probabilities.");
    }

    // dump the coincidence probabilities to the file
    for (int i = 0; i < seq1->size(); ++i) {
        for (const auto &item : coinc_prob[i]) {
            const int j = item.first;
            const double prob = item.second;
            if (prob < threshold) continue;

            // output i, j, and the probability to the file
            file << i << " " << j << " " << std::fixed << std::setprecision(4) << prob << std::endl;
        }
    }
};

// legacy methods below
// ------------------------------------------------------------------------------------------------------------------------
// void LinearAlign::run_backward_phase(bool verbose_output) {
//     auto start_time = std::chrono::high_resolution_clock::now();
//     if (verbose_output) {
//         std::cout << "[LinearAlign] Running Outside Algorithm:" << std::endl;
//     }
//     bool use_match_score = (pm1 != nullptr && pm2 != nullptr);
//     for (int s = seq_len_sum; s >= 0; --s) {
//         if (verbose_output) {
//             Utility::showProgressBar(seq_len_sum - s, seq_len_sum);
//         }
//         for (const HStateType h : hstate_types) {
//             std::unordered_map<std::pair<int, int>, HState, PairHash> *beam = get_beam(h);
//             for (const auto &item : beam[s]) {
//                 int i = item.first.first;
//                 int j = item.first.second;
//                 HState &state = beam[s][{i, j}];

//                 // INS1
//                 if (i < seq1->size() && j <= seq2->size()) {
//                     double prob = get_trans_emit_prob(i + 1, j, HStateType::INS1, h);
//                     double score = LOG_MUL(prob, bestINS1[s + 1][{i + 1, j}].beta);
//                     update_state_beta(state, score);
//                 }

//                 // INS2
//                 if (i <= seq1->size() && j < seq2->size()) {
//                     double prob = get_trans_emit_prob(i, j + 1, HStateType::INS2, h);
//                     double score = LOG_MUL(prob, bestINS2[s + 1][{i, j + 1}].beta);
//                     update_state_beta(state, score);
//                 }

//                 // ALN
//                 const bool end_check = (i == seq1->size() && j == seq2->size());
//                 if ((i < seq1->size() && j < seq2->size()) || end_check) {
//                     double prob = get_trans_emit_prob(i + 1, j + 1, HStateType::ALN, h);
//                     double score = LOG_MUL(prob, bestALN[s + 2][{i + 1, j + 1}].beta);
//                     if (use_match_score) {
//                         double match_score = get_match_score(i, j);
//                         score = LOG_MUL(score, match_score);
//                     }
//                     update_state_beta(state, score);
//                 }
//             }
//         }
//     }
//     // update/print time stats
//     auto end_time = std::chrono::high_resolution_clock::now();
//     outside_execution_time = std::chrono::duration_cast<std::chrono::milliseconds>(end_time - start_time).count();
//     if (verbose_output) {
//         printf("  - Execution Time: %.2f ms (%.2f%% of inside time)\n", outside_execution_time,
//                100.0 * outside_execution_time / std::max(inside_execution_time, 1.0));
//     }
// }

// MultiSeq LinearAlign::old_traceback() {
//     HStateType h = bestALN[seq_len_sum + 2][{seq1->size() + 1, seq2->size() + 1}].pre;

//     int i = seq1->size();
//     int j = seq2->size();

//     std::string aln1 = "";
//     std::string aln2 = "";

//     while (i != 0 || j != 0) {
//         switch (h) {
//         case ALN:
//             h = bestALN[i + j][{i, j}].pre;
//             i -= 1;
//             j -= 1;
//             aln1 += std::to_string(seq1->at(i));
//             aln2 += std::to_string(seq2->at(j));
//             break;

//         case INS1:
//             h = bestINS1[i + j][{i, j}].pre;
//             i -= 1;
//             aln1 += std::to_string(seq1->at(i));
//             aln2 += "-";
//             break;

//         case INS2:
//             h = bestINS2[i + j][{i, j}].pre;
//             j -= 1;
//             aln1 += "-";
//             aln2 += std::to_string(seq2->at(j));
//             break;
//         }
//     }

//     std::reverse(aln1.begin(), aln1.end());
//     std::reverse(aln2.begin(), aln2.end());

//     MultiSeq alignment;
//     alignment.add_sequence(Seq(this->sequence1->id, aln1));
//     alignment.add_sequence(Seq(this->sequence2->id, aln2));
//     return alignment;
// }
