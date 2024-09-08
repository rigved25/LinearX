#include <unordered_map>
#include <vector>

#include "./utility/log_math.hpp"
#include "Seq.hpp"
#include "constants.hpp"
#include "utility.hpp"
#include <iostream> // NOTE: for debugging, remove later

enum HStateType {
    ALN,  // 0
    INS1, // 1
    INS2, // 2
};

struct HState {

    HStateType type;
    HStateType pre;
    float alpha;
    float beta;

    HState() : type(ALN), alpha(VALUE_MIN), beta(VALUE_MIN) {};
};

class LinearAlign {
    inline const static std::vector<HStateType> hstate_types = {ALN, INS1, INS2};

  public:
    LinearAlign(const Seq &seq1, const Seq &seq2, const double (*transition_probs)[3][3] = &trans_probs,
                const double (*emission_probs)[27][3] = &emit_probs)
        : seq1(Seq::encode(seq1, nuc_encoding_scheme)), seq2(Seq::encode(seq2, nuc_encoding_scheme)),
          trans_probs_(transition_probs), emit_probs_(emission_probs) {
        len_sum = this->seq1.size() + this->seq2.size();

        bestALN = new std::unordered_map<std::pair<int, int>, HState, PairHash>[len_sum + 3];
        bestINS1 = new std::unordered_map<std::pair<int, int>, HState, PairHash>[len_sum + 1];
        bestINS2 = new std::unordered_map<std::pair<int, int>, HState, PairHash>[len_sum + 1];

        // print seq1 and seq2
        std::cout << "Seq1: ";
        for (auto &nuc : this->seq1) {
            std::cout << nuc << " ";
        }
        std::cout << std::endl;

        std::cout << "Seq2: ";
        for (auto &nuc : this->seq2) {
            std::cout << nuc << " ";
        }
        std::cout << std::endl;


        // std::cout << "\nStarting Inside: " << std::endl;
        run_forward_phase();
        // debug_beam();

        // std::cout << "\nStarting Outside: " << std::endl;
        run_backward_phase();
        debug_beam();
        // traceback();
    }

    void debug_beam() {

        for (int s = 0; s < len_sum + 3; ++s) {
            for (auto &item : bestALN[s]) {
                std::cout << "ALN: (" << item.first.first << ", " << item.first.second << ") : " << item.second.alpha
                          << " " << item.second.beta << std::endl;
            }
        }
        std::cout << "--------------------------------" << std::endl;
        for (int s = 0; s < len_sum + 1; ++s) {
            for (auto &item : bestINS1[s]) {
                std::cout << "INS1: (" << item.first.first << ", " << item.first.second << ") : " << item.second.alpha
                          << " " << item.second.beta << std::endl;
            }
        }
        std::cout << "--------------------------------" << std::endl;
        for (int s = 0; s < len_sum + 1; ++s) {
            for (auto &item : bestINS2[s]) {
                std::cout << "INS2: (" << item.first.first << ", " << item.first.second << ") : " << item.second.alpha
                          << " " << item.second.beta << std::endl;
            }
        }
    }

    const std::vector<int> seq1;
    const std::vector<int> seq2;
    const double (*trans_probs_)[3][3];
    const double (*emit_probs_)[27][3];

    int len_sum;

    std::unordered_map<std::pair<int, int>, HState, PairHash> *bestALN;
    std::unordered_map<std::pair<int, int>, HState, PairHash> *bestINS1;
    std::unordered_map<std::pair<int, int>, HState, PairHash> *bestINS2;
    std::unordered_map<std::pair<int, int>, float, PairHash> coinc_prob;

    double get_trans_emit_prob(const int i, const int j, const HStateType h, const HStateType h_prev) {
        int emit_idx;

        if (i == seq1.size() + 1 && j == seq2.size() + 1) {
            emit_idx = 26;
        } else {
            // seq encoding is (0:N, 1:A, 2:C, 3:G, 4:U, 5:T)
            // lookup table encoding is (0:A, 1:C, 2:G, 3:U, 3:T, 4:.)
            // gap (-) is encoded as 4 in the lookup table
            int nuci = 4, nucj = 4;
            if (h == HStateType::INS1) {
                nuci = seq1.at(i - 1);
            } else if (h == HStateType::INS2) {
                nucj = seq2.at(j - 1);
            } else {
                nuci = seq1.at(i - 1);
                nucj = seq2.at(j - 1);
            }
            emit_idx = nuci * 5 + nucj;
        }

        int prev_h = static_cast<int>(h_prev);
        int curr_h = static_cast<int>(h);

        return xlog_mul((*trans_probs_)[prev_h][curr_h], (*emit_probs_)[emit_idx][curr_h]);
    }

    void update_state_alpha(HState &state, const double new_score, const HStateType h, const HStateType pre,
                const bool best_only = false) {
        if (best_only) {
            if (new_score > state.alpha) {
                state.alpha = new_score;
                state.type = h;
                state.pre = pre;
            }
        } else {
            state.alpha = xlog_sum(state.alpha, new_score);
        }
    }

    void update_state_beta(HState &state, const double new_score) {
        state.beta = xlog_sum(state.beta, new_score);
    }

    void update_coinc_prob(const std::pair<int, int> &pos, const double prob) {
        coinc_prob[pos] = xlog_sum(coinc_prob[pos], prob);
    }

    void run_forward_phase() {

        bestALN[0][{0, 0}] = HState();
        bestALN[0][{0, 0}].alpha = 1.0;

        for (int s = 0; s <= len_sum; ++s) {
            for (const HStateType h : hstate_types) {
                std::unordered_map<std::pair<int, int>, HState, PairHash> *beam;
                if (h == HStateType::ALN) {
                    beam = bestALN;
                } else if (h == HStateType::INS1) {
                    beam = bestINS1;
                } else {
                    beam = bestINS2;
                }

                for (const auto &item : beam[s]) {
                    int i = item.first.first;
                    int j = item.first.second;
                    HState &state = beam[s][{i, j}];

                    // ALN
                    const bool end_check = (i == seq1.size() && j == seq2.size());
                    if ((i < seq1.size() && j < seq2.size()) || end_check) {
                        double prob = get_trans_emit_prob(i + 1, j + 1, HStateType::ALN, h);
                        double new_score = xlog_mul(state.alpha, prob);
                        update_state_alpha(bestALN[s + 2][{i + 1, j + 1}], new_score, HStateType::ALN, h);
                        // if (end_check) {
                        //     std::cout << "ALN (END): " << h << " : " << i << " " << j << " " << " : " << state.alpha
                        //               << ", " << bestALN[s + 2][{i + 1, j + 1}].alpha << std::endl;
                        // }
                    }

                    // INS1
                    if (i < seq1.size() && j <= seq2.size()) {
                        double prob = get_trans_emit_prob(i + 1, j, HStateType::INS1, h);
                        double new_score = xlog_mul(state.alpha, prob);
                        update_state_alpha(bestINS1[s + 1][{i + 1, j}], new_score, HStateType::INS1, h);
                    }

                    // INS2
                    if (i <= seq1.size() && j < seq2.size()) {
                        double prob = get_trans_emit_prob(i, j + 1, HStateType::INS2, h);
                        double new_score = xlog_mul(state.alpha, prob);
                        update_state_alpha(bestINS2[s + 1][{i, j + 1}], new_score, HStateType::INS2, h);
                    }
                }
            }
        }
    }


    void run_backward_phase(bool compute_coinc_probs = true) {
        bestALN[len_sum + 2][{seq1.size() + 1, seq2.size() + 1}].beta = 1.0;


        float p_xy = bestALN[len_sum + 2][{seq1.size() + 1, seq2.size() + 1}].alpha;
        if (compute_coinc_probs) {
            coinc_prob[{0, 0}] = p_xy;
        }

        for (int s = len_sum; s >= 0; --s) {
            for (const HStateType h : hstate_types) {
                std::unordered_map<std::pair<int, int>, HState, PairHash> *beam;
                if (h == HStateType::ALN) {
                    beam = bestALN;
                } else if (h == HStateType::INS1) {
                    beam = bestINS1;
                } else {
                    beam = bestINS2;
                }

                for (const auto &item : beam[s]) {
                    int i = item.first.first;
                    int j = item.first.second;

                    HState &state = beam[s][{i, j}];

                    const bool end_check = (i == seq1.size() && j == seq2.size());
                    // ALN
                    if ((i < seq1.size() && j < seq2.size()) || end_check) {
                        double prob = get_trans_emit_prob(i + 1, j + 1, HStateType::ALN, h);
                        double score = xlog_mul(prob, bestALN[s + 2][{i + 1, j + 1}].beta);
                        update_state_beta(state, score);
                    }

                    // INS1
                    if (i < seq1.size() && j <= seq2.size()) {
                        double prob = get_trans_emit_prob(i + 1, j, HStateType::INS1, h);
                        double score = xlog_mul(prob, bestINS1[s + 1][{i + 1, j}].beta);
                        update_state_beta(state, score);
                    }

                    // INS2
                    if (i <= seq1.size() && j < seq2.size()) {
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
    }

    void traceback() {
        // std::cout << "\nSteps:" << std::endl;

        HStateType h = bestALN[len_sum + 2][{seq1.size() + 1, seq2.size() + 1}].pre;

        int i = seq1.size();
        int j = seq2.size();

        std::string aln1 = "";
        std::string aln2 = "";

        while (i != 0 || j != 0) {
            switch (h) {
            case ALN:
                // std::cout << "ALN" << " " << i << " " << j << std::endl;
                h = bestALN[i + j][{i, j}].pre;
                i -= 1;
                j -= 1;
                aln1 += std::to_string(seq1.at(i));
                aln2 += std::to_string(seq2.at(j));
                break;

            case INS1:
                // std::cout << "INS1" << " " << i << " " << j << std::endl;
                h = bestINS1[i + j][{i, j}].pre;
                i -= 1;
                aln1 += std::to_string(seq1.at(i));
                aln2 += "-";
                break;

            case INS2:
                // std::cout << "INS2" << " " << i << " " << j << std::endl;
                h = bestINS2[i + j][{i, j}].pre;
                j -= 1;
                aln1 += "-";
                aln2 += std::to_string(seq2.at(j));
                break;
            }
        }

        std::reverse(aln1.begin(), aln1.end());
        std::reverse(aln2.begin(), aln2.end());

        std::cout << "\nAlignment: " << std::endl;
        std::cout << aln1 << std::endl;
        std::cout << aln2 << std::endl;
    }
};
