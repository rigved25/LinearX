// src/alignment/linear_align.cpp
#include <fstream>
#include <linearx/alignment/linear_align.hpp>

using namespace linearx::utils;
using namespace linearx::constants::math;
using namespace std;

LinearAlignment::LinearAlignment(Sequence &seq1, Sequence &seq2, bool verbose, double alpha1, double alpha2,
                                 double alpha3)
    : seq1(seq1),
      seq2(seq2),
      seq_len_sum(seq1.length() + seq2.length()),
      alpha1(alpha1),
      alpha2(alpha2),
      alpha3(alpha3) {
    reset_beams();
    seq1.randomize_N();
    seq2.randomize_N();

    if (verbose) {
        cout << "LinearAlignment initialized with sequences of lengths: " << seq1.length() << " and " << seq2.length()
             << endl;
    }
}

void LinearAlignment::use_prob_set1() {
    delete phmm;
    phmm = new Phmm(Phmm::EMIT_PROBS, Phmm::TRANS_PROBS);
}

void LinearAlignment::use_prob_set2(float similarity) {
    delete phmm;
    string phmm_pars_fp = PARAM_FILE_PATH;
    phmm = new Phmm(phmm_pars_fp.c_str());
    phmm->set_parameters_by_sim(similarity);
}

void LinearAlignment::reset_beams(bool freeMemory) {
    if (freeMemory) {
        bestALN.clear();
        bestINS1.clear();
        bestINS2.clear();
    }

    bestALN = vector<unordered_map<pair<int, int>, HState, PairHash>>(seq_len_sum + 3);
    bestINS1 = vector<unordered_map<pair<int, int>, HState, PairHash>>(seq_len_sum + 1);
    bestINS2 = vector<unordered_map<pair<int, int>, HState, PairHash>>(seq_len_sum + 1);

    bestALN[0][{0, 0}].alpha = LOG(1.0);
    bestALN[seq_len_sum + 2][{seq1.length() + 1, seq2.length() + 1}].beta = LOG(1.0);
}

double LinearAlignment::get_trans_emit_prob(const int i, const int j, const HStateType h,
                                            const HStateType h_prev) const {
    int emit_idx;

    if (i == seq1.length() + 1 && j == seq2.length() + 1) {
        emit_idx = 26;
    } else {
        // seq encoding is: N->0, A->1, C->2, G->3, U->4
        // lookup table encoding is: A->0, C->1, G->2, U->3, .->4
        // gap (-) is encoded as 4 in the lookup table
        int nuci = 4, nucj = 4;
        if (h == HStateType::INS1) {
            nuci = seq1.enc[i - 1] - 1;
        } else if (h == HStateType::INS2) {
            nucj = seq2.enc[j - 1] - 1;
        } else {
            nuci = seq1.enc[i - 1] - 1;
            nucj = seq2.enc[j - 1] - 1;
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

double LinearAlignment::get_match_score(const int i, const int j) const {
    if (i >= seq1.length() || j >= seq2.length()) {
        return 0;
    }

    const double t1 = sqrt(pm1->upstrm.at(i) * pm2->upstrm.at(j));
    const double t2 = sqrt(pm1->dwnstrm.at(i) * pm2->dwnstrm.at(j));
    const double t3 = sqrt(max(1 - pm1->upstrm.at(i) - pm1->dwnstrm.at(i), 0.0) *
                           max(1 - pm2->upstrm.at(j) - pm2->dwnstrm.at(j), 0.0));

    const double output = ((t1 + t2) * alpha1) + (t3 * alpha2) + (alpha3);
    return LOG(output);
}

void LinearAlignment::set_prob_accm(ProbAccm &prob_accm1, ProbAccm &prob_accm2) {
    pm1 = &prob_accm1;
    pm2 = &prob_accm2;
}

void LinearAlignment::compute_coincidence_probabilities(bool verbose_output) {
    // clear the previous matrix
    coinc_prob.clear();
    prob_rev_idx.clear();
    coinc_prob.resize(seq1.length());
    prob_rev_idx.resize(seq2.length());

    double p_xy = bestALN[seq_len_sum + 2][{seq1.length() + 1, seq2.length() + 1}].alpha;
    for (int s = 0; s <= seq_len_sum; ++s) {
        for (const HStateType h : hstate_types) {
            vector<unordered_map<pair<int, int>, HState, PairHash>> &beam = get_beam(h);
            for (const auto &item : beam[s]) {
                const int i = item.first.first;
                const int j = item.first.second;
                HState &state = beam[s][{i, j}];

                const double prob = LOG_DIV(LOG_MUL(state.alpha, state.beta), p_xy);
                if (prob > -linearx::constants::limits::DEVIATION_THRESHOLD && i > 0 && j > 0) {
                    auto [ptr_cprob_ij, inserted] = coinc_prob[i - 1].try_emplace(j - 1, LOG(0.0));
                    ptr_cprob_ij->second = LOG_SUM(ptr_cprob_ij->second, prob);
                }
            }
        }
    }

    unsigned long num_pruned = 0;  // for keeping track of pruned P(i,j)s
    unsigned long num_saved = 0;   // for keeping track of saved P(i,j)s
    for (int i = 0; i < seq1.length(); ++i) {
        for (auto it = coinc_prob[i].begin(); it != coinc_prob[i].end();) {
            const int j = it->first;
            double &prob = it->second;

            if (prob < phmm->get_fam_threshold()) {
                it = coinc_prob[i].erase(it);  // erase and get the next valid iterator
                ++num_pruned;
            } else {
                prob = EXP(prob);
                if (prob > 1.001) {
                    fprintf(stderr,
                            "[LinearAlignment: Warning] BPP value too high, something is wrong! bpp(%d, %d): %.5f\n", i,
                            j, prob);
                }
                prob = min(prob, 1.0);
                prob_rev_idx[j].push_back(i);
                ++num_saved;
                ++it;  // move to the next element if not erased
            }
        }
    }

    if (verbose_output) {
        fprintf(stderr, "[LinearAlignment] Coincidence Probabilities Computed: %lu (saved) + %lu (pruned)\n", num_saved,
                num_pruned);
    }
}

void LinearAlignment::dump_coinc_probs(const string &filepath, const float threshold) const {
    if (coinc_prob.empty()) {
        throw runtime_error(
            "[LinearAlignment Error] Coincidence probabilities not computed yet! You must run "
            "compute_coincidence_probabilities() first.");
    }

    // open the file for writing
    ofstream file(filepath);
    if (!file) {
        cerr << "[Hint] The directory for the output file may not exist. Please create it before running the method."
             << endl;
        throw runtime_error("[LinearAlignment Error] Unable to open the file " + filepath +
                            " for writing coincidence probabilities.");
    }

    // dump the coincidence probabilities to the file
    for (int i = 0; i < seq1.length(); ++i) {
        for (const auto &item : coinc_prob[i]) {
            const int j = item.first;
            const double prob = item.second;
            if (prob < threshold) continue;

            // output i, j, and the probability to the file
            file << i << " " << j << " " << fixed << setprecision(4) << prob << endl;
        }
    }
};

void LinearAlignment::print_alpha_beta() const {
    cerr << "Alpha(ALN, n1 + 1, n2 + 1): " << bestALN[seq_len_sum + 2].at({seq1.length() + 1, seq2.length() + 1}).alpha
         << endl;
    cerr << "Beta(ALN, 0, 0): " << bestALN[0].at({0, 0}).beta << endl << endl;
}

void LinearAlignment::print_seqs() const {
    seq1.print(cerr);
    seq2.print(cerr);
    cerr << endl;
}

void LinearAlignment::print_beams() const {
    for (int s = 0; s < seq_len_sum + 3; ++s) {
        for (auto &item : bestALN[s]) {
            std::cout << "ALN: (" << item.first.first << ", " << item.first.second << ") : " << item.second.alpha << " "
                      << item.second.beta << std::endl;
        }
    }
    std::cout << "--------------------------------" << std::endl;
    for (int s = 0; s < seq_len_sum + 1; ++s) {
        for (auto &item : bestINS1[s]) {
            std::cout << "INS1: (" << item.first.first << ", " << item.first.second << ") : " << item.second.alpha
                      << " " << item.second.beta << std::endl;
        }
    }
    std::cout << "--------------------------------" << std::endl;
    for (int s = 0; s < seq_len_sum + 1; ++s) {
        for (auto &item : bestINS2[s]) {
            std::cout << "INS2: (" << item.first.first << ", " << item.first.second << ") : " << item.second.alpha
                      << " " << item.second.beta << std::endl;
        }
    }
}

// legacy methods below
// ------------------------------------------------------------------------------------------------------------------------
// MultiSeq LinearAlignment::old_traceback() {
//     HStateType h = bestALN[seq_len_sum + 2][{seq1.length() + 1, seq2.length() + 1}].pre;

//     int i = seq1.length();
//     int j = seq2.length();

//     string aln1 = "";
//     string aln2 = "";

//     while (i != 0 || j != 0) {
//         switch (h) {
//         case ALN:
//             h = bestALN[i + j][{i, j}].pre;
//             i -= 1;
//             j -= 1;
//             aln1 += to_string(seq1->at(i));
//             aln2 += to_string(seq2->at(j));
//             break;

//         case INS1:
//             h = bestINS1[i + j][{i, j}].pre;
//             i -= 1;
//             aln1 += to_string(seq1->at(i));
//             aln2 += "-";
//             break;

//         case INS2:
//             h = bestINS2[i + j][{i, j}].pre;
//             j -= 1;
//             aln1 += "-";
//             aln2 += to_string(seq2->at(j));
//             break;
//         }
//     }

//     reverse(aln1.begin(), aln1.end());
//     reverse(aln2.begin(), aln2.end());

//     MultiSeq alignment;
//     alignment.add_sequence(Sequence(this->sequence1->id, aln1));
//     alignment.add_sequence(Sequence(this->sequence2->id, aln2));
//     return alignment;
// }
