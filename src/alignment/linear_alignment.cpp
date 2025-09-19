// src/alignment/linear_align.cpp
#include <fstream>
#include <linearx/alignment/config.hpp>
#include <linearx/alignment/linear_align.hpp>
#include <set>

using namespace linearx::utils;
using namespace linearx::constants::math;
using namespace linearx::math;
using namespace std;

template <typename T>
LinearAlignmentInterface<T>::LinearAlignmentInterface(Sequence& seq1, Sequence& seq2, float alpha1, float alpha2,
                                                      float alpha3)
    : seq1(seq1),
      seq2(seq2),
      seq_len_sum(seq1.length() + seq2.length()),
      alpha1(alpha1),
      alpha2(alpha2),
      alpha3(alpha3) {
    seq1.randomize_N();
    seq2.randomize_N();
}

template <typename T>
void LinearAlignmentInterface<T>::use_prob_set1() {
    delete phmm;
    phmm = new Phmm(Phmm::EMIT_PROBS, Phmm::TRANS_PROBS);
}

template <typename T>
void LinearAlignmentInterface<T>::use_prob_set2(float similarity) {
    delete phmm;
    string phmm_pars_fp = PARAM_FILE_PATH;
    phmm = new Phmm(phmm_pars_fp.c_str());
    phmm->set_parameters_by_sim(similarity);
}

template <typename T>
void LinearAlignmentInterface<T>::reset_beams(unsigned beam_size) {
    reset_beam_vector(bestALN, seq_len_sum + 3);
    reset_beam_vector(bestINS1, seq_len_sum + 1);
    reset_beam_vector(bestINS2, seq_len_sum + 1);

    bestALN[0][{0, 0}].alpha = LOG_ONE;
    bestALN[seq_len_sum + 2][{seq1.length() + 1, seq2.length() + 1}].beta = LOG_ONE;
}

template <typename T>
void LinearAlignmentInterface<T>::set_prob_accm(ProbAccm& prob_accm1, ProbAccm& prob_accm2) {
    pm1 = &prob_accm1;
    pm2 = &prob_accm2;
}

template <typename T>
value_type LinearAlignmentInterface<T>::get_trans_emit_prob(const int i, const int j, const HStateType h,
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

    const value_type tp_val = phmm->get_trans_prob(prev_h, curr_h);
    const value_type ep_val = phmm->get_emit_prob(emit_idx, curr_h);
    const value_type score = LOG_MUL(tp_val, ep_val);
    return score;
}

template <typename T>
void LinearAlignmentInterface<T>::compute_coincidence_probabilities(bool verbose_output) {
    // clear the previous matrix
    reset_beam_vector(coinc_prob, seq1.length());
    reset_beam_vector(prob_rev_idx, seq2.length());

    const value_type p_xy = bestALN[seq_len_sum + 2][{seq1.length() + 1, seq2.length() + 1}].alpha;
    for (int s = 0; s <= seq_len_sum; ++s) {
        for (const HStateType h : hstate_types) {
            vector<unordered_map<pair<int, int>, HState, PairHash>>& beam = get_beams(h);
            for (auto& item : beam[s]) {
                const int i = item.first.first;
                const int j = item.first.second;
                HState& state = item.second;

                const value_type prob = LOG_DIV(LOG_MUL(state.alpha, state.beta), p_xy);
                if (prob > -linearx::constants::limits::DEVIATION_THRESHOLD && i > 0 && j > 0) {
                    auto [ptr_cprob_ij, inserted] = coinc_prob[i - 1].try_emplace(j - 1, LOG_ZERO);
                    ptr_cprob_ij->second = LOG_SUM(ptr_cprob_ij->second, prob);
                }
            }
        }
    }

    unsigned long num_pruned = 0;  // for keeping track of pruned P(i,j)s
    unsigned long num_saved = 0;   // for keeping track of saved P(i,j)s
    const value_type fam_threshold = phmm->get_fam_threshold();
    for (int i = 0; i < seq1.length(); ++i) {
        for (auto it = coinc_prob[i].begin(); it != coinc_prob[i].end();) {
            const int j = it->first;
            value_type& prob = it->second;
            if (prob < fam_threshold) {
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

template <typename T>
void LinearAlignmentInterface<T>::dump_coinc_probs(const std::string& out_dir) const {
    if (coinc_prob.empty()) {
        throw std::runtime_error(
            "[LinearAlignment Error] Coincidence probabilities not computed yet! You must run "
            "compute_coincidence_probabilities() first.");
    }

    // create output directory if it doesn't exist
    std::filesystem::create_directories(out_dir);

    // construct output file path: cp_{seq1ID}_{seq2ID}.txt
    std::string filename = out_dir + "/cp_" + std::to_string(seq1.id) + "_" + std::to_string(seq2.id) + ".txt";

    // open file for writing
    std::ofstream file(filename);
    if (!file) {
        throw std::runtime_error("[LinearAlignment Error] Unable to open file " + filename +
                                 " for writing coincidence probabilities.");
    }

    // write all coincidence probabilities to the file
    for (int i = 0; i < seq1.length(); ++i) {
        for (const auto& item : coinc_prob[i]) {
            const int j = item.first;
            const value_type prob = item.second;
            file << i + 1 << " " << j + 1 << " " << std::fixed << std::setprecision(4) << prob << "\n";
        }
    }
}

template <typename T>
void LinearAlignmentInterface<T>::print_alpha_beta() const {
    cerr << "Alpha(ALN, n1 + 1, n2 + 1): " << bestALN[seq_len_sum + 2].at({seq1.length() + 1, seq2.length() + 1}).alpha
         << endl;
    cerr << "Beta(ALN, 0, 0): " << bestALN[0].at({0, 0}).beta << endl << endl;
}

template <typename T>
void LinearAlignmentInterface<T>::print_seqs() const {
    seq1.print(cerr);
    seq2.print(cerr);
    cerr << endl;
}

template <typename T>
void LinearAlignmentInterface<T>::print_beams() const {
    for (int s = 0; s < seq_len_sum + 3; ++s) {
        for (auto& item : bestALN[s]) {
            std::cout << "ALN: (" << item.first.first << ", " << item.first.second << ") : " << item.second.alpha << " "
                      << item.second.beta << std::endl;
        }
    }
    std::cout << "--------------------------------" << std::endl;
    for (int s = 0; s < seq_len_sum + 1; ++s) {
        for (auto& item : bestINS1[s]) {
            std::cout << "INS1: (" << item.first.first << ", " << item.first.second << ") : " << item.second.alpha
                      << " " << item.second.beta << std::endl;
        }
    }
    std::cout << "--------------------------------" << std::endl;
    for (int s = 0; s < seq_len_sum + 1; ++s) {
        for (auto& item : bestINS2[s]) {
            std::cout << "INS2: (" << item.first.first << ", " << item.first.second << ") : " << item.second.alpha
                      << " " << item.second.beta << std::endl;
        }
    }
}

// Instantiate templates for LinearAlignmentInterface with desired types
#define X(TYPE) template class LinearAlignmentInterface<TYPE>;
LA_TEMPLATE_TYPES
#undef X
