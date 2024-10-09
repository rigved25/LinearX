#include "turbofold.hpp"

int LinearTurboFold::get_seq_pair_index(const int k1, const int k2) {
    // Ensure the indices are within bounds
    if (k1 >= multi_seq->size() || k2 >= multi_seq->size()) {
        throw std::out_of_range("Seq index k out of range");
    }

    // Ensure k1 != k2 and k1 < k2 by swapping if necessary
    int a = k1, b = k2;
    if (k1 == k2) {
        throw std::invalid_argument("k1 and k2 must be different");
    }
    if (k1 > k2) {
        std::swap(a, b);
    }

    // Number of sequences
    int n = multi_seq->size();

    // Compute the index using the triangular indexing formula
    int index = a * (2 * n - a - 1) / 2 + (b - a - 1);

    return index;
}

float LinearTurboFold::get_extrinsic_info(const Seq &x, const int i, const int j) {
    float output = 0.0;
    if (this->itr == 0)
        return output;

    if (ext_info_cache.exists(x.k_id, i, j))
        return ext_info_cache.get(x.k_id, i, j);

    for (const Seq &y : *multi_seq) {
        if (x.k_id == y.k_id)
            continue;

        const int seq_pair_idx = get_seq_pair_index(x.k_id, y.k_id);
        const float seq_idnty = seq_identities[seq_pair_idx];
        const LinearAlign &aln = alns[seq_pair_idx];
        const TurboPartition &y_pf = pfs.at(y.k_id);

        if (x.k_id == aln.sequence1->k_id) {
            for (const auto &itr1 : aln.coinc_prob[i]) {
                for (const auto &itr2 : aln.coinc_prob[j]) {
                    const int k = itr1.first;
                    const int l = itr2.first;
                    const double aln_prob_ik = itr1.second;
                    const double aln_prob_jl = itr2.second;
                    const double y_bpp_kl = y_pf.get_bpp(k, l);

                    output += y_bpp_kl * aln_prob_ik * aln_prob_jl;
                }
            }
        } else {
            for (const int k : aln.prob_rev_idx[i]) {
                for (const int l : aln.prob_rev_idx[j]) {
                    const double aln_prob_ik = aln.get_bpp(k, i);
                    const double aln_prob_jl = aln.get_bpp(l, j);
                    const double y_bpp_kl = y_pf.get_bpp(k, l);

                    output += y_bpp_kl * aln_prob_ik * aln_prob_jl;
                }
            }
        }

        output *= (1 - seq_idnty);
    }

    output = xlog(output);
    ext_info_cache.insert(x.k_id, i, j, output);
    return output;
}

void LinearTurboFold::run() {
    int max_itr = 4;
    for (itr = 0; itr <= max_itr; ++itr) {
        // Utility::showProgressBar(itr, max_itr);
        // align step
        if (itr > 0) {
            for (LinearAlign &aln : alns) {
                const int k1 = aln.sequence1->k_id;
                const int k2 = aln.sequence2->k_id;
                const int aln_pair_index = get_seq_pair_index(k1, k2);
                const float seq_idnty = seq_identities[aln_pair_index];

                // get the alignments
                aln.reset_beams();
                // aln.prob_set1();
                itr <= 1 ? aln.prob_set1() : aln.prob_set2(seq_idnty);
                aln.set_prob_accm(pfs[k1].prob_accm, pfs[k2].prob_accm);
                aln.compute_inside(true, 100, verbose_state == VerboseState::DEBUG);
                MultiSeq alignment = aln.get_alignment();
                seq_identities[aln_pair_index] = alignment.get_seq_identity();

                if (verbose_state == VerboseState::DEBUG) {
                    std::cout << "Alignment: " << k1 << " " << k2 << std::endl;
                    std::cout << alignment[0].sequence << std::endl;
                    std::cout << alignment[1].sequence << std::endl;
                    aln.print_alpha_beta();
                    std::cout << alignment.get_seq_identity() << std::endl;
                }

                // compute partition function
                aln.reset_beams();
                // aln.prob_set1();
                aln.prob_set2(seq_idnty);
                aln.set_prob_accm(pfs[k1].prob_accm, pfs[k2].prob_accm);
                aln.compute_inside(false, 100, verbose_state == VerboseState::DEBUG);
                aln.compute_outside(verbose_state == VerboseState::DEBUG);
                aln.compute_coincidence_probabilities(verbose_state == VerboseState::DEBUG);
                if (verbose_state == VerboseState::DEBUG) {
                    aln.print_alpha_beta();
                }
            }
        }

        if (itr < max_itr) {
            // fold step
            for (TurboPartition &pf : pfs) {
                pf.reset_beams();
                pf.compute_inside();
                pf.compute_outside();
                // std::cout << "Ensemble Energy: " << pf.get_ensemble_energy() << std::endl;
                // pf.print_alpha_beta();
            }
            for (TurboPartition &pf : pfs) {
                pf.compute_bpp_matrix();
                pf.calc_prob_accm();
            }
            this->ext_info_cache.clear();
        } else {
            // std::cout << std::endl << std::endl;
            for (TurboPartition &pf : pfs) {
                std::cout << ">" << pf.sequence->id << std::endl;
                std::cout << pf.get_threshknot_structure() << std::endl;
            }
        }
    }
}
