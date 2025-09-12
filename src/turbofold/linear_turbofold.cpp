#include <linearx/turbofold/linear_turbofold.hpp>

using namespace std;
using namespace linearx::utils::io;
using namespace linearx::utils;

LinearTurbofold::LinearTurbofold(MultiSeq &msa, EnergyModel &energy_model, float lambda,
                                 float alignment_pruning_threshold, float folding_pruning_threshold,
                                 float threshknot_threshold, float min_helix_size, bool allow_sharp_turn, float alpha1,
                                 float alpha2, float alpha3)
    : msa(msa),
      energy_model(energy_model),
      lambda(lambda),
      alignment_pruning_threshold(alignment_pruning_threshold),
      folding_pruning_threshold(folding_pruning_threshold),
      threshknot_threshold(threshknot_threshold),
      min_helix_size(min_helix_size) {
    unsigned num_pairs = msa.size() * (msa.size() - 1) / 2;
    // reserve space for sequence pairs and sequence identities
    alns.reserve(num_pairs);
    pfs.reserve(msa.size());
    seq_identities.reserve(num_pairs);

    for (int i = 0; i < msa.size(); i++) {
        // create partition function object for each sequence
        pfs.emplace_back(*this, msa[i], energy_model, allow_sharp_turn);

        // enumerate all possible k^2 sequence pairs and create LinearAlign objects
        for (int j = i + 1; j < msa.size(); j++) {
            alns.emplace_back(*this, msa[i], msa[j], alpha1, alpha2, alpha3);
        }
    }
}

int LinearTurbofold::get_seq_pair_index(const int k1, const int k2) const {
    // ensure the indices are within bounds
    if (k1 >= msa.size() || k2 >= msa.size()) {
        throw std::out_of_range("Seq index k out of range");
    }

    // ensure k1 != k2 and k1 < k2 by swapping if necessary
    int a = k1, b = k2;
    if (k1 == k2) {
        throw std::invalid_argument("k1 and k2 must be different");
    }
    if (k1 > k2) {
        std::swap(a, b);
    }

    // number of sequences
    int n = msa.size();

    // compute the index using the triangular indexing formula
    int index = a * (2 * n - a - 1) / 2 + (b - a - 1);

    return index;
}

value_type LinearTurbofold::get_extrinsic_info(const Sequence &x, const int i, const int j) const {
    value_type result = 0.0;
    for (const Sequence &y : msa) {
        if (x.id == y.id) {
            continue;
        }
        const int seq_pair_idx = get_seq_pair_index(x.id, y.id);
        const double seq_idnty = seq_identities[seq_pair_idx];
        const TurboAlignment &aln = alns[seq_pair_idx];
        const TurboPartition &y_pf = pfs[y.id];
        if (x.id == aln.seq1.id) {  // x is the first sequence in the alignment object
            for (const auto itr1 : aln.coinc_prob[i]) {
                for (const auto itr2 : aln.coinc_prob[j]) {
                    const int k = itr1.first;
                    const int l = itr2.first;
                    const value_type aln_prob_ik = itr1.second;
                    const value_type aln_prob_jl = itr2.second;
                    const value_type y_bpp_kl = y_pf.get_bpp(k, l);
                    result += (1.0 - seq_idnty) * y_bpp_kl * aln_prob_ik * aln_prob_jl;
                }
            }
        } else {  // x is the second sequence in the alignment object
            for (const int k : aln.prob_rev_idx[i]) {
                for (const int l : aln.prob_rev_idx[j]) {
                    const value_type aln_prob_ik = aln.get_bpp(k, i);
                    const value_type aln_prob_jl = aln.get_bpp(l, j);
                    const value_type y_bpp_kl = y_pf.get_bpp(k, l);
                    result += (1.0 - seq_idnty) * y_bpp_kl * aln_prob_ik * aln_prob_jl;
                }
            }
        }
    }
    result = LOG(2 * result);
    return result;
}

void LinearTurbofold::run(const int num_itr, const bool use_lazy_outside, const bool use_prev_itr_beta,
                          const bool restrict_search, const bool verbose_output, const bool save_logs,
                          const bool save_probs, const std::string out_dir) {
    unsigned beam_size = 100;
    use_lazy_outside_ = use_lazy_outside;
    use_prev_itr_beta_ = use_prev_itr_beta;
    restrict_search_ = restrict_search;
    TurboFoldLog log(num_itr, use_lazy_outside, use_prev_itr_beta, restrict_search, alignment_pruning_threshold,
                     folding_pruning_threshold);
    std::cout << "Running LinearTurboFold with " << msa.size() << " sequences\n";
    for (curr_itr = 0; curr_itr <= num_itr; ++curr_itr) {
        showProgressBar(curr_itr, num_itr, 1);
        if (verbose_output) {
            fprintf(stderr, "\n------------------------------ LTF Iteration %d / %d ------------------------------\n\n",
                    curr_itr, num_itr);
        }

        // align step
        if (curr_itr > 0) {
            auto start_time = std::chrono::high_resolution_clock::now();
            for (TurboAlignment &aln : alns) {
                // BEST MODE
                const int k1 = aln.seq1.id;
                const int k2 = aln.seq2.id;
                const int aln_pair_index = get_seq_pair_index(k1, k2);
                aln.reset_beams(beam_size);
                aln.use_prob_set1();
                aln.set_prob_accm(pfs[k1].prob_accm, pfs[k2].prob_accm);
                aln.compute_inside<BEST>(beam_size, verbose_output);
                seq_identities[aln_pair_index] = aln.get_alignment().average_seq_identity();

                // PARTITION MODE
                aln.reset_beams(beam_size);
                aln.use_prob_set2(seq_identities[aln_pair_index]);
                aln.set_prob_accm(pfs[k1].prob_accm, pfs[k2].prob_accm);
                aln.compute_inside<PARTITION_INSIDE>(beam_size, verbose_output);
                AlignmentLog aln_log =
                    aln.compute_outside(use_lazy_outside, alignment_pruning_threshold, verbose_output);
                aln_log.seq_identity = seq_identities[aln_pair_index];
                aln_log.effective_beam_size = beam_size;
                log.aln_logs[curr_itr].emplace_back(std::move(aln_log));
                aln.compute_coincidence_probabilities(verbose_output);
                if (restrict_search && curr_itr < num_itr) {
                    aln.save_partition_function(true);
                }
            }
            auto end_time = std::chrono::high_resolution_clock::now();
            log.aln_itr_exec_times[curr_itr] =
                chrono::duration_cast<chrono::milliseconds>(end_time - start_time).count();
            if (save_probs && curr_itr == num_itr) {
                for (TurboAlignment &aln : alns) {
                    aln.dump_coinc_probs(out_dir + "/aln_cp_itr_" + std::to_string(curr_itr));
                }
            }
        }

        // fold step
        auto start_time = std::chrono::high_resolution_clock::now();
        for (TurboPartition &pf : pfs) {
            pf.reset_beams(beam_size);
            pf.compute_inside(Mode::PARTITION_INSIDE, beam_size, verbose_output);
            PartitionLog pf_log = pf.compute_outside(use_lazy_outside, folding_pruning_threshold, verbose_output);
            pf_log.effective_beam_size = beam_size;
            log.pf_logs[curr_itr].emplace_back(std::move(pf_log));
            reset_beam_vector(pf.extinfo_cache, pf.seq_length, beam_size);  // clear extrinsic info cache
        }
        for (TurboPartition &pf : pfs) {
            pf.compute_bpp_matrix(beam_size);
            pf.calc_prob_accm();
            if (restrict_search && curr_itr < num_itr) {
                pf.save_partition_function(true, beam_size);
            }
        }
        auto end_time = std::chrono::high_resolution_clock::now();
        log.pf_itr_exec_times[curr_itr] = chrono::duration_cast<chrono::milliseconds>(end_time - start_time).count();
        if (save_probs && curr_itr == num_itr) {
            for (TurboPartition &pf : pfs) {
                pf.dump_bpp(out_dir + "/pf_bpp_itr_" + std::to_string(curr_itr));
            }
        }
    }  // iterations loop

    // calculate and print total inside and outside times for partition and alignment
    log.itrs_exec_time = 0.0;
    for (const auto &t : log.pf_itr_exec_times) log.itrs_exec_time += t;
    for (const auto &t : log.aln_itr_exec_times) log.itrs_exec_time += t;
    fprintf(stderr, "Total iterations time: %.2fms\n", log.itrs_exec_time);

    if (save_logs || save_probs) {
        std::cout << "\nSaving output to " << out_dir << "\n";
        if (save_logs) log.save_logs(out_dir);
    }
}
