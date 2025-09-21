#include <linearx/turbofold/linear_turbofold.hpp>

using namespace std;
using namespace linearx::utils::io;
using namespace linearx::utils;

bool DEBUG = false;

LinearTurbofold::LinearTurbofold(MultiSeq& msa, EnergyModel& energy_model, const unsigned alignment_beam_size,
                                 const unsigned folding_beam_size, const float alignment_pruning_threshold,
                                 const float folding_pruning_threshold, const float lambda,
                                 const float threshknot_threshold, const float min_helix_size,
                                 const bool allow_sharp_turn, const float alpha1, const float alpha2,
                                 const float alpha3)
    : msa(msa),
      energy_model(energy_model),
      alignment_beam_size(alignment_beam_size),
      folding_beam_size(folding_beam_size),
      alignment_pruning_threshold(alignment_pruning_threshold),
      folding_pruning_threshold(folding_pruning_threshold),
      lambda(lambda),
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

int LinearTurbofold::get_seq_pair_index(int k1, int k2) const {
    int n = msa.size();
    // ensure the indices are within bounds
    if (k1 >= n || k2 >= n) {
        return -1;  // invalid indices
    }
    // ensure k1 != k2 and k1 < k2 by swapping if necessary
    if (k1 == k2) return -1;  // invalid as both indices are the same
    if (k1 > k2) std::swap(k1, k2);
    // compute the index using the triangular indexing formula
    int index = k1 * (2 * n - k1 - 1) / 2 + (k2 - k1 - 1);
    return index;
}

value_type LinearTurbofold::get_extrinsic_info(const Sequence& x, const int i, const int j) const {
    value_type result = 0.0;
    for (const Sequence& y : msa) {
        if (x.id == y.id) {
            continue;
        }
        const int seq_pair_idx = get_seq_pair_index(x.id, y.id);
        const double seq_idnty = seq_identities[seq_pair_idx];
        const TurboAlignment& aln = alns[seq_pair_idx];
        const TurboPartition& y_pf = pfs[y.id];
        if (x.id == aln.seq1.id) {  // x is the first sequence in the alignment object
            for (const auto itr1 : aln.cp[i]) {
                for (const auto itr2 : aln.cp[j]) {
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
                    const value_type aln_prob_ik = aln.get_cp(k, i);
                    const value_type aln_prob_jl = aln.get_cp(l, j);
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
            for (TurboAlignment& aln : alns) {
                // BEST MODE
                const int k1 = aln.seq1.id;
                const int k2 = aln.seq2.id;
                const int aln_pair_index = get_seq_pair_index(k1, k2);
                aln.use_prob_set1();
                aln.set_prob_accm(pfs[k1].prob_accm, pfs[k2].prob_accm);
                aln.compute_inside<BEST>(alignment_beam_size, verbose_output);
                seq_identities[aln_pair_index] = aln.get_alignment().average_seq_identity();
                aln.log.seq_identity = seq_identities[aln_pair_index];
                aln.reset_beams(alignment_beam_size);

                // PARTITION MODE
                aln.use_prob_set2(seq_identities[aln_pair_index]);
                aln.set_prob_accm(pfs[k1].prob_accm, pfs[k2].prob_accm);
                aln.compute_inside<PARTITION_INSIDE>(alignment_beam_size, verbose_output);
                aln.compute_outside(use_lazy_outside, alignment_pruning_threshold, verbose_output);
                aln.compute_coincidence_probabilities(verbose_output);
                if (restrict_search && curr_itr < num_itr) {
                    aln.save_partition_function(true);
                }
                aln.reset_beams(alignment_beam_size);
                log.aln_logs[curr_itr].emplace_back(aln.log);
            }
            auto end_time = std::chrono::high_resolution_clock::now();
            log.aln_itr_exec_times[curr_itr] =
                chrono::duration_cast<chrono::milliseconds>(end_time - start_time).count();
            if (out_dir.length() > 0 && save_probs && curr_itr == num_itr) {
                for (TurboAlignment& aln : alns) {
                    aln.dump_coinc_probs(out_dir + "/cps/itr_" + std::to_string(curr_itr));
                }
            }
        }

        // fold step
        auto start_time = std::chrono::high_resolution_clock::now();
        for (TurboPartition& pf : pfs) {
            pf.compute_inside(Mode::PARTITION_INSIDE, folding_beam_size, verbose_output);
            pf.compute_outside(use_lazy_outside, folding_pruning_threshold, verbose_output);
            reset_beam_vector(pf.extinfo_cache, pf.seq_length, folding_beam_size);  // clear extrinsic info cache
        }
        for (TurboPartition& pf : pfs) {
            pf.compute_bpp_matrix(folding_beam_size);
            pf.calc_prob_accm();
            if (restrict_search && curr_itr > -1 && curr_itr < num_itr) {
                pf.save_partition_function(true, folding_beam_size);
            }
            pf.reset_beams(folding_beam_size);
            log.pf_logs[curr_itr].emplace_back(pf.log);
        }
        auto end_time = std::chrono::high_resolution_clock::now();
        log.pf_itr_exec_times[curr_itr] = chrono::duration_cast<chrono::milliseconds>(end_time - start_time).count();
        if (out_dir.length() > 0 && save_probs && curr_itr == num_itr) {
            for (TurboPartition& pf : pfs) {
                pf.dump_bpp(out_dir + "/bpps/itr_" + std::to_string(curr_itr));
            }
        }
    }  // iterations loop

    // calculate and print total inside and outside times for partition and alignment
    log.itrs_exec_time = 0.0;
    for (const auto& t : log.pf_itr_exec_times) log.itrs_exec_time += t;
    for (const auto& t : log.aln_itr_exec_times) log.itrs_exec_time += t;
    fprintf(stdout, "Total iterations time: %.2fms\n", log.itrs_exec_time);

    // get threshknot structures for each sequence
    auto start_time = std::chrono::high_resolution_clock::now();
    for (TurboPartition& pf : pfs) {
        Structure* struc_ptr = new Structure(pf.get_threshknot_structure(threshknot_threshold, min_helix_size));
        log.seqs_strucs.emplace_back((&(pf.seq)), struc_ptr);
    }
    if (out_dir.length() > 0) {
        std::cout << "Saving output to " << out_dir;
        log.save_threshknot_db(out_dir);
        if (save_logs) log.save_logs(out_dir);
    } else {
        // print to console
        std::cout << "\nThreshknot Structures:\n";
        for (const auto& [seq_ptr, struc_ptr] : log.seqs_strucs) {
            if (!seq_ptr || !struc_ptr) continue;  // safety
            std::cout << ">" << seq_ptr->name << "\n";
            std::cout << struc_ptr->getDotBracket() << "\n";
        }
    }
    auto end_time = std::chrono::high_resolution_clock::now();
    value_type exec_time = chrono::duration_cast<chrono::milliseconds>(end_time - start_time).count();
    fprintf(stdout, "\nOutput writing time: %.2fms\n", exec_time);
}
