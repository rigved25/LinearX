#include <iomanip>
#include <linearx/turbofold/linear_turbofold.hpp>
#include <chrono>
#ifdef _OPENMP
#include <omp.h>
#endif

using namespace std;
using namespace linearx::utils::io;
using namespace linearx::utils;

bool DEBUG = false;

LinearTurbofold::LinearTurbofold(MultiSeq& multi_seq, EnergyModel& energy_model, const unsigned alignment_beam_size,
                                 const unsigned folding_beam_size, const float alignment_pruning_threshold,
                                 const float folding_pruning_threshold, const float lambda,
                                 const float threshknot_threshold, const float min_helix_size,
                                 const bool allow_sharp_turn, const float alpha1, const float alpha2,
                                 const float alpha3, const int num_itr, const bool verbose_output,
                                 const bool save_logs, const bool save_probs, const std::string out_dir)
    : multi_seq(multi_seq),
      energy_model(energy_model),
      alignment_beam_size(alignment_beam_size),
      folding_beam_size(folding_beam_size),
      alignment_pruning_threshold(alignment_pruning_threshold),
      folding_pruning_threshold(folding_pruning_threshold),
      lambda(lambda),
      threshknot_threshold(threshknot_threshold),
      min_helix_size(min_helix_size),
      num_itr_(num_itr),
      verbose_output_(verbose_output),
      save_logs_(save_logs),
      save_probs_(save_probs),
      out_dir_(out_dir) {
    unsigned num_pairs = multi_seq.size() * (multi_seq.size() - 1) / 2;
    // reserve space for sequence pairs and sequence identities
    alns.reserve(num_pairs);
    pfs.reserve(multi_seq.size());
    seq_identities.reserve(num_pairs);
    posterior.resize(multi_seq.size());
    
    for (int i = 0; i < multi_seq.size(); i++) {
        // create partition function object for each sequence
        pfs.emplace_back(*this, multi_seq[i], energy_model, allow_sharp_turn);
        // enumerate all possible k^2 sequence pairs and create LinearAlign objects
        for (int j = i + 1; j < multi_seq.size(); j++) {
            alns.emplace_back(*this, multi_seq[i], multi_seq[j], alpha1, alpha2, alpha3);
        }
        posterior[i].resize(multi_seq.size(), nullptr);
    }
}

LinearTurbofold::~LinearTurbofold() {
    for (int i = 0; i < multi_seq.size(); ++i) {
        for (int j = 0; j < multi_seq.size(); ++j) {
            if (posterior[i][j] != nullptr) {
                delete[] posterior[i][j];
                posterior[i][j] = nullptr;
            }
        }
    }
}

void LinearTurbofold::dump_coinc_probs(const std::string& out_dir, std::unordered_map<int, value_type>* posterior, int seq_len, int k1, int k2) {
    if (posterior == nullptr) {
        throw std::runtime_error(
            "[LinearTurbofold Error] Coincidence probabilities not computed yet! You must run "
            "compute_posterior() first.");
    }

    // create output directory if it doesn't exist
    std::filesystem::create_directories(out_dir);

    std::string filename = out_dir + "/" + std::to_string(k1) + "_" + std::to_string(k2)+ "_posterior.txt";
    // open file for writing
    std::ofstream file(filename);
    if (!file) {
        throw std::runtime_error("[LinearTurbofold Error] Unable to open file " + filename +
                                 " for writing coincidence probabilities.");
    }

    // write all coincidence probabilities to the file
    for (int i = 0; i < seq_len; ++i) {
        for (const auto& item : posterior[i]) {
            const int j = item.first;
            const value_type prob = item.second;
            file << i + 1 << " " << j + 1 << " " << std::fixed << std::setprecision(4) << prob << "\n";
        }
    }
}

int LinearTurbofold::get_seq_pair_index(int k1, int k2) const {
    int n = multi_seq.size();
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
    for (const Sequence& y : multi_seq) {
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

void LinearTurbofold::run_phmm_alignment(TurboFoldLog& log){
    auto start_time = std::chrono::high_resolution_clock::now();
    for (TurboAlignment& aln : alns) {
        // BEST MODE
        const int k1 = aln.seq1.id;
        const int k2 = aln.seq2.id;
        const int aln_pair_index = get_seq_pair_index(k1, k2);
        aln.use_prob_set1();
        aln.set_prob_accm(pfs[k1].prob_accm, pfs[k2].prob_accm);
        aln.compute_inside<BEST>(alignment_beam_size, verbose_output_);
        seq_identities[aln_pair_index] = aln.get_alignment().average_seq_identity();
        aln.log.seq_identity = seq_identities[aln_pair_index];
        aln.reset_beams(alignment_beam_size);

        // PARTITION MODE
        aln.use_prob_set2(seq_identities[aln_pair_index]);
        aln.set_prob_accm(pfs[k1].prob_accm, pfs[k2].prob_accm);
        aln.compute_inside<PARTITION_INSIDE>(alignment_beam_size, verbose_output_);
        aln.compute_outside(use_lazy_outside_, alignment_pruning_threshold, verbose_output_);
        if (curr_itr == num_itr_ + 1) {
            aln.compute_posterior(posterior, verbose_output_);
            dump_coinc_probs(out_dir_ + "/cps/itr_" + std::to_string(curr_itr), posterior[k1][k2], aln.seq1.length(), k1, k2);
        } else {
            aln.compute_coincidence_probabilities(verbose_output_);
        }
        if (restrict_search_ && curr_itr <= num_itr_) {
            aln.save_partition_function(true);
        }
        aln.reset_beams(alignment_beam_size);
        log.aln_logs[curr_itr].emplace_back(aln.log);
    }
    auto end_time = std::chrono::high_resolution_clock::now();
    log.aln_itr_exec_times[curr_itr] = chrono::duration_cast<chrono::milliseconds>(end_time - start_time).count();
    if (out_dir_.length() > 0 && save_probs_) {
        for (TurboAlignment& aln : alns) {
            aln.dump_coinc_probs(out_dir_ + "/cps/itr_" + std::to_string(curr_itr));
        }        
    }
}

void LinearTurbofold::multiple_sequence_alignment(TurboFoldLog& log, unsigned int beam_size){
    unsigned int seq_len = multi_seq.size();

    vector<vector<value_type>> distances (seq_len, vector<value_type> (seq_len, 0));
    ProbabilisticModel model(out_dir_);
    
    fprintf(stderr, "[Multi Seq Align] Starting the Max Exp Accuracy calculation for all pairs\n");

    auto mea_start_time = std::chrono::high_resolution_clock::now();

    // Parallelize the MEA calculations using OpenMP
    #pragma omp parallel for schedule(dynamic) if(seq_len > 1)
    for (int i = 0; i < seq_len; i++) {
        // // Progress reporting (only from master thread to avoid interleaved output)
        // #pragma omp critical
        // {
        //     if (omp_get_thread_num() == 0) {
        //         fprintf(stderr, "[MSA] Processing sequence %d/%d\n", i + 1, seq_len);
        //     }
        // }

        for (int j = i + 1; j < seq_len; j++) {
            size_t seq1length = multi_seq.at(i).length();
            size_t seq2length = multi_seq.at(j).length();

            pair<string*, value_type> pair_alignment = model.LinearComputeAlignment(beam_size, seq1length, seq2length, posterior[i][j]);

            // Guard against inflated MEA by capping per-position reward at 1.0
            value_type distance = pair_alignment.second / static_cast<value_type>(std::min(seq1length, seq2length));

            // Each thread writes to independent locations in the distance matrix
            distances[i][j] = distances[j][i] = distance;

            // Thread-safe print result
            // #pragma omp critical
            // {
            //     fprintf(stderr, "[MSA] D[%d][%d] = %f\n", i, j, distance);
            // }

            delete pair_alignment.first;
        }
    }

    fprintf(stderr, "[Multi Seq Align] Distance matrix computation completed\n");
    auto mea_end_time = std::chrono::high_resolution_clock::now();
    log.msa_mea_calc_time = std::chrono::duration_cast<std::chrono::milliseconds>(mea_end_time - mea_start_time).count();

    // For now, set total MSA time to MEA calculation time
    // When full MSA is implemented, this will include process tree and iterative refinement
    log.msa_total_time = log.msa_mea_calc_time;

    fprintf(stderr, "[Multi Seq Align] Max Exp Accuracy calculation for all pairs completed (%.2f ms)\n", log.msa_mea_calc_time);

    for (int r = 0; r<model.num_consistency_reps_; r++ ) {
        posterior = model.LinearMultiConsistencyTransform(multi_seq, posterior);
        // for(int x = 0; x < multi_seq.size(); x++){
        //     for(int y = 0; y < multi_seq.size(); y++){
        //         if(x == y) continue;
        //         dump_coinc_probs(out_dir_ + "/cps/ct_" + std::to_string(r), posterior[x][y], multi_seq.at(x).length(), x, y);
        //     }
        // }

    }
    
    TreeNode *tree = TreeNode::ComputeTree(distances);
    
    multi_alignment = model.LinearComputeFinalAlignment(tree, &multi_seq, posterior, beam_size);

    // delete tree;

    // return alignment;

}

void LinearTurbofold::run() {
    TurboFoldLog log(num_itr_, use_lazy_outside_, use_prev_itr_beta_, restrict_search_, alignment_pruning_threshold,
                     folding_pruning_threshold);
    std::cout << "Running LinearTurboFold with " << multi_seq.size() << " sequences\n";
    for (curr_itr = 0; curr_itr <= num_itr_; ++curr_itr) {
        showProgressBar(curr_itr, num_itr_+1, 1);
        if (verbose_output_) {
            fprintf(stderr, "\n------------------------------ LTF Iteration %d / %d ------------------------------\n\n",
                    curr_itr, num_itr_+1);
        }

        // align step
        if (curr_itr > 0) {
            run_phmm_alignment(log);
        }

        // fold step
        auto start_time = std::chrono::high_resolution_clock::now();
        for (TurboPartition& pf : pfs) {
            pf.compute_inside(Mode::PARTITION_INSIDE, folding_beam_size, verbose_output_);
            pf.compute_outside(use_lazy_outside_, folding_pruning_threshold, verbose_output_);
            reset_beam_vector(pf.extinfo_cache, pf.seq_length, folding_beam_size);  // clear extrinsic info cache
        }
        for (TurboPartition& pf : pfs) {
            pf.compute_bpp_matrix(folding_beam_size);
            pf.calc_prob_accm();
            if (restrict_search_ && curr_itr > -1 && curr_itr < num_itr_) {
                pf.save_partition_function(true, folding_beam_size);
            }
            pf.reset_beams(folding_beam_size);
            log.pf_logs[curr_itr].emplace_back(pf.log);
        }
        auto end_time = std::chrono::high_resolution_clock::now();
        log.pf_itr_exec_times[curr_itr] = chrono::duration_cast<chrono::milliseconds>(end_time - start_time).count();
        if (out_dir_.length() > 0 && save_probs_ && curr_itr == num_itr_) {
            for (TurboPartition& pf : pfs) {
                pf.dump_bpp(out_dir_ + "/bpps/itr_" + std::to_string(curr_itr));
            }
        }

        // Multi Seq Alignment
        if(curr_itr == num_itr_)
        {   
            curr_itr++;
            showProgressBar(curr_itr, num_itr_+1, 1);

            if (verbose_output_) {
                fprintf(stderr, "\n------------------------------ LTF Iteration %d / %d ------------------------------\n\n",
                        curr_itr, num_itr_+1);
            }
            
            run_phmm_alignment(log);

            auto msa_start_time = std::chrono::high_resolution_clock::now();
            multiple_sequence_alignment(log);
            auto msa_end_time = std::chrono::high_resolution_clock::now();

            (void)msa_start_time; (void)msa_end_time; // silence unused if not used elsewhere

            multi_alignment->write_fasta(out_dir_ + "/probcons_msa.fasta");
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
    if (out_dir_.length() > 0) {
        std::cout << "Saving output to " << out_dir_;
        log.save_threshknot_db(out_dir_);
        if (save_logs_) log.save_logs(out_dir_);
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
