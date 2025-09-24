#include <chrono>

#include "turbofold.hpp"

void LinearTurboFold::dump_coinc_probs(const std::string &filepath, const float threshold, std::unordered_map<int, AlnProbs>* coinc_prob, int seqlen) {
    if (coinc_prob == nullptr) {
        throw std::runtime_error(
            "[LinearTurboFold Error] Coincidence probabilities not computed yet! You must run "
            "compute_coincidence_probabilities() first.");
    }

    // open the file for writing
    std::ofstream file(filepath);
    if (!file) {
        std::cerr
            << "[Hint] The directory for the output file may not exist. Please create it before running the method."
            << std::endl;
        throw std::runtime_error("[LinearTurboFold Error] Unable to open the file " + filepath +
                                 " for writing coincidence probabilities.");
    }

    // dump the coincidence probabilities to the file
    for (int i = 0; i < seqlen; ++i) {
        for (const auto &item : coinc_prob[i]) {
            const int j = item.first;
            const double prob = item.second.prob;
            //if (prob < threshold) continue;

            // output i, j, and the probability to the file
            file << "i=" << i 
                    << ", j=" << j 
                    << ", probs=" << std::scientific << std::setprecision(6) << prob 
                    << std::endl;

        }
    }
}

void LinearTurboFold::dump_coinc_aln_probs(const std::string &filepath, const float threshold, std::unordered_map<int, AlnProbs>* coinc_prob, int seqlen) {
    if (coinc_prob == nullptr) {
        throw std::runtime_error(
            "[LinearTurboFold Error] Coincidence probabilities not computed yet! You must run "
            "compute_coincidence_probabilities() first.");
    }

    // open the file for writing
    std::ofstream file(filepath);
    if (!file) {
        std::cerr
            << "[Hint] The directory for the output file may not exist. Please create it before running the method."
            << std::endl;
        throw std::runtime_error("[LinearTurboFold Error] Unable to open the file " + filepath +
                                 " for writing coincidence probabilities.");
    }

    // dump the coincidence probabilities to the file
    for (int i = 0; i < seqlen; ++i) {
        for (const auto &item : coinc_prob[i]) {
            const int j = item.first;
            const double prob = item.second.aln_prob;
            //if (prob < threshold) continue;

            // output i, j, and the probability to the file
            file << "i=" << i 
                    << ", j=" << j 
                    << ", probs=" << std::scientific << std::setprecision(6) << prob 
                    << std::endl;

        }
    }
}

int LinearTurboFold::get_seq_pair_index(const int k1, const int k2) {
    // ensure the indices are within bounds
    if (k1 >= multi_seq->size() || k2 >= multi_seq->size()) {
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
    int n = multi_seq->size();

    // compute the index using the triangular indexing formula
    int index = a * (2 * n - a - 1) / 2 + (b - a - 1);

    return index;
}

double LinearTurboFold::get_extrinsic_info(const Seq &x, const int i, const int j) {
    if (this->itr == 0) return 0.0;

    const auto it = extinf_cache[x.k_id][j].find(i);
    if (it != extinf_cache[x.k_id][j].end()) {
        return it->second;  // cache lookup
    }

    double output = 0;
    for (const Seq &y : *multi_seq) {
        if (x.k_id == y.k_id) continue;

        const int seq_pair_idx = get_seq_pair_index(x.k_id, y.k_id);
        const double seq_idnty = seq_identities[seq_pair_idx];
        const LinearAlign &aln = alns[seq_pair_idx];
        const TurboPartition &y_pf = pfs.at(y.k_id);

        if (x.k_id == aln.sequence1->k_id) {
            for (const auto &itr1 : aln.coinc_prob1[i]) {
                for (const auto &itr2 : aln.coinc_prob1[j]) {
                    const int k = itr1.first;
                    const int l = itr2.first;
                    const double aln_prob_ik = itr1.second.prob;
                    const double aln_prob_jl = itr2.second.prob;
                    const double y_bpp_kl = y_pf.get_bpp(k, l);

                    output += (1 - seq_idnty) * y_bpp_kl * aln_prob_ik * aln_prob_jl;
                }
            }
        } else {
            for (const int k : aln.prob_rev_idx[i]) {
                for (const int l : aln.prob_rev_idx[j]) {
                    const double aln_prob_ik = aln.get_bpp(k, i);
                    const double aln_prob_jl = aln.get_bpp(l, j);
                    const double y_bpp_kl = y_pf.get_bpp(k, l);

                    output += (1 - seq_idnty) * y_bpp_kl * aln_prob_ik * aln_prob_jl;
                }
            }
        }
    }

    output = xlog(2 * output);
    extinf_cache[x.k_id][j][i] = output;  // cache the result
    return output;
}

void LinearTurboFold::reset_extinf_cache() {
    for (int k = 0; k < multi_seq->size(); k++) {
        for (int j = 0; j < multi_seq->at(k).length(); j++) {
            extinf_cache[k][j].clear();
        }
    }
}

void LinearTurboFold::run_phmm_alignment(){
    auto align_start_time = std::chrono::high_resolution_clock::now();
    int align_total_inside_time = 0;
    int align_total_outside_time = 0;
    for (TurboAlign &aln : alns) {
        const int k1 = aln.sequence1->k_id;
        const int k2 = aln.sequence2->k_id;
        const int aln_pair_index = get_seq_pair_index(k1, k2);
        float seq_idnty = seq_identities[aln_pair_index];

        // get the alignments
        aln.reset_beams(use_prev_outside_score || restrict_search ? false : true);
        aln.prob_set1();
        // itr <= 1 ? aln.prob_set1() : aln.prob_set2(seq_idnty);
        if (itr > 0) aln.set_prob_accm(pfs[k1].prob_accm, pfs[k2].prob_accm);
        aln.compute_inside(true, beam_size, verbose_state == VerboseState::DEBUG);
        MultiSeq alignment = aln.get_alignment();
        seq_idnty = alignment.get_seq_identity();    // get the new sequence identity using the new alignment
        seq_identities[aln_pair_index] = seq_idnty;  // store the updated sequence identity

        // if (verbose_state == VerboseState::DEBUG) {
        //     std::cerr << "Alignment: " << k1 << " " << k2 << std::endl;
        //     std::cerr << alignment[0].sequence << std::endl;
        //     std::cerr << alignment[1].sequence << std::endl;
        //     // aln.print_alpha_beta();
        //     std::cerr << "Sequence Identity: " << seq_idnty << std::endl;
        // }

        // compute partition function
        aln.reset_beams(true);
        aln.prob_set2(seq_idnty);
        if (itr > 0) aln.set_prob_accm(pfs[k1].prob_accm, pfs[k2].prob_accm);

        // equivalent to ml_alignment
        auto align_inside_start_time = std::chrono::high_resolution_clock::now();
        aln.compute_inside(false, beam_size, verbose_state == VerboseState::DEBUG);
        auto align_inside_end_time = std::chrono::high_resolution_clock::now();

        auto align_outside_start_time = std::chrono::high_resolution_clock::now();
        aln.compute_outside(use_lazy_outside, alignment_pruning_threshold,
                            verbose_state == VerboseState::DEBUG);
        auto align_outside_end_time = std::chrono::high_resolution_clock::now();

        align_total_inside_time += std::chrono::duration_cast<std::chrono::milliseconds>(
                                        align_inside_end_time - align_inside_start_time)
                                        .count();
        align_total_outside_time += std::chrono::duration_cast<std::chrono::milliseconds>(
                                        align_outside_end_time - align_outside_start_time)
                                        .count();

        // equivalent to cal_align_prob
        aln.compute_coincidence_probabilities(verbose_state == VerboseState::DEBUG);

        if(itr == num_itr + 1) {
            consistency_transform[k1][k2] = aln.coinc_prob1;
            consistency_transform[k2][k1] = aln.coinc_prob2;
            // initialize Probability Consistency Transform with the posterior matrix (also called coincidence prob matrix)
            // Deep copy so CT owns memory, avoiding aliasing/dangling pointers
            // const int len1 = aln.sequence1->length();
            // const int len2 = aln.sequence2->length();

            // if (consistency_transform[k1][k2]) { delete[] consistency_transform[k1][k2]; }
            // if (consistency_transform[k2][k1]) { delete[] consistency_transform[k2][k1]; }

            // std::unordered_map<int, AlnProbs>* copy12 = new std::unordered_map<int, AlnProbs>[len1];
            // for (int i = 0; i < len1; ++i) {
            //     for (const auto &kv : aln.coinc_prob1[i]) {
            //         const int j = kv.first;
            //         if (j >= 0 && j < len2) {
            //             copy12[i][j] = kv.second;
            //         }
            //     }
            // }

            // std::unordered_map<int, AlnProbs>* copy21 = new std::unordered_map<int, AlnProbs>[len2];
            // for (int j = 0; j < len2; ++j) {
            //     for (const auto &kv : aln.coinc_prob2[j]) {
            //         const int i = kv.first;
            //         if (i >= 0 && i < len1) {
            //             copy21[j][i] = kv.second;
            //         }
            //     }
            // }

            // consistency_transform[k1][k2] = copy12;
            // consistency_transform[k2][k1] = copy21;
        }

        // dump_coinc_probs(("./vb_info/" + std::to_string(itr) + "_aln_" + std::to_string(k1) + "_" +
        //                             std::to_string(k2) + ".bpp.txt"), 0.0, aln.coinc_prob1, aln.sequence1->length());
        
        // dump_coinc_aln_probs(("./vb_info/" + std::to_string(itr) + "_aln_prob_" + std::to_string(k1) + "_" +
        //                             std::to_string(k2) + ".bpp.txt"), 0.0, aln.coinc_prob1, aln.sequence1->length());
        

        if (verbose_state == VerboseState::DEBUG) {
            aln.print_alpha_beta();
        } else if (verbose_state == VerboseState::DETAIL) {
            cerr << "dumping consistency Matrix" << endl;
            aln.dump_coinc_probs("./vb_info/" + std::to_string(itr) + "_aln_" + std::to_string(k1) + "_" +
                                std::to_string(k2) + ".bpp.txt", 0.0);
        }

        // save alignment beams for the next iteration
        if (use_prev_outside_score || restrict_search) {
            aln.ab.free();
            aln.ab.save(aln);
        }
    }
    auto align_end_time = std::chrono::high_resolution_clock::now();

    if (verbose_state == VerboseState::DEBUG) {
        std::cerr
            << "[ALIGNMENT] Total Time taken for iteration " << itr << ": "
            << std::chrono::duration_cast<std::chrono::milliseconds>(align_end_time - align_start_time).count()
            << "ms" << std::endl;
        std::cerr << "[ALIGNMENT] Total inside time for iteration " << itr << ": " << align_total_inside_time
                    << "ms" << std::endl;
        std::cerr << "[ALIGNMENT] Total outside time for iteration " << itr << ": " << align_total_outside_time
                    << "ms\n"
                    << std::endl;
    }
}

void LinearTurboFold::run() {
    for (itr = 0; itr <= num_itr; ++itr) {
        // Utility::showProgressBar(itr, num_itr);

        // if (shrink_beam && itr > 1) {
        //     beam_size -= 25;
        //     beam_size = std::max(min_beam_size, beam_size);
        // }

        // align step
        std::cerr << "-------------------------CURRENT ITERATION: " << itr << "-------------------------\n"
                  << "BEAM SIZE: " << beam_size << std::endl;
        if (itr > 0) {
            // Alignment step
            run_phmm_alignment();
        }

        // fold step
        auto fold_start_time = std::chrono::high_resolution_clock::now();
        int fold_total_inside_time = 0;
        int fold_total_outside_time = 0;
        for (TurboPartition &pf : pfs) {
            pf.reset_beams(use_prev_outside_score || restrict_search ? false : true);

            auto fold_inside_start_time = std::chrono::high_resolution_clock::now();
            pf.compute_inside(beam_size);
            auto fold_inside_end_time = std::chrono::high_resolution_clock::now();

            auto fold_outside_start_time = std::chrono::high_resolution_clock::now();
            pf.compute_outside(use_lazy_outside ? folding_pruning_threshold : NEG_INF);
            auto fold_outside_end_time = std::chrono::high_resolution_clock::now();

            fold_total_inside_time +=
                std::chrono::duration_cast<std::chrono::milliseconds>(fold_inside_end_time - fold_inside_start_time)
                    .count();
            fold_total_outside_time +=
                std::chrono::duration_cast<std::chrono::milliseconds>(fold_outside_end_time - fold_outside_start_time)
                    .count();

            // // save partition function beams for the next iteration
            if (use_prev_outside_score || restrict_search) {
                pf.pfb.free();
                pf.pfb.save(pf);
            }

            if (verbose_state == VerboseState::DEBUG) {
                pf.print_alpha_beta();
            }
        }

        for (TurboPartition &pf : pfs) {
            pf.compute_bpp_matrix();
            pf.calc_prob_accm();
            if (verbose_state == VerboseState::DETAIL) {
                pf.dump_bpp("./vb_info/" + std::to_string(itr) + "_pf_" + pf.sequence->id + ".bpp.txt");
            }

            if (itr == num_itr) {
                std::cout << ">" << pf.sequence->id << std::endl;
                std::cout << pf.get_threshknot_structure() << std::endl;
            }
        }
        auto fold_end_time = std::chrono::high_resolution_clock::now();

        if (VerboseState::DEBUG) {
            std::cerr << "\n[FOLDING] Total Time taken for iteration " << itr << ": "
                      << std::chrono::duration_cast<std::chrono::milliseconds>(fold_end_time - fold_start_time).count()
                      << "ms" << std::endl;
            std::cerr << "[FOLDING] Total inside time for iteration " << itr << ": " << fold_total_inside_time << "ms"
                      << std::endl;
            std::cerr << "[FOLDING] Total outside time for iteration " << itr << ": " << fold_total_outside_time
                      << "ms\n"
                      << std::endl;
        }

        // Multi Seq Alignment
        if(itr == num_itr)
        {   
            itr++;
            run_phmm_alignment();

            auto msa_start_time = std::chrono::high_resolution_clock::now();

            std::cerr << "Starting the Multi Sequence Alignment process" << std::endl;
            multiple_sequence_alignment();
            
            auto msa_end_time = std::chrono::high_resolution_clock::now();
            if (verbose_state == VerboseState::DEBUG) {
                std::cerr
                    << "[Multi Sequence Alignment] Total Time taken in iteration " << itr << ": "
                    << std::chrono::duration_cast<std::chrono::milliseconds>(msa_end_time - msa_start_time).count()
                    << "ms" << std::endl;
            }

            // save to file
            std::cout
                << std::endl
                << "[Multi Sequence Alignment] "
                << std::endl;
            multi_alignment->write_fasta(cout);
        }

        reset_extinf_cache();
    }

}

int LinearTurboFold::multiple_sequence_alignment()
{
    unsigned int seq_len = this->multi_seq->size();
    unsigned int hmm_beam = 100; //naukarkr, make this a param
    unsigned int num_consistency_reps = 2; //naukarkr, make this a param

    vector<vector<float>> distances (seq_len, vector<float> (seq_len, 0));
    ProbabilisticModel model;

    // TOCHECK if removing this mul_aln_results works in this func. Directly use cons_trans
    // vector<vector<unordered_map<int, AlnProbs>*>> mul_aln_results;
    // mul_aln_results.resize(seq_len);
    // for(unsigned int i_seq1 = 0; i_seq1 < seq_len; i_seq1++){
    //     mul_aln_results[i_seq1].resize(seq_len);
    //     for(unsigned int i_seq2 = 0; i_seq2 < seq_len; i_seq2++) {
    //         if(i_seq1 == i_seq2) continue;
    //         mul_aln_results[i_seq1][i_seq2] = new unordered_map<int, AlnProbs>[multi_seq->at(i_seq1).length()];
    //     }
    // }

    std::cerr << "[Multi Seq Align] Starting the Max Exp Accuracy calculation for all pairs" << std::endl;
    auto mea_start_time = std::chrono::high_resolution_clock::now();
    for(unsigned int i_seq1 = 0; i_seq1 < seq_len; i_seq1++)
    {
        for(unsigned int i_seq2 = i_seq1+1; i_seq2 < seq_len; i_seq2++)
        {
            if(i_seq1 != i_seq2)
            {   
                size_t seq1length = multi_seq->at(i_seq1).length();
                size_t seq2length = multi_seq->at(i_seq2).length();

                // Maximum Expected Accuracy calculation
                for (int i = 0; i < seq1length; ++i) {
                    auto& cons_trans_i = consistency_transform[i_seq1][i_seq2][i];

                    for (auto it = cons_trans_i.begin(); it != cons_trans_i.end(); ) {
                        if (it->second.aln_prob < 0.01) {
                            int k = it->first;
                            consistency_transform[i_seq2][i_seq1][k].erase(i);
                            it = cons_trans_i.erase(it); // returns iterator to next
                        } else {
                            ++it;
                        }
                    }
                }

                // for(int i = 0; i < seq1length; i++) {
                //     for(auto &item : consistency_transform[i_seq1][i_seq2][i]){
                //         int k = item.first;
                //         double value = item.second.aln_prob;
                //         if (value >= 0.01) {
                //             mul_aln_results[i_seq1][i_seq2][i][k].aln_prob = consistency_transform[i_seq1][i_seq2][i][k].aln_prob;
                //             mul_aln_results[i_seq2][i_seq1][k][i].aln_prob = consistency_transform[i_seq1][i_seq2][i][k].aln_prob;
                //         }
                //     }
                // }

                pair<string*, float> pair_alignment = model.LinearComputeAlignment(hmm_beam, seq1length, 
                    seq2length, consistency_transform[i_seq1][i_seq2]);

                std::cerr << "pair: ("<< i_seq1 <<","<< i_seq2 <<") "
                    << " lengths: "<< seq1length
                    << ","<< seq2length
                    << "  transform-size: "<< consistency_transform.size()
                    << "x"<< consistency_transform[i_seq1].size()
                    //<< " posterior size: "<< consistency_transform[i_seq1][i_seq2]->size()
                    << " Alignment: " << (*pair_alignment.first)
                    << " MEA: " << pair_alignment.second
                    << std::endl;

                std::cerr << "completed running MEA" << std::endl;

                // Guard against inflated MEA by capping per-position reward at 1.0
                float distance = pair_alignment.second / min (seq1length, seq2length);

                // float normalized_mea = std::min(pair_alignment.second, (float)std::min(seq1length, seq2length));
                // float distance = normalized_mea / std::min (seq1length, seq2length);

                distances[i_seq1][i_seq2] = distances[i_seq2][i_seq1] = distance;
                delete pair_alignment.first;
            }
        } // i_seq2 loop.
    } // i_seq1 loop.
    auto mea_end_time = std::chrono::high_resolution_clock::now();
    if (verbose_state == VerboseState::DEBUG) {
        std::cerr << "[Multi Sequence Alignment] Total Time taken for Max expected accuracy calculation: "
                  << std::chrono::duration_cast<std::chrono::milliseconds>(mea_end_time - mea_start_time).count()
                  << "ms" << std::endl;
    }

    std::cerr << "[Multi Seq Align] Probabilistic Consistency Transformation" << std::endl;

    // Probabilistic consistency transformation
    // for (int r = 0; r<num_consistency_reps; r++ ) {
    //     model.LinearMultiConsistencyTransform(multi_seq, consistency_transform);
    // }

    auto pct_start_time = std::chrono::high_resolution_clock::now();
    for (int r = 0; r<num_consistency_reps; r++ ) {
        consistency_transform = model.LinearMultiConsistencyTransform(multi_seq, consistency_transform);
        // now replace the old posterior probs.
        // for (int i = 0; i < seq_len; i++){
        //     for (int j = 0; j < seq_len; j++){
        //         if (i == j) continue;
        //         delete[] mul_aln_results[i][j];
        //         mul_aln_results[i][j] = new_aln_results[i][j];
                
        //         int seq1Length = multi_seq->at(i).length();
        //         dump_coinc_aln_probs(("./vb_info/" + std::to_string(r) + "_msa_pct_alnprobs_" + std::to_string(i) + "_" +
        //         std::to_string(j) + ".bpp.txt"), 0.0, mul_aln_results[i][j], seq1Length);

        //     }
        // }
    }
    auto pct_end_time = std::chrono::high_resolution_clock::now();
    if (verbose_state == VerboseState::DEBUG) {
        std::cerr << "[Multi Sequence Alignment] Total Time taken for Multi consistency transformation: "
                  << std::chrono::duration_cast<std::chrono::milliseconds>(pct_end_time - pct_start_time).count()
                  << "ms" << std::endl;
    }

    std::cerr << "[Multi Seq Align] Computing the Guide Tree" << std::endl;
    std::cerr << "[MSA] Distance matrix:" << std::endl;
    for (unsigned int i = 0; i < seq_len; ++i) {
        std::cerr << "[MSA] D[" << i << "]";
        for (unsigned int j = 0; j < seq_len; ++j) {
            std::cerr << " " << distances[i][j];
        }
        std::cerr << std::endl;
    }
    this->multi_alignment=NULL;
    auto tree_start_time = std::chrono::high_resolution_clock::now();
    TreeNode *tree = TreeNode::ComputeTree(distances); // lisiz, guide tree 
    auto tree_end_time = std::chrono::high_resolution_clock::now();

    if (verbose_state == VerboseState::DEBUG) {
        long long compute_tree_ms = std::chrono::duration_cast<std::chrono::milliseconds>(tree_end_time - tree_start_time).count();
        std::cerr << "[Multi Sequence Alignment] Total Time taken for Compute tree: " << compute_tree_ms << "ms" << std::endl;
    }

    std::cerr << "[Multi Seq Align] Computed the Guide Tree" << std::endl;

    // make the final alignment
    // 1. Guide tree computation
    // 2. Progressive alignment
    // 3. Iterative Refinement
    std::cerr << "[Multi Seq Align] Initiating Final Alignment" << std::endl;

    this->multi_alignment = model.LinearComputeFinalAlignment(tree, this->multi_seq, consistency_transform, model, hmm_beam);
    
    std::cerr << "[Multi Seq Align] Completed Final Alignment" << std::endl;

    // int numSeqs = this->multi_alignment->size();
    // for (int i = 0; i < numSeqs; i++){
    //     mul_aln_results[i].clear();
    // }
    // mul_aln_results.clear();

    delete tree;

    return(0);
}