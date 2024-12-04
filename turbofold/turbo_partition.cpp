#include "turbofold.hpp"

void TurboPartition::compute_inside() {
    auto start_time = std::chrono::high_resolution_clock::now();
    if (verbose_output) {
        fprintf(stderr, "[LinearPartition] (Seq k_id: %d) Running Inside Algorithm:\n", sequence->k_id);
    }
    for (int j = 0; j < seq->size(); j++) {
        if (verbose_output) {
            Utility::showProgressBar(j, seq->size() - 1);
        }

        // beam of H
        beam_prune(bestH[j]);
        beamstep_H(j, next_pair);
        if (j == 0) continue;
        // beam of Multi
        beam_prune(bestMulti[j]);
        beamstep_Multi(j, next_pair);
        // beam of P
        auto it = bestP[j].begin();
        while (it != bestP[j].end()) {
            int i = it->first;
            State &state = it->second;
            const double ext_info = turbofold->get_extrinsic_info(*(this->sequence), i, j);
            if (ext_info <= LOG_OF_ZERO) {
                it = bestP[j].erase(it);  // erase the element and update the iterator
            } else {
                state.alpha = xlog_mul(state.alpha, ext_info * 0.3);  // adjust the weight as needed
                ++it;                                                 // only increment if not erased
            }
        }
        beam_prune(bestP[j]);
        beamstep_P(j, next_pair);
        // beam of M2
        beam_prune(bestM2[j]);
        beamstep_M2(j, next_pair);
        // beam of M
        beam_prune(bestM[j]);
        beamstep_M(j);
        // beam of C
        beamstep_C(j);
    }

    // set viterbi pointer
    this->viterbi = &bestC[seq->size() - 1];

    // update/print time stats
    auto end_time = std::chrono::high_resolution_clock::now();
    inside_execution_time = std::chrono::duration_cast<std::chrono::milliseconds>(end_time - start_time).count();
    if (verbose_output) {
        std::cerr << "  - Execution Time: " << inside_execution_time << " ms" << std::endl;
        if (mode == InsideMode::MFE) {
            std::cerr << "  - MFE (Minimum Free Energy): " << bestC[seq->size() - 1].alpha / -100.0 << " kcal/mol"
                      << std::endl;
            printf("\n%s (%.2f)\n", get_mfe_structure().c_str(), bestC[seq->size() - 1].alpha / -100.0);
        } else {
            std::cerr << "  - Free Energy of the Ensemble: " << get_ensemble_energy() << " kcal/mol" << std::endl;
        }
        std::cerr << std::endl;
    }
}

void TurboPartition::get_incoming_edges_state(const int i, const int j, const StateType type,
                                              std::vector<HEdge> &incoming_hedges) {
    Partition::get_incoming_edges_state(i, j, type, incoming_hedges);

    if (type == StateType::P) {
        for (auto &hedge : incoming_hedges) {
            const double ext_info = turbofold->get_extrinsic_info(*(this->sequence), i, j);
            hedge.weight = xlog_mul(hedge.weight, ext_info * 0.3);
        }
    }
}

void TurboPartition::calc_prob_accm() {
    prob_accm.upstrm.assign(seq->size(), 0);
    prob_accm.dwnstrm.assign(seq->size(), 0);

    for (int j = 0; j < seq->size(); ++j) {
        for (const auto &item : this->bpp[j]) {
            const int i = item.first;
            double prob = item.second;

            prob_accm.upstrm[j] += prob;
            prob_accm.dwnstrm[i] += prob;
        }
    }

    for (int j = 0; j < seq->size(); ++j) {
        if (prob_accm.upstrm[j] > 1.0) {
            prob_accm.upstrm[j] = 1.0;
            prob_accm.dwnstrm[j] = 0.0;
        }
        if (prob_accm.dwnstrm[j] > 1.0) {
            prob_accm.upstrm[j] = 0.0;
            prob_accm.dwnstrm[j] = 1.0;
        }
    }
}
