#include "turbofold.hpp"

void TurboPartition::compute_inside() {
    auto start_time = std::chrono::high_resolution_clock::now();
    if (verbose_output) {
        std::cout << "\n[LinearPartition] Running Inside Algorithm:" << std::endl;
    }
    for (int j = 0; j < seq->size(); j++) {
        if (verbose_output) {
            Utility::showProgressBar(j, seq->size() - 1);
        }

        // beam of H
        beam_prune(bestH[j]);
        beamstep_H(j, next_pair);
        if (j == 0)
            continue;
        // beam of Multi
        beam_prune(bestMulti[j]);
        beamstep_Multi(j, next_pair);
        // beam of P
        for (auto &item : bestP[j]) {
            int i = item.first;
            State &state = item.second;
            const float ext_info = turbofold->get_extrinsic_info(*(this->sequence), i, j);
            state.alpha += ext_info * 0.3; // adjust the weight as needed
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
        std::cout << "  - Execution Time: " << inside_execution_time << " ms" << std::endl;
        if (mode == InsideMode::MFE) {
            std::cout << "  - MFE (Minimum Free Energy): " << bestC[seq->size() - 1].alpha / -100.0 << " kcal/mol"
                      << std::endl;
            printf("\n%s (%.2f)\n", get_mfe_structure().c_str(), bestC[seq->size() - 1].alpha / -100.0);
        } else {
            std::cout << "  - Free Energy of the Ensemble: " << get_ensemble_energy() << " kcal/mol" << std::endl;
        }
    }
}

void TurboPartition::get_incoming_edges_state(const int i, const int j, const StateType type,
                                              std::vector<HEdge> &incoming_hedges) {
    switch (type) {
    case H:
        break;
    case Multi:
        get_incoming_hedges_Multi(i, j, &incoming_hedges);
        break;
    case P:
        get_incoming_hedges_P(i, j, &incoming_hedges);
        for (auto &hedge : incoming_hedges) {
            const float ext_info = turbofold->get_extrinsic_info(*(this->sequence), i, j);
            hedge.weight += ext_info * 0.3;
        }
        break;
    case M2:
        get_incoming_hedges_M2(i, j, &incoming_hedges);
        break;
    case M:
        get_incoming_hedges_M(i, j, &incoming_hedges);
        break;
    case C:
        get_incoming_hedges_C(j, &incoming_hedges);
        break;
    }
}