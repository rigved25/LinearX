#include "turbofold.hpp"

float TurboPartition::turbo_beam_prune(std::unordered_map<int, State> &beamstep, const int j) {
    if (beam_size == 0 || beamstep.size() <= beam_size) {
        return VALUE_MIN;
    }
    std::vector<std::pair<double, int>> scores;
    for (auto &item : beamstep) {
        int i = item.first;
        State &cand = item.second;
        int k = i - 1;

        const float ext_info = turbofold->get_extrinsic_info(this->sequence, i, j);
        float new_alpha = (k >= 0 ? bestC[k].alpha : 0.0) + (cand.alpha * std::pow(ext_info, 0.3));
        scores.push_back(std::make_pair(new_alpha, i));
    }
    if (scores.size() <= beam_size) {
        return VALUE_MIN;
    }
    double threshold = quickselect(scores, 0, scores.size() - 1, scores.size() - beam_size);
    for (auto &p : scores) {
        if (p.first < threshold)
            beamstep.erase(p.second);
    }
    return threshold;
}

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
        turbo_beam_prune(bestH[j], j);
        beamstep_H(j, next_pair);
        if (j == 0)
            continue;
        // beam of Multi
        turbo_beam_prune(bestMulti[j], j);
        beamstep_Multi(j, next_pair);
        // beam of P
        turbo_beam_prune(bestP[j], j);
        beamstep_P(j, next_pair);
        // beam of M2
        turbo_beam_prune(bestM2[j], j);
        beamstep_M2(j, next_pair);
        // beam of M
        turbo_beam_prune(bestM[j], j);
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