#include <chrono>

#include "turbofold.hpp"

TurboPartition::~TurboPartition() { turbofold = nullptr; }

double TurboPartition::beam_prune(StateType type, int j, int beam_size) {
    std::unordered_map<int, State> *beamstep;
    std::unordered_map<int, State> *prev_beamstep;

    switch (type) {
        case StateType::H:
            beamstep = &bestH[j];
            prev_beamstep = &pfb.bestH[j];
            break;
        case StateType::P:
            beamstep = &bestP[j];
            prev_beamstep = &pfb.bestP[j];
            break;
        case StateType::Multi:
            beamstep = &bestMulti[j];
            prev_beamstep = &pfb.bestMulti[j];
            break;
        case StateType::M2:
            beamstep = &bestM2[j];
            prev_beamstep = &pfb.bestM2[j];
            break;
        case StateType::M:
            beamstep = &bestM[j];
            prev_beamstep = &pfb.bestM[j];
            break;
        default:
            return LOG_OF_ZERO;
    }

    if (beam_size == 0 || beamstep->size() <= beam_size) {
        return LOG_OF_ZERO;
    }

    std::vector<std::pair<double, int>> scores;
    for (auto it = beamstep->begin(); it != beamstep->end();) {
        int i = it->first;
        State &cand = it->second;
        int k = i - 1;

        double offset_alpha = 0.0;
        if (!use_prev_outside_score || turbofold->itr == 0) {
            offset_alpha = (k >= 0 ? bestC[k].alpha : 0.0);
        } else {
            auto prev_it = prev_beamstep->find(i);
            if (prev_it == prev_beamstep->end()) {
                it = beamstep->erase(it);
                continue;
            }
            offset_alpha = prev_it->second.beta;
        }

        double newalpha = offset_alpha + cand.alpha;
        scores.push_back(std::make_pair(newalpha, i));

        ++it;  // Move to the next iterator
    }

    if (scores.size() <= beam_size) {
        return VALUE_MIN;
    }

    double threshold = Utility::quickselect(scores, 0, scores.size() - 1, scores.size() - beam_size);
    for (auto &p : scores) {
        if (p.first <= threshold) (*beamstep).erase(p.second);
    }

    return threshold;
}

bool TurboPartition::check_state(StateType type, int i, int j) const {
    if (use_prev_outside_score && turbofold->itr > 0) {
        std::unordered_map<int, State> *prev_beamstep;

        switch (type) {
            case StateType::H:
                prev_beamstep = &pfb.bestH[j];
                break;
            case StateType::P:
                prev_beamstep = &pfb.bestP[j];
                break;
            case StateType::Multi:
                prev_beamstep = &pfb.bestMulti[j];
                break;
            case StateType::M2:
                prev_beamstep = &pfb.bestM2[j];
                break;
            case StateType::M:
                prev_beamstep = &pfb.bestM[j];
                break;
            default:
                return true;
        }

        const auto it = prev_beamstep->find(i);
        bool keep_state = false;
        if (it != prev_beamstep->end()) {
            double score = xlog_div(xlog_mul(it->second.alpha, it->second.beta), pfb.total_alpha);
            keep_state = score > turbofold->folding_pruning_threshold;
        }

        return keep_state;
    }

    return true;
}

void TurboPartition::compute_inside(int beam_size) {
    auto start_time = std::chrono::high_resolution_clock::now();
    if (verbose_output) {
        fprintf(stderr, "[LinearPartition] (Seq k_id: %d) Running Inside Algorithm:\n", sequence->k_id);
    }
    for (int j = 0; j < seq->size(); j++) {
        if (verbose_output) {
            Utility::showProgressBar(j, seq->size() - 1);
        }

        // beam of H
        beam_prune(StateType::H, j, beam_size);
        beamstep_H(j, next_pair);
        if (j == 0) continue;
        // beam of Multi
        beam_prune(StateType::Multi, j, beam_size);
        beamstep_Multi(j, next_pair);
        // beam of P
        if (turbofold->itr > 0) {
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
        }
        beam_prune(StateType::P, j, beam_size);
        beamstep_P(j, next_pair);
        // beam of M2
        beam_prune(StateType::M2, j, beam_size);
        beamstep_M2(j, next_pair);
        // beam of M
        beam_prune(StateType::M, j, beam_size);
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
