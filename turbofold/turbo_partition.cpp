#include "turbofold.hpp"

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
            return xlog(0.0);
    }

    if (beam_size == 0 || beamstep->size() <= beam_size) {
        return xlog(0.0);
    }
    std::vector<std::pair<double, int>> scores;
    for (auto &item : *beamstep) {
        int i = item.first;
        State &cand = item.second;
        int k = i - 1;

        double offset_alpha = 0.0;
        if (!use_prev_outside_score || turbofold->itr < 2) {
            offset_alpha = (k >= 0 ? bestC[k].alpha : double(0.0));
        } else if (prev_beamstep->find(i) != prev_beamstep->end()) {
            offset_alpha = prev_beamstep->at(i).beta;
        } else {
            offset_alpha = NEG_INF;
        }

        double newalpha = offset_alpha + cand.alpha;
        scores.push_back(std::make_pair(newalpha, i));
    }
    if (scores.size() <= beam_size) {
        return xlog(0.0);
    }
    double threshold = Utility::quickselect(scores, 0, scores.size() - 1, scores.size() - beam_size);
    for (auto &p : scores) {
        if (p.first < threshold) (*beamstep).erase(p.second);
    }
    return threshold;
}

bool TurboPartition::check_state(const StateType type, const int i, const int j) {
    // return true;

    // if (!pfb.has_data || j < 2) {
    //     return true;
    // }

    // switch (type) {
    //     case StateType::H:
    //         return true;
    //         // if (pfb.bestH[j].find(i) == pfb.bestH[j].end()) {
    //         //     // for (auto item : pfb.bestH[j]) {
    //         //     //     printf("H[%d][%d]: %.5f\n", item.first, j, item.second.alpha);
    //         //     // }
    //         //     return false;
    //         // }
    //     case StateType::P:
    //         if (pfb.bestP[j].find(i) == pfb.bestP[j].end()) {
    //             // std::cout << "Skipped P" << std::endl;
    //             return false;
    //         }
    //     case StateType::Multi:
    //         if (pfb.bestMulti[j].find(i) == pfb.bestMulti[j].end()) {
    //             // std::cout << "Skipped Multi" << std::endl;
    //             return false;
    //         }
    //     case StateType::M2:
    //         if (pfb.bestM2[j].find(i) == pfb.bestM2[j].end()) {
    //             // std::cout << "Skipped M2" << std::endl;
    //             return false;
    //         }
    //     case StateType::M:
    //         if (pfb.bestM[j].find(i) == pfb.bestM[j].end()) {
    //             // std::cout << "Skipped M" << std::endl;
    //             return false;
    //         }
    //     default:
    //         return true;
    // }
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

void TurboPartition::compute_outside(bool use_lazy_outside) {
    double deviation_threshold = use_lazy_outside ? DEVIATION_THRESHOLD : POS_INF;
    double global_threshold = bestC[seq->size() - 1].alpha - deviation_threshold;

    unsigned long total_states = 0, states_visited = 0;
    unsigned long edges_saved = 0, edges_pruned = 0;

    auto process_beam = [&](const int j, std::unordered_map<int, State> &beam, const StateType type) {
        for (auto &item : beam) {
            const int i = item.first;
            State &state = item.second;

            // if (state.beta > -deviation_threshold) {    // Major Bug Here (fixed!)
            if (state.alpha + state.beta > global_threshold) {
                double edge_threshold = global_threshold - state.beta;
                std::pair<int, int> local_edges_info = backward_update(i, j, state, type, edge_threshold);
                edges_saved += local_edges_info.first;
                edges_pruned += local_edges_info.second;
                states_visited += 1;
            }
            total_states += 1;
        }
    };

    if (verbose_output) {
        fprintf(stderr, "[LinearPartition] (Seq k_id: %d) Running Outside Algorithm:\n", sequence->k_id);
    }
    auto start_time = std::chrono::high_resolution_clock::now();
    for (int j = seq->size() - 1; j >= 0; --j) {
        // reverse topological order: C->M->M2->P->Multi
        if (verbose_output) {
            Utility::showProgressBar(seq->size() - 1 - j, seq->size() - 1);
        }
        // if (bestC[j].beta > -deviation_threshold) {    // Major Bug Here (fixed!)
        if (bestC[j].alpha + bestC[j].beta > global_threshold) {
            double edge_threshold = global_threshold - bestC[j].beta;
            backward_update(0, j, bestC[j], StateType::C, edge_threshold);
        }
        process_beam(j, bestM[j], StateType::M);
        process_beam(j, bestM2[j], StateType::M2);
        process_beam(j, bestP[j], StateType::P);
        process_beam(j, bestMulti[j], StateType::Multi);
    }
    auto end_time = std::chrono::high_resolution_clock::now();
    outside_execution_time = std::chrono::duration_cast<std::chrono::milliseconds>(end_time - start_time).count();

    if (verbose_output) {
        fprintf(stderr, "  - Execution Time: %.2f ms (%.2f%% of inside time)\n", outside_execution_time,
                100.0 * outside_execution_time / std::max(inside_execution_time, 1.0));
        fprintf(stderr, "  - Visited Edges: %lu (saved) + %lu (pruned)\n", edges_saved, edges_pruned);
        fprintf(stderr, "  - Visited Nodes (%.2f%%): %lu (visited) / %lu (total)\n",
                100.0 * states_visited / total_states, states_visited, total_states);
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

std::pair<int, int> TurboPartition::backward_update(const int i, const int j, State &state, const StateType type,
                                                    const double edge_threshold) {
    std::vector<HEdge> incoming_hedges;
    get_incoming_edges_state(i, j, type, incoming_hedges);
    if (incoming_hedges.empty()) {
        return std::make_pair(0, 0);
    }

    std::vector<HEdge *> saved_hedges;
    HEdge *best_hedge = nullptr;

    double best_inside = xlog(0.0);
    double saved_inside = xlog(0.0);

    int num_local_edges_pruned = 0;
    int num_local_edges_saved = 0;

    for (auto &hedge : incoming_hedges) {
        double edge_inside = hedge.weight + hedge.left->alpha + (hedge.right ? hedge.right->alpha : 0);

        if (edge_inside > edge_threshold) {  // keep the edge
            Fast_LogPlusEquals(saved_inside, edge_inside);
            saved_hedges.push_back(&hedge);
        } else {  // prune the edge
            num_local_edges_pruned++;
            if (saved_hedges.empty() && edge_inside > best_inside) {
                best_inside = edge_inside;
                best_hedge = &hedge;
            }
        }
    }

    double delta;  // scaling factor to compensate for edge pruning
    if (!saved_hedges.empty()) {
        delta = state.alpha - saved_inside;
    } else {
        delta = state.alpha - best_inside;
        saved_hedges.push_back(best_hedge);
        num_local_edges_pruned -= 1;  // one more edge recovered
    }

    for (auto &hedge : saved_hedges) {
        State *left = hedge->left, *right = hedge->right;

        // double ext_info = 0.0;
        // if (type == StateType::P) {
        //     ext_info = turbofold->get_extrinsic_info(*(this->sequence), i, j);
        //     // delta += ext_info * 0.3;
        // }

        double new_score = state.beta + hedge->weight;

        if (!right) {
            Fast_LogPlusEquals(left->beta, new_score + delta);
        } else {
            Fast_LogPlusEquals(left->beta, right->alpha + new_score + delta);
            Fast_LogPlusEquals(right->beta, left->alpha + new_score + delta);
        }
    }

    num_local_edges_saved += saved_hedges.size();
    return std::make_pair(num_local_edges_saved, num_local_edges_pruned);
}

void TurboPartition::compute_bpp_matrix() {
    // clear the existing bpp matrix
    delete[] this->bpp;

    // compute the new bpp matrix
    this->bpp = new std::unordered_map<int, double>[seq->size()];  // reallocate memory
    for (int j = 0; j < seq->size(); ++j) {
        for (auto &item : bestP[j]) {
            const int i = item.first;
            State &state = item.second;

            // const double ext_info = turbofold->get_extrinsic_info(*(this->sequence), i, j);
            // double offset = ext_info * 0.3;
            // state.beta = xlog_div(state.beta, offset);

            double prob = xlog_div(xlog_mul(state.alpha, state.beta), viterbi->alpha);
            if (prob > -DEVIATION_THRESHOLD) {
                prob = xexp(prob);  // Convert log prob to regular prob
                if (prob > 1.001) {
                    fprintf(stderr,
                            "[LinearPartition Warning] BPP value too high, something is wrong! bpp(%d, %d): %.5f\n", i,
                            j, prob);
                }
                prob = std::min(prob, 1.0);  // Clamp the probability to [0, 1]
                this->bpp[j][i] = prob;      // Set the bpp value
            }
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
