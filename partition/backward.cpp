#include "./partition.hpp"

void Partition::lazy_outside() {
    bestC[seq.size() - 1].beta = 0;
    bestC[-1].alpha = 0;
    float deviation_threshold = VALUE_MIN;
    // float global_threshold = bestC[seq.size() - 1].alpha - deviation_threshold;

    auto process_beam = [&](const int j, std::unordered_map<int, State> &beam, const StateType type) {
        for (auto &item : beam) {
            const int i = item.first;
            State &state = item.second;
            // if (state.beta > deviation_threshold) {
            backward_update(i, j, state, type, deviation_threshold);
            // }
        }
    };

    for (int j = seq.size() - 1; j > 0; --j) {
        // reverse topological order: C->M->M2->P->Multi
        // if (bestC[j].beta > deviation_threshold) {
        backward_update(0, j, bestC[j], StateType::C, deviation_threshold);
        // }
        process_beam(j, bestM[j], StateType::M);
        process_beam(j, bestM2[j], StateType::M2);
        process_beam(j, bestP[j], StateType::P);
        process_beam(j, bestMulti[j], StateType::Multi);
    }
}

std::pair<int, int> Partition::backward_update(const int i, const int j, State &state, const StateType type,
                                               const float edge_threshold) {

    std::vector<HEdge> incoming_hedges;
    switch (type) {
    case H:
        // get_incoming_hedges_H(j, incoming_hedges);
        break;
    case Multi:
        get_incoming_hedges_Multi(i, j, &incoming_hedges);
        break;
    case P:
        get_incoming_hedges_P(i, j, &incoming_hedges);
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
    if (incoming_hedges.empty())
        return std::make_pair(0, 0);

    std::vector<HEdge *> saved_hedges;
    HEdge *best_hedge = nullptr;

    double best_inside = VALUE_MIN;
    double saved_inside = VALUE_MIN;

    int local_pruned = 0;
    int local_saved = 0;

    for (auto &hedge : incoming_hedges) {
        hedge.weight *= INV_KT;
        double edge_inside = hedge.weight + hedge.left->alpha + (hedge.right ? hedge.right->alpha : 0);
        // if (edge_inside >= edge_threshold) { // keep the edge
        if (true) {
            Fast_LogPlusEquals(saved_inside, edge_inside);
            saved_hedges.push_back(&hedge);
        } else { // prune the edge
            local_pruned++;
            if (saved_hedges.empty() && edge_inside > best_inside) {
                best_inside = edge_inside;
                best_hedge = &hedge;
            }
        }
    }

    float delta; // scaling factor to compensate for edge pruning
    if (!saved_hedges.empty()) {
        delta = state.alpha - saved_inside;
    } else {
        delta = state.alpha - best_inside;
        saved_hedges.push_back(best_hedge);
        local_pruned -= 1; // one more edge recovered
    }

    // delta = 0;
    for (auto &hedge : saved_hedges) {
        State *left = hedge->left, *right = hedge->right;
        if (!right) {
            Fast_LogPlusEquals(left->beta, hedge->weight + state.beta);
            num_beta_updates += 1;
        } else {
            Fast_LogPlusEquals(left->beta, right->alpha + hedge->weight + state.beta);
            Fast_LogPlusEquals(right->beta, left->alpha + hedge->weight + state.beta);
            num_beta_updates += 2;
        }
    }

    local_saved += saved_hedges.size();

    // std::cout << "Backward Update: " << i << " " << j << " : " << type << ", Pruned: " << local_pruned
    //           << ", Delta: " << delta << "(" << state.alpha << ", " << saved_inside << ")" << std::endl;
    // std::cout << local_pruned << " " << local_saved << std::endl;

    incoming_hedges.clear();
    return std::make_pair(local_pruned, local_saved);
}

TraceInfo Partition::get_incoming_hedges_C(const int j, std::vector<HEdge> *incoming_hedges) {
    // no need to check condition
    // if (j < 1) return;
    TraceInfo best_trace, new_trace;
    HEdge best_hedge, new_hedge;

    // C = C + U
    float new_score = 0;
    // incoming_hedges->push_back(HEdge(new_score, &bestC[j - 1], nullptr));
    new_hedge.set(new_score, &bestC[j - 1], nullptr);
    new_trace.set(0, j - 1, j - 1, StateType::C, StateType::C);
    hedges_trace_helper(incoming_hedges, best_hedge, new_hedge, best_trace, new_trace);

    // C = C + P
    for (auto &item : bestP[j]) {
        int i = item.first;
        State &state = item.second;

        int new_score = -energy_model->score_external_paired(i, j, (i > 0 ? seq.at(i - 1) : -1), seq.at(i), seq.at(j),
                                                             (j + 1 < seq.size() ? seq.at(j + 1) : -1), seq.size());
        new_hedge.set(new_score, &state, &bestC[i - 1]);
        new_trace.set(0, j, i - 1, StateType::C, StateType::P);
        hedges_trace_helper(incoming_hedges, best_hedge, new_hedge, best_trace, new_trace);

        // int h = i - 1;
        // if (h >= 0) {
        //     int new_score = -energy_model->score_external_paired(i, j, seq.at(h), seq.at(h + 1), seq.at(j),
        //                                                          (j + 1 < seq.size() ? seq.at(j + 1) : -1),
        //                                                          seq.size());
        //     // incoming_hedges->push_back(HEdge(new_score, &state, &bestC[h]));
        //     new_hedge.set(new_score, &state, &bestC[h]);
        //     new_trace.set(0, j, h, StateType::C, StateType::P);
        // } else {
        //     int new_score = -energy_model->score_external_paired(0, j, -1, seq.at(0), seq.at(j),
        //                                                          (j + 1 < seq.size() ? seq.at(j + 1) : -1),
        //                                                          seq.size());
        //     // incoming_hedges->push_back(HEdge(new_score, &state, nullptr));
        //     new_hedge.set(new_score, &state, nullptr);
        //     new_trace.set(0, j, j, StateType::P, StateType::P);
        // }
        // hedges_trace_helper(incoming_hedges, best_hedge, new_hedge, best_trace, new_trace);
    }

    return best_trace;
}

TraceInfo Partition::get_incoming_hedges_P(const int i, const int j, std::vector<HEdge> *incoming_hedges) {
    TraceInfo best_trace, new_trace;
    HEdge best_hedge, new_hedge;

    // P = H
    auto itr = bestH[j].find(i);
    if (itr != bestH[j].end()) {
        float new_score = 0;
        // incoming_hedges->push_back(HEdge(new_score, &itr->second, nullptr));
        new_hedge.set(new_score, &itr->second, nullptr);
        new_trace.set(i, j, j, StateType::H, StateType::H);
        hedges_trace_helper(incoming_hedges, best_hedge, new_hedge, best_trace, new_trace);
    }

    // P = P (scan left & jump right)
    for (int p = i + 1; (p < j - 1) && (p - i <= MAXLOOPSIZE); ++p) {
        int q = prev_pair[seq[p]][j];
        while ((q != -1) && (p < q) && ((p - i) + (j - q) - 2 <= MAXLOOPSIZE)) {
            itr = bestP[q].find(p);
            if (itr != bestP[q].end()) {
                // std::cout << "P: " << i << " " << j << " " << p << " " << q << std::endl;
                // current shape is: i...p (pair) q...j
                float new_score = -energy_model->score_single_loop(i, j, p, q, seq[i], seq[i + 1], seq[j - 1], seq[j],
                                                                   seq[p - 1], seq[p], seq[q], seq[q + 1]);
                // incoming_hedges->push_back(HEdge(new_score, &itr->second, nullptr));
                new_hedge.set(new_score, &itr->second, nullptr);
                new_trace.set(p, q, q, StateType::P, StateType::P);
                hedges_trace_helper(incoming_hedges, best_hedge, new_hedge, best_trace, new_trace);
            }
            q = prev_pair[seq[p]][q];
        }
    }

    // P = Multi
    itr = bestMulti[j].find(i);
    if (itr != bestMulti[j].end()) {
        float new_score = -energy_model->score_multi(i, j, seq[i], seq[i + 1], seq[j - 1], seq[j], seq.size());
        // incoming_hedges->push_back(HEdge(new_score, &itr->second, nullptr));
        new_hedge.set(new_score, &itr->second, nullptr);
        new_trace.set(i, j, j, StateType::Multi, StateType::Multi);
        hedges_trace_helper(incoming_hedges, best_hedge, new_hedge, best_trace, new_trace);
    }

    return best_trace;
}

TraceInfo Partition::get_incoming_hedges_M(const int i, const int j, std::vector<HEdge> *incoming_hedges) {
    TraceInfo best_trace, new_trace;
    HEdge best_hedge, new_hedge;

    // M = M + U
    if (j > 0) {
        auto itr = bestM[j - 1].find(i);
        if (itr != bestM[j - 1].end()) {
            float new_score = -energy_model->score_multi_unpaired(j - 1, j);
            // incoming_hedges->push_back(HEdge(new_score, &itr->second, nullptr));
            new_hedge.set(new_score, &itr->second, nullptr);
            new_trace.set(i, j - 1, j - 1, StateType::M, StateType::M);
            hedges_trace_helper(incoming_hedges, best_hedge, new_hedge, best_trace, new_trace);
        }
    }

    // M = P
    auto itr = bestP[j].find(i);
    if (itr != bestP[j].end()) {
        float new_score = -energy_model->score_M1(i, j, j, seq[i - 1], seq[i], seq[j],
                                                  (j + 1 < seq.size() ? seq.at(j + 1) : -1), seq.size());
        // incoming_hedges->push_back(HEdge(new_score, &itr->second, nullptr));
        new_hedge.set(new_score, &itr->second, nullptr);
        new_trace.set(i, j, j, StateType::P, StateType::P);
        hedges_trace_helper(incoming_hedges, best_hedge, new_hedge, best_trace, new_trace);
    }

    // M = M2
    itr = bestM2[j].find(i);
    if (itr != bestM2[j].end()) {
        float new_score = 0;
        // incoming_hedges->push_back(HEdge(new_score, &itr->second, nullptr));
        new_hedge.set(new_score, &itr->second, nullptr);
        new_trace.set(i, j, j, StateType::M2, StateType::M2);
        hedges_trace_helper(incoming_hedges, best_hedge, new_hedge, best_trace, new_trace);
    }

    return best_trace;
}

TraceInfo Partition::get_incoming_hedges_M2(const int i, const int j, std::vector<HEdge> *incoming_hedges) {
    TraceInfo best_trace, new_trace;
    HEdge best_hedge, new_hedge;

    // [TODO] sort P?

    // M2 = M + P
    for (auto &item : bestP[j]) {
        int t = item.first;
        State &state = item.second;

        if (t > i) {
            auto itr = bestM[t - 1].find(i);
            if (itr != bestM[t - 1].end()) {
                float new_score = -energy_model->score_M1(t, j, j, seq[t - 1], seq[t], seq[j],
                                                          (j + 1 < seq.size() ? seq.at(j + 1) : -1), seq.size());
                // incoming_hedges->push_back(HEdge(new_score, &itr->second, &state));
                new_hedge.set(new_score, &itr->second, &state);
                new_trace.set(i, j, t - 1, StateType::M, StateType::P);
                hedges_trace_helper(incoming_hedges, best_hedge, new_hedge, best_trace, new_trace);
            }
        }
    }

    return best_trace;
}

TraceInfo Partition::get_incoming_hedges_Multi(const int i, const int j, std::vector<HEdge> *incoming_hedges) {
    TraceInfo best_trace, new_trace;
    HEdge best_hedge, new_hedge;

    int jprev = prev_pair[seq[i]][j];

    // Multi = Multi (jump right)
    auto itr = bestMulti[jprev].find(i);
    if (itr != bestMulti[jprev].end()) {
        float new_score = -energy_model->score_multi_unpaired(jprev, j);
        // incoming_hedges->push_back(HEdge(new_score, &itr->second, nullptr));
        new_hedge.set(new_score, &itr->second, nullptr);
        new_trace.set(i, jprev, jprev, StateType::Multi, StateType::Multi);
        hedges_trace_helper(incoming_hedges, best_hedge, new_hedge, best_trace, new_trace);
    }

    // Multi = M2 (scan left & jump right)
    for (int q = j - 1; q >= jprev; --q) {
        for (int p = i + 1; (p < q) && (p - i <= MAXLOOPSIZE); ++p) {
            auto itr = bestM2[q].find(p);
            if (itr != bestM2[q].end()) {
                // the current shape is i..p M2 q..j
                float new_score =
                    -(energy_model->score_multi_unpaired(i, p - 1) + energy_model->score_multi_unpaired(q, j - 1));
                // incoming_hedges->push_back(HEdge(new_score, &itr->second, nullptr));
                new_hedge.set(new_score, &itr->second, nullptr);
                new_trace.set(p, q, q, StateType::M2, StateType::M2);
                hedges_trace_helper(incoming_hedges, best_hedge, new_hedge, best_trace, new_trace);
            }
        }
    }
    // std::cout << "MULTI!!" << i << " "  << j << " :" << best_trace.i << " " << best_trace.j << " : " <<
    // best_hedge.weight << std::endl;
    return best_trace;
}