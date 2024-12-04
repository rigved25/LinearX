#include <sys/time.h>

#include "./partition.hpp"

void Partition::hedges_trace_update_helper(std::vector<HEdge> *incoming_hedges, HEdge &best_hedge,
                                           const HEdge &new_hedge, TraceInfo &best_trace, const TraceInfo &new_trace) {
    if (incoming_hedges) {
        incoming_hedges->push_back(new_hedge);
    }
    const double best_value = best_hedge.weight + (best_hedge.left ? best_hedge.left->alpha : 0) +
                              (best_hedge.right ? best_hedge.right->alpha : 0);
    const double new_value = new_hedge.weight + (new_hedge.left ? new_hedge.left->alpha : 0) +
                             (new_hedge.right ? new_hedge.right->alpha : 0);
    if (new_value >= best_value) {
        best_hedge = new_hedge;
        best_trace = new_trace;
    }
}

void Partition::mfe_backtrack(int i, int j, StateType type, std::string &structure) {
    if (i >= j) return;

    TraceInfo trace;
    switch (type) {
        case H:
            return;
        case Multi:
            trace = get_incoming_hedges_Multi(i, j, nullptr);
            break;
        case P:
            structure[i] = '(';
            structure[j] = ')';
            trace = get_incoming_hedges_P(i, j, nullptr);
            break;
        case M2:
            trace = get_incoming_hedges_M2(i, j, nullptr);
            break;
        case M:
            trace = get_incoming_hedges_M(i, j, nullptr);
            break;
        case C:
            trace = get_incoming_hedges_C(j, nullptr);
            break;
    }

    mfe_backtrack(trace.i, trace.t, trace.type_left, structure);
    mfe_backtrack(trace.t + 1, trace.j, trace.type_right, structure);
}

void Partition::compute_outside(bool use_lazy_outside) {
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
        std::cerr << "\n[LinearPartition] Running Outside Algorithm:" << std::endl;
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

void Partition::get_incoming_edges_state(const int i, const int j, const StateType type,
                                         std::vector<HEdge> &incoming_hedges) {
    switch (type) {
        case H:
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
}

std::pair<int, int> Partition::backward_update(const int i, const int j, State &state, const StateType type,
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
        // hedge.weight *= INV_KT;         [TODO] FIX THIS
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
        if (!right) {
            Fast_LogPlusEquals(left->beta, state.beta + hedge->weight + delta);
        } else {
            Fast_LogPlusEquals(left->beta, state.beta + right->alpha + hedge->weight + delta);
            Fast_LogPlusEquals(right->beta, state.beta + left->alpha + hedge->weight + delta);
        }
    }

    num_local_edges_saved += saved_hedges.size();
    return std::make_pair(num_local_edges_saved, num_local_edges_pruned);
}

TraceInfo Partition::get_incoming_hedges_C(const int j, std::vector<HEdge> *incoming_hedges) {
    TraceInfo best_trace, new_trace;
    HEdge best_hedge, new_hedge;

    // C = C + U
    double new_score = 0;
    new_hedge.set(new_score, &bestC[j - 1], nullptr);
    new_trace.set(0, j - 1, j - 1, StateType::C, StateType::C);
    hedges_trace_update_helper(incoming_hedges, best_hedge, new_hedge, best_trace, new_trace);

    // C = C + P
    for (auto &item : bestP[j]) {
        int i = item.first;
        State &state = item.second;

        int new_score =
            -energy_model->score_external_paired(i, j, (i > 0 ? seq->at(i - 1) : -1), seq->at(i), seq->at(j),
                                                 (j + 1 < seq->size() ? seq->at(j + 1) : -1), seq->size());
        new_hedge.set(new_score, &bestC[i - 1], &state);
        new_trace.set(0, j, i - 1, StateType::C, StateType::P);
        hedges_trace_update_helper(incoming_hedges, best_hedge, new_hedge, best_trace, new_trace);
    }

    return best_trace;
}

TraceInfo Partition::get_incoming_hedges_P(const int i, const int j, std::vector<HEdge> *incoming_hedges) {
    TraceInfo best_trace, new_trace;
    HEdge best_hedge, new_hedge;

    // P = H
    auto itr = bestH[j].find(i);
    if (itr != bestH[j].end()) {
        double new_score = 0;
        new_hedge.set(new_score, &itr->second, nullptr);
        new_trace.set(i, j, j, StateType::H, StateType::H);
        hedges_trace_update_helper(incoming_hedges, best_hedge, new_hedge, best_trace, new_trace);
    }

    // P = P (scan left & jump right)
    for (int p = i + 1; (p < j - 1) && (p - i <= MAXLOOPSIZE); ++p) {
        int q = prev_pair[seq->at(p)][j];
        while ((q != -1) && (p < q) && ((p - i) + (j - q) - 2 <= MAXLOOPSIZE)) {
            itr = bestP[q].find(p);
            if (itr != bestP[q].end()) {
                // current shape is: i...p (pair) q...j
                double new_score =
                    -energy_model->score_single_loop(i, j, p, q, seq->at(i), seq->at(i + 1), seq->at(j - 1), seq->at(j),
                                                     seq->at(p - 1), seq->at(p), seq->at(q), seq->at(q + 1));
                new_hedge.set(new_score, &itr->second, nullptr);
                new_trace.set(p, q, q, StateType::P, StateType::P);
                hedges_trace_update_helper(incoming_hedges, best_hedge, new_hedge, best_trace, new_trace);
            }
            q = prev_pair[seq->at(p)][q];
        }
    }

    // P = Multi
    itr = bestMulti[j].find(i);
    if (itr != bestMulti[j].end()) {
        double new_score =
            -energy_model->score_multi(i, j, seq->at(i), seq->at(i + 1), seq->at(j - 1), seq->at(j), seq->size());
        new_hedge.set(new_score, &itr->second, nullptr);
        new_trace.set(i, j, j, StateType::Multi, StateType::Multi);
        hedges_trace_update_helper(incoming_hedges, best_hedge, new_hedge, best_trace, new_trace);
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
            double new_score = -energy_model->score_multi_unpaired(j - 1, j);
            new_hedge.set(new_score, &itr->second, nullptr);
            new_trace.set(i, j - 1, j - 1, StateType::M, StateType::M);
            hedges_trace_update_helper(incoming_hedges, best_hedge, new_hedge, best_trace, new_trace);
        }
    }

    // M = P
    auto itr = bestP[j].find(i);
    if (itr != bestP[j].end()) {
        double new_score = -energy_model->score_M1(i, j, j, seq->at(i - 1), seq->at(i), seq->at(j),
                                                   (j + 1 < seq->size() ? seq->at(j + 1) : -1), seq->size());
        new_hedge.set(new_score, &itr->second, nullptr);
        new_trace.set(i, j, j, StateType::P, StateType::P);
        hedges_trace_update_helper(incoming_hedges, best_hedge, new_hedge, best_trace, new_trace);
    }

    // M = M2
    itr = bestM2[j].find(i);
    if (itr != bestM2[j].end()) {
        double new_score = 0;
        new_hedge.set(new_score, &itr->second, nullptr);
        new_trace.set(i, j, j, StateType::M2, StateType::M2);
        hedges_trace_update_helper(incoming_hedges, best_hedge, new_hedge, best_trace, new_trace);
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
                double new_score = -energy_model->score_M1(t, j, j, seq->at(t - 1), seq->at(t), seq->at(j),
                                                           (j + 1 < seq->size() ? seq->at(j + 1) : -1), seq->size());
                // incoming_hedges->push_back(HEdge(new_score, &itr->second, &state));
                new_hedge.set(new_score, &itr->second, &state);
                new_trace.set(i, j, t - 1, StateType::M, StateType::P);
                hedges_trace_update_helper(incoming_hedges, best_hedge, new_hedge, best_trace, new_trace);
            }
        }
    }

    return best_trace;
}

TraceInfo Partition::get_incoming_hedges_Multi(const int i, const int j, std::vector<HEdge> *incoming_hedges) {
    TraceInfo best_trace, new_trace;
    HEdge best_hedge, new_hedge;

    int jprev = prev_pair[seq->at(i)][j];
    if (jprev == -1) {
        return best_trace;
    }

    // Multi = Multi (jump right)
    auto itr = bestMulti[jprev].find(i);
    if (itr != bestMulti[jprev].end()) {
        double new_score = -energy_model->score_multi_unpaired(jprev, j);
        new_hedge.set(new_score, &itr->second, nullptr);
        new_trace.set(i, jprev, jprev, StateType::Multi, StateType::Multi);
        hedges_trace_update_helper(incoming_hedges, best_hedge, new_hedge, best_trace, new_trace);
    }

    // Multi = M2 (scan left & jump right)
    for (int q = j - 1; q >= jprev; --q) {
        for (int p = i + 1; (p < q) && (p - i <= MAXLOOPSIZE); ++p) {
            auto itr = bestM2[q].find(p);
            if (itr != bestM2[q].end()) {
                // the current shape is i..p M2 q..j
                double new_score =
                    -(energy_model->score_multi_unpaired(i, p - 1) + energy_model->score_multi_unpaired(q, j - 1));
                new_hedge.set(new_score, &itr->second, nullptr);
                new_trace.set(p, q, q, StateType::M2, StateType::M2);
                hedges_trace_update_helper(incoming_hedges, best_hedge, new_hedge, best_trace, new_trace);
            }
        }
    }

    return best_trace;
}
