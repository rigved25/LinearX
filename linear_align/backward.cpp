#include "./linear_align.hpp"

void LinearAlign::edges_helper(std::vector<AlnEdge> *incoming_edges, const AlnEdge &new_edge) {
    incoming_edges->push_back(new_edge);
}

void LinearAlign::compute_outside() {
    float deviation_threshold = DEVIATION_THRESHOLD;
    float global_threshold = bestALN[seq_len_sum + 2][{seq1->size() + 1, seq2->size() + 1}].alpha - deviation_threshold;

    unsigned long total_states = 0, states_visited = 0;
    unsigned long edges_saved = 0, edges_pruned = 0;

    auto process_beam = [&](const int s, std::unordered_map<std::pair<int, int>, HState, PairHash> &beam,
                            const HStateType type) {
        for (auto &item : beam) {
            const int i = item.first.first;
            const int j = item.first.second;
            HState &state = item.second;
            if (state.beta > -deviation_threshold) {
                float edge_threshold = global_threshold - state.beta;
                std::pair<int, int> local_edges_info = backward_update(i, j, state, type, edge_threshold);
                edges_saved += local_edges_info.first;
                edges_pruned += local_edges_info.second;
                states_visited += 1;
            }
            total_states += 1;
        }
    };

    if (verbose_output) {
        std::cout << "\n[LinearAlign] Running Outside Algorithm:" << std::endl;
    }
    auto start_time = std::chrono::high_resolution_clock::now();
    process_beam(seq_len_sum + 2, bestALN[seq_len_sum + 2], HStateType::ALN);
    for (int s = seq_len_sum; s > 0; --s) {
        if (verbose_output) {
            Utility::showProgressBar(seq_len_sum - s, seq_len_sum - 1);
        }

        // reverse topological order: ALN->INS2->INS1
        process_beam(s, bestALN[s], HStateType::ALN);
        process_beam(s, bestINS2[s], HStateType::INS2);
        process_beam(s, bestINS1[s], HStateType::INS1);
    }
    auto end_time = std::chrono::high_resolution_clock::now();
    outside_execution_time = std::chrono::duration_cast<std::chrono::milliseconds>(end_time - start_time).count();

    if (verbose_output) {
        printf("  - Execution Time: %.2f ms (%.2f%% of inside time)\n", outside_execution_time,
               100.0 * outside_execution_time / std::max(inside_execution_time, 1.0f));
        printf("  - Visited Edges: %lu (saved) + %lu (pruned)\n", edges_saved, edges_pruned);
        printf("  - Visited Nodes (%.2f%%): %lu (visited) / %lu (total)\n", 100.0 * states_visited / total_states,
               states_visited, total_states);
    }
}

std::pair<int, int> LinearAlign::backward_update(const int i, const int j, const HState &state, const HStateType type,
                                                 const float edge_threshold) {

    if ((i == 0 || j == 0) && type == HStateType::ALN)
        return std::make_pair(0, 0);
    std::vector<AlnEdge> incoming_hedges;
    get_incoming_edges(i, j, state, type, &incoming_hedges);
    if (incoming_hedges.empty())
        return std::make_pair(0, 0);

    std::vector<AlnEdge *> saved_edges;
    AlnEdge *best_edge = nullptr;

    double best_inside = VALUE_MIN;
    double saved_inside = VALUE_MIN;

    int num_local_edges_pruned = 0;
    int num_local_edges_saved = 0;

    for (auto &edge : incoming_hedges) {
        double edge_inside = xlog_mul(edge.prev->alpha, edge.weight); // xlog_mul(a, b) -> a + b
        if (edge_inside > edge_threshold) {                           // keep the edge
            saved_inside = xlog_sum(saved_inside, edge_inside);       // Fast_LogPlusEquals(saved_inside, edge_inside);
            saved_edges.push_back(&edge);
        } else { // prune the edge
            num_local_edges_pruned++;
            if (saved_edges.empty() && edge_inside > best_inside) {
                best_inside = edge_inside;
                best_edge = &edge;
            }
        }
    }

    float delta; // scaling factor to compensate for edge pruning
    if (!saved_edges.empty()) {
        delta = xlog_div(state.alpha, saved_inside); // xlog_div(a, b) -> a - b
    } else {
        delta = xlog_div(state.alpha, best_inside); // state.alpha - best_inside
        saved_edges.push_back(best_edge);
        num_local_edges_pruned -= 1; // one more edge recovered
    }

    for (auto &edge : saved_edges) {
        edge->prev->beta = xlog_sum(edge->prev->beta, state.beta + edge->weight + delta);
    }

    num_local_edges_saved += saved_edges.size();
    return std::make_pair(num_local_edges_saved, num_local_edges_pruned);
}

void LinearAlign::get_incoming_edges(const int i, const int j, const HState &state, const HStateType type,
                                     std::vector<AlnEdge> *incoming_edges) {
    AlnEdge new_edge;

    int p = i, q = j;
    switch (type) {
    case HStateType::ALN:
        p = i - 1;
        q = j - 1;
        break;
    case HStateType::INS1:
        p = i - 1;
        break;
    case HStateType::INS2:
        q = j - 1;
        break;
    }

    bool use_match_score = (pm1 != nullptr && pm2 != nullptr);
    for (const HStateType h_prev : hstate_types) {
        std::unordered_map<std::pair<int, int>, HState, PairHash> *beam = get_beam(h_prev);
        if (beam[p + q].find({p, q}) != beam[p + q].end()) {
            double edge_weight = get_trans_emit_prob(i, j, type, h_prev);
            if (type == HStateType::ALN && use_match_score) {
                float match_score = get_match_score(i, j);
                edge_weight = xlog_mul(edge_weight, match_score);
            }
            new_edge.set(&beam[p + q][{p, q}], edge_weight);
            edges_helper(incoming_edges, new_edge);
        }
    }
}

float LinearAlign::get_bpp(int i, int j) {
    float pij = 0;
    for (const HStateType h_prev : hstate_types) {
        std::unordered_map<std::pair<int, int>, HState, PairHash> *beam = get_beam(h_prev);
        if (beam->find({i, j}) != beam->end()) {
            pij += beam[i + j][{i, j}].alpha + beam[i + j][{i, j}].beta;
        }
    }
    pij /= bestALN[seq_len_sum + 2][{seq1->size() + 1, seq2->size() + 1}].alpha;
    return pij;
}