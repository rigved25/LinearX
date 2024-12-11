#include "./linear_align.hpp"

MultiSeq LinearAlign::get_alignment() {
    int i = seq1->size();
    int j = seq2->size();
    HStateType h = get_incoming_edges(i + 1, j + 1, HStateType::ALN, nullptr);

    std::string aln1 = "";
    std::string aln2 = "";

    while (i > 0 || j > 0) {
        HStateType h_prev = get_incoming_edges(i, j, h, nullptr);
        switch (h) {
            case HStateType::ALN:
                i -= 1;
                j -= 1;
                aln1 += std::to_string(seq1->at(i));
                aln2 += std::to_string(seq2->at(j));
                break;

            case HStateType::INS1:
                i -= 1;
                aln1 += std::to_string(seq1->at(i));
                aln2 += "-";
                break;

            case HStateType::INS2:
                j -= 1;
                aln1 += "-";
                aln2 += std::to_string(seq2->at(j));
                break;
        }
        h = h_prev;
    }

    std::reverse(aln1.begin(), aln1.end());
    std::reverse(aln2.begin(), aln2.end());

    MultiSeq alignment;
    alignment.add_sequence(Seq(this->sequence1->id, aln1));
    alignment.add_sequence(Seq(this->sequence2->id, aln2));
    return alignment;
}

void LinearAlign::edges_trace_update_helper(std::vector<AlnEdge> *incoming_edges, const AlnEdge &new_edge,
                                            AlnEdge &best_edge, const HStateType &new_trace, HStateType &best_trace) {
    if (incoming_edges) {
        incoming_edges->push_back(new_edge);
    }
    const double new_score = new_edge.weight + (new_edge.prev ? new_edge.prev->alpha : 0);
    const double best_score = best_edge.weight + (best_edge.prev ? best_edge.prev->alpha : 0);
    if (new_score >= best_score) {
        best_edge = new_edge;
        best_trace = new_trace;
    }
}

void LinearAlign::compute_outside(bool use_lazy_outside, bool verbose_output) {
    double deviation_threshold = use_lazy_outside ? DEVIATION_THRESHOLD : POS_INF;
    double global_threshold =
        bestALN[seq_len_sum + 2][{seq1->size() + 1, seq2->size() + 1}].alpha - deviation_threshold;

    unsigned long total_states = 0, states_visited = 0;
    unsigned long edges_saved = 0, edges_pruned = 0;

    auto process_beam = [&](const int s, std::unordered_map<std::pair<int, int>, HState, PairHash> &beam,
                            const HStateType type) {
        for (auto &item : beam) {
            const int i = item.first.first;
            const int j = item.first.second;
            HState &state = item.second;
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
        std::cerr << "[LinearAlign] Running Outside Algorithm:" << std::endl;
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
        fprintf(stderr, "  - Execution Time: %.2f ms (%.2f%% of inside time)\n", outside_execution_time,
                100.0 * outside_execution_time / std::max(inside_execution_time, 1.0));
        fprintf(stderr, "  - Visited Edges: %lu (saved) + %lu (pruned)\n", edges_saved, edges_pruned);
        fprintf(stderr, "  - Visited Nodes (%.2f%%): %lu (visited) / %lu (total)\n\n",
                100.0 * states_visited / total_states, states_visited, total_states);
    }
}

std::pair<int, int> LinearAlign::backward_update(const int i, const int j, const HState &state, const HStateType type,
                                                 const double edge_threshold) {
    if ((i == 0 || j == 0) && type == HStateType::ALN) {
        return std::make_pair(0, 0);
    }
    std::vector<AlnEdge> incoming_hedges;
    get_incoming_edges(i, j, type, &incoming_hedges);
    if (incoming_hedges.empty()) {
        return std::make_pair(0, 0);
    }

    std::vector<AlnEdge *> saved_edges;
    saved_edges.reserve(incoming_hedges.size());
    AlnEdge *best_edge = nullptr;

    double best_inside = LOG(0.0);
    double saved_inside = LOG(0.0);

    int num_local_edges_pruned = 0;
    int num_local_edges_saved = 0;

    for (auto &edge : incoming_hedges) {
        double edge_inside = LOG_MUL(edge.prev->alpha, edge.weight);  // LOG_MUL(a, b) -> a + b
        if (edge_inside > edge_threshold) {                           // keep the edge
            saved_inside = LOG_SUM(saved_inside, edge_inside);        // Fast_LogPlusEquals(saved_inside, edge_inside);
            saved_edges.push_back(&edge);
        } else {  // prune the edge
            num_local_edges_pruned++;
            if (saved_edges.empty() && edge_inside > best_inside) {
                best_inside = edge_inside;
                best_edge = &edge;
            }
        }
    }

    double delta;  // scaling factor to compensate for edge pruning
    if (!saved_edges.empty()) {
        delta = LOG_DIV(state.alpha, saved_inside);  // LOG_DIV(a, b) -> a - b
    } else {
        delta = LOG_DIV(state.alpha, best_inside);  // state.alpha - best_inside
        saved_edges.push_back(best_edge);
        num_local_edges_pruned -= 1;  // one more edge recovered
    }

    for (auto &edge : saved_edges) {
        edge->prev->beta = LOG_SUM(edge->prev->beta, state.beta + edge->weight + delta);
    }

    num_local_edges_saved += saved_edges.size();
    return std::make_pair(num_local_edges_saved, num_local_edges_pruned);
}

HStateType LinearAlign::get_incoming_edges(const int i, const int j, const HStateType type,
                                           std::vector<AlnEdge> *incoming_edges) {
    AlnEdge new_edge, best_edge;
    HStateType new_trace, best_trace;

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
        const auto it = beam[p + q].find({p, q});
        if (it != beam[p + q].end()) {
            double edge_weight = get_trans_emit_prob(i, j, type, h_prev);
            if (use_match_score && type == HStateType::ALN) {
                double match_score = get_match_score(i - 1, j - 1);
                edge_weight = LOG_MUL(edge_weight, match_score);
            }
            new_edge.set(&(it->second), edge_weight);
            new_trace = h_prev;
            edges_trace_update_helper(incoming_edges, new_edge, best_edge, new_trace, best_trace);
        }
    }

    return best_trace;
}

void LinearAlign::compute_coincidence_probabilities(bool verbose_output) {
    // clear the previous matrix
    delete[] coinc_prob;
    delete[] prob_rev_idx;
    coinc_prob = new std::unordered_map<int, double>[seq1->size()];  // reallocate memory
    prob_rev_idx = new std::vector<int>[seq2->size()];               // reallocate memory

    double p_xy = bestALN[seq_len_sum + 2][{seq1->size() + 1, seq2->size() + 1}].alpha;
    for (int s = 0; s <= seq_len_sum; ++s) {
        for (const HStateType h : hstate_types) {
            std::unordered_map<std::pair<int, int>, HState, PairHash> *beam = get_beam(h);
            for (const auto &item : beam[s]) {
                const int i = item.first.first;
                const int j = item.first.second;
                HState &state = beam[s][{i, j}];

                const double prob = LOG_DIV(LOG_MUL(state.alpha, state.beta), p_xy);
                if (prob > -DEVIATION_THRESHOLD && i > 0 && j > 0) {
                    auto [ptr_cprob_ij, inserted] = coinc_prob[i - 1].try_emplace(j - 1, LOG(0.0));
                    ptr_cprob_ij->second = LOG_SUM(ptr_cprob_ij->second, prob);
                }
            }
        }
    }

    unsigned long num_pruned = 0;  // for keeping track of pruned P(i,j)s
    unsigned long num_saved = 0;   // for keeping track of saved P(i,j)s
    for (int i = 0; i < seq1->size(); ++i) {
        for (auto it = coinc_prob[i].begin(); it != coinc_prob[i].end();) {
            const int j = it->first;
            double &prob = it->second;

            if (prob < phmm->get_fam_threshold()) {
                it = coinc_prob[i].erase(it);  // erase and get the next valid iterator
                ++num_pruned;
            } else {
                prob = EXP(prob);
                if (prob > 1.00001) {
                    fprintf(stderr,
                            "[LinearAlign: Warning] BPP value too high, something is wrong! bpp(%d, %d): %.5f\n", i, j,
                            prob);
                }
                prob = std::min(prob, 1.0);
                prob_rev_idx[j].push_back(i);
                ++num_saved;
                ++it;  // move to the next element if not erased
            }
        }
    }

    if (verbose_output) {
        fprintf(stderr, "[LinearAlign] Coincidence Probabilities Computed: %lu (saved) + %lu (pruned)\n", num_saved,
                num_pruned);
    }
}

double LinearAlign::get_bpp(const int i, const int j) const {
    if (coinc_prob == nullptr || prob_rev_idx == nullptr)
        throw std::runtime_error("[LinearAlign: Error] Coincidence probability matrix not computed yet!");

    const auto it = coinc_prob[i].find(j);
    if (it == coinc_prob[i].end()) {
        return 0.0;
    }

    return it->second;  // access the value directly from the iterator
}
