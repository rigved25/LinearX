// src/alignment/forward.cpp
#include <linearx/alignment/config.hpp>
#include <linearx/alignment/linear_align.hpp>

#include <boost/heap/fibonacci_heap.hpp>
#include <algorithm>
#include <cmath>
#include <limits>
#include <unordered_map>
#include <unordered_set>

using namespace std;
using namespace linearx::utils;

template <typename T>
void LinearAlignmentInterface<T>::compute_inside_BEST(const unsigned beam_size, bool verbose_output) {
    run_beam_size_ = beam_size;
    const auto start_time = chrono::high_resolution_clock::now();
    if (verbose_output) {
        fprintf(stderr, "[LinearAlignment] Running Inside Algorithm:\n");
    }

    struct NodeKey {
        int i;
        int j;
        HStateType h;

        bool operator==(const NodeKey& other) const noexcept {
            return i == other.i && j == other.j && h == other.h;
        }
    };

    struct NodeKeyHash {
        size_t operator()(const NodeKey& k) const noexcept {
            size_t h1 = std::hash<int>{}(k.i);
            size_t h2 = std::hash<int>{}(k.j);
            size_t h3 = std::hash<int>{}(static_cast<int>(k.h));
            size_t seed = h1;
            seed ^= h2 + 0x9e3779b9 + (seed << 6) + (seed >> 2);
            seed ^= h3 + 0x9e3779b9 + (seed << 6) + (seed >> 2);
            return seed;
        }
    };

    struct HeapEntry {
        value_type dist;
        int i;
        int j;
        HStateType h;
    };

    struct EntryCompare {
        bool operator()(const HeapEntry& lhs, const HeapEntry& rhs) const {
            return lhs.dist > rhs.dist;
        }
    };

    using Heap = boost::heap::fibonacci_heap<HeapEntry, boost::heap::compare<EntryCompare>>;
    using Handle = typename Heap::handle_type;

    const int goal_i = seq1.length() + 1;
    const int goal_j = seq2.length() + 1;

    unordered_map<NodeKey, value_type, NodeKeyHash> best_alpha;
    unordered_map<NodeKey, value_type, NodeKeyHash> best_cost;
    unordered_map<NodeKey, Handle, NodeKeyHash> handle;
    unordered_set<NodeKey, NodeKeyHash> visited;
    Heap heap;

    const value_type EPS = 1e-9;

    auto relax = [&](int ni, int nj, HStateType nh, value_type new_cost, value_type new_alpha) {
        NodeKey nk{ni, nj, nh};
        auto it_cost = best_cost.find(nk);
        if (it_cost == best_cost.end()) {
            best_cost.emplace(nk, new_cost);
            best_alpha[nk] = new_alpha;
            HeapEntry entry{new_cost, ni, nj, nh};
            Handle hdl = heap.push(entry);
            handle.emplace(nk, hdl);
        } else {
            value_type& stored_cost = it_cost->second;
            value_type& stored_alpha = best_alpha[nk];
            if (new_cost + EPS < stored_cost || (std::fabs(new_cost - stored_cost) <= EPS && new_alpha > stored_alpha)) {
                stored_cost = new_cost;
                stored_alpha = new_alpha;
                Handle& hdl = handle[nk];
                (*hdl).dist = new_cost;
                heap.update(hdl);
            } else if (new_alpha > stored_alpha) {
                stored_alpha = new_alpha;
            }
        }
    };

    auto push_start = [&](HStateType h) {
        NodeKey nk{0, 0, h};
        best_cost.emplace(nk, 0.0);
        best_alpha.emplace(nk, linearx::math::LOG_ONE);
        HeapEntry entry{0.0, 0, 0, h};
        Handle hdl = heap.push(entry);
        handle.emplace(nk, hdl);
    };

    push_start(HStateType::ALN);

    auto restoration = [&](const NodeKey& nk, value_type best_score) {
        auto& beam = get_beams(nk.h)[nk.i + nk.j];
        auto it = beam.find({nk.i, nk.j});
        if (it == beam.end()) {
            auto [new_it, inserted] = beam.emplace(make_pair(nk.i, nk.j), HState());
            it = new_it;
        }
        it->second.alpha = max(it->second.alpha, best_score);
    };

    value_type best_score = linearx::math::LOG_ZERO;
    bool goal_reached = false;
    value_type best_cost_to_goal = std::numeric_limits<value_type>::infinity();

    while (!heap.empty()) {
        HeapEntry top = heap.top();
        heap.pop();
        NodeKey current{top.i, top.j, top.h};

        if (visited.count(current)) {
            continue;
        }
        visited.insert(current);

        const value_type cur_cost = best_cost[current];
        const value_type cur_alpha = best_alpha[current];
        restoration(current, cur_alpha);

        if (current.i == goal_i && current.j == goal_j && current.h == HStateType::ALN) {
            best_score = cur_alpha;
            best_cost_to_goal = cur_cost;
            goal_reached = true;
            break;
        }

        const int s = current.i + current.j;
        const HState& state = get_beams(current.h)[s][{current.i, current.j}];

        auto process_edge = [&](int ni, int nj, HStateType nh, bool is_match) {
            if (ni < 0 || nj < 0) {
                return;
            }
            if (ni > static_cast<int>(seq1.length() + 1) || nj > static_cast<int>(seq2.length() + 1)) {
                return;
            }

            const value_type trans_emit = get_trans_emit_prob(ni, nj, nh, current.h);
            value_type edge_weight = trans_emit;
            if (is_match) {
                edge_weight = LOG_MUL(edge_weight, get_match_score(current.i, current.j));
            }

            const value_type edge_prob = EXP(edge_weight);
            const value_type edge_cost = value_type(1.0) - edge_prob;
            const value_type new_alpha = LOG_MUL(state.alpha, edge_weight);
            if (IS_LOG_ZERO(new_alpha)) {
                return;
            }

            const value_type new_cost = cur_cost + edge_cost;
            if (goal_reached && new_cost >= best_cost_to_goal) {
                return;
            }

            relax(ni, nj, nh, new_cost, new_alpha);
        };

        if (current.i < seq1.length()) {
            process_edge(current.i + 1, current.j, HStateType::INS1, false);
        }
        if (current.j < seq2.length()) {
            process_edge(current.i, current.j + 1, HStateType::INS2, false);
        }
        bool end_check = (current.i == seq1.length() && current.j == seq2.length());
        if ((current.i < seq1.length() && current.j < seq2.length()) || end_check) {
            process_edge(current.i + 1, current.j + 1, HStateType::ALN, true);
        }
    }

    if (!goal_reached) {
        best_score = linearx::math::LOG_ZERO;
    }

    bestALN[goal_i + goal_j][{goal_i, goal_j}].alpha = best_score;

    const auto end_time = chrono::high_resolution_clock::now();
    const value_type execution_time = chrono::duration_cast<chrono::milliseconds>(end_time - start_time).count();
    if (verbose_output) {
        fprintf(stderr, "  - Execution Time: %.3f ms\n", execution_time);
        fprintf(stderr, "  - States Pruned: %d\n", 0);
        fprintf(stderr, "  - Score: %.4f\n", best_score);
    }
    beam_scores.clear();

    log.total_score.first = best_score;
    log.best_exec_time = execution_time;
    log.states_pruned.first = 0;
    log.effective_beam_size.first = beam_size;
}

// template <typename T>
// void LinearAlignmentInterface<T>::compute_inside_BEST2(const unsigned beam_size, bool verbose_output) {
//     run_beam_size_ = beam_size;
//     const auto start_time = chrono::high_resolution_clock::now();
//     if (verbose_output) {
//         fprintf(stderr, "[LinearAlignment] Running Inside Algorithm:\n");
//     }

//     struct HeapEntry {
//         value_type priority;  // -(alpha + beta)
//         int i;
//         int j;
//     };

//     struct EntryCompare {
//         bool operator()(const HeapEntry& lhs, const HeapEntry& rhs) const {
//             return lhs.priority > rhs.priority;  // min-heap
//         }
//     };

//     using Heap = boost::heap::fibonacci_heap<HeapEntry, boost::heap::compare<EntryCompare>>;
//     Heap frontier;

//     const int dest_seq1 = seq1.length() + 1;
//     const int dest_seq2 = seq2.length() + 1;

//     // Push starting state (0, 0)
//     HState& start_state = get_beams(HStateType::ALN)[0][{0, 0}];
//     value_type start_priority = -(start_state.alpha + start_state.beta);
//     frontier.push({start_priority, 0, 0});

//     value_type best_score = linearx::math::LOG_ZERO;

//     while (!frontier.empty()) {
//         HeapEntry top = frontier.top();
//         frontier.pop();

//         const int i = top.i;
//         const int j = top.j;
//         const int s = i + j;

//         // Check goal
//         if (i == dest_seq1 && j == dest_seq2) {
//             auto& dest_beam = get_beams(HStateType::ALN);
//             best_score = dest_beam[s][{i, j}].alpha;
//             break;
//         }

//         // Process all three state types at (i, j)
//         for (const HStateType h : hstate_types) {
//             auto& beam = get_beams(h);
//             auto it = beam[s].find({i, j});
//             if (it == beam[s].end()) {
//                 continue;  // This state type not present at (i,j)
//             }

//             HState& state = it->second;

//             // INS1
//             if (i < seq1.length() && j <= seq2.length()) {
//                 const value_type new_score = get_trans_emit_prob(i + 1, j, HStateType::INS1, h);
//                 HState& next_state = bestINS1[s + 1][{i + 1, j}];
//                 next_state.alpha = max(next_state.alpha, LOG_MUL(state.alpha, new_score));

//                 value_type next_priority = -next_state.alpha;
//                 HState* saved_state = get_saved_state(HStateType::INS1, i + 1, j);
//                 if (saved_state && saved_state->beta > LOG_ZERO) {
//                     next_priority = -LOG_MUL(next_state.alpha, saved_state->beta);
//                 }
//                 frontier.push({next_priority, i + 1, j});

//             }

//             // INS2
//             if (i <= seq1.length() && j < seq2.length()) {
//                 const value_type new_score = get_trans_emit_prob(i, j + 1, HStateType::INS2, h);
//                 HState& next_state = bestINS2[s + 1][{i, j + 1}];
//                 next_state.alpha = max(next_state.alpha, LOG_MUL(state.alpha, new_score));

//                 value_type next_priority = -next_state.alpha;
//                 HState* saved_state = get_saved_state(HStateType::INS2, i, j + 1);
//                 if (saved_state && saved_state->beta > LOG_ZERO) {
//                     next_priority = -LOG_MUL(next_state.alpha, saved_state->beta);
//                 }
//                 frontier.push({next_priority, i, j + 1});
//             }

//             // ALN
//             const bool end_check = (i == seq1.length() && j == seq2.length());
//             if ((i < seq1.length() && j < seq2.length()) || end_check) {
//                 value_type new_score = get_trans_emit_prob(i + 1, j + 1, HStateType::ALN, h);
//                 new_score = LOG_MUL(new_score, get_match_score(i, j));
//                 HState& next_state = bestALN[s + 2][{i + 1, j + 1}];
//                 next_state.alpha = max(next_state.alpha, LOG_MUL(state.alpha, new_score));

//                 value_type next_priority = -next_state.alpha;
//                 HState* saved_state = get_saved_state(HStateType::ALN, i + 1, j + 1);
//                 if (saved_state && saved_state->beta > LOG_ZERO) {
//                     next_priority = -LOG_MUL(next_state.alpha, saved_state->beta);
//                 }
//                 frontier.push({next_priority, i + 1, j + 1});
//             }
//         }
//     }

//     const auto end_time = chrono::high_resolution_clock::now();
//     const value_type execution_time = chrono::duration_cast<chrono::milliseconds>(end_time - start_time).count();
//     if (verbose_output) {
//         fprintf(stderr, "  - Execution Time: %.3f ms\n", execution_time);
//         fprintf(stderr, "  - States Pruned: %d\n", 0);
//         fprintf(stderr, "  - Score: %.4f\n", best_score);
//     }

//     log.total_score.first = best_score;
//     log.best_exec_time = execution_time;
//     log.states_pruned.first = 0;
//     log.effective_beam_size.first = beam_size;
// }

template <typename T>
void LinearAlignmentInterface<T>::compute_inside_BEST3(const unsigned beam_size, bool verbose_output) {
    run_beam_size_ = beam_size;
    const auto start_time = chrono::high_resolution_clock::now();
    if (verbose_output) {
        fprintf(stderr, "[LinearAlignment] Running Inside Algorithm:\n");
    }

    struct HeapEntry {
        value_type priority;  // -(alpha + beta)
        int i;
        int j;
    };

    struct EntryCompare {
        bool operator()(const HeapEntry& lhs, const HeapEntry& rhs) const {
            return lhs.priority > rhs.priority;  // min-heap
        }
    };

    using Heap = boost::heap::fibonacci_heap<HeapEntry, boost::heap::compare<EntryCompare>>;
    using Handle = typename Heap::handle_type;
    using Key = std::pair<int, int>;

    Heap frontier;
    std::unordered_map<Key, Handle, PairHash> handle_map;

    const int dest_seq1 = seq1.length() + 1;
    const int dest_seq2 = seq2.length() + 1;

    auto get_priority = [&](const HStateType type, const int i, const int j, const value_type alpha) -> value_type {
        value_type combined = alpha;
        if (HState* saved = get_saved_state(type, i, j); saved && saved->beta > linearx::math::LOG_ZERO) {
            combined = LOG_MUL(alpha, saved->beta);
        }
        return -combined;
    };

    auto push_or_update = [&](int i, int j, HStateType h, value_type alpha) {
        const Key key{i, j};
        const value_type new_priority = get_priority(h, i, j, alpha);
        auto it = handle_map.find(key);
        if (it == handle_map.end()) {
            Handle hdl = frontier.push({new_priority, i, j});
            handle_map.emplace(key, hdl);
        } else {
            Handle hdl = it->second;
            if (new_priority < (*hdl).priority) {
                (*hdl).priority = new_priority;
                (*hdl).i = i;
                (*hdl).j = j;
                frontier.update(hdl);
            }
        }
    };

    // Push starting state (0, 0)
    HState& start_state = get_beams(HStateType::ALN)[0][{0, 0}];
    push_or_update(0, 0, HStateType::ALN, start_state.alpha);

    value_type best_score = linearx::math::LOG_ZERO;

    while (!frontier.empty()) {
        HeapEntry top = frontier.top();
        frontier.pop();

        const int i = top.i;
        const int j = top.j;
        const int s = i + j;

        handle_map.erase({i, j});

        // Check goal
        if (i == dest_seq1 && j == dest_seq2) {
            best_score = get_beams(HStateType::ALN)[s][{i, j}].alpha;
            break;
        }

        // Process all three state types at (i, j)
        for (const HStateType h : hstate_types) {
            auto& beam = get_beams(h);
            auto it = beam[s].find({i, j});
            if (it == beam[s].end()) {
                continue;
            }

            HState& state = it->second;

            auto propagate = [&](int ni, int nj, HStateType nh) {
                if (ni < 0 || nj < 0 || ni > dest_seq1 || nj > dest_seq2) {
                    return;
                }

                value_type edge_weight = get_trans_emit_prob(ni, nj, nh, h);
                if (nh == HStateType::ALN) {
                    edge_weight = LOG_MUL(edge_weight, get_match_score(i, j));
                }

                const value_type candidate_alpha = LOG_MUL(state.alpha, edge_weight);
                auto& next_beam = get_beams(nh);
                HState& next_state = next_beam[ni + nj][{ni, nj}];
                if (candidate_alpha <= next_state.alpha) {  // alpha is already better and beta stays the same
                    return;
                }
                next_state.alpha = candidate_alpha;
                push_or_update(ni, nj, nh, candidate_alpha);
            };

            if (i < seq1.length() && j <= seq2.length()) {                
                propagate(i + 1, j, HStateType::INS1);
            }

            if (i <= seq1.length() && j < seq2.length()) {
                propagate(i, j + 1, HStateType::INS2);
            }

            const bool end_check = (i == seq1.length() && j == seq2.length());
            if ((i < seq1.length() && j < seq2.length()) || end_check) {
                propagate(i + 1, j + 1, HStateType::ALN);
            }
        }
    }

    const auto end_time = chrono::high_resolution_clock::now();
    const value_type execution_time = chrono::duration_cast<chrono::milliseconds>(end_time - start_time).count();
    if (verbose_output) {
        fprintf(stderr, "  - Execution Time: %.3f ms\n", execution_time);
        fprintf(stderr, "  - States Pruned: %d\n", 0);
        fprintf(stderr, "  - Score: %.4f\n", best_score);
    }

    log.total_score.first = best_score;
    log.best_exec_time = execution_time;
    log.states_pruned.first = 0;
    log.effective_beam_size.first = beam_size;
}

template <typename T>
template <Mode mode>
void LinearAlignmentInterface<T>::compute_inside(const unsigned beam_size, bool verbose_output) {
    run_beam_size_ = beam_size;
    const auto start_time = chrono::high_resolution_clock::now();
    if (verbose_output) {
        fprintf(stderr, "[LinearAlignment] Running Inside Algorithm:\n");

    }
    unsigned long state_pruned = 0;
    for (int s = 0; s <= seq_len_sum; ++s) {
        if (verbose_output) {
            linearx::utils::io::showProgressBar(s, seq_len_sum);
        }
        for (const HStateType h : hstate_types) {
            vector<unordered_map<pair<int, int>, HState, PairHash>>& beam = get_beams(h);
            state_pruned += beam_prune(beam[s]);
            for (auto& item : beam[s]) {
                const int i = item.first.first;
                const int j = item.first.second;
                HState& state = item.second;

                // INS1
                if (i < seq1.length() && j <= seq2.length()) {
                    if constexpr (mode == Mode::BEST) {
                        const value_type new_score = get_trans_emit_prob(i + 1, j, HStateType::INS1, h);
                        HState& next_state = bestINS1[s + 1][{i + 1, j}];
                        next_state.alpha = max(next_state.alpha, LOG_MUL(state.alpha, new_score));
                    } else {
                        HState* next_state = check_state(HStateType::INS1, i + 1, j);
                        if (next_state) {
                            const value_type new_score = get_trans_emit_prob(i + 1, j, HStateType::INS1, h);
                            AlnEdge(&state, new_score).update_state_alpha(*next_state);
                        }
                    }
                }

                // INS2
                if (i <= seq1.length() && j < seq2.length()) {
                    if constexpr (mode == Mode::BEST) {
                        const value_type new_score = get_trans_emit_prob(i, j + 1, HStateType::INS2, h);
                        HState& next_state = bestINS2[s + 1][{i, j + 1}];
                        next_state.alpha = max(next_state.alpha, LOG_MUL(state.alpha, new_score));
                    } else {
                        HState* next_state = check_state(HStateType::INS2, i, j + 1);
                        if (next_state) {
                            const value_type new_score = get_trans_emit_prob(i, j + 1, HStateType::INS2, h);
                            AlnEdge(&state, new_score).update_state_alpha(*next_state);
                        }
                    }
                }

                // ALN
                const bool end_check = (i == seq1.length() && j == seq2.length());
                if ((i < seq1.length() && j < seq2.length() || end_check)) {
                    if constexpr (mode == Mode::BEST) {
                        value_type new_score = get_trans_emit_prob(i + 1, j + 1, HStateType::ALN, h);
                        new_score = LOG_MUL(new_score, get_match_score(i, j));
                        HState& next_state = bestALN[s + 2][{i + 1, j + 1}];
                        next_state.alpha = max(next_state.alpha, LOG_MUL(state.alpha, new_score));
                    } else {
                        HState* next_state = check_state(HStateType::ALN, i + 1, j + 1);
                        if (next_state) {
                            value_type new_score = get_trans_emit_prob(i + 1, j + 1, HStateType::ALN, h);
                            new_score = LOG_MUL(new_score, get_match_score(i, j));
                            AlnEdge(&state, new_score).update_state_alpha(*next_state);
                        }
                    }
                }
            }
        }
    }
    // update/print time stats
    const auto end_time = chrono::high_resolution_clock::now();
    const value_type execution_time = chrono::duration_cast<chrono::milliseconds>(end_time - start_time).count();
    if (verbose_output) {
        fprintf(stderr, "  - Execution Time: %.3f ms\n", execution_time);
        fprintf(stderr, "  - States Pruned: %lu\n", state_pruned);
        fprintf(stderr, "  - Score: %.4f\n", bestALN[seq_len_sum + 2][{seq1.length() + 1, seq2.length() + 1}].alpha);
    }

    // update logs
    log.total_score.first = bestALN[seq_len_sum + 2][{seq1.length() + 1, seq2.length() + 1}].alpha;
    if constexpr (mode == Mode::BEST) {
        log.best_exec_time = execution_time;
    } else {
        log.exec_time.first = execution_time;
    }
    log.states_pruned.first = state_pruned;
    log.effective_beam_size.first = beam_size;
}

// instantiate templates for LinearAlignmentInterface with desired types
#define DECLARE_FUNCS(TYPE)                                                                           \
    template class LinearAlignmentInterface<TYPE>;                                                    \
    template void LinearAlignmentInterface<TYPE>::compute_inside<Mode::BEST>(unsigned, bool); \
    template void LinearAlignmentInterface<TYPE>::compute_inside<Mode::PARTITION_INSIDE>(unsigned, bool);

#define X(TYPE) DECLARE_FUNCS(TYPE)
LA_TEMPLATE_TYPES
#undef X
#undef DECLARE_FUNCS
