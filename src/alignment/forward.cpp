// src/alignment/forward.cpp
#include <linearx/alignment/config.hpp>
#include <linearx/alignment/linear_align.hpp>

// #include <boost/heap/fibonacci_heap.hpp>
#include <boost/heap/d_ary_heap.hpp>
#include <algorithm>
#include <cmath>
#include <limits>
#include <queue>
#include <unordered_map>
#include <unordered_set>

using namespace std;
using namespace linearx::utils;

template <typename T>
void LinearAlignmentInterface<T>::compute_inside_Astar(const bool use_lazy_outside, const unsigned beam_size, bool verbose_output) {
    const auto start_time = chrono::high_resolution_clock::now();
    if (verbose_output) {
        fprintf(stderr, "[LinearAlignment] Running Inside Algorithm (A*):\n");
    }

    struct HeapEntry {
        value_type priority;  // -LOG_MUL(alpha, beta) for A* heuristic
        int i;
        int j;
    };

    struct EntryCompare {
        bool operator()(const HeapEntry& lhs, const HeapEntry& rhs) const {
            return lhs.priority > rhs.priority;  // min-heap
        }
    };

    // using Heap = boost::heap::fibonacci_heap<HeapEntry, boost::heap::compare<EntryCompare>>;
    using Heap = boost::heap::d_ary_heap<HeapEntry, boost::heap::compare<EntryCompare>, boost::heap::mutable_<true>, boost::heap::arity<2>>;
    using Handle = typename Heap::handle_type;
    using Key = std::pair<int, int>;

    Heap frontier;
    std::unordered_map<Key, Handle, PairHash> handle_map;

    const int dest_seq1 = seq1.length() + 1;
    const int dest_seq2 = seq2.length() + 1;

    auto try_push_state = [&](int i, int j, HStateType h, value_type alpha) -> bool {
        HState* saved = get_saved_state(h, i, j);
        
        // Compute priority
        value_type priority;
        if (saved && saved->beta > linearx::math::LOG_ZERO) {
            priority = -LOG_MUL(alpha, saved->beta);  // A* heuristic
        } else {
            priority = -alpha;  // No heuristic available
        }
        
        // Update heap and handle map
        const Key key{i, j};
        auto it = handle_map.find(key);
        if (it == handle_map.end()) {
            Handle hdl = frontier.push({priority, i, j});
            handle_map.emplace(key, hdl);
        } else {
            Handle hdl = it->second;
            if (priority < (*hdl).priority) {
                (*hdl).priority = priority;
                (*hdl).i = i;
                (*hdl).j = j;
                frontier.update(hdl);
            }
        }
        return true;
    };

    // Push starting state (0, 0)
    HState& start_state = get_beams(HStateType::ALN)[0][{0, 0}];
    try_push_state(0, 0, HStateType::ALN, start_state.alpha);

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
                
                // Check if we should explore
                HState* next_state = check_state_AStar(nh, ni, nj);
                if (!next_state) {
                    return;
                }

                value_type edge_weight = get_trans_emit_prob(ni, nj, nh, h);
                if (nh == HStateType::ALN) {
                    edge_weight = LOG_MUL(edge_weight, get_match_score(i, j));
                }

                const value_type candidate_alpha = LOG_MUL(state.alpha, edge_weight);
                if (candidate_alpha <= next_state->alpha) {
                    return;  // alpha is already better and beta stays the same
                }

                next_state->alpha = candidate_alpha;
                try_push_state(ni, nj, nh, candidate_alpha);  // Check + push in one step
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
        fprintf(stderr, "  - Score: %.4f\n", best_score);
    }

    log.total_score.first = best_score;
    log.best_exec_time.first = execution_time;
    log.states_pruned.first = 0;
    log.effective_beam_size.first = beam_size;
}

template <typename T>
void LinearAlignmentInterface<T>::compute_inside_Astar_lazy(const bool use_lazy_outside, const unsigned beam_size, bool verbose_output) {
    run_beam_size_ = beam_size;
    const auto start_time = chrono::high_resolution_clock::now();
    if (verbose_output) {
        fprintf(stderr, "[LinearAlignment] Running Inside Algorithm (A* lazy):\n");
    }
    (void)use_lazy_outside;  // kept for API parity

    struct HeapEntry {
        value_type priority;  // -LOG_MUL(alpha, beta) for A* heuristic
        int i;
        int j;
    };

    struct EntryCompare {
        bool operator()(const HeapEntry& lhs, const HeapEntry& rhs) const {
            return lhs.priority > rhs.priority;  // min-heap
        }
    };

    using Heap = std::priority_queue<HeapEntry, std::vector<HeapEntry>, EntryCompare>;
    using Key = std::pair<int, int>;

    Heap frontier;
    // Track the best (lowest) priority seen so far for each (i, j).
    // This emulates decrease-key behavior: only the best entry for a key is considered;
    // worse duplicates become stale and are skipped when popped.
    std::unordered_map<Key, value_type, PairHash> best_priority;

    const int dest_seq1 = seq1.length() + 1;
    const int dest_seq2 = seq2.length() + 1;

    auto push_state = [&](int i, int j, HStateType h, value_type alpha) {
        HState* saved = get_saved_state(h, i, j);

        value_type priority;
        if (saved && saved->beta > linearx::math::LOG_ZERO) {
            priority = -LOG_MUL(alpha, saved->beta);  // A* heuristic
        } else {
            priority = -alpha;  // No heuristic available
        }

        const Key key{i, j};
        auto it = best_priority.find(key);
        if (it == best_priority.end() || priority < it->second) {
            best_priority[key] = priority;
            frontier.push({priority, i, j});
        }
    };

    // Push starting state (0, 0)
    HState& start_state = get_beams(HStateType::ALN)[0][{0, 0}];
    push_state(0, 0, HStateType::ALN, start_state.alpha);

    value_type best_score = linearx::math::LOG_ZERO;

    while (!frontier.empty()) {
        HeapEntry top = frontier.top();
        frontier.pop();

        const int i = top.i;
        const int j = top.j;
        const Key key{i, j};
        // Skip stale entries that are not the current best for this (i, j)
        auto best_it = best_priority.find(key);
        if (best_it != best_priority.end() && top.priority > best_it->second) {
            continue;
        }

        const int s = i + j;

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

                HState* next_state = check_state(nh, ni, nj);
                if (!next_state) {
                    return;
                }

                value_type edge_weight = get_trans_emit_prob(ni, nj, nh, h);
                if (nh == HStateType::ALN) {
                    edge_weight = LOG_MUL(edge_weight, get_match_score(i, j));
                }

                const value_type candidate_alpha = LOG_MUL(state.alpha, edge_weight);
                if (candidate_alpha <= next_state->alpha) {
                    return;  // alpha is already better and beta stays the same
                }

                next_state->alpha = candidate_alpha;
                push_state(ni, nj, nh, candidate_alpha);
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
        fprintf(stderr, "  - Score: %.4f\n", best_score);
    }

    log.total_score.first = best_score;
    log.best_exec_time.first = execution_time;
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
        log.best_exec_time.first = execution_time;
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
