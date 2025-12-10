#include "linearx/turbofold/utils/ProbabilisticModel.hpp"
#include <linearx/turbofold/utils/random.h>
#include <algorithm>
#include <vector>
#include <chrono>
#include <boost/heap/d_ary_heap.hpp>
#include <limits>


/////////////////////////////////////////////////////////////////
// ProbabilisticModel::LinearComputeAlignment()
//
// Computes an alignment based on the given posterior matrix using
// a push-based DP approach (like Viterbi). This finds the maximum
// summing path through the posterior matrix. The final alignment
// is returned as a pair consisting of:
//    (1) a string (e.g., XXXBBXXXBBBBBBYYYYBBB) where X's and
//        denote insertions in one of the two sequences and
//        B's denote that both sequences are present (i.e.
//        matches).
//    (2) a value_type indicating the sum achieved
/////////////////////////////////////////////////////////////////

pair<string *, value_type> ProbabilisticModel::LinearComputeAlignment(int hmmBeam, int seq1Length, int seq2Length, const unordered_map<int, value_type>* posterior) {
    vector<unordered_map<pair<int, int>, MEAState, PairHash>> beam(seq1Length + seq2Length + 3);
    size_t states_explored = 0;

    // Initialize
    beam[0][{0, 0}] = MEAState(0.0, '\0');

    for (int s = 0; s <= seq1Length + seq2Length; ++s) {
        // Beam pruning
        if (beam[s].size() > hmmBeam + skip_beam_prune_) beam_prune(beam[s], hmmBeam);
        // if (beam[s].size() > hmmBeam) beam_prune(beam[s], hmmBeam);


        // Process each state in current step
        for (auto& item : beam[s]) {
            int i = item.first.first;
            int j = item.first.second;
            value_type score = item.second.score;
            states_explored++;

            // ALN
            if (i + 1 <= seq1Length && j + 1 <= seq2Length) {
                int next_i = i + 1;
                int next_j = j + 1;
                auto row_it = posterior[i].find(j); // 0-indexed, we do next_i-1 == i
                if (row_it != posterior[i].end()) {
                    value_type new_score = score + row_it->second;
                    MEAState& next_state = beam[s + 2][{next_i, next_j}];
                    // Update always: it is always the best
                    next_state.score = new_score;
                    next_state.traceback = 'D'; // Diagonal (match)
                }
            }

            // INS1
            if (i + 1 <= seq1Length) {
                MEAState& next_state = beam[s + 1][{i + 1, j}];

                // Update if better
                if (score > next_state.score) {
                    next_state.score = score;
                    next_state.traceback = 'U'; // Up (ins1)
                }
            }

            // INS2
            if (j + 1 <= seq2Length) {
                MEAState& next_state = beam[s + 1][{i, j + 1}];

                // Update if better
                if (score > next_state.score) {
                    next_state.score = score;
                    next_state.traceback = 'L'; // Left (ins2)
                }
            }
        }
    }

    // Find the best score at the final position
    value_type total = 0.0;
    unsigned final_step = seq1Length + seq2Length;

    if (beam[final_step].find({seq1Length, seq2Length}) != beam[final_step].end()) {
        total = beam[final_step][{seq1Length, seq2Length}].score;
    }

    // Traceback to reconstruct alignment
    string* alignment = new string();
    unsigned current_i = seq1Length;
    unsigned current_j = seq2Length;
    unsigned current_step = final_step;

    while (current_i > 0 || current_j > 0) {
        auto state_it = beam[current_step].find({current_i, current_j});
        if (state_it == beam[current_step].end()) break;

        char ch = state_it->second.traceback;
        switch (ch) {
        case 'D': // Diagonal (match)
            alignment->push_back('B');
            current_i--;
            current_j--;
            current_step -= 2;
            break;
        case 'U': // Up (ins1)
            alignment->push_back('X');
            current_i--;
            current_step -= 1;
            break;
        case 'L': // Left (ins2)
            alignment->push_back('Y');
            current_j--;
            current_step -= 1;
            break;
        default:
            current_step = 0; // safety
            current_i = current_j = 0;
            break;
        }
    }

    // Reverse the alignment string
    reverse(alignment->begin(), alignment->end());

    fprintf(stderr, "[LinearComputeAlignment] States explored: %zu\n", states_explored);

    return make_pair(alignment, total);
}

pair<string *, value_type> ProbabilisticModel::LinearComputeAlignmentDijkstraCost(int hmmBeam, int seq1Length, int seq2Length, const unordered_map<int, value_type>* posterior){
    struct HeapEntry {
        value_type cost;
        int i;
        int j;

        bool operator>(const HeapEntry& other) const {
            return cost > other.cost;
        }
    };

    vector<unordered_map<pair<int, int>, DijkstraState, PairHash>> beam(seq1Length + seq2Length + 3);
    priority_queue<HeapEntry, vector<HeapEntry>, greater<HeapEntry>> heap;

    auto& start_state = beam[0][{0, 0}];
    start_state.cost = 0.0;
    start_state.traceback = 'S';
    heap.push({0.0, 0, 0});

    const size_t max_states_per_step =
        hmmBeam > 0 ? static_cast<size_t>(hmmBeam) : std::numeric_limits<size_t>::max();

    size_t states_explored = 0;

    while (!heap.empty()) {
        HeapEntry top = heap.top();
        heap.pop();

        int i = top.i;
        int j = top.j;
        int step = i + j;

        auto it = beam[step].find({i, j});
        if (it == beam[step].end() || top.cost > it->second.cost) {
            continue;
        }

        states_explored++;

        if (i == seq1Length && j == seq2Length) {
            break;
        }

        const value_type curr_cost = it->second.cost;

        auto relax = [&](int next_i, int next_j, value_type move_cost, char move) {
            if (next_i > seq1Length || next_j > seq2Length) return;

            int next_step = next_i + next_j;
            auto& next_map = beam[next_step];
            auto key = make_pair(next_i, next_j);
            if (next_map.size() >= max_states_per_step && next_map.find(key) == next_map.end()) {
                return;
            }

            auto insert_result = next_map.emplace(key, DijkstraState());
            DijkstraState& next_state = insert_result.first->second;

            value_type new_cost = curr_cost + move_cost;
            if (new_cost < next_state.cost) {
                next_state.cost = new_cost;
                next_state.traceback = move;
                heap.push({new_cost, next_i, next_j});
            }
        };

        if (i + 1 <= seq1Length && j + 1 <= seq2Length) {
            value_type prob = 0.0;
            auto row_it = posterior[i].find(j);
            if (row_it != posterior[i].end()) {
                prob = row_it->second;
            }
            value_type diag_cost = (prob > 0.0) ? std::max<value_type>(0.0, 1.0 - prob) : 1.0;
            relax(i + 1, j + 1, diag_cost, 'D');
        }

        if (i + 1 <= seq1Length) {
            relax(i + 1, j, 0.5, 'U');
        }

        if (j + 1 <= seq2Length) {
            relax(i, j + 1, 0.5, 'L');
        }
    }

    unsigned final_step = seq1Length + seq2Length;
    auto final_it = beam[final_step].find({seq1Length, seq2Length});
    if (final_it == beam[final_step].end()) {
        return make_pair(new string(), 0.0);
    }

    string* alignment = new string();
    unsigned current_i = seq1Length;
    unsigned current_j = seq2Length;
    unsigned current_step = final_step;

    size_t diag_count = 0;
    while (current_i > 0 || current_j > 0) {
        auto state_it = beam[current_step].find({current_i, current_j});
        if (state_it == beam[current_step].end()) break;

        char ch = state_it->second.traceback;
        switch (ch) {
        case 'D':
            alignment->push_back('B');
            diag_count++;
            current_i--;
            current_j--;
            current_step -= 2;
            break;
        case 'U':
            alignment->push_back('X');
            current_i--;
            current_step -= 1;
            break;
        case 'L':
            alignment->push_back('Y');
            current_j--;
            current_step -= 1;
            break;
        default:
            current_i = current_j = 0;
            current_step = 0;
            break;
        }
    }

    reverse(alignment->begin(), alignment->end());

    value_type total_cost = final_it->second.cost;
    const size_t aln_length = alignment->size();
    const size_t gap_count = aln_length - diag_count;
    // Recover MEA-equivalent reward: |B| + 0.5 * gaps - accumulated cost
    value_type reconstructed_score =
        static_cast<value_type>(diag_count) + 0.5 * static_cast<value_type>(gap_count) - total_cost;

    fprintf(stderr, "[LinearComputeAlignmentDijkstraCost] States explored: %zu\n", states_explored);

    return make_pair(alignment, reconstructed_score);
}

pair<string *, value_type> ProbabilisticModel::LinearComputeAlignmentDijkstraScore(int hmmBeam, int seq1Length, int seq2Length, const unordered_map<int, value_type>* posterior){
    struct HeapEntry {
        value_type cost;
        int i;
        int j;

        bool operator>(const HeapEntry& other) const {
            return cost > other.cost;
        }
    };

    vector<unordered_map<pair<int, int>, DijkstraLazyState, PairHash>> beam(seq1Length + seq2Length + 3);
    priority_queue<HeapEntry, vector<HeapEntry>, greater<HeapEntry>> heap;

    auto& start_state = beam[0][{0, 0}];
    start_state.cost = 0.0;
    start_state.score = 0.0;
    start_state.traceback = '\0';
    heap.push({0.0, 0, 0});

    const size_t max_states_per_step =
        hmmBeam > 0 ? static_cast<size_t>(hmmBeam) : std::numeric_limits<size_t>::max();

    while (!heap.empty()) {
        HeapEntry top = heap.top();
        heap.pop();

        int i = top.i;
        int j = top.j;
        int step = i + j;

        auto it = beam[step].find({i, j});
        if (it == beam[step].end() || top.cost > it->second.cost) {
            continue;
        }

        if (i == seq1Length && j == seq2Length) {
            break;
        }

        const value_type curr_cost = it->second.cost;
        const value_type curr_score = it->second.score;

        auto relax = [&](int next_i, int next_j, value_type move_cost, value_type score_delta, char move) {
            if (next_i > seq1Length || next_j > seq2Length) return;

            int next_step = next_i + next_j;
            auto& next_map = beam[next_step];
            auto key = make_pair(next_i, next_j);
            // if (next_map.size() >= max_states_per_step && next_map.find(key) == next_map.end()) {
            //     return;
            // }

            auto insert_result = next_map.emplace(key, DijkstraLazyState());
            DijkstraLazyState& next_state = insert_result.first->second;

            value_type new_cost = curr_cost + move_cost;
            value_type new_score = curr_score + score_delta;
            if (new_cost < next_state.cost) {
                next_state.cost = new_cost;
                next_state.score = new_score;
                next_state.traceback = move;
                heap.push({new_cost, next_i, next_j});
            }
        };

        if (i + 1 <= seq1Length && j + 1 <= seq2Length) {
            value_type prob = 0.0;
            auto row_it = posterior[i].find(j);
            if (row_it != posterior[i].end()) {
                prob = row_it->second;
            }
            value_type diag_cost = (prob > 0.0) ? std::max<value_type>(0.0, 1.0 - prob) : 1.0;
            value_type diag_reward = (prob > 0.0) ? prob : 0.0;
            relax(i + 1, j + 1, diag_cost, diag_reward, 'D');
        }

        if (i + 1 <= seq1Length) {
            relax(i + 1, j, 0.5, 0.0, 'U');
        }

        if (j + 1 <= seq2Length) {
            relax(i, j + 1, 0.5, 0.0, 'L');
        }
    }

    unsigned final_step = seq1Length + seq2Length;
    auto final_it = beam[final_step].find({seq1Length, seq2Length});
    if (final_it == beam[final_step].end()) {
        return make_pair(new string(), 0.0);
    }

    string* alignment = new string();
    unsigned current_i = seq1Length;
    unsigned current_j = seq2Length;
    unsigned current_step = final_step;

    while (current_i > 0 || current_j > 0) {
        auto state_it = beam[current_step].find({current_i, current_j});
        if (state_it == beam[current_step].end()) break;

        char ch = state_it->second.traceback;
        switch (ch) {
        case 'D':
            alignment->push_back('B');
            current_i--;
            current_j--;
            current_step -= 2;
            break;
        case 'U':
            alignment->push_back('X');
            current_i--;
            current_step -= 1;
            break;
        case 'L':
            alignment->push_back('Y');
            current_j--;
            current_step -= 1;
            break;
        default:
            current_i = current_j = 0;
            current_step = 0;
            break;
        }
    }

    reverse(alignment->begin(), alignment->end());

    value_type total_score = final_it->second.score;
    return make_pair(alignment, total_score);
}

pair<string *, value_type> ProbabilisticModel::LinearComputeAlignmentDijkstraLogs(int hmmBeam, int seq1Length, int seq2Length, const unordered_map<int, value_type>* posterior){
    /* 
        Lazy Dijkstra implementation:
        - Uses min-heap where multiple entries for same node can exist (lazy deletion)
        - When we pop a node, check if it's stale by comparing with beam's cost
        - Beam stores the best known state for each (i,j) including traceback
    */

    struct HeapEntry {
        value_type priority;  // Dijkstra cost (sum of 1-p)
        int i, j;
        
        bool operator>(const HeapEntry& other) const {
            return priority > other.priority;  // min-heap
        }
    };

    // Beam indexed by step (i+j), stores finalized/best-known states
    vector<unordered_map<pair<int, int>, DijkstraLazyState, PairHash>> beam(seq1Length + seq2Length + 3);
    priority_queue<HeapEntry, vector<HeapEntry>, greater<HeapEntry>> heap;
    size_t states_explored = 0;
    size_t max_heap_size = 0;
    const size_t kHeapLogStart = 5000;
    size_t next_heap_log = kHeapLogStart;
    const size_t kBeamLogStart = 5000;
    size_t next_beam_log = kBeamLogStart;
    size_t total_beam_entries = 1; // initial (0,0)
    constexpr size_t kMaxStatesPerStep = 100;

    // Initialize: start at (0,0)
    beam[0][{0, 0}] = DijkstraLazyState(0.0, 0.0, '\0');
    heap.push({0.0, 0, 0});

    auto log_heap_growth = [&](const char* move, int from_i, int from_j, int to_i, int to_j, int target_step) {
        size_t current_size = heap.size();
        if (current_size >= next_heap_log) {
            fprintf(stderr,
                    "[LinearComputeAlignmentDijkstra-Lazy] Heap size %zu reached at step %d via %s (from %d,%d -> %d,%d)\n",
                    current_size, target_step, move, from_i, from_j, to_i, to_j);
            next_heap_log *= 2;
        }
    };

    auto log_beam_growth = [&](int target_step, int to_i, int to_j, const char* move) {
        if (total_beam_entries >= next_beam_log) {
            fprintf(stderr,
                    "[LinearComputeAlignmentDijkstra-Lazy] Beam states %zu reached at step %d via %s (state %d,%d)\n",
                    total_beam_entries, target_step, move, to_i, to_j);
            next_beam_log *= 2;
        }
    };

    auto record_new_state = [&](int target_step, int to_i, int to_j, const char* move) {
        total_beam_entries++;
        log_beam_growth(target_step, to_i, to_j, move);
    };

    // Dijkstra search with lazy deletion
    while (!heap.empty()) {
        // Track maximum heap size
        if (heap.size() > max_heap_size) {
            max_heap_size = heap.size();
        }
        
        HeapEntry top = heap.top();
        heap.pop();
        
        int i = top.i;
        int j = top.j;
        int step = i + j;

        // Check if this is a stale entry
        auto it = beam[step].find({i, j});
        if (it == beam[step].end() || top.priority > it->second.cost) {
            continue;  // Stale entry, skip
        }

        states_explored++;

        // Goal reached
        if (i == seq1Length && j == seq2Length) break;

        value_type curr_cost = top.priority;
        value_type curr_score = it->second.score;

        // ALN: diagonal move (reward from posterior[i][j], which is 0-indexed)
        if (i + 1 <= seq1Length && j + 1 <= seq2Length) {
            int next_i = i + 1;
            int next_j = j + 1;
            int next_step = next_i + next_j;

            auto& next_map = beam[next_step];
            // if (next_map.size() < kMaxStatesPerStep) {
                auto row_it = posterior[i].find(j); // posterior is 0-indexed
                value_type prob = (row_it != posterior[i].end()) ? row_it->second : 0.0;
                value_type diag_cost = (prob > 0.0) ? std::max<value_type>(0.0, 1.0 - prob) : 1.0;
                value_type new_cost = curr_cost + diag_cost;
                value_type new_score = curr_score + ((prob > 0.0) ? prob : 0.0);

                auto key = make_pair(next_i, next_j);
                auto insert_result = next_map.emplace(key, DijkstraLazyState());
                if (insert_result.second) {
                    record_new_state(next_step, next_i, next_j, "ALN");
                }
                auto& next_state = insert_result.first->second;
                if (new_cost < next_state.cost) {
                    next_state.cost = new_cost;
                    next_state.score = new_score;
                    next_state.traceback = 'D';
                    heap.push({new_cost, next_i, next_j});
                    log_heap_growth("ALN", i, j, next_i, next_j, next_step);
                }
            // }
        }

        // Explore neighbors
        // INS1: move down in seq1 (no reward, cost increases by 0.5)
        if (i + 1 <= seq1Length) {
            int next_i = i + 1;
            int next_j = j;
            int next_step = next_i + next_j;
            value_type new_cost = curr_cost + 0.5;

            auto& next_map = beam[next_step];
            // if (next_map.size() < kMaxStatesPerStep) {
                auto key = make_pair(next_i, next_j);
                auto insert_result = next_map.emplace(key, DijkstraLazyState());
                if (insert_result.second) {
                    record_new_state(next_step, next_i, next_j, "INS1");
                }
                auto& next_state = insert_result.first->second;
                if (new_cost < next_state.cost) {
                    next_state.cost = new_cost;
                    next_state.score = curr_score;
                    next_state.traceback = 'U';
                    heap.push({new_cost, next_i, next_j});
                    log_heap_growth("INS1", i, j, next_i, next_j, next_step);
                }
            // }
        }

        // INS2: move right in seq2 (no reward, cost increases by 0.5)
        if (j + 1 <= seq2Length) {
            int next_i = i;
            int next_j = j + 1;
            int next_step = next_i + next_j;
            value_type new_cost = curr_cost + 0.5;

            auto& next_map = beam[next_step];
            // if (next_map.size() < kMaxStatesPerStep) {
                auto key = make_pair(next_i, next_j);
                auto insert_result = next_map.emplace(key, DijkstraLazyState());
                if (insert_result.second) {
                    record_new_state(next_step, next_i, next_j, "INS2");
                }
                auto& next_state = insert_result.first->second;
                if (new_cost < next_state.cost) {
                    next_state.cost = new_cost;
                    next_state.score = curr_score;
                    next_state.traceback = 'L';
                    heap.push({new_cost, next_i, next_j});
                    log_heap_growth("INS2", i, j, next_i, next_j, next_step);
                }
            // }
        }
    }

    // Traceback
    value_type total = 0.0;
    unsigned final_step = seq1Length + seq2Length;
    
    auto final_it = beam[final_step].find({seq1Length, seq2Length});
    if (final_it != beam[final_step].end()) {
        total = final_it->second.score;
    }

    string* alignment = new string();
    unsigned current_i = seq1Length;
    unsigned current_j = seq2Length;
    unsigned current_step = final_step;

    while (current_i > 0 || current_j > 0) {
        auto state_it = beam[current_step].find({current_i, current_j});
        if (state_it == beam[current_step].end()) break;

        char ch = state_it->second.traceback;
        switch (ch) {
        case 'D': // Diagonal (match)
            alignment->push_back('B');
            current_i--;
            current_j--;
            current_step -= 2;
            break;
        case 'U': // Up (ins1)
            alignment->push_back('X');
            current_i--;
            current_step -= 1;
            break;
        case 'L': // Left (ins2)
            alignment->push_back('Y');
            current_j--;
            current_step -= 1;
            break;
        default:
            current_step = 0;
            current_i = current_j = 0;
            break;
        }
    }

    reverse(alignment->begin(), alignment->end());

    fprintf(stderr, "[LinearComputeAlignmentDijkstra-Lazy] States explored: %zu, Max heap size: %zu\n", states_explored, max_heap_size);

    return make_pair(alignment, total);
}

pair<string *, value_type> ProbabilisticModel::LinearComputeAlignmentDijkstra(int hmmBeam, int seq1Length, int seq2Length, const unordered_map<int, value_type>* posterior){
    struct HeapEntry {
        value_type priority;
        int i;
        int j;

        bool operator>(const HeapEntry& other) const {
            return priority > other.priority;
        }
    };

    const int rows = seq1Length;
    const int cols = seq2Length;
    const size_t max_states_per_step =
        hmmBeam > 0 ? static_cast<size_t>(hmmBeam) : std::numeric_limits<size_t>::max();

    vector<unordered_map<pair<int, int>, DijkstraState, PairHash>> beam(rows + cols + 3);
    priority_queue<HeapEntry, vector<HeapEntry>, greater<HeapEntry>> heap;

    beam[0][{0, 0}] = DijkstraState(0.0, 'S');
    heap.push({0.0, 0, 0});

    auto acquire_state = [&](int step, int next_i, int next_j) -> DijkstraState* {
        auto& next_map = beam[step];
        auto key = make_pair(next_i, next_j);
        auto it = next_map.find(key);
        if (it != next_map.end()) {
            return &it->second;
        }
        if (next_map.size() >= max_states_per_step) {
            return nullptr;
        }
        return &next_map.emplace(key, DijkstraState()).first->second;
    };

    auto acquire_state_ALN = [&](int step, int next_i, int next_j) -> DijkstraState* {
        auto& next_map = beam[step];
        auto key = make_pair(next_i, next_j);
        auto it = next_map.find(key);
        if (it != next_map.end()) {
            return &it->second;
        }

        return &next_map.emplace(key, DijkstraState()).first->second;
    };

    size_t states_explored = 0;

    while (!heap.empty()) {
        HeapEntry top = heap.top();
        heap.pop();

        int i = top.i;
        int j = top.j;
        int step = i + j;

        auto it = beam[step].find({i, j});
        if (it == beam[step].end() || top.priority > it->second.cost) {
            continue;
        }

        states_explored++;

        if (i == rows && j == cols) {
            break;
        }

        const value_type curr_cost = it->second.cost;

        // ALN (match)
        if (i + 1 <= rows && j + 1 <= cols) {
            int next_i = i + 1;
            int next_j = j + 1;
            int next_step = next_i + next_j;

            const auto& row_map = posterior[i];
            auto row_it = row_map.find(j);
            value_type prob = (row_it != row_map.end()) ? row_it->second : 0.0;
            value_type diag_cost = (prob > 0.0) ? std::max<value_type>(0.0, 1.0 - prob) : 1.0;
            value_type new_cost = curr_cost + diag_cost;

            if (DijkstraState* next_state = acquire_state(next_step, next_i, next_j)) {
                if (new_cost < next_state->cost) {
                    next_state->cost = new_cost;
                    next_state->traceback = 'D';
                    heap.push({new_cost, next_i, next_j});
                }
            }
        }

        int next_step = i + j + 1;
        value_type new_cost = curr_cost + 0.5;

        // INS1 (gap in seq2)
        if (i + 1 <= rows) {
            int next_i = i + 1;
            int next_j = j;
            
            if (DijkstraState* next_state = acquire_state(next_step, next_i, next_j)) {
                if (new_cost < next_state->cost) {
                    next_state->cost = new_cost;
                    next_state->traceback = 'U';
                    heap.push({new_cost, next_i, next_j});
                }
            }
        }

        // INS2 (gap in seq1)
        if (j + 1 <= cols) {
            int next_i = i;
            int next_j = j + 1;
            if (DijkstraState* next_state = acquire_state(next_step, next_i, next_j)) {
                if (new_cost < next_state->cost) {
                    next_state->cost = new_cost;
                    next_state->traceback = 'L';
                    heap.push({new_cost, next_i, next_j});
                }
            }
        }
    }

    unsigned final_step = seq1Length + seq2Length;
    auto final_it = beam[final_step].find({seq1Length, seq2Length});
    if (final_it == beam[final_step].end()) {
        return make_pair(new string(), 0.0);
    }

    string* alignment = new string();
    unsigned current_i = seq1Length;
    unsigned current_j = seq2Length;
    unsigned current_step = final_step;

    size_t diag_count = 0;
    while (current_i > 0 || current_j > 0) {
        auto state_it = beam[current_step].find({current_i, current_j});
        if (state_it == beam[current_step].end()) break;

        char ch = state_it->second.traceback;
        switch (ch) {
        case 'D':
            alignment->push_back('B');
            diag_count++;
            current_i--;
            current_j--;
            current_step -= 2;
            break;
        case 'U':
            alignment->push_back('X');
            current_i--;
            current_step -= 1;
            break;
        case 'L':
            alignment->push_back('Y');
            current_j--;
            current_step -= 1;
            break;
        default:
            current_i = current_j = 0;
            current_step = 0;
            break;
        }
    }

    reverse(alignment->begin(), alignment->end());

    value_type total_cost = final_it->second.cost;
    const size_t gap_count = alignment->size() - diag_count;
    // Recover MEA-equivalent reward: |B| + 0.5 * gaps - accumulated cost
    value_type reconstructed_score =
        static_cast<value_type>(diag_count) + (0.5 * static_cast<value_type>(gap_count)) - total_cost;

    fprintf(stderr, "[LinearComputeAlignmentDijkstra] States explored: %zu\n", states_explored);

    return make_pair(alignment, reconstructed_score);
}

// pair<string *, value_type> ProbabilisticModel::LinearComputeAlignmentDijkstra5(int hmmBeam, int seq1Length, int seq2Length, const unordered_map<int, value_type>* posterior){
//     struct HeapEntry {
//         value_type priority;
//         int i;
//         int j;

//         bool operator>(const HeapEntry& other) const {
//             return priority > other.priority;
//         }
//     };

//     const int rows = seq1Length;
//     const int cols = seq2Length;
//     const size_t max_states_per_step =
//         hmmBeam > 0 ? static_cast<size_t>(hmmBeam) : std::numeric_limits<size_t>::max();

//     vector<unordered_map<pair<int, int>, DijkstraState, PairHash>> beam(rows + cols + 3);
//     priority_queue<HeapEntry, vector<HeapEntry>, greater<HeapEntry>> heap;

//     beam[0][{0, 0}] = DijkstraState(0.0, 'S');
//     heap.push({0.0, 0, 0});

//     auto acquire_state = [&](int step, int next_i, int next_j) -> DijkstraState* {
//         auto& next_map = beam[step];
//         auto key = make_pair(next_i, next_j);
//         auto it = next_map.find(key);
//         if (it != next_map.end()) {
//             return &it->second;
//         }
//         // if (next_map.size() >= max_states_per_step) {
//         //     return nullptr;
//         // }
//         return &next_map.emplace(key, DijkstraState()).first->second;
//     };

//     auto acquire_state_ALN = [&](int step, int next_i, int next_j) -> DijkstraState* {
//         auto& next_map = beam[step];
//         auto key = make_pair(next_i, next_j);
//         auto it = next_map.find(key);
//         if (it != next_map.end()) {
//             return &it->second;
//         }

//         return &next_map.emplace(key, DijkstraState()).first->second;
//     };

//     size_t states_explored = 0;

//     while (!heap.empty()) {
//         HeapEntry top = heap.top();
//         heap.pop();

//         int i = top.i;
//         int j = top.j;
//         int step = i + j;

//         auto it = beam[step].find({i, j});
//         if (it == beam[step].end() || top.priority > it->second.cost) {
//             continue;
//         }

//         states_explored++;

//         if (i == rows && j == cols) {
//             break;
//         }

//         const value_type curr_cost = it->second.cost;
        
//         // set common values for INS beams
//         int next_step = i + j + 1;
//         value_type new_cost = curr_cost + 0.5;
//         const auto& next_row_map = posterior[i]; // next row -- 0-indexed posterior; (i+1) - 1 == i
//         const auto& curr_row_map = posterior[i-1]; // curr row

//         // INS1 (gap in seq2)
//         if (i + 1 <= rows) {
//             int next_i = i + 1;
//             int next_j = j;
            
//             auto row_it = next_row_map.find(j-1);
//             if (row_it != next_row_map.end()){
//                 if (DijkstraState* next_state = acquire_state(next_step, next_i, next_j)) {
//                     if (new_cost < next_state->cost) {
//                         next_state->cost = new_cost;
//                         next_state->traceback = 'U';
//                         heap.push({new_cost, next_i, next_j});
//                     }
//                 }
//             }
//         }

//         // INS2 (gap in seq1)
//         if (j + 1 <= cols) {
//             int next_i = i;
//             int next_j = j + 1;

//             auto row_it = curr_row_map.find(j);
//             if (row_it != curr_row_map.end()){
//                 if (DijkstraState* next_state = acquire_state(next_step, next_i, next_j)) {
//                     if (new_cost < next_state->cost) {
//                         next_state->cost = new_cost;
//                         next_state->traceback = 'L';
//                         heap.push({new_cost, next_i, next_j});
//                     }
//                 }
//             }
//         }

//         // ALN (match)
//         if (i + 1 <= rows && j + 1 <= cols) {
//             int next_i = i + 1;
//             int next_j = j + 1;
//             int next_step = next_i + next_j;

//             auto row_it = row_map.find(j);
//             value_type prob = (row_it != row_map.end()) ? row_it->second : 0.0;
//             value_type diag_cost = std::max<value_type>(0.0, 1.0 - prob);
//             value_type new_cost = curr_cost + diag_cost;

//             if (DijkstraState* next_state = acquire_state_ALN(next_step, next_i, next_j)) {
//                 if (new_cost < next_state->cost) {
//                     next_state->cost = new_cost;
//                     next_state->traceback = 'D';
//                     heap.push({new_cost, next_i, next_j});
//                 }
//             }
//         }
//     }

//     unsigned final_step = seq1Length + seq2Length;
//     auto final_it = beam[final_step].find({seq1Length, seq2Length});
//     if (final_it == beam[final_step].end()) {
//         return make_pair(new string(), 0.0);
//     }

//     string* alignment = new string();
//     unsigned current_i = seq1Length;
//     unsigned current_j = seq2Length;
//     unsigned current_step = final_step;

//     size_t diag_count = 0;
//     while (current_i > 0 || current_j > 0) {
//         auto state_it = beam[current_step].find({current_i, current_j});
//         if (state_it == beam[current_step].end()) break;

//         char ch = state_it->second.traceback;
//         switch (ch) {
//         case 'D':
//             alignment->push_back('B');
//             diag_count++;
//             current_i--;
//             current_j--;
//             current_step -= 2;
//             break;
//         case 'U':
//             alignment->push_back('X');
//             current_i--;
//             current_step -= 1;
//             break;
//         case 'L':
//             alignment->push_back('Y');
//             current_j--;
//             current_step -= 1;
//             break;
//         default:
//             current_i = current_j = 0;
//             current_step = 0;
//             break;
//         }
//     }

//     reverse(alignment->begin(), alignment->end());

//     value_type total_cost = final_it->second.cost;
//     const size_t gap_count = alignment->size() - diag_count;
//     // Recover MEA-equivalent reward: |B| + 0.5 * gaps - accumulated cost
//     value_type reconstructed_score =
//         static_cast<value_type>(diag_count) + (0.5 * static_cast<value_type>(gap_count)) - total_cost;

//     fprintf(stderr, "[LinearComputeAlignmentDijkstra] States explored: %zu\n", states_explored);

//     return make_pair(alignment, reconstructed_score);
// }

// pair<string *, value_type> ProbabilisticModel::LinearComputeAlignmentDijkstraEfficient(int hmmBeam, int seq1Length, int seq2Length, const unordered_map<int, value_type>* posterior){
//     /* 
//         Efficient Dijkstra implementation with decrease-key:
//         - Uses mutable boost heap with decrease_key operation
//         - No stale entries - when we find a better path, we decrease the key
//         - Store traceback info in heap entries, save to beam when finalized (popped)
//     */

//     struct HeapEntry {
//         value_type priority;  // Dijkstra cost (sum of 1-p)
//         value_type score;     // Actual reward (sum of probabilities)
//         int i, j;
//         char traceback;
        
//         bool operator>(const HeapEntry& other) const {
//             return priority > other.priority;  // min-heap
//         }
//     };

//     using Handle = typename boost::heap::d_ary_heap<
//         HeapEntry,
//         boost::heap::arity<2>,
//         boost::heap::compare<greater<HeapEntry>>,
//         boost::heap::mutable_<true>
//     >::handle_type;

//     boost::heap::d_ary_heap<HeapEntry, boost::heap::arity<2>, boost::heap::compare<greater<HeapEntry>>, boost::heap::mutable_<true>> heap;
    
//     // Map from (i,j) to heap handle for decrease-key operation
//     unordered_map<pair<int, int>, Handle, PairHash> handles;
    
//     // Beam stores finalized states (popped from heap) for traceback
//     vector<unordered_map<pair<int, int>, MEAState, PairHash>> beam(seq1Length + seq2Length + 3);
    
//     size_t states_explored = 0;
//     size_t max_heap_size = 0;
//     const size_t kHeapLogStart = 5000;
//     size_t next_heap_log = kHeapLogStart;
//     const size_t kBeamLogStart = 5000;
//     size_t next_beam_log = kBeamLogStart;
//     size_t total_beam_entries = 1; // initial (0,0)
    
//     // Initialize: start at (0,0)
//     Handle h0 = heap.push({0.0, 0.0, 0, 0, '\0'});
//     handles[{0, 0}] = h0;

//     auto log_heap_growth = [&](const char* move, int from_i, int from_j, int to_i, int to_j, int target_step) {
//         size_t current_size = heap.size();
//         if (current_size >= next_heap_log) {
//             fprintf(stderr,
//                     "[LinearComputeAlignmentDijkstra-Efficient] Heap size %zu reached at step %d via %s (from %d,%d -> %d,%d)\n",
//                     current_size, target_step, move, from_i, from_j, to_i, to_j);
//             next_heap_log *= 2;
//         }
//     };

//     // Dijkstra search with decrease-key
//     while (!heap.empty()) {
//         // Track maximum heap size
//         if (handles.size() > max_heap_size) {
//             max_heap_size = handles.size();
//         }
        
//         HeapEntry top = heap.top();
//         heap.pop();
        
//         int i = top.i;
//         int j = top.j;
//         int step = i + j;
        
//         // Remove from handles (this node is now finalized)
//         handles.erase({i, j});
        
//         // Store finalized state in beam for traceback
//         beam[step][{i, j}] = MEAState(top.score, top.priority, top.traceback);
        
//         states_explored++;

//         // Goal reached
//         if (i == seq1Length && j == seq2Length) break;

//         value_type curr_cost = top.priority;
//         value_type curr_score = top.score;

//         // Helper lambda to update or insert neighbor
//         auto update_neighbor = [&](int next_i, int next_j, value_type new_cost, value_type new_score, char trace, const char* move_label) {
//             auto key = make_pair(next_i, next_j);
//             auto it = handles.find(key);
            
//             if (it != handles.end()) {
//                 // Node already in heap - check if we found better path
//                 HeapEntry& entry = (*it->second);
//                 if (new_cost < entry.priority) {
//                     // Update with better path
//                     entry.priority = new_cost;
//                     entry.score = new_score;
//                     entry.traceback = trace;
//                     heap.increase(it->second);  // Actually decrease in min-heap (boost naming)
//                 }
//             } else {
//                 // Check if already finalized
//                 int next_step = next_i + next_j;
//                 if (beam[next_step].find(key) == beam[next_step].end()) {
//                     // Not in handles and not finalized - insert into heap
//                     Handle h = heap.push({new_cost, new_score, next_i, next_j, trace});
//                     handles[key] = h;
//                     log_heap_growth(move_label, i, j, next_i, next_j, next_step);
//                 }
//             }
//         };

//         // Explore neighbors
//         // INS1: move down in seq1 (no reward, cost increases by 0.5)
//         if (i + 1 <= seq1Length) {
//             update_neighbor(i + 1, j, curr_cost + 0.5, curr_score, 'U', "INS1");
//         }

//         // INS2: move right in seq2 (no reward, cost increases by 0.5)
//         if (j + 1 <= seq2Length) {
//             update_neighbor(i, j + 1, curr_cost + 0.5, curr_score, 'L', "INS2");
//         }

//         // ALN: diagonal move (posterior stored 0-indexed at [i][j])
//         if (i + 1 <= seq1Length && j + 1 <= seq2Length) {
//             int next_i = i + 1;
//             int next_j = j + 1;
//             auto row_it = posterior[i].find(j);
//             if (row_it != posterior[i].end() && row_it->second > 0) {
//                 value_type prob = row_it->second;
//                 value_type new_cost = curr_cost + max(0.0, 1.0 - prob);
//                 value_type new_score = curr_score + prob;
//                 update_neighbor(next_i, next_j, new_cost, new_score, 'D', "ALN");
//             }
//         }
//     }

//     // Traceback
//     value_type total = 0.0;
//     unsigned final_step = seq1Length + seq2Length;
    
//     auto final_it = beam[final_step].find({seq1Length, seq2Length});
//     if (final_it != beam[final_step].end()) {
//         total = final_it->second.score;
//     }

//     string* alignment = new string();
//     unsigned current_i = seq1Length;
//     unsigned current_j = seq2Length;
//     unsigned current_step = final_step;

//     while (current_i > 0 || current_j > 0) {
//         auto state_it = beam[current_step].find({current_i, current_j});
//         if (state_it == beam[current_step].end()) break;

//         char ch = state_it->second.traceback;
//         switch (ch) {
//         case 'D': // Diagonal (match)
//             alignment->push_back('B');
//             current_i--;
//             current_j--;
//             current_step -= 2;
//             break;
//         case 'U': // Up (ins1)
//             alignment->push_back('X');
//             current_i--;
//             current_step -= 1;
//             break;
//         case 'L': // Left (ins2)
//             alignment->push_back('Y');
//             current_j--;
//             current_step -= 1;
//             break;
//         default:
//             current_step = 0;
//             current_i = current_j = 0;
//             break;
//         }
//     }

//     reverse(alignment->begin(), alignment->end());

//     fprintf(stderr, "[LinearComputeAlignmentDijkstra-Efficient] States explored: %zu, Max heap size: %zu\n", states_explored, max_heap_size);

//     return make_pair(alignment, total);
// }

void ProbabilisticModel::LinearConsistencyTransform(int seq1Length, unordered_map<int, value_type>* &xz_posterior, unordered_map<int, value_type>* &zy_posterior, unordered_map<int, value_type>* &new_xy_posterior){
    for(int i = 0; i < seq1Length; i++){
        for(auto &xz_cand : xz_posterior[i]){
            int k = xz_cand.first;

            for(auto &zy_cand : zy_posterior[k]){
                int j = zy_cand.first;
                new_xy_posterior[i][j] += xz_cand.second * zy_cand.second;
            }
        }
    }
}

vector<vector<unordered_map<int, value_type>*>> ProbabilisticModel::LinearMultiConsistencyTransform(MultiSeq &sequences, vector<vector<unordered_map<int, value_type>*>> &posterior){
    const int numSeqs = sequences.size();
    
    // allocate space for the new posterior matrix
    vector<vector<unordered_map<int, value_type>*>> new_posterior(numSeqs);
    for (int i = 0; i < numSeqs; i++){
        new_posterior[i].resize(numSeqs);
        for (int j = 0; j < numSeqs; j++){
            if (i == j) continue;
            new_posterior[i][j] = new unordered_map<int, value_type>[sequences.at(i).length()];
        }
    }

    // For every pair of sequences
    for (int x = 0; x < numSeqs; x++){
        for (int y = x+1; y < numSeqs; y++){
            Sequence seq1 = sequences.at(x);
            Sequence seq2 = sequences.at(y);

            const int seq1Length = seq1.length();
            const int seq2Length = seq2.length();

            // allocate space for temporary results
            unordered_map<int, value_type>* transformation = new unordered_map<int, value_type>[seq1Length];

            // Get the original alignment result
            unordered_map<int, value_type>* &original = posterior[x][y];

            // Contribution from the summation where z = x and z = y
            for (int i = 0; i < seq1Length; i++){
                for(auto &item : original[i]){
                    int j = item.first;
                    transformation[i][j] = 2 * item.second;
                }
            }

            // Contribution from all other sequences
            for (int z = 0; z < numSeqs; z++) {
                if (z == x || z == y) continue;
                LinearConsistencyTransform(seq1Length, posterior[x][z], posterior[z][y], transformation);
            }

            // Renormalization
            for (int i = 0; i < seq1Length; i++){
                for(auto &item : transformation[i]){
                    int j = item.first;
                    transformation[i][j] /= numSeqs;
                }
            }

            // // Renormalization and mask out smaller values
            // for (int i = 0; i < seq1Length; ++i) {
            //     auto& cons_trans_i = transformation[i];

            //     for (auto it = cons_trans_i.begin(); it != cons_trans_i.end(); ) {
            //         int j = it->first;
            //         value_type prob = it->second;
            //         prob /= numSeqs;

            //         if (prob < 0.01) {
            //             it = cons_trans_i.erase(it); // returns iterator to next
            //         } else {
            //             posterior[x][y][i][j] = prob;
            //             posterior[y][x][j][i] = prob;
            //             ++it;
            //         }
            //     }
            // }

            // Mask out positions not originally in the posterior matrix
            for (int i = 0; i < seq1Length; i++){
                for(auto &item : original[i]){
                    int j = item.first;
                    if (transformation[i].find(j) == transformation[i].end()) continue;
                    if (transformation[i][j] >= 0.01){
                        new_posterior[x][y][i][j] = transformation[i][j];
                        new_posterior[y][x][j][i] = transformation[i][j];
                        // cout << i << " " << j << " " << k << " " << l << " " << temp_pair_CT[k][l].value << endl;
                    }
                }
            }
            delete[] transformation;
        }
    }

    return new_posterior;
    // return posterior;
}

/// Equivalent to BuildPosterior from the Probcons 
/////////////////////////////////////////////////////////////////
// ProbabilisticModel::BuildPosterior()
//
// Builds a posterior probability matrix needed to align a pair
// of alignments.  Mathematically, the returned matrix M is
// defined as follows:
//    M[i,j] =     sum          sum      f(s,t,i,j)
//             s in align1  t in align2
// where
//                  [  P(s[i'] <--> t[j'])
//                  [       if s[i'] is a letter in the ith column of align1 and
//                  [          t[j'] it a letter in the jth column of align2
//    f(s,t,i,j) =  [
//                  [  0    otherwise
//
/////////////////////////////////////////////////////////////////
unordered_map<int, value_type>* ProbabilisticModel::LinearMultiAlnResults(MultiSeq* align1, MultiSeq* align2, const vector<vector<unordered_map<int, value_type>*>> &posterior, value_type cutoff) const {
    const int seq1Length = align1->at(0).length();
    const int seq2Length = align2->at(0).length();
    
    unordered_map<int, value_type>* transformation = new unordered_map<int, value_type>[seq1Length];
    vector<unordered_map<int, size_t>> contribution_counts(seq1Length);
    for (int i = 0; i < align1->size(); i++){
        int first = align1->at(i).id;
        vector<int>* mapping1 = align1->at(i).get_mapping();

        // Loops through align2
        for (int j = 0; j < align2->size(); j++){
            int second = align2->at(j).id;
            vector<int>* mapping2 = align2->at(j).get_mapping();

            if(first < second){
                unordered_map<int, value_type>* original = posterior[first][second];

                int seq1len = mapping1->size();
                for (int ii = 0; ii < seq1len; ii++){
                    int ibase = (*mapping1)[ii];

                    for (auto &item : original[ii]) {
                        int jbase = (*mapping2)[item.first];
                        if (item.second < cutoff) continue;

                        transformation[ibase][jbase] += item.second;
                        contribution_counts[ibase][jbase] += 1;
                    }
                }
            } else {
                unordered_map<int, value_type>* original = posterior[second][first];

                int seq2len = mapping2->size();
                for (int jj = 0; jj < seq2len; jj++){
                    int jbase = (*mapping2)[jj];

                    for (auto &item : original[jj]) {
                        int ibase = (*mapping1)[item.first];
                        if (item.second < cutoff) continue;

                        transformation[ibase][jbase] += item.second;
                        contribution_counts[ibase][jbase] += 1;
                    }
                }
            }

            delete mapping2;
        }

        delete mapping1;
    }
    for (int i = 0; i < seq1Length; ++i) {
        auto& row = transformation[i];
        auto& count_row = contribution_counts[i];
        for (auto& kv : row) {
            auto count_it = count_row.find(kv.first);
            if (count_it != count_row.end() && count_it->second > 0) {
                kv.second /= static_cast<value_type>(count_it->second);
            }
        }
    }
    return transformation;
}

pair<MultiSeq*, value_type> ProbabilisticModel::LinearAlignAlignments (MultiSeq* align1, MultiSeq* align2, const vector<vector<unordered_map<int, value_type>*>> &posterior, int hmmBeam){

    // Choose the alignment routine depending on the "cosmetic" gap penalties used
    unordered_map<int, value_type>* posterior_matrix = LinearMultiAlnResults(align1, align2, posterior, 0.01f);

    pair<string*, value_type> alignment = LinearComputeAlignmentDijkstra(hmmBeam, align1->at(0).length(), align2->at(0).length(), posterior_matrix);
    delete[] posterior_matrix;

    fprintf(stderr, "Alignment: %s\n", alignment.first->c_str());
    fprintf(stderr, "MEA score: %f\n", alignment.second);
    
    // Build final alignment
    MultiSeq *result = new MultiSeq();
    for (int i = 0; i < align1->size(); i++)
        result->add_sequence (  *(align1->at(i).add_gaps(alignment.first, 'X')) );
    for (int i = 0; i < align2->size(); i++)
        result->add_sequence (  *(align2->at(i).add_gaps(alignment.first, 'Y')) );

    value_type mea_score = alignment.second;
    delete alignment.first;

    return make_pair(result, mea_score);
}

pair<MultiSeq*, value_type> ProbabilisticModel::LinearProcessTree (const TreeNode *tree, MultiSeq* sequences, const vector<vector<unordered_map<int, value_type>*>> &posterior, int hmmBeam){
    MultiSeq *result;
    value_type mea_score = 0.0;

    // Check if this is an internal node of the alignment tree
    if (tree->GetSequenceLabel() == -1){
        auto left_result = LinearProcessTree (tree->GetLeftChild(), sequences, posterior, hmmBeam);
        auto right_result = LinearProcessTree (tree->GetRightChild(), sequences, posterior, hmmBeam);

        MultiSeq *alignLeft = left_result.first;
        MultiSeq *alignRight = right_result.first;

        assert (alignLeft);
        assert (alignRight);

        auto alignment_result = LinearAlignAlignments (alignLeft, alignRight, posterior, hmmBeam);
        result = alignment_result.first;
        mea_score = alignment_result.second;
        assert (result);

        delete alignLeft;
        delete alignRight;
    }
    // Otherwise, this is a leaf of the alignment tree
    else {
        result = new MultiSeq(); assert (result);
        result->add_sequence( *(sequences->at(tree->GetSequenceLabel()).clone()) );
        mea_score = 0.0; // No alignment for single sequences
    }

    return make_pair(result, mea_score);
}

pair<MultiSeq*, value_type> ProbabilisticModel::LinearDoIterativeRefinement (const vector<vector<unordered_map<int, value_type>*>> &posterior, MultiSeq* alignment, int i, int hmmBeam){
    set<int> groupOne, groupTwo;

    randomnumber rn;
    rn.seed(1234+i);

    // Create two separate groups
    for (int j = 0; j < alignment->size(); j++){
        int x = rn.roll_int(1,10);

        if (x % 2)
            groupOne.insert (j);
        else
            groupTwo.insert (j);
    }

    if (groupOne.empty() || groupTwo.empty()) {
        // Return the original alignment with a default MEA score of 0
        return make_pair(new MultiSeq(*alignment), 0.0);
    }

    // Project into the two groups
    MultiSeq *groupOneSeqs = alignment->Project (groupOne); assert (groupOneSeqs);
    MultiSeq *groupTwoSeqs = alignment->Project (groupTwo); assert (groupTwoSeqs);

    // Realign
    auto alignment_result = LinearAlignAlignments (groupOneSeqs, groupTwoSeqs, posterior, hmmBeam);

    delete groupOneSeqs;
    delete groupTwoSeqs;

    return alignment_result;
}

MultiSeq* ProbabilisticModel::LinearComputeFinalAlignment (const TreeNode *tree, MultiSeq* sequences, const vector<vector<unordered_map<int, value_type>*>> &posterior, int hmmBeam, TurboFoldLog* log){
    auto process_tree_start = std::chrono::high_resolution_clock::now();
    auto tree_result = LinearProcessTree (tree, sequences, posterior, hmmBeam);
    auto process_tree_end = std::chrono::high_resolution_clock::now();
    
    if (log) {
        log->msa_process_tree_time = std::chrono::duration_cast<std::chrono::milliseconds>(process_tree_end - process_tree_start).count();
        fprintf(stderr, "[MSA] Tree processing completed (%.2f ms)\n", log->msa_process_tree_time);
    }

    // Initialize best alignment with the tree result and its actual MEA score
    MultiSeq* best_alignment = tree_result.first;
    value_type best_mea = tree_result.second;

    fprintf(stderr, "Initial tree alignment MEA: %f\n", best_mea);
    fprintf(stderr, "\nStarting iterative refinement\n");

    auto iterative_refine_start = std::chrono::high_resolution_clock::now();
    
    // Iterative refinement - track the best alignment by MEA score
    for (int i = 0; i < num_iterative_refinement_reps_; i++) {
        auto refinement_result = LinearDoIterativeRefinement (posterior, best_alignment, i, hmmBeam);
        MultiSeq* new_alignment = refinement_result.first;
        value_type new_mea = refinement_result.second;

        fprintf(stderr, "Iteration %d - New MEA: %f, Current best: %f\n", i, new_mea, best_mea);

        // Keep the alignment with higher MEA score
        if (new_mea > best_mea) {
            fprintf(stderr, "  -> Updating best alignment (new MEA %f > old %f)\n", new_mea, best_mea);
            delete best_alignment;
            best_alignment = new_alignment;
            best_mea = new_mea;
        } else {
            fprintf(stderr, "  -> Keeping current alignment\n");
            delete new_alignment;
        }
    }
    
    auto iterative_refine_end = std::chrono::high_resolution_clock::now();
    
    if (log) {
        log->msa_iterative_refine_time = std::chrono::duration_cast<std::chrono::milliseconds>(iterative_refine_end - iterative_refine_start).count();
        fprintf(stderr, "[MSA] Iterative refinement completed (%.2f ms)\n", log->msa_iterative_refine_time);
    }

    fprintf(stderr, "Final best MEA: %f\n", best_mea);
    return best_alignment;
}

