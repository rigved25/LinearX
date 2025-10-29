// linearx/linear_align/utility.hpp
#pragma once
#include <linearx/math/log_math.hpp>

enum HStateType {
    INS1,  // 0
    INS2,  // 1
    ALN,   // 2
};

struct AlignmentLog {
    value_type seq_identity;
    bool lazy_outside;
    value_type best_exec_time;
    value_type cp_exec_time;
    std::pair<value_type, value_type> exec_time;             // (inside, outside)
    std::pair<value_type, value_type> total_score;           // (inside, outside)
    std::pair<float, float> effective_beam_size;             // (inside, outside)
    std::pair<unsigned long, unsigned long> states_visited;  // (inside, outside)
    std::pair<unsigned long, unsigned long> states_pruned;   // (inside, outside)
    unsigned long edges_saved;
    unsigned long edges_pruned;
};

struct HState {
    value_type alpha;
    value_type beta;

    HState() : alpha(linearx::math::LOG_ZERO), beta(linearx::math::LOG_ZERO) {};
};

struct AlnEdge {
    value_type weight;
    HState* prev;

    AlnEdge() : weight(linearx::math::LOG_ZERO), prev(nullptr) {};                 // Default constructor
    AlnEdge(HState* prev, value_type weight = 0) : prev(prev), weight(weight) {};  // Parameterized constructor

    void reset() {
        weight = linearx::math::LOG_ZERO;
        prev = nullptr;
    }

    void set(HState* prev, const value_type weight = 0) {
        this->prev = prev;
        this->weight = weight;
    }

    void update_state_alpha(HState& next_state) {
        next_state.alpha = LOG_SUM(next_state.alpha, LOG_MUL(prev->alpha, weight));
    }

    void update_state_beta(HState& next_state) {
        prev->beta = LOG_SUM(prev->beta, LOG_MUL(next_state.beta, weight));
    }
};
