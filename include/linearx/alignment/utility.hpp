// linearx/linear_align/utility.hpp
#pragma once
#include <linearx/math/log_math.hpp>

enum HStateType {
    INS1,  // 0
    INS2,  // 1
    ALN,   // 2
};

struct AlignmentLog {
    value_type seq_identity = 0.0;
    bool lazy_outside = false;
    std::pair<value_type, value_type> best_exec_time = {0.0, 0.0};        // (inside, outside)
    value_type cp_exec_time = 0.0;
    std::pair<value_type, value_type> exec_time = {0.0, 0.0};             // (inside, outside)
    std::pair<value_type, value_type> total_score = {0.0, 0.0};           // (inside, outside)
    std::pair<float, float> effective_beam_size = {0.0f, 0.0f};           // (inside, outside)
    std::pair<unsigned long, unsigned long> states_visited = {0, 0};      // (inside, outside)
    std::pair<unsigned long, unsigned long> states_pruned = {0, 0};       // (inside, outside)
    unsigned long edges_saved = 0;
    unsigned long edges_pruned = 0;

    // Overhead timing fields (ms) for profiling non-core alignment steps
    value_type best_prob_setup_time = 0.0;       // use_prob_set1() + set_prob_accm() before best mode
    value_type save_best_pf_time = 0.0;          // save_partition_function(false)
    value_type traceback_time = 0.0;             // get_alignment().average_seq_identity()
    value_type best_reset_beams_time = 0.0;      // reset_beams() after best mode
    value_type partition_prob_setup_time = 0.0;   // use_prob_set2() + set_prob_accm() before partition mode
    value_type save_partition_pf_time = 0.0;     // save_partition_function(true)
    value_type partition_reset_beams_time = 0.0; // reset_beams() after partition mode
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
