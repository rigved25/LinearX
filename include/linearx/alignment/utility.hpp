// linearx/linear_align/utility.hpp
#pragma once
#include <linearx/math/log_math.hpp>

enum HStateType {
    INS1,  // 0
    INS2,  // 1
    ALN,   // 2
};

struct AlignmentInsideLog {
    value_type score;  // corresponds to bestALN[seq_len_sum+2][{seq1.length()+1, seq2.length()+1}].alpha
    value_type execution_time;
    unsigned beam_size;
    unsigned long nodes_pruned;
};

struct AlignmentOutsideLog {
    value_type score;  // corresponds to bestALN[0][{0, 0}].beta
    value_type execution_time;
    value_type deviation_threshold;
    float effective_beam_size;
    unsigned long nodes_visited;
    unsigned long nodes_pruned;
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
    HState *prev;

    AlnEdge() : weight(linearx::math::LOG_ZERO), prev(nullptr) {};                // Default constructor
    AlnEdge(HState *prev, value_type weight = 0) : prev(prev), weight(weight) {};  // Parameterized constructor

    inline void reset() {
        weight = linearx::math::LOG_ZERO;
        prev = nullptr;
    }

    inline void set(HState *prev, const value_type weight = 0) {
        this->prev = prev;
        this->weight = weight;
    }
};
