// linearx/linear_align/utility.hpp
#pragma once
#include <linearx/math/log_math.hpp>

inline double LOG(double x) { return linearx::math::xlog(x); }

inline double EXP(double x) { return linearx::math::xexp(x); }

inline double LOG_SUM(double a, double b) {
    return linearx::math::xlog_sum(a, b);
    // linearx::math::Fast_LogPlusEquals(a, b);
    // return a;
}

inline double LOG_MUL(double a, double b) { return linearx::math::xlog_mul(a, b); }

inline double LOG_DIV(double a, double b) { return linearx::math::xlog_div(a, b); }

enum HStateType {
    INS1,  // 0
    INS2,  // 1
    ALN,   // 2
};

struct AlignmentInsideLog {
    double score;  // corresponds to bestALN[seq_len_sum+2][{seq1.length()+1, seq2.length()+1}].alpha
    double execution_time;
    unsigned beam_size;
    unsigned long nodes_pruned;
};

struct AlignmentOutsideLog {
    double score;  // corresponds to bestALN[0][{0, 0}].beta
    double execution_time;
    double deviation_threshold;
    double effective_beam_size;
    unsigned long nodes_visited;
    unsigned long nodes_pruned;
    unsigned long edges_saved;
    unsigned long edges_pruned;
};

struct HState {
    double alpha;
    double beta;

    HState() : alpha(linearx::math::xlog(0.0)), beta(linearx::math::xlog(0.0)) {};
};

struct AlnEdge {
    double weight;
    HState *prev;

    AlnEdge() : weight(linearx::math::xlog(0.0)), prev(nullptr) {};            // Default constructor
    AlnEdge(HState *prev, double weight = 0) : prev(prev), weight(weight) {};  // Parameterized constructor

    inline void reset() {
        weight = linearx::math::xlog(0.0);
        prev = nullptr;
    }

    inline void set(HState *prev, const double weight = 0) {
        this->prev = prev;
        this->weight = weight;
    }
};
