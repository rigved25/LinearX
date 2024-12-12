#ifndef LA_UTILS_HPP
#define LA_UTILS_HPP

#include "./../utility/log_math.hpp"


inline double LOG(double x) {
    return xlog(x);
}

inline double EXP(double x) {
    return xexp(x);
}

inline double LOG_SUM(double a, double b) {
    return xlog_sum(a, b);
    // Fast_LogPlusEquals(a, b);
    // return a;
}

inline double LOG_MUL(double a, double b) {
    return xlog_mul(a, b);
}

inline double LOG_DIV(double a, double b) {
    return xlog_div(a, b);
}

enum HStateType {
    INS1, // 0
    INS2, // 1
    ALN,  // 2
};

struct HState {
    double alpha;
    double beta;

    HState() : alpha(xlog(0.0)), beta(xlog(0.0)) {};
};

struct AlnEdge {
    double weight;
    HState *prev;

    AlnEdge() : weight(xlog(0.0)), prev(nullptr) {};

    void set(HState *prev, const double weight = 0) {
        this->prev = prev;
        this->weight = weight;
    }
};

#endif // LA_UTILS_HPP
