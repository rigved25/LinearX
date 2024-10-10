#ifndef LFP_UTILITIES_HPP
#define LFP_UTILITIES_HPP

#include "./../utility/constants.hpp"

enum InsideMode {
    MFE,
    PARTITION,
};

enum StateType {
    H,
    Multi,
    P,
    M2,
    M,
    C,
};

struct State {
    double alpha;
    double beta;
    State() : alpha(VALUE_MIN), beta(VALUE_MIN) {};
};

struct TraceInfo {
    int i;
    int j;
    int t; // split point
    StateType type_left;
    StateType type_right;

    TraceInfo() : i(-1), j(-1), t(-1), type_left(H), type_right(H) {};

    void set(const int i, const int j, const int t, const StateType type_left, const StateType type_right) {
        this->i = i;
        this->j = j;
        this->t = t;
        this->type_left = type_left;
        this->type_right = type_right;
    }
};

struct HEdge {
    double weight;
    State *left;
    State *right; // right=null <=> Edge

    HEdge() : weight(VALUE_MIN), left(nullptr), right(nullptr) {};
    HEdge(double weight, State *left, State *right) : weight(weight), left(left), right(right) {};

    void set(double weight, State *left, State *right) {
        this->weight = weight;
        this->left = left;
        this->right = right;
    }
};

#endif // LFP_UTILITIES_HPP
