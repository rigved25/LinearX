// linearx/partition/utility.hpp
#pragma once
#include <linearx/math/log_math.hpp>

enum StateType {
    H,
    Multi,
    P,
    M2,
    M,
    C,
};

struct PartitionLog {
    value_type free_energy_of_ensemble;  // should be equal to -energy
    value_type total_inside_energy;      // corresponds to bestC[seq_length - 1].alpha
    value_type total_outside_energy;     // corresponds to bestC[-1].beta
    value_type inside_exec_time;
    value_type outside_exec_time;
    float effective_beam_size;
    std::string run_mode;
    unsigned long nodes_visited;
    unsigned long nodes_pruned;
    unsigned long edges_saved;
    unsigned long edges_pruned;
};

struct State {
    value_type alpha;
    value_type beta;
    State() : alpha(linearx::math::LOG_ZERO), beta(linearx::math::LOG_ZERO) {};
};

struct TraceInfo {
    int i;
    int j;
    int t;  // split point
    StateType type_left;
    StateType type_right;

    constexpr TraceInfo() : i(-1), j(-1), t(-1), type_left(H), type_right(H) {}  // default constructor
    constexpr TraceInfo(int i, int j, int t, StateType type_left, StateType type_right)
        : i(i), j(j), t(t), type_left(type_left), type_right(type_right) {}  // parameterized constructor

    inline void reset() {
        i = -1;
        j = -1;
        t = -1;
        type_left = H;
        type_right = H;
    }

    inline void set(int i, int j, int t, StateType type_left, StateType type_right) {
        this->i = i;
        this->j = j;
        this->t = t;
        this->type_left = type_left;
        this->type_right = type_right;
    }

    // TraceInfo(const TraceInfo &) = default;             // copy constructor
    // TraceInfo &operator=(const TraceInfo &) = default;  // copy assignment
    // TraceInfo(TraceInfo &&) = default;                  // move constructor
    // TraceInfo &operator=(TraceInfo &&) = default;       // move assignment
};

struct HEdge {
    value_type weight;
    State *left;
    State *right;  // right == nullptr <=> unary edge

    HEdge() : weight(linearx::math::LOG_ZERO), left(nullptr), right(nullptr) {}  // default constructor

    HEdge(value_type weight, State *left, State *right)
        : weight(weight), left(left), right(right) {}  // parameterized constructor

    // HEdge(const HEdge &) = default;             // copy constructor
    // HEdge &operator=(const HEdge &) = default;  // copy assignment
    // HEdge(HEdge &&) = default;                  // move constructor
    // HEdge &operator=(HEdge &&) = default;       // move assignment

    inline void reset() {
        weight = linearx::math::LOG_ZERO;
        left = nullptr;
        right = nullptr;
    }

    inline void set(value_type weight, State *left, State *right) {
        this->weight = weight;
        this->left = left;
        this->right = right;
    }

    // set next_state alpha using current states and the edge
    inline void update_state_alpha(State &next_state) {
        value_type prev_score = (left ? left->alpha : 0) + (right ? right->alpha : 0);
        // std::cout << prev_score << " " << weight << std::endl;
        next_state.alpha = LOG_SUM(next_state.alpha, prev_score + (weight * linearx::constants::energy::INV_KT));
    }

    // set current states beta using next_state and the edge
    inline void update_state_beta(State &next_state) {
        if (!right) {
            left->beta = LOG_SUM(left->beta, weight + next_state.beta);
        } else {
            left->beta = LOG_SUM(left->beta, right->alpha + weight + next_state.beta);
            right->beta = LOG_SUM(right->beta, left->alpha + weight + next_state.beta);
        }
    }
};

template <typename T>
class VectorWithNegOneIndex {
   private:
    T specialCase;        // For index -1
    std::vector<T> data;  // For non-negative indices

   public:
    // Constructor to initialize the vector with a given size and optional default value
    VectorWithNegOneIndex(size_t size, const T &defaultValue = T())
        : specialCase(defaultValue), data(size, defaultValue) {}

    // Clear the vector and reset the special case
    void clear(const T &defaultValue = T()) {
        specialCase = defaultValue;  // Reset special case to the default value
        data.clear();                // Clear the vector
    }

    void reset(const T &defaultValue = T()) {
        specialCase = defaultValue;                         // Reset special case
        std::fill(data.begin(), data.end(), defaultValue);  // Reset all vector elements
    }

    // Access operator (supports -1 and non-negative indices)
    T &operator[](int index) {
        if (index == -1) return specialCase;
        if (index < 0 || static_cast<size_t>(index) >= data.size()) throw std::out_of_range("Index out of range");
        return data[index];
    }

    // Const access operator
    const T &operator[](int index) const {
        if (index == -1) return specialCase;
        if (index < 0 || static_cast<size_t>(index) >= data.size()) throw std::out_of_range("Index out of range");
        return data[index];
    }

    // Get the size of the vector (excluding the special case)
    size_t size() const { return data.size(); }

    // Resize the vector (does not affect the special case)
    void resize(size_t newSize, const T &defaultValue = T()) { data.resize(newSize, defaultValue); }
};
