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
    State() : alpha(xlog(0.0)), beta(xlog(0.0)) {};
};

struct TraceInfo {
    int i;
    int j;
    int t;  // split point
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
    State *right;  // right=null <=> Edge

    HEdge() : weight(xlog(0.0)), left(nullptr), right(nullptr) {};
    HEdge(double weight, State *left, State *right) : weight(weight), left(left), right(right) {};

    void set(double weight, State *left, State *right) {
        this->weight = weight;
        this->left = left;
        this->right = right;
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

#endif  // LFP_UTILITIES_HPP
