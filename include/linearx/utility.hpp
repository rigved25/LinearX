// linearx/utility.hpp
#pragma once

#include <unistd.h>

#include <linearx/config.hpp>
#include <linearx/constants.hpp>
#include <sstream>
#include <vector>
#include <unordered_map>

enum Mode { BEST, PARTITION_INSIDE, PARTITION_OUTSIDE };

namespace linearx::utils {

struct ProbAccm {
    std::vector<value_type> upstrm;
    std::vector<value_type> dwnstrm;
};

inline std::unordered_map<char, int> VIENNA_NUC_ENCODING_SCHEME = {{'N', 0}, {'A', 1}, {'C', 2}, {'G', 3}, {'U', 4},
                                                                   {'-', 5}, {'n', 0}, {'a', 1}, {'c', 2}, {'g', 3},
                                                                   {'u', 4}, {'.', 5}};  // Vienna encoding scheme

// used for hashing a pair of integers
struct PairHash {
    std::size_t operator()(const std::pair<int, int>& p) const noexcept {
        return (static_cast<std::size_t>(p.first) << 32) | static_cast<std::size_t>(p.second);
    }
};

// A cache-friendly 2D view backed by a single contiguous 1D array.
// Provides row access via operator[] so callers can use a[i][j].
template <typename T>
class Flat2D {
   public:
    Flat2D() : rows_(0), cols_(0) {}
    Flat2D(const int rows, const int cols, const T& init = T()) { resize(rows, cols, init); }

    inline void resize(const int rows, const int cols, const T& init = T()) {
        rows_ = rows;
        cols_ = cols;
        data_.assign(static_cast<size_t>(rows) * static_cast<size_t>(cols), init);
    }

    inline int rows() const noexcept { return rows_; }
    inline int cols() const noexcept { return cols_; }
    inline size_t size() const noexcept { return data_.size(); }

    inline T* data() noexcept { return data_.data(); }
    inline const T* data() const noexcept { return data_.data(); }

    // Row proxy for a[i][j] indexing
    class Row {
       public:
        explicit Row(T* ptr) : ptr_(ptr) {}
        inline T& operator[](const int j) noexcept { return ptr_[j]; }
        inline const T& operator[](const int j) const noexcept { return ptr_[j]; }

       private:
        T* ptr_;
    };

    class ConstRow {
       public:
        explicit ConstRow(const T* ptr) : ptr_(ptr) {}
        inline const T& operator[](const int j) const noexcept { return ptr_[j]; }

       private:
        const T* ptr_;
    };

    inline Row operator[](const int i) noexcept { return Row(&data_[static_cast<size_t>(i) * cols_]); }
    inline ConstRow operator[](const int i) const noexcept {
        return ConstRow(&data_[static_cast<size_t>(i) * cols_]);
    }

    inline T& operator()(const int i, const int j) noexcept { return data_[static_cast<size_t>(i) * cols_ + j]; }
    inline const T& operator()(const int i, const int j) const noexcept {
        return data_[static_cast<size_t>(i) * cols_ + j];
    }

   private:
    int rows_;
    int cols_;
    std::vector<T> data_;
};

template <typename T>
inline void reset_beam_vector(std::vector<T>& vec, const int outer_size, const int inner_size = 0) {
    vec.resize(outer_size);
    if (inner_size > 0) {
        for (auto& map : vec) {
            map.clear();              // reuses internal memory
            map.reserve(inner_size);  // reserve space for efficiency (if needed)
        }
    } else {
        std::fill(vec.begin(), vec.end(), T());
    }
}

inline bool check_valid_pair(const int nuc1, const int nuc2) { return (nuc1 + nuc2) > 3 && (nuc1 + nuc2) % 2 != 0; }

template <typename T>
inline int quickselect_partition(std::vector<std::pair<value_type, T>>& scores, int lower, int upper) {
    value_type pivot = scores[upper].first;
    while (lower < upper) {
        while (scores[lower].first < pivot) ++lower;
        while (scores[upper].first > pivot) --upper;
        if (scores[lower].first == scores[upper].first)
            ++lower;
        else if (lower < upper)
            std::swap(scores[lower], scores[upper]);
    }
    return upper;
}

template <typename T>
inline value_type quickselect(std::vector<std::pair<value_type, T>>& scores, int lower, int upper, int k) {
    if (lower == upper) return scores[lower].first;
    int split = quickselect_partition(scores, lower, upper);
    int length = split - lower + 1;
    if (length == k)
        return scores[split].first;
    else if (k < length)
        return quickselect(scores, lower, split - 1, k);
    else
        return quickselect(scores, split + 1, upper, k - length);
}
}  // namespace linearx::utils

namespace linearx::utils::io {

/// Returns a writable file stream or throws on failure
/// @param fp File path to open
/// @param mode Mode in which to open the file (e.g., "r" for reading, "w" for writing)
/// @throws std::invalid_argument if fp or mode is nullptr
/// @throws std::runtime_error if the file cannot be opened
/// @return FILE pointer to the opened file
inline FILE* open_f(const char* fp, const char* mode) {
    if (!fp || !mode) throw std::invalid_argument("open_f received nullptr");

    FILE* f = fopen(fp, mode);
    if (!f) {
        std::string action = (mode[0] == 'r') ? "reading" : "writing";
        throw std::runtime_error("Could not open " + std::string(fp) + " for " + action);
    }

    return f;
}

/// Displays a progress bar in the terminal
/// @param cur_itr Current iteration number
/// @param total_itr Total number of iterations
/// @param step_percent Percentage step to update the progress bar (default is 5% increments)
/// @note It writes directly to stderr (file descriptor 2)
inline void showProgressBar(int cur_itr, int total_itr, int step_percent = 10) {
    static int last_percent_shown = -1;
    static char buffer[128];

    if (total_itr <= 0 || step_percent <= 0) return;

    int percent = (cur_itr * 100) / total_itr;
    if (percent == last_percent_shown || percent % step_percent != 0) return;
    last_percent_shown = percent;

    constexpr int barWidth = 50;
    int pos = (percent * barWidth) / 100;
    int offset = 0;

    // Start bar
    buffer[offset++] = '[';
    for (int i = 0; i < barWidth; ++i) {
        buffer[offset++] = (i < pos) ? '=' : (i == pos) ? '>' : ' ';
    }
    buffer[offset++] = ']';
    buffer[offset++] = ' ';

    // Add percentage as text: "100 %" or "  7 %"
    if (percent < 10) {
        buffer[offset++] = ' ';
        buffer[offset++] = ' ';
        buffer[offset++] = '0' + percent;
    } else if (percent < 100) {
        buffer[offset++] = ' ';
        buffer[offset++] = '0' + (percent / 10);
        buffer[offset++] = '0' + (percent % 10);
    } else {
        buffer[offset++] = '1';
        buffer[offset++] = '0';
        buffer[offset++] = '0';
    }
    buffer[offset++] = ' ';
    buffer[offset++] = '%';
    buffer[offset++] = '\r';

    // Write to stderr (fd = 2)
    write(2, buffer, offset);

    if (percent == 100) {
        write(2, "\n", 1);
        last_percent_shown = -1;  // Reset
    }
}

}  // namespace linearx::utils::io
