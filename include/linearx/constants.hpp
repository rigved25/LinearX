// linearx/constants.hpp
#pragma once

#include <limits>
#include <linearx/config.hpp>

namespace linearx::constants::math {

inline constexpr int INT_POS_INF = std::numeric_limits<int>::max();
inline constexpr int INT_NEG_INF = std::numeric_limits<int>::lowest();
inline constexpr value_type VALUE_MIN = std::numeric_limits<value_type>::lowest();
inline constexpr value_type VALUE_MAX = std::numeric_limits<value_type>::max();
inline constexpr float GOLDEN_RATIO = 0x9e3779b1;

}  // namespace linearx::constants::math

namespace linearx::constants::energy {

inline constexpr int MIN_HELIX_SIZE = 3;
inline constexpr int MAXLOOPSIZE = 30;
inline constexpr value_type LXC37 = 107.856;
inline constexpr value_type kT = 61.63207755;
inline constexpr value_type INV_KT = 0.016225317;

}  // namespace linearx::constants::energy

namespace linearx::constants::limits {

inline constexpr value_type DEVIATION_THRESHOLD = 9.91152;
inline constexpr value_type EPSILON = 1e-300;

}  // namespace linearx::constants::limits

namespace linearx::constants::alignment {

/// Above this seq_len_sum (seq1.len + seq2.len) we use swap in save_partition_function; below it we use copy/move.
inline constexpr int SAVE_SWAP_THRESHOLD = 1500;

}  // namespace linearx::constants::alignment
