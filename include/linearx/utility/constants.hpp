// linearx/utility/constants.hpp
#pragma once
#include <limits>

namespace linearx::constants::math {

inline constexpr int INT_POS_INF = std::numeric_limits<int>::max();
inline constexpr int INT_NEG_INF = std::numeric_limits<int>::lowest();
inline constexpr float VALUE_MIN = std::numeric_limits<float>::lowest();
inline constexpr float VALUE_MAX = std::numeric_limits<float>::max();
inline constexpr float GOLDEN_RATIO = 0x9e3779b1;

}  // namespace linearx::constants

namespace linearx::constants::energy {

inline constexpr int MIN_HELIX_SIZE = 3;
inline constexpr int MAXLOOPSIZE = 30;
inline constexpr double LXC37 = 107.856;
inline constexpr double kT = 61.63207755;
inline constexpr double INV_KT = 0.016225317;

}  // namespace linearx::constants::energy

namespace linearx::constants::limits {

inline constexpr double DEVIATION_THRESHOLD = 9.91152;
inline constexpr double EPSILON = 1e-300;

}  // namespace linearx::constants::limits
