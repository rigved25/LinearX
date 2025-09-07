// linearx/config.hpp
#pragma once

// #define FAST_MATH  // comment this line to use value_type as double
#define FAST_LOG   // comment this line to use approximate log from log_math.hpp

#ifdef FAST_MATH
using value_type = float;
#else
using value_type = double;
#endif
