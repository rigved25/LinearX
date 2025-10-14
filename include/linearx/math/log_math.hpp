// linearx/math/log_math.hpp
#pragma once

#include <cassert>
#include <cfloat>
#include <cmath>
#include <iostream>
#include <limits>
#include <linearx/constants.hpp>
#include <stdexcept>

namespace linearx::math {

static const value_type LOG_ZERO = -log(DBL_MAX) * 1000;
static const value_type LOG_ONE = 0.0;
#define IS_LOG_ZERO(x) (x <= LOG_ZERO)

// Calculates log1p(exp(x))  -- used for calculating xlog_sum.
inline value_type log1pexp(value_type x) noexcept {
    if (IS_LOG_ZERO(x)) return 0.0;  // if x <= LOG_ZERO, then exp(x) is 0, and log(1+0)=log(1)=0, so return 0.
#ifdef USE_LOGP1_LOOKUP_SUM
    return LogTable.exp_lu(x);
#else
    return log1p(exp(x));
#endif
}

// calculates log1p(-exp(x))  -- used for calculating xlog_sub.
inline value_type log1pNexp(value_type x) noexcept {
    if (IS_LOG_ZERO(x)) return 0.0;  // if x <= LOG_ZERO, then exp(x) is 0, and log(1-0)=log(1)=0, so return 0.
#ifdef USE_LOGP1_LOOKUP_DIFF
    return LogTableDiff.exp_lu(x);
#else
    return log1p(-exp(x));
#endif
}

inline value_type xlog(value_type value) noexcept {
    if (value == 0.0) return LOG_ZERO;  // the log function itself will throw an exception if value < 0.
    return log(value);                  // the log function itself will throw an exception if value < 0.
}

inline value_type xexp(value_type xlog) noexcept {
    if (xlog <= LOG_ZERO) return LOG_ONE;
    return exp(xlog);
}

// Computes log(exp(a)+exp(b))
inline value_type xlog_sum(value_type a, value_type b) {
    if (a < b) std::swap(a, b);  // now: a >= b (a is max)
    if (IS_LOG_ZERO(b)) return a;
    return a + log1p(exp(b - a));
}

inline value_type xlog_sum2(value_type a, value_type b) {
    if (a < b) std::swap(a, b);  // now: a >= b (a is max)
    if (IS_LOG_ZERO(b)) return a;
    return a + log1pexp(b - a);  // b - a <= 0
}

// Computes log(exp(a)-exp(b))
inline value_type xlog_sub(value_type a, value_type b) {
    // if b is effectively 0 (i.e., LOG_ZERO), then result is a
    if (IS_LOG_ZERO(b)) return a;
    // undefined: can't subtract larger from smaller in log-space
    if (a < b) {
        throw std::runtime_error("xlog_sub: log(exp(a) - exp(b)) is undefined because a < b (in " __FILE__ ":" +
                                 std::to_string(__LINE__) + ")");
    }
    // if a == b, result is log(0)
    return (a == b) ? LOG_ZERO : a + log1pNexp(b - a);
}

inline value_type xlog_mul(value_type log1, value_type log2) {
    return (IS_LOG_ZERO(log1) || IS_LOG_ZERO(log2)) ? LOG_ZERO : log1 + log2;
}

// Returns 0 if log1 is 0 no matter what log2 is.
inline value_type xlog_div(value_type log1, value_type log2) {
    if (IS_LOG_ZERO(log1)) return LOG_ZERO;
    if (IS_LOG_ZERO(log2)) return -100 * LOG_ZERO;
    return log1 - log2;
}

inline value_type xlog_pow(value_type log_value, value_type pow) {
    return IS_LOG_ZERO(log_value) ? LOG_ZERO : log_value * pow;
}

inline value_type Fast_LogExpPlusOne(value_type x) {
    assert(value_type(0.0000000000) <= x && x <= value_type(11.8624794162) && "Argument out-of-range.");
    if (x < value_type(3.3792499610)) {
        if (x < value_type(1.6320158198)) {
            if (x < value_type(0.6615367791))
                return ((value_type(-0.0065591595) * x + value_type(0.1276442762)) * x + value_type(0.4996554598)) * x +
                       value_type(0.6931542306);
            return ((value_type(-0.0155157557) * x + value_type(0.1446775699)) * x + value_type(0.4882939746)) * x +
                   value_type(0.6958092989);
        }
        if (x < value_type(2.4912588184))
            return ((value_type(-0.0128909247) * x + value_type(0.1301028251)) * x + value_type(0.5150398748)) * x +
                   value_type(0.6795585882);
        return ((value_type(-0.0072142647) * x + value_type(0.0877540853)) * x + value_type(0.6208708362)) * x +
               value_type(0.5909675829);
    }
    if (x < value_type(5.7890710412)) {
        if (x < value_type(4.4261691294))
            return ((value_type(-0.0031455354) * x + value_type(0.0467229449)) * x + value_type(0.7592532310)) * x +
                   value_type(0.4348794399);
        return ((value_type(-0.0010110698) * x + value_type(0.0185943421)) * x + value_type(0.8831730747)) * x +
               value_type(0.2523695427);
    }
    if (x < value_type(7.8162726752))
        return ((value_type(-0.0001962780) * x + value_type(0.0046084408)) * x + value_type(0.9634431978)) * x +
               value_type(0.0983148903);
    return ((value_type(-0.0000113994) * x + value_type(0.0003734731)) * x + value_type(0.9959107193)) * x +
           value_type(0.0149855051);
}

inline value_type Fast_LogPlusEquals(value_type x, value_type y) {
    if (x < y) {
        std::swap(x, y);
    }  // now: x >= y (x is max)
    if (y > value_type(LOG_ZERO) && x - y < value_type(11.8624794162)) {
        x = Fast_LogExpPlusOne(x - y) + y;
    }
    return x;
}

// Approximating interval: (-9.91152, -5.86228) --> ((T(0.0000803850)*x+T(0.0021627428))*x+T(0.0194708555))*x+T(0.0588080014);
// Approximating interval: (-5.86228, -3.83966) --> ((T(0.0013889414)*x+T(0.0244676474))*x+T(0.1471290604))*x+T(0.3042757740);
// Approximating interval: (-3.83966, -2.4915) --> ((T(0.0072335607)*x+T(0.0906002677))*x+T(0.3983111356))*x+T(0.6245959221);
// Approximating interval: (-2.4915, -1.48054) --> ((T(0.0232410351)*x+T(0.2085645908))*x+T(0.6906367911))*x+T(0.8682322329);
// Approximating interval: (-1.48054, -0.672505) --> ((T(0.0573782771)*x+T(0.3580258429))*x+T(0.9121133217))*x+T(0.9793091728);
// Approximating interval: (-0.672505, -3.9145e-11) --> ((T(0.1199175927)*x+T(0.4815668234))*x+T(0.9975991939))*x+T(0.9999505077);
// 6 polynomials needed.
inline value_type Fast_Exp(const value_type x) {
    if (x < value_type(-2.4915033807)) {
        if (x < value_type(-5.8622823336)) {
            if (x < value_type(-9.91152)) return value_type(0);
            return ((value_type(0.0000803850) * x + value_type(0.0021627428)) * x + value_type(0.0194708555)) * x +
                   value_type(0.0588080014);
        }
        if (x < value_type(-3.8396630909))
            return ((value_type(0.0013889414) * x + value_type(0.0244676474)) * x + value_type(0.1471290604)) * x +
                   value_type(0.3042757740);
        return ((value_type(0.0072335607) * x + value_type(0.0906002677)) * x + value_type(0.3983111356)) * x +
               value_type(0.6245959221);
    }
    if (x < value_type(-0.6725053211)) {
        if (x < value_type(-1.4805375919))
            return ((value_type(0.0232410351) * x + value_type(0.2085645908)) * x + value_type(0.6906367911)) * x +
                   value_type(0.8682322329);
        return ((value_type(0.0573782771) * x + value_type(0.3580258429)) * x + value_type(0.9121133217)) * x +
               value_type(0.9793091728);
    }
    if (x < value_type(0))
        return ((value_type(0.1199175927) * x + value_type(0.4815668234)) * x + value_type(0.9975991939)) * x +
               value_type(0.9999505077);
    return (x > value_type(46.052) ? value_type(1e20) : expf(x));
}
}  // namespace linearx::math

// Wrapper functions for easier use in the rest of the codebase
FORCE_INLINE value_type LOG(value_type x) { return linearx::math::xlog(x); }
FORCE_INLINE value_type EXP(value_type x) {
#ifdef FAST_LOG
    return linearx::math::Fast_Exp(x);
#else
    return linearx::math::xexp(x);
#endif
}
FORCE_INLINE value_type LOG_SUM(value_type a, value_type b) {
#ifdef FAST_LOG
    return linearx::math::Fast_LogPlusEquals(a, b);
#else
    return linearx::math::xlog_sum(a, b);
#endif
}
FORCE_INLINE value_type LOG_MUL(value_type a, value_type b) { return linearx::math::xlog_mul(a, b); }
FORCE_INLINE value_type LOG_DIV(value_type a, value_type b) { return linearx::math::xlog_div(a, b); }
