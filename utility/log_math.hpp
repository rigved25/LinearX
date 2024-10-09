#ifndef LOG_MATH_HPP
#define LOG_MATH_HPP

#include <cassert>
#include <cfloat>
#include <cmath>
#include <iostream>
#include <limits>
#include <stdexcept>

#include "./constants.hpp"

// This header contains extended logarithmic/exponential functions derived from
// methods described in TP Mann's 2006 paper:
//   "Numerically stable hidden Markov model implementation." (Tobias P. Mann, 2006)
//   http://bozeman.genome.washington.edu/compbio/mbt599_2006/hmm_scaling_revised.pdf
// Our functions are more optimized than those Mann described, and they use a
// different value to represent LOG_OF_ZERO.
// The functions here are used to implement PPF_Math, which is used in the PHMM forward-backward algorithm.
// They were later adapted to implement the log-scale calculations used in the partition function (pfunction).

#define USE_XLOG_ZERO

#ifdef USE_XLOG_ZERO
static const double LOG_OF_ZERO = -log(DBL_MAX) * 1000;
#define IS_LOG_ZERO(x) (x <= LOG_OF_ZERO)
#else
// Log-scale values less than or equal to LOG_OF_ZERO are assumed to represent zero.
static const double LOG_OF_ZERO = -std::numeric_limits<double>::infinity();
#define IS_LOG_ZERO(x) false
// #define LOG_OF_ZERO -1073741824 // -2^(30), can be represented exactly.
#endif
static const double LOG_OF_ONE = 0.0;

// Calculates log1p(exp(x))  -- used for calculating xlog_sum.
inline double log1pexp(const double &x) noexcept {
    if (IS_LOG_ZERO(x))
        return 0.0; // if x <= LOG_OF_ZERO, then exp(x) is 0, and log(1+0)=log(1)=0, so return 0.
#ifdef USE_LOGP1_LOOKUP_SUM
    return LogTable.exp_lu(x);
#else
    return log1p(exp(x));
#endif
}

// calculates log1p(-exp(x))  -- used for calculating xlog_sub.
inline double log1pNexp(const double &x) noexcept {
    if (IS_LOG_ZERO(x))
        return 0.0; // if x <= LOG_OF_ZERO, then exp(x) is 0, and log(1-0)=log(1)=0, so return 0.
#ifdef USE_LOGP1_LOOKUP_DIFF
    return LogTableDiff.exp_lu(x);
#else
    return log1p(-exp(x));
#endif
}

// Convert a number (in standard, linear scale) to a log-scale value (i.e. stored as the natural logarithm of the
// original number).
double xlog(const double &linear_value) noexcept;
// Convert a log-scale value back to its original linear-scale number (by taking the exponent of it, base e).
double xexp(const double &xlog) noexcept;
// Calculates the product of two log-scale values.
double xlog_mul(const double &xlog1, const double &xlog2);
// Calculates the quotient of two log-scale values.
double xlog_div(const double &xlog1, const double &xlog2);
// Calculates the sum of two log-scale values.
double xlog_sum(const double &xlog1, const double &xlog2);
// Calculates the difference between two log-scale values.
double xlog_sub(const double &xlog1, const double &xlog2);
// Calculates the result of raising a log-scale value to a power. (The power should be in linear scale.)
double xlog_pow(const double &xlog, const double &power);

inline double xlog(const double &value) noexcept {
#ifdef USE_XLOG_ZERO
    if (value == 0.0)
        return LOG_OF_ZERO; // the log function itself will throw an exception if value < 0.
#endif
    return log(value); // the log function itself will throw an exception if value < 0.
}

inline double xexp(const double &xlog) noexcept {
#ifdef USE_XLOG_ZERO
    if (xlog <= LOG_OF_ZERO)
        return 0.0;
#endif
    return exp(xlog);
}

// Computes log(exp(a)+exp(b))
inline double xlog_sum(const double &a, const double &b) {
    // Derivation:   log(exp(a)-exp(b))
    //             = log(exp(a) * (1-exp(b-a)) )
    //             = a+log(1-exp(b-a))
    //             = a+log1p(-exp(b-a))
#ifdef USE_XLOG_ZERO
    if (IS_LOG_ZERO(a))
        return b;
    if (IS_LOG_ZERO(b))
        return a;
#endif
    // Note: The test of a>b is important when A or B is greater than MAX_DOUBLE.
    // e.g. A=1E+500 and B=1 --- exp(a-b) will overflow, while exp(b-a) will not.
    // As long as a>=b, the value of exp(b-a) will always be in the range (0,1],
    // so it will not overflow (even when one or both numbers are greater than
    // MAX_DOUBLE in the linear scale).
    // In fact, the result will always be max(a,b)+Q where Q is in the range [0..ln(2)].
    return a > b ? a + log1p(exp(b - a)) : b + log1p(exp(a - b));
    /*
    // -----Possible alternatives (todo: test for speed)-----
// -----alt.1----- (temporary vars for min,max)
    const double mx = std::max(a,b), mn=std::min(a,b);
return mx+log1pexp(mn-mx);
    // -----alt.2----- (implicit compiler-defined temporary vars for min,max)
return std::max(a,b)+log1pexp(std::min(a,b)-std::max(a,b));
    // -----alt.3----- (max and -abs)
    return std::max(a,b)+log1pexp(-std::abs(a-b));
    */
}

inline double xlog_sum2(const double &a, const double &b) {
    // Derivation:   log(exp(a)-exp(b))
    //             = log(exp(a) * (1-exp(b-a)) )
    //             = a+log(1-exp(b-a))
    //             = a+log1p(-exp(b-a))
#ifdef USE_XLOG_ZERO
    if (IS_LOG_ZERO(a))
        return b;
    if (IS_LOG_ZERO(b))
        return a;
#endif
    // Note: The test of a>b is important when A or B is greater than MAX_DOUBLE.
    // e.g. A=1E+500 and B=1 --- exp(a-b) will overflow, while exp(b-a) will not.
    // As long as a>=b, the value of exp(b-a) will always be in the range (0,1],
    // so it will not overflow (even when one or both numbers are greater than
    // MAX_DOUBLE in the linear scale).
    // In fact, the result will always be max(a,b)+Q where Q is in the range [0..ln(2)].
    return a > b ? a + log1pexp(b - a) : b + log1pexp(a - b);
    /*
    // -----Possible alternatives (todo: test for speed)-----
// -----alt.1----- (temporary vars for min,max)
    const double mx = std::max(a,b), mn=std::min(a,b);
return mx+log1pexp(mn-mx);
    // -----alt.2----- (implicit compiler-defined temporary vars for min,max)
return std::max(a,b)+log1pexp(std::min(a,b)-std::max(a,b));
    // -----alt.3----- (max and -abs)
    return std::max(a,b)+log1pexp(-std::abs(a-b));
    */
}
inline double xlog_sum3(const double &a, const double &b) {
    double c, d, x;
    if (a < b) {
        c = b;
        d = a;
    } else {
        c = a;
        d = b;
    }
    x = c - d;

    if (x > 11.8624794162)
        return d + x;
    if (x < double(3.3792499610)) {
        if (x < double(1.6320158198)) {
            if (x < double(0.6615367791))
                return d + ((double(-0.0065591595) * x + double(0.1276442762)) * x + double(0.4996554598)) * x +
                       double(0.6931542306);
            return d + ((double(-0.0155157557) * x + double(0.1446775699)) * x + double(0.4882939746)) * x +
                   double(0.6958092989);
        }
        if (x < double(2.4912588184))
            return d + ((double(-0.0128909247) * x + double(0.1301028251)) * x + double(0.5150398748)) * x +
                   double(0.6795585882);
        return d + ((double(-0.0072142647) * x + double(0.0877540853)) * x + double(0.6208708362)) * x +
               double(0.5909675829);
    }
    if (x < double(5.7890710412)) {
        if (x < double(4.4261691294))
            return d + ((double(-0.0031455354) * x + double(0.0467229449)) * x + double(0.7592532310)) * x +
                   double(0.4348794399);
        return d + ((double(-0.0010110698) * x + double(0.0185943421)) * x + double(0.8831730747)) * x +
               double(0.2523695427);
    }
    if (x < double(7.8162726752))
        return d + ((double(-0.0001962780) * x + double(0.0046084408)) * x + double(0.9634431978)) * x +
               double(0.0983148903);
    return d + ((double(-0.0000113994) * x + double(0.0003734731)) * x + double(0.9959107193)) * x +
           double(0.0149855051);
}

// Computes log(exp(a)-exp(b))
inline double xlog_sub(const double &a, const double &b) { // (previously xlog_sub)
// Derivation:   log(exp(a)-exp(b))
//             = log(exp(a) * (1-exp(b-a)) )
//             = a+log(1-exp(b-a))
//             = a+log1p(-exp(b-a))
#ifdef USE_XLOG_ZERO
    // `a` must be >= `b` so that exp(a) >= exp(b), otherwise the value inside the log would be negative.
    if (b <= LOG_OF_ZERO)
        return a; // Important to do this test, because (a<b) is NOT an error if (b <= LOG_OF_ZERO)
    if (a < b)
        throw std::runtime_error(
            "Subtraction of xlog values resulted in an unrepresentable negative number. (in " __FILE__ ")");
    return a == b ? LOG_OF_ZERO : a + log1pNexp(b - a); // note that (b-a) is always < 0, which is required because we
                                                        // are computing log(1-exp(b-a))
#else
    return a + log1pNexp(
                   b - a); // note that (b-a) is always < 0, which is required because we are computing log(1-exp(b-a))
#endif
}

inline double xlog_mul(const double &log1, const double &log2) {
#ifdef USE_XLOG_ZERO
    return (log1 <= LOG_OF_ZERO || log2 <= LOG_OF_ZERO) ? LOG_OF_ZERO : log1 + log2;
#else
    return log1 + log2;
#endif
}

// Returns 0 if log1 is 0 no matter what log2 is.
inline double xlog_div(const double &log1, const double &log2) {
#ifdef USE_XLOG_ZERO
    if (log1 <= LOG_OF_ZERO)
        return LOG_OF_ZERO;
    if (log2 <= LOG_OF_ZERO)
        throw std::runtime_error("Division by xlog zero-value (in " __FILE__ ")");
#endif
    return log1 - log2;
}

inline double xlog_pow(const double &log_value, const double &pow) {
#ifdef USE_XLOG_ZERO
    return log_value <= LOG_OF_ZERO ? LOG_OF_ZERO : log_value * pow;
#else
    return log_value * pow;
#endif
}

inline double Fast_LogExpPlusOne(double x) {

    // Bounds for tolerance of 7.05e-06: (0, 11.8625)
    // Approximating interval: (0, 0.661537) -->
    // ((T(-0.0065591595)*x+T(0.1276442762))*x+T(0.4996554598))*x+T(0.6931542306); Approximating interval:
    // (0.661537, 1.63202) --> ((T(-0.0155157557)*x+T(0.1446775699))*x+T(0.4882939746))*x+T(0.6958092989); Approximating
    // interval: (1.63202, 2.49126) --> ((T(-0.0128909247)*x+T(0.1301028251))*x+T(0.5150398748))*x+T(0.6795585882);
    // Approximating interval: (2.49126, 3.37925) -->
    // ((T(-0.0072142647)*x+T(0.0877540853))*x+T(0.6208708362))*x+T(0.5909675829); Approximating interval:
    // (3.37925, 4.42617) --> ((T(-0.0031455354)*x+T(0.0467229449))*x+T(0.7592532310))*x+T(0.4348794399); Approximating
    // interval: (4.42617, 5.78907) --> ((T(-0.0010110698)*x+T(0.0185943421))*x+T(0.8831730747))*x+T(0.2523695427);
    // Approximating interval: (5.78907, 7.81627) -->
    // ((T(-0.0001962780)*x+T(0.0046084408))*x+T(0.9634431978))*x+T(0.0983148903); Approximating interval:
    // (7.81627, 11.8625) --> ((T(-0.0000113994)*x+T(0.0003734731))*x+T(0.9959107193))*x+T(0.0149855051); 8 polynomials
    // needed.

    assert(double(0.0000000000) <= x && x <= double(11.8624794162) && "Argument out-of-range.");
    if (x < double(3.3792499610)) {
        if (x < double(1.6320158198)) {
            if (x < double(0.6615367791))
                return ((double(-0.0065591595) * x + double(0.1276442762)) * x + double(0.4996554598)) * x +
                       double(0.6931542306);
            return ((double(-0.0155157557) * x + double(0.1446775699)) * x + double(0.4882939746)) * x +
                   double(0.6958092989);
        }
        if (x < double(2.4912588184))
            return ((double(-0.0128909247) * x + double(0.1301028251)) * x + double(0.5150398748)) * x +
                   double(0.6795585882);
        return ((double(-0.0072142647) * x + double(0.0877540853)) * x + double(0.6208708362)) * x +
               double(0.5909675829);
    }
    if (x < double(5.7890710412)) {
        if (x < double(4.4261691294))
            return ((double(-0.0031455354) * x + double(0.0467229449)) * x + double(0.7592532310)) * x +
                   double(0.4348794399);
        return ((double(-0.0010110698) * x + double(0.0185943421)) * x + double(0.8831730747)) * x +
               double(0.2523695427);
    }
    if (x < double(7.8162726752))
        return ((double(-0.0001962780) * x + double(0.0046084408)) * x + double(0.9634431978)) * x +
               double(0.0983148903);
    return ((double(-0.0000113994) * x + double(0.0003734731)) * x + double(0.9959107193)) * x + double(0.0149855051);

    /*
    // Bounds for tolerance of 9.99e-05: (0, 9.21129)
    // Approximating interval: (0, 1.40131) -->
    ((T(-0.0118287252)*x+T(0.1342168806))*x+T(0.4976005362))*x+T(0.6932470806);
    // Approximating interval: (1.40131, 3.06792) -->
    ((T(-0.0117040733)*x+T(0.1232945547))*x+T(0.5276092444))*x+T(0.6721240615);
    // Approximating interval: (3.06792, 5.15409) -->
    ((T(-0.0027005983)*x+T(0.0419040665))*x+T(0.7762991688))*x+T(0.4152395732);
    // Approximating interval: (5.15409, 9.21129) -->
    ((T(-0.0001617326)*x+T(0.0040111354))*x+T(0.9666890441))*x+T(0.0929363811);
    // 4 polynomials needed.

    assert(double (0.0000000000) <= x && x <= double (9.2112909219), "Argument out-of-range.");
    if (x < double (3.0679202382))
    {
        if (x < double (1.4013117629))
            return ((double (-0.0118287252)*x+double (0.1342168806))*x+double (0.4976005362))*x+double (0.6932470806);
        return ((double (-0.0117040733)*x+double (0.1232945547))*x+double (0.5276092444))*x+double (0.6721240615);
    }
    if (x < double (5.1540922927))
        return ((double (-0.0027005983)*x+double (0.0419040665))*x+double (0.7762991688))*x+double (0.4152395732);
    return ((double (-0.0001617326)*x+double (0.0040111354))*x+double (0.9666890441))*x+double (0.0929363811);
    */
}

inline void Fast_LogPlusEquals(double &x, double y) {
    if (x < y)
        std::swap(x, y);
    if (y > double(NEG_INF / 2) && x - y < double(11.8624794162))
        x = Fast_LogExpPlusOne(x - y) + y;
}

inline double Fast_Exp(double x) {
    // Bounds for tolerance of 4.96e-05: (-9.91152, 0)
    // Approximating interval: (-9.91152, -5.86228) -->
    // ((T(0.0000803850)*x+T(0.0021627428))*x+T(0.0194708555))*x+T(0.0588080014); Approximating interval: (-5.86228,
    // -3.83966) --> ((T(0.0013889414)*x+T(0.0244676474))*x+T(0.1471290604))*x+T(0.3042757740); Approximating interval:
    // (-3.83966, -2.4915) --> ((T(0.0072335607)*x+T(0.0906002677))*x+T(0.3983111356))*x+T(0.6245959221); Approximating
    // interval: (-2.4915, -1.48054) --> ((T(0.0232410351)*x+T(0.2085645908))*x+T(0.6906367911))*x+T(0.8682322329);
    // Approximating interval: (-1.48054, -0.672505) -->
    // ((T(0.0573782771)*x+T(0.3580258429))*x+T(0.9121133217))*x+T(0.9793091728); Approximating interval: (-0.672505,
    // -3.9145e-11) --> ((T(0.1199175927)*x+T(0.4815668234))*x+T(0.9975991939))*x+T(0.9999505077); 6 polynomials needed.

    if (x < double(-2.4915033807)) {
        if (x < double(-5.8622823336)) {
            if (x < double(-9.91152))
                return double(0);
            return ((double(0.0000803850) * x + double(0.0021627428)) * x + double(0.0194708555)) * x +
                   double(0.0588080014);
        }
        if (x < double(-3.8396630909))
            return ((double(0.0013889414) * x + double(0.0244676474)) * x + double(0.1471290604)) * x +
                   double(0.3042757740);
        return ((double(0.0072335607) * x + double(0.0906002677)) * x + double(0.3983111356)) * x +
               double(0.6245959221);
    }
    if (x < double(-0.6725053211)) {
        if (x < double(-1.4805375919))
            return ((double(0.0232410351) * x + double(0.2085645908)) * x + double(0.6906367911)) * x +
                   double(0.8682322329);
        return ((double(0.0573782771) * x + double(0.3580258429)) * x + double(0.9121133217)) * x +
               double(0.9793091728);
    }
    if (x < double(0))
        return ((double(0.1199175927) * x + double(0.4815668234)) * x + double(0.9975991939)) * x +
               double(0.9999505077);
    return (x > double(46.052) ? double(1e20) : expf(x));
}

#endif // LOG_MATH_HPP