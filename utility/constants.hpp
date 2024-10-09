#ifndef CONSTANTS_HPP
#define CONSTANTS_HPP

#include <limits>
#include <unordered_map>

// macro definitions
#define VALUE_MIN std::numeric_limits<float>::lowest()
#define VALUE_MAX std::numeric_limits<float>::max()
#define POS_INF std::numeric_limits<int>::max()
#define NEG_INF std::numeric_limits<int>::lowest()

#define MAXLOOPSIZE 30
#define LXC37 107.856
#define kT 61.63207755
#define INV_KT 0.016225317
// #define INV_KT 1
#define GOLDEN_RATIO 0x9e3779b1
#define DEVIATION_THRESHOLD 9.91152

inline std::unordered_map<char, int> nuc_encoding_scheme = {{'N', 0}, {'A', 1}, {'C', 2}, {'G', 3}, {'U', 4},
                                                            {'-', 5}, {'n', 0}, {'a', 1}, {'c', 2}, {'g', 3},
                                                            {'u', 4}, {'.', 5}}; // Vienna encoding scheme

#endif // CONSTANTS_HPP

// AC 3 (invalid)
// AG 4 (invalid)
// AU 5 (valid)
// CG 5 (valid)
// CU 6 (invalid)
// GU 7 (valid)

// pairs: 0:NP 1:CG 2:GC 3:GU 4:UG 5:AU 6:UA 7:NN
// nucleotides: CONTRAfold: 0:A 1:C 2:G 3:U 4:N ; Vienna: 0:N 1:A 2:C 3:G 4:U