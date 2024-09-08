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
// #define INV_KT 0.016225317
#define INV_KT 1
#define GOLDEN_RATIO 0x9e3779b1
#define DEVIATION_THRESHOLD 9.91152


inline const std::unordered_map<char, int> nuc_encoding_scheme = {{'N', 0}, {'A', 1}, {'C', 2},
                                                              {'G', 3}, {'U', 4}, {'-', 5}, 
															  {'n', 0}, {'a', 1}, {'c', 2},
															  {'g', 3}, {'u', 4}, {'.', 5}}; // Vienna encoding scheme


const static double emit_probs[27][3] = 
{
     // ALN    // INS1   // INS2
	{0.134009, 0.000000, 0.000000}, // AA
	{0.027164, 0.000000, 0.000000}, // AC
	{0.049659, 0.000000, 0.000000}, // AG
	{0.028825, 0.000000, 0.000000}, // AU
	{0.000000, 0.211509, 0.000000}, // A.
	{0.027164, 0.000000, 0.000000}, // CA
	{0.140242, 0.000000, 0.000000}, // CC
	{0.037862, 0.000000, 0.000000}, // CG
	{0.047735, 0.000000, 0.000000}, // CU
	{0.000000, 0.257349, 0.000000}, // C.
	{0.049659, 0.000000, 0.000000}, // GA
	{0.037862, 0.000000, 0.000000}, // GC
	{0.178863, 0.000000, 0.000000}, // GG
	{0.032351, 0.000000, 0.000000}, // GU
	{0.000000, 0.271398, 0.000000}, // G.
	{0.028825, 0.000000, 0.000000}, // UA
	{0.047735, 0.000000, 0.000000}, // UC
	{0.032351, 0.000000, 0.000000}, // UG
	{0.099694, 0.000000, 0.000000}, // UU
	{0.000000, 0.259744, 0.000000}, // U.
	{0.000000, 0.000000, 0.211509}, // .A
	{0.000000, 0.000000, 0.257349}, // .C
	{0.000000, 0.000000, 0.271398}, // .G
	{0.000000, 0.000000, 0.259744}, // .U
	{0.000000, 0.000000, 0.000000}, // ..
	{1.000000, 0.000000, 0.000000}, // START
	{1.000000, 0.000000, 0.000000}, // END
};

const static double trans_probs[3][3] = 
{   // ALN     // INS1   // INS2
    {0.954668, 0.022666, 0.022666}, // ALN
	{0.292242, 0.666439, 0.041319}, // INS1
	{0.292242, 0.041319, 0.666439}, // INS2
};

#endif // CONSTANTS_HPP

// AC 3 (invalid)
// AG 4 (invalid)
// AU 5 (valid)
// CG 5 (valid)
// CU 6 (invalid)
// GU 7 (valid)

// pairs: 0:NP 1:CG 2:GC 3:GU 4:UG 5:AU 6:UA 7:NN
// nucleotides: CONTRAfold: 0:A 1:C 2:G 3:U 4:N ; Vienna: 0:N 1:A 2:C 3:G 4:U