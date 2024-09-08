#include "./../../Seq.hpp"
#include "./../../LinearAlign.hpp"
#include <iostream>

using namespace std;

int main() {
    cout << "Hello, World!" << endl;


    Seq seq1("seq1", "gcggggguggcccagccugguacggcgugggacugcuaaucccauugggcaacgcccagcccggguucaagucccggcccccgcg");
    Seq seq2("seq2", "gcccggguggcggaauggcagacgcgcuagcuugaggugcuaguguccuauuaacggacguggggguucaaguccccccccgggca");

    LinearAlign linear_align(seq1, seq2);

    int l1 = seq1.length();
    int l2 = seq2.length();
    
    // alpha
    double alpha_start = linear_align.beamALN[0][{0, 0}].alpha;
    double alpha_end = linear_align.beamALN[l1 + l2 + 2][{l1 + 1, l2 + 1}].alpha;

    // beta
    double beta_start = linear_align.beamALN[l1 + l2 + 2][{l1 + 1, l2 + 1}].beta;
    double beta_end = linear_align.beamALN[0][{0, 0}].beta;


    cout << "Alpha Start: " << alpha_start << endl;
    cout << "Alpha End: " << alpha_end << endl;

    cout << "Beta Start: " << beta_start << endl;
    cout << "Beta End: " << beta_end << endl;

    return 0;
}

// "--"
// "AC"


