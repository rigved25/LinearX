#include "./../../sequence/multi_seq.hpp"
#include "./../../turbofold/turbofold.hpp"
#include <iostream>

using namespace std;

int main() {
    MultiSeq mseq;
    // mseq.add_sequence(Seq("seq1", "CCCAAAGGG", 0, &nuc_encoding_scheme));
    // mseq.add_sequence(Seq("seq2", "GGGAAUACCC", 0, &nuc_encoding_scheme));

    mseq.add_sequence(Seq("seq1", "ugccuggcggccguagcgcgguggucccaccugaccccaugccgaacucagaagugaaacgccguagcgccgaugguaguguggggucuccccaugcgagaguagggaacugccaggcau", 0, &nuc_encoding_scheme));
    mseq.add_sequence(Seq("seq2", "ugguggcgaugcgcuuaggggaaacacccguucccaucucgaacacgaugguuaagacuuaagcggccgaugguacuuugcuggagacggcacgggagaguagguggcugccagauu", 0, &nuc_encoding_scheme));

    LinearTurboFold ltf(&mseq, 1);
    ltf.run();

    return 0;
}