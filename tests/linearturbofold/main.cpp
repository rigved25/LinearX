#include "./../../sequence/multi_seq.hpp"
#include "./../../turbofold/turbofold.hpp"
#include <iostream>
#include <string>

using namespace std;

int main(int argc, char* argv[]) {
    // Check if the number of arguments is correct (must be an even number for name-sequence pairs)
    if (argc < 3 || (argc - 1) % 2 != 0) {
        cerr << "Usage: " << argv[0] << " seq_name1 sequence1 seq_name2 sequence2 ..." << endl;
        return 1;
    }

    // Initialize MultiSeq object
    MultiSeq mseq;

    // Number of sequences
    int num_sequences = (argc - 1) / 2;

    // Add each sequence and its corresponding name from the command-line arguments
    for (int i = 0; i < num_sequences; ++i) {
        string name = argv[2 * i + 1];       // Sequence name
        string sequence = argv[2 * i + 2];   // Sequence string

        // Add the sequence to the MultiSeq object
        mseq.add_sequence(Seq(name, sequence, i, &nuc_encoding_scheme));
    }

    // Initialize LinearTurboFold object with the input sequences
    LinearTurboFold ltf(&mseq, 1, VerboseState::SILENT);

    // Run the LinearTurboFold algorithm
    ltf.run();

    return 0;


    // // std::cout << Fast_Exp(VALUE_MIN) << std::endl;
    // MultiSeq mseq;
    // // mseq.add_sequence(Seq("seq1", "CCCAAAGGG", 0, &nuc_encoding_scheme));
    // // mseq.add_sequence(Seq("seq2", "GGGAAUACCC", 0, &nuc_encoding_scheme));


    // // mseq.add_sequence(Seq("seq1", "ccgggugccuauaccggaggggccacacccguucccauuccgaacacggucguuaagcccuccagggccgaugguacuggggcguuaccgcccugggagaguaggucggugcccggg", 0, &nuc_encoding_scheme));
    // mseq.add_sequence(Seq("seq1", "ugcuuggcgaccauagcguuauggacccaccugaucccaugccgaacucaguagugaaacguaauagcgccgaugguaguguggggucuccccaugugagaguaggacaucgccaggcau", 0, &nuc_encoding_scheme));
    // mseq.add_sequence(Seq("seq2", "ccgggugccuauaccggaggggccacacccguucccauuccgaacacggucguuaagcccuccagggccgaugguacuggggcguuaccgcccugggagaguaggucggugcccggg", 1, &nuc_encoding_scheme));
    // // mseq.add_sequence(Seq("seq3", "ccgggugccuauaccggaggggccacacccguucccauuccgaacacggucguuaagcccuccagggccgaugguacuggggcguuaccgcccugggagaguaggucggugccggg", 2, &nuc_encoding_scheme));
    
    

    // LinearTurboFold ltf(&mseq, 1, VerboseState::SILENT);
    // ltf.run();

    // return 0;
}