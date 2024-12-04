#include <iostream>
#include <string>

#include "./../../sequence/multi_seq.hpp"
#include "./../../turbofold/turbofold.hpp"

using namespace std;

MultiSeq read_msa_file(const std::string& filePath, std::unordered_map<char, int>* encoding_scheme = nullptr) {
    std::ifstream file(filePath);
    if (!file.is_open()) {
        throw std::runtime_error("Could not open file: " + filePath);
    }

    MultiSeq mseq;
    std::string line, currentName, currentDotBracket;

    while (std::getline(file, line)) {
        if (line.empty()) {
            continue;  // Skip empty lines
        }
        if (line[0] == '>') {  // Header line
            if (!currentName.empty() && !currentDotBracket.empty()) {
                mseq.add_sequence(Seq(currentName, currentDotBracket, mseq.size(), encoding_scheme));
                currentDotBracket.clear();
            }
            currentName = line.substr(1);  // Remove '>'
        } else {
            currentDotBracket = line;  // Dot-bracket line
        }
    }

    // Add the last sequence
    if (!currentName.empty() && !currentDotBracket.empty()) {
        mseq.add_sequence(Seq(currentName, currentDotBracket, mseq.size(), encoding_scheme));
    }

    return mseq;
}

int main(int argc, char* argv[]) {
    // Check if the number of arguments is correct (must be an even number for name-sequence pairs)
    if (argc != 5) {
        std::cerr << "Usage: " << argv[0] << " <energy_params> <msa_file_path> <num_iterations> <use_lazy_outside>\n";
        return EXIT_FAILURE;
    }

    std::string msaFilePath = argv[1];
    EnergyParamsType energy_params = std::stoi(argv[2]) == 0 ? EnergyParamsType::VIENNA : EnergyParamsType::BL_STAR;
    int num_itr = std::stoi(argv[3]);
    bool use_lazy_outside = std::stoi(argv[4]);

    try {
        // Read MSA file and generate MultiSeq
        MultiSeq mseq = read_msa_file(msaFilePath, &VIENNA_NUC_ENCODING_SCHEME);

        LinearTurboFold ltf(&mseq, energy_params, num_itr, use_lazy_outside, VerboseState::DEBUG);
        ltf.run();

    } catch (const std::exception& e) {
        std::cerr << "Error: " << e.what() << "\n";
        return EXIT_FAILURE;
    }

    return EXIT_SUCCESS;

    // MultiSeq mseq;
    // mseq.add_sequence(Seq("seq1", "CCCAAAGGG", 0, &nuc_encoding_scheme));
    // mseq.add_sequence(Seq("seq2", "GGGAAUACCC", 0, &nuc_encoding_scheme));

    // mseq.add_sequence(Seq("seq1",
    // "ccgggugccuauaccggaggggccacacccguucccauuccgaacacggucguuaagcccuccagggccgaugguacuggggcguuaccgcccugggagaguaggucggugcccggg",
    // 0, &nuc_encoding_scheme)); mseq.add_sequence(Seq("seq1",
    // "ugcuuggcgaccauagcguuauggacccaccugaucccaugccgaacucaguagugaaacguaauagcgccgaugguaguguggggucuccccaugugagaguaggacaucgccaggcau",
    // 0, &nuc_encoding_scheme)); mseq.add_sequence(Seq("seq2",
    // "ccgggugccuauaccggaggggccacacccguucccauuccgaacacggucguuaagcccuccagggccgaugguacuggggcguuaccgcccugggagaguaggucggugcccggg",
    // 1, &nuc_encoding_scheme)); mseq.add_sequence(Seq("seq3",
    // "ccgggugccuauaccggaggggccacacccguucccauuccgaacacggucguuaagcccuccagggccgaugguacuggggcguuaccgcccugggagaguaggucggugccggg",
    // 2, &nuc_encoding_scheme));

    // mseq.add_sequence(
    //     Seq("seq1",
    //         "acaugcaagucgagcgccccgcaaggggagcggcagacgggugaguaacgcgugggaaucuacccaucucuacggaacaacuccgggaaacuggagcuaauaccgu"
    //         "auacguccuucgggagaaagauuuaucggagauggaugagcccgcguuggauuagcuaguuggugggguaauggccuaccaaggcgacgauccauagcuggucuga"
    //         "gaggaugaucagccacacugggacugagacacggcccagacuccuacgggaggcagcaguggggaauauuggacaaugggcgcaagccugauccagccaugccgcg"
    //         "ugagugaugaaggcccuaggguuguaaagcucuuucaacggugaagauaaugacgguaaccguagaagaagccccggcuaacuucgugccagcagccgcgguaaua"
    //         "cgaagggggcuagcguuguucggaauuacugggcguaaagcgcacguaggcggauauuuaagucaggggugaaaucccggggcucaaccccggaacugccuuugau"
    //         "acuggguaucucgaguccggaagaggugaguggaauuccgaguguagaggugaaauucguagauauucggaggaacaccaguggcgaaggcggcucacugguccgg"
    //         "uacugacgcugaggugcgaaagcguggggagcaaacaggauuagauacccugguaguccacgccguaaacgauggaagcuagccguuggcaaguuuacuugucggu"
    //         "ggcgcagcuaacgcauuaagcuucccgccuggggaguacggucgcaagauuaaaacucaaaggaauugacgggggcccgcacaagcgguggagcaugugguuuaau"
    //         "ucgaagcaacgcgcagaaccuuaccagcccuugacaucccggucgcgguuuccagagauggaaaccuucaguucggcuggaccggugacaggugcugcauggcugu"
    //         "cgucagcucgugucgugagauguuggguuaagucccgcaacgagcgcaacccucgcccuuaguugccagcauucaguugggcacucuaaggggacugccggugaua"
    //         "agccgagaggaagguggggaugacgucaaguccucauggcccuuacgggcugggcuacacacgugcuacaaugguggugacagugggcagcgagaccgcgaggucg"
    //         "agcuaaucuccaaaagccaucucaguucggauugcacucugcaacucgagugcaugaaguuggaaucgcuaguaaucgcggaucagcaugccgcggugaauacguu"
    //         "cccgggccuuguacacaccgcccgucacaccaugggaguugguuuuacccgaaggcgcugugcuaaccgcaa",
    //         0, &nuc_encoding_scheme));
    // mseq.add_sequence(Seq(
    //     "seq2",
    //     "gcgaacgcuggcggcaggcuuaacacaugcaagucgagcgggcauagcaauaugucagcggcagacgggugaguaacgcgugggaacguaccuuuugguucggaacaaca"
    //     "cagggaaacuugugcuaauaccggauaagcccuuacggggaaagauuuaucgccgaaagaucggcccgcgucugauuagcuaguuggugagguaauggcucaccaaggcg"
    //     "acgaucaguagcuggucugagaggaugaucagccacauugggacugagacacggcccaaacuccuacgggaggcagcaguggggaauauuggacaaugggcgcaagccug"
    //     "auccagccaugccgcgugagugaugaaggcccuaggguuguaaagcucuuuugugcgggaagauaaugacgguaccgcaagaauaagccccggcuaacuucgugccagca"
    //     "gccgcgguaauacgaagggggcuagcguugcucggaaucacugggcguaaagggugcguaggcgggucuuuaagucaggggugaaauccuggagcucaacuccagaacug"
    //     "ccuuugauacugaagaucuugaguucgggagaggugaguggaacugcgaguguagaggugaaauucguagauauucgcaagaacaccaguggcgaaggcggcucacuggc"
    //     "ccgauacugacgcugaggcacgaaagcguggggagcaaacaggauuagauacccugguaguccacgccguaaacgaugaaugccagccguuaguggguuuacucacuagu"
    //     "ggcgcagcuaacgcuuuaagcauuccgccuggggaguacggucgcaagauuaaaacucaaaggaauugacgggggcccgcacaagcgguggagcaugugguuuaauucga"
    //     "cgcaacgcgcagaaccuuaccagcccuugacaucccggucgcggacuccagagacggaguucuucaguucggcuggaccggagacaggugcugcauggcugucgucagcu"
    //     "cgugucgugagauguuggguuaagucccgcaacgagcgcaacccccguccuuaguugcuaccauuuaguugagcacucuaaggagacugccggugauaagccgcgaggaa"
    //     "gguggggaugacgucaaguccucauggcccuuacgggcugggcuacacacgugcuacaauggcggugacaaugggaugcuaaggggcgacccuucgcaaaucucaaaaag"
    //     "ccgucucaguucggauugggcucugcaacucgagcccaugaaguuggaaucgcuaguaaucguggaucagcacgccacggugaauacguucccgggccuuguacacaccg"
    //     "cccgucacaccaugggaguugguuuuaccugaagacggugcgcuaacccgcaagggaggcagccggccacgguagggucagcgacuggggug",
    //     0, &nuc_encoding_scheme));

    // LinearTurboFold ltf(&mseq, 1, VerboseState::DEBUG);
    // ltf.run();

    // return 0;
}