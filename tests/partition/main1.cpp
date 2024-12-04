#include <iostream>

#include "./../../partition/partition.hpp"
#include "./../../sequence/multi_seq.hpp"

using namespace std;

// example run: // ./main ./../../eval/rnastralign/data/v1/no_aln/5S.RNAstralign_LZ_format_k_30_1.auto_mafft.fasta 1

// Function to read MSA file and populate MultiSeq
MultiSeq read_msa_file(const std::string& filePath) {
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
                mseq.add_sequence(Seq(currentName, currentDotBracket));
                currentDotBracket.clear();
            }
            currentName = line.substr(1);  // Remove '>'
        } else {
            currentDotBracket = line;  // Dot-bracket line
        }
    }

    // Add the last sequence
    if (!currentName.empty() && !currentDotBracket.empty()) {
        mseq.add_sequence(Seq(currentName, currentDotBracket));
    }

    return mseq;
}

int main(int argc, char* argv[]) {
    if (argc != 3) {
        std::cerr << "Usage: " << argv[0] << " <msa_file_path> <use_lazy_outside>\n";
        return EXIT_FAILURE;
    }

    std::string msaFilePath = argv[1];
    bool use_lazy_outside = std::stoi(argv[2]);

    try {
        // Read MSA file and generate MultiSeq
        MultiSeq mseq = read_msa_file(msaFilePath);

        // Display the sequences
        std::cerr << "MultiSeq contains " << mseq.size() << " sequences:\n";

        // Print the sequences
        for (int k = 0; k < mseq.size(); ++k) {
            std::cerr << ">" << mseq[k].id << std::endl;
            std::cerr << mseq[k].sequence << std::endl;
        }

        for (int k = 0; k < mseq.size(); ++k) {
            std::cerr << "\n---------- Sequence " << mseq[k].id << " ----------\n";
            Seq& seq = mseq[k];

            seq.set_encoding(VIENNA_NUC_ENCODING_SCHEME);
            InsideMode mode = InsideMode::PARTITION;
            Partition partition(&seq, *(new EnergyModel(EnergyParamsType::VIENNA)), mode);
            partition.compute_inside();
            if (mode == InsideMode::PARTITION) {
                partition.compute_outside(use_lazy_outside);
                partition.print_alpha_beta();
                partition.compute_bpp_matrix();
                std::cout << ">" << partition.sequence->id << std::endl;
                std::cout << partition.get_threshknot_structure() << std::endl;
            }
        }

    } catch (const std::exception& e) {
        std::cerr << "Error: " << e.what() << "\n";
        return EXIT_FAILURE;
    }

    return EXIT_SUCCESS;
}
