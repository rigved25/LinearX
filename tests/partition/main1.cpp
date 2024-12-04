#include <fstream>
#include <iostream>
#include <stdexcept>
#include <string>

#include "./../../partition/partition.hpp"
#include "./../../sequence/multi_seq.hpp"

using namespace std;

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
    if (argc < 2 || argc > 4) {
        std::cerr << "Usage: " << argv[0] << " <msa_file_path> [energy_params=0|1] [use_lazy_outside=1|0]\n";
        std::cerr << "  [energy_params=0|1]: 0 for VIENNA, 1 for BL_STAR (default: 0)\n";
        std::cerr << "  [use_lazy_outside=1|0]: 1 for on, 0 for off (default: 1)\n";
        return EXIT_FAILURE;
    }

    // Parse input arguments
    std::string msaFilePath = argv[1];
    EnergyParamsType energy_params =
        (argc >= 3 && std::stoi(argv[2]) == 1) ? EnergyParamsType::BL_STAR : EnergyParamsType::VIENNA;
    bool use_lazy_outside = (argc == 4) ? std::stoi(argv[3]) : true;

    try {
        // Read MSA file and generate MultiSeq
        MultiSeq mseq = read_msa_file(msaFilePath);

        // Display the sequences
        std::cerr << "MultiSeq contains " << mseq.size() << " sequences:\n";
        for (int k = 0; k < mseq.size(); ++k) {
            std::cerr << ">" << mseq[k].id << std::endl;
            std::cerr << mseq[k].sequence << std::endl;
        }

        for (int k = 0; k < mseq.size(); ++k) {
            std::cerr << "\n---------- Sequence " << mseq[k].id << " ----------\n";
            Seq& seq = mseq[k];

            seq.set_encoding(VIENNA_NUC_ENCODING_SCHEME);
            InsideMode mode = InsideMode::PARTITION;
            Partition partition(&seq, *(new EnergyModel(energy_params)), mode);
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
