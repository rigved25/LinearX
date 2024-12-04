#include <iostream>
#include <string>

#include "./../../partition/partition.hpp"
#include "./../../sequence/seq.hpp"

using namespace std;

// example run: awk "NR==2" ./seqs/1.txt | ./main2 0

int main(int argc, char* argv[]) {
    // Default lazy outside option
    bool use_lazy_outside = true;

    // Parse optional argument for lazy outside option
    if (argc > 2) {
        std::cerr << "Usage: echo <sequence> | " << argv[0] << " [0|1]\n";
        return EXIT_FAILURE;
    }

    if (argc == 2) {
        use_lazy_outside = std::stoi(argv[1]);
    }

    // Read sequence from standard input
    std::string sequence;
    if (!std::getline(std::cin, sequence)) {
        std::cerr << "Error: Failed to read sequence from standard input.\n";
        return EXIT_FAILURE;
    }

    if (sequence.empty()) {
        std::cerr << "Error: Input sequence is empty.\n";
        return EXIT_FAILURE;
    }

    std::cerr << "Sequence: " << sequence << std::endl;
    std::cout << "Use Lazy Outside " << use_lazy_outside << std::endl;

    try {
        // Create a Seq object from the input sequence
        Seq seq("0", sequence);

        seq.set_encoding(VIENNA_NUC_ENCODING_SCHEME);
        InsideMode mode = InsideMode::PARTITION;

        // Create Partition object and compute values
        Partition partition(&seq, *(new EnergyModel(EnergyParamsType::VIENNA)), mode);
        partition.compute_inside();
        if (mode == InsideMode::PARTITION) {
            partition.compute_outside(use_lazy_outside);
            partition.print_alpha_beta();
            partition.compute_bpp_matrix();

            // Print results
            std::cout << ">" << partition.sequence->id << std::endl;
            std::cout << partition.get_threshknot_structure() << std::endl;
        }

    } catch (const std::exception& e) {
        std::cerr << "Error: " << e.what() << "\n";
        return EXIT_FAILURE;
    }

    return EXIT_SUCCESS;
}
