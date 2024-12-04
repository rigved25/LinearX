#include <iostream>
#include <string>

#include "./../../partition/partition.hpp"
#include "./../../sequence/seq.hpp"

using namespace std;

// example run: awk "NR==2" ./test_seqs/4.fasta | ./main2 0 1 3>&1

int main(int argc, char* argv[]) {
    // Default options
    EnergyParamsType energy_params = EnergyParamsType::VIENNA;
    bool use_lazy_outside = true;

    // Parse optional arguments
    if (argc > 3) {
        std::cerr << "Usage: echo <sequence> | " << argv[0] << " [0|1] [0|1]\n";
        std::cerr << "  [0|1]: Energy Parameters (0 = VIENNA, 1 = BL_STAR)\n";
        std::cerr << "  [0|1]: Use Lazy Outside (0 = off, 1 = on)\n";
        return EXIT_FAILURE;
    }

    if (argc >= 2) {
        int energy_param_arg = std::stoi(argv[1]);
        if (energy_param_arg == 0) {
            energy_params = EnergyParamsType::VIENNA;
        } else if (energy_param_arg == 1) {
            energy_params = EnergyParamsType::BL_STAR;
        } else {
            std::cerr << "Error: Invalid value for 'Energy Parameters'. Use 0 (VIENNA) or 1 (BL_STAR).\n";
            return EXIT_FAILURE;
        }
    }

    if (argc == 3) {
        use_lazy_outside = std::stoi(argv[2]);
        if (use_lazy_outside != 0 && use_lazy_outside != 1) {
            std::cerr << "Error: Invalid value for 'Use Lazy Outside'. Use 0 or 1.\n";
            return EXIT_FAILURE;
        }
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
    std::cout << "Energy Parameters: " << (energy_params == EnergyParamsType::VIENNA ? "VIENNA" : "BL_STAR")
              << std::endl;
    std::cout << "Use Lazy Outside: " << use_lazy_outside << std::endl;

    try {
        // Create a Seq object from the input sequence
        Seq seq("0", sequence);

        seq.set_encoding(VIENNA_NUC_ENCODING_SCHEME);
        InsideMode mode = InsideMode::PARTITION;

        // Create Partition object and compute values
        Partition partition(&seq, *(new EnergyModel(energy_params)), mode);
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
