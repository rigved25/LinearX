// src/turbofold/main.cpp
#include <iostream>
#include <linearx/constants.hpp>
#include <linearx/energy/energy_model.hpp>
#include <linearx/sequence/multi_sequence.hpp>
#include <linearx/turbofold/linear_turbofold.hpp>

int main(int argc, char* argv[]) {
    if (argc != 15) {
        std::cerr << "Usage: " << argv[0]
                  << " <msa_file_path> <out_dir> <energy_model> <num_iterations> "
                     "<use_lazy_outside> <restrict_search> <astar_viterbi> <max_marginal> "
                     "<verbose> <save_logs> <save_probs> "
                     "<alignment_threshold> <folding_threshold> <max_marginal_pruning_threshold>\n";
        return EXIT_FAILURE;
    }

    std::string msa_file_path = argv[1];
    std::string out_dir = argv[2];
    int energy_model_choice = std::stoi(argv[3]);
    int num_itr = std::stoi(argv[4]);
    bool use_lazy_outside = std::stoi(argv[5]) != 0;
    bool restrict_search = std::stoi(argv[6]) != 0;
    bool astar_viterbi = std::stoi(argv[7]) != 0;
    bool max_marginal = std::stoi(argv[8]) != 0;
    bool verbose = std::stoi(argv[9]) != 0;
    bool save_logs = std::stoi(argv[10]) != 0;
    bool save_probs = std::stoi(argv[11]) != 0;
    float alignment_threshold = std::stof(argv[12]);
    float folding_threshold = std::stof(argv[13]);
    float max_marginal_pruning_threshold = std::stof(argv[14]);

    std::cerr << "Arguments:\n";
    std::cerr << "  msa_file_path: " << msa_file_path << "\n";
    std::cerr << "  out_dir: " << (out_dir.empty() ? "[not provided]" : out_dir) << "\n";
    std::cerr << "  energy_model: " << energy_model_choice << "\n";
    std::cerr << "  num_iterations: " << num_itr << "\n";
    std::cerr << "  use_lazy_outside: " << use_lazy_outside << "\n";
    std::cerr << "  restrict_search: " << restrict_search << "\n";
    std::cerr << "  astar_viterbi: " << astar_viterbi << "\n";
    std::cerr << "  max_marginal: " << max_marginal << "\n";
    std::cerr << "  verbose: " << verbose << "\n";
    std::cerr << "  save_logs: " << save_logs << "\n";
    std::cerr << "  save_probs: " << save_probs << "\n";
    std::cerr << "  alignment_threshold: " << alignment_threshold << "\n";
    std::cerr << "  folding_threshold: " << folding_threshold << "\n";
    std::cerr << "  max_marginal_pruning_threshold: " << max_marginal_pruning_threshold << "\n";

    try {
        MultiSeq ms;
        ms.read_fasta(msa_file_path, linearx::utils::VIENNA_NUC_ENCODING_SCHEME);
        EnergyModel energy_model(energy_model_choice == 0 ? EnergyParamsType::VIENNA : EnergyParamsType::BL_STAR);
        LinearTurbofold ltf(ms, energy_model, 100, 100,
                           alignment_threshold,
                           folding_threshold,
                           max_marginal_pruning_threshold,
                           0.3, 0.3, 3, false, 1.0, 0.8, 0.5,
                           num_itr, verbose, save_logs, save_probs, out_dir);
        ltf.run(use_lazy_outside, restrict_search, astar_viterbi, max_marginal);
    } catch (const std::exception& e) {
        std::cerr << "Error: " << e.what() << "\n";
        return EXIT_FAILURE;
    }
    return EXIT_SUCCESS;
}
