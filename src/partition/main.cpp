// src/partition/main.cpp
#include <fstream>
#include <iostream>
#include <linearx/constants.hpp>
#include <linearx/energy/energy_model.hpp>
#include <linearx/partition/linear_partition.hpp>
#include <linearx/sequence/sequence.hpp>
#include <sstream>

bool file_exists(const std::string& path) {
    std::ifstream f(path);
    return f.good();
}

int main(int argc, char* argv[]) {
    if (argc != 9) {
        std::cerr << "Usage: " << argv[0]
                  << " <sequence_or_fasta> <energy_model> <use_lazy_outside> <mfe_mode> <verbose> "
                     "<use_threshknot> <compute_bpp> <bpp_path>\n";
        std::cerr << "  energy_model: 0 = Vienna, 1 = BL*\n";
        std::cerr << "  use_lazy_outside: 0 = false, 1 = true\n";
        std::cerr << "  mfe_mode: 0 = partition function, 1 = MFE (best mode)\n";
        std::cerr << "  verbose: 0 = false, 1 = true\n";
        std::cerr << "  use_threshknot: 0 = false, 1 = true\n";
        std::cerr << "  compute_bpp: 0 = false, 1 = true\n";
        std::cerr << "  bpp_path: directory path (ignored if compute_bpp=0)\n";
        return EXIT_FAILURE;
    }

    std::string seq_or_file = argv[1];
    int energy_model_choice = std::stoi(argv[2]);
    bool use_lazy_outside = std::stoi(argv[3]) != 0;
    bool use_mfe = std::stoi(argv[4]) != 0;
    bool verbose = std::stoi(argv[5]) != 0;
    bool use_threshknot = std::stoi(argv[6]) != 0;
    bool compute_bpp = std::stoi(argv[7]) != 0;
    std::string bpp_path = argv[8];

    if (bpp_path == "none") {
        bpp_path.clear();
    }

    std::cerr << "Arguments:\n";
    std::cerr << "  input: " << seq_or_file << "\n";
    std::cerr << "  energy_model: " << energy_model_choice << "\n";
    std::cerr << "  use_lazy_outside: " << use_lazy_outside << "\n";
    std::cerr << "  mfe_mode: " << use_mfe << "\n";
    std::cerr << "  verbose: " << verbose << "\n";
    std::cerr << "  threshknot: " << use_threshknot << "\n";
    std::cerr << "  bpp: " << compute_bpp << "\n";
    if (compute_bpp) {
        std::cerr << "  bpp_path: " << (bpp_path.empty() ? "[not provided]" : bpp_path) << "\n";
    }

    try {
        Sequence seq;
        if (file_exists(seq_or_file)) {
            seq.read_fasta(seq_or_file);
        } else if (!seq_or_file.empty() && seq_or_file[0] == '>') {
            std::istringstream fasta_stream(seq_or_file);
            seq.read_fasta_stream(fasta_stream);
        } else {
            seq.set_seq(seq_or_file);
        }
        seq.set_encoding(linearx::utils::VIENNA_NUC_ENCODING_SCHEME);

        EnergyModel energy_model(energy_model_choice == 0 ? EnergyParamsType::VIENNA : EnergyParamsType::BL_STAR);

        LinearPartition lp(seq, energy_model);
        lp.reset_beams(100);

        if (use_mfe) {
            lp.compute_inside(Mode::BEST, 100, verbose);
            lp.compute_outside(use_lazy_outside, linearx::constants::limits::DEVIATION_THRESHOLD, verbose);
            auto structure = lp.get_mfe_structure();
            std::cout << "MFE Structure: " << structure.getDotBracket() << "\n";
        } else {
            lp.compute_inside(Mode::PARTITION_INSIDE, 100, verbose);
            lp.compute_outside(use_lazy_outside, linearx::constants::limits::DEVIATION_THRESHOLD, verbose);
            if (use_threshknot || compute_bpp) {
                lp.compute_bpp_matrix(100, verbose);
            }
            fprintf(stderr, "Free Energy of the Ensemble: %.5f kcal/mol\n", lp.get_ensemble_energy());
            if (use_threshknot) {
                Structure tk_struct = lp.get_threshknot_structure();
                std::cout << "ThreshKnot Structure: " << tk_struct.getDotBracket() << "\n";
            }
            if (compute_bpp) {
                if (!bpp_path.empty()) {
                    lp.dump_bpp(bpp_path);
                    std::cout << "BPP matrix written to " << bpp_path << "\n";
                } else {
                    std::cout << "BPP matrix computed (not saved, no bpp_path provided)\n";
                }
            }
        }

    } catch (const std::exception& e) {
        std::cerr << "Error: " << e.what() << "\n";
        return EXIT_FAILURE;
    }
    return EXIT_SUCCESS;
}
