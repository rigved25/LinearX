// src/partition/linear_partition.cpp
#include <fstream>
#include <iomanip>
#include <linearx/partition/config.hpp>
#include <linearx/partition/linear_partition.hpp>
#include <set>

using namespace linearx::math;
using namespace std;
using namespace linearx::utils;

template <typename T>
LinearPartitionInterface<T>::LinearPartitionInterface(const Sequence& seq, const EnergyModel& energy_model,
                                                      bool allow_sharp_turn)
    : seq(seq),
      seq_length(seq.length()),
      energy_model(energy_model),
      allow_sharp_turn(allow_sharp_turn),
      bestC(seq.length()) {
    for (int nuc = 0; nuc < 5; ++nuc) {
        prev_pair[nuc].resize(seq_length, -1);
        next_pair[nuc].resize(seq_length, seq_length);
        if (nuc == 0) {
            continue;  // skip N
        }
        int prev = -1;
        int next = seq_length;
        for (int j = 0; j < seq_length; ++j) {
            prev_pair[nuc][j] = prev;
            if (check_valid_pair(nuc, seq.enc[j])) {
                prev = j;
            }
        }
        for (int j = seq_length - 1; j >= 0; --j) {
            next_pair[nuc][j] = next;
            if (check_valid_pair(nuc, seq.enc[j])) {
                next = j;
            }
        }
    }
    energy_model.init_tetra_hex_tri(seq.seq, seq_length, if_tetraloops, if_hexaloops, if_triloops);
}

template <typename T>
State& LinearPartitionInterface<T>::get_viterbi() {
    return bestC[seq_length - 1];
}

template <typename T>
void LinearPartitionInterface<T>::reset_beams(const unsigned beam_size) {
    reset_beam_vector(bestH, seq_length, beam_size);
    reset_beam_vector(bestP, seq_length, beam_size);
    reset_beam_vector(bestM, seq_length, beam_size);
    reset_beam_vector(bestM2, seq_length, beam_size);
    reset_beam_vector(bestMulti, seq_length, beam_size);
    bestC.reset();

    bestC[-1].alpha = LOG_ONE;
    bestC[seq_length - 1].beta = LOG_ONE;
    if (seq_length > 0) {
        bestC[0].alpha = LOG_ONE;
        if (seq_length > 1) {
            bestC[1].alpha = LOG_ONE;
        }
    }
}

template <typename T>
void LinearPartitionInterface<T>::compute_bpp_matrix(const unsigned beam_size, const bool verbose_output) {
    auto start_time = chrono::high_resolution_clock::now();
    // clear the existing bpp matrix
    reset_beam_vector(bpp, seq_length, beam_size);

    State& viterbi = get_viterbi();
    for (int j = 0; j < seq_length; ++j) {
        for (const auto& item : bestP[j]) {
            const int i = item.first;
            const State& state = item.second;

            value_type prob = LOG_DIV(LOG_MUL(state.alpha, state.beta), viterbi.alpha);
            if (prob > -linearx::constants::limits::DEVIATION_THRESHOLD) {
                prob = xexp(prob);  // Convert log prob to regular prob
                if (prob > 1.001) {
                    fprintf(stderr,
                            "[LinearPartition Warning] BPP value too high, something is wrong! bpp(%d, %d): %.5f\n", i,
                            j, prob);
                }
                this->bpp[j][i] = min(prob, (value_type)1.0);  // Clamp the probability to [0, 1]
            }
        }
    }
    auto end_time = chrono::high_resolution_clock::now();
    const value_type execution_time = chrono::duration_cast<chrono::milliseconds>(end_time - start_time).count();
    if (verbose_output) {
        fprintf(stderr, "[LinearPartition] BPP Matrix Computation Time: %.2f ms\n", execution_time);
    }
    log.bpp_exec_time = execution_time;
}

template <typename T>
value_type LinearPartitionInterface<T>::get_ensemble_energy() const {
    return -linearx::constants::energy::kT * (bestC[seq_length - 1].alpha) / 100.0;
}  // -kT log(Q(x))

template <typename T>
void LinearPartitionInterface<T>::print_alpha_beta() const {
    const int j = seq.length() - 1;
    printf("C[%d]: %.5f\t%.5f\n", j, bestC[j].alpha, bestC[j].beta);
}

template <typename T>
Structure LinearPartitionInterface<T>::get_mfe_structure() {
    Structure structure(seq_length);
    mfe_backtrack(0, seq_length - 1, C, structure);
    return structure;
}

template <typename T>
Structure LinearPartitionInterface<T>::get_threshknot_structure(float threshknot_threshold, int min_helix_size) const {
    vector<value_type> best_prob(seq_length, 0.0);
    Structure structure(seq_length);
    set<int> visited;

    for (int j = 0; j < seq_length; j++) {
        for (const auto& item : bpp[j]) {
            const int i = item.first;
            if (j - i < min_helix_size + 1) {
                continue;
            }
            const value_type prob = item.second;
            best_prob[i] = max(best_prob[i], prob);
            best_prob[j] = max(best_prob[j], prob);
        }
    }

    for (int j = 0; j < seq_length; j++) {
        for (const auto& item : bpp[j]) {
            const int i = item.first;
            if (j - i < min_helix_size + 1) {
                continue;
            }

            const value_type prob = item.second;
            if (prob >= threshknot_threshold && prob == best_prob[i] && prob == best_prob[j]) {
                if (visited.find(i) != visited.end() || visited.find(j) != visited.end()) {
                    continue;
                }
                structure.addPair(i, j);
                visited.insert(i);
                visited.insert(j);
            }
        }
    }

    structure.removeShortHelices(min_helix_size);
    return structure;
}

template <typename T>
void LinearPartitionInterface<T>::dump_bpp(const std::string& out_dir) const {
    if (bpp.empty()) {
        throw std::runtime_error(
            "[LinearPartition: Error] BPP matrix is empty! You must run compute_bpp_matrix() first.");
    }

    // create directory if it doesn't exist
    std::filesystem::create_directories(out_dir);

    // construct output filename
    std::string filename = out_dir + "/bpp_" + std::to_string(seq.id) + ".txt";

    // open the file for writing
    std::ofstream file(filename);
    if (!file) {
        throw std::runtime_error("[LinearPartition: Error] Unable to open the file " + filename +
                                 " for writing BPP matrix.");
    }

    // write bpp matrix
    for (int j = 0; j < seq_length; ++j) {
        for (const auto& item : bpp[j]) {
            const int i = item.first;
            const value_type prob = item.second;
            file << i + 1 << " " << j + 1 << " " << std::fixed << std::setprecision(4) << prob << "\n";
        }
    }
}

// Instantiate templates for LinearPartitionInterface with desired types
#define X(TYPE) template class LinearPartitionInterface<TYPE>;
LP_TEMPLATE_TYPES
#undef X
