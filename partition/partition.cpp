#include "partition.hpp"

PartitionFunctionBeam::~PartitionFunctionBeam() { free(); }

void PartitionFunctionBeam::free() {
    delete[] bestH;
    delete[] bestP;
    delete[] bestM;
    delete[] bestM2;
    delete[] bestMulti;

    bestH = nullptr;
    bestP = nullptr;
    bestM = nullptr;
    bestM2 = nullptr;
    bestMulti = nullptr;
    total_alpha = 0.0;
}

void PartitionFunctionBeam::save(std::unordered_map<int, State> *bestH, std::unordered_map<int, State> *bestP,
                                 std::unordered_map<int, State> *bestM, std::unordered_map<int, State> *bestM2,
                                 std::unordered_map<int, State> *bestMulti, double total_alpha) {
    this->bestH = bestH;
    this->bestP = bestP;
    this->bestM = bestM;
    this->bestM2 = bestM2;
    this->bestMulti = bestMulti;
    this->total_alpha = total_alpha;
}

void PartitionFunctionBeam::save(Partition &pf) {
    save(pf.bestH, pf.bestP, pf.bestM, pf.bestM2, pf.bestMulti, pf.bestC[pf.seq->size() - 1].alpha);
}

void Partition::reset_beams(bool freeMemory) {
    if (freeMemory) {
        delete[] bestH;
        delete[] bestP;
        delete[] bestM;
        delete[] bestM2;
        delete[] bestMulti;
    }

    bestH = new std::unordered_map<int, State>[seq->size()];
    bestP = new std::unordered_map<int, State>[seq->size()];
    bestM = new std::unordered_map<int, State>[seq->size()];
    bestM2 = new std::unordered_map<int, State>[seq->size()];
    bestMulti = new std::unordered_map<int, State>[seq->size()];

    // bestC is not a pointer but unordered_map, [TODO] need to fix this!
    bestC.reset();

    bestC[-1].alpha = 0;
    bestC[seq->size() - 1].beta = 0;
    if (seq->size() > 0) bestC[0].alpha = 0.0;
    if (seq->size() > 1) bestC[1].alpha = 0.0;
}

void Partition::compute_bpp_matrix() {
    // clear the existing bpp matrix
    delete[] this->bpp;

    // compute the new bpp matrix
    this->bpp = new std::unordered_map<int, double>[seq->size()];  // reallocate memory
    for (int j = 0; j < seq->size(); ++j) {
        for (const auto &item : bestP[j]) {
            const int i = item.first;
            const State &state = item.second;

            double prob = xlog_div(xlog_mul(state.alpha, state.beta), viterbi->alpha);
            if (prob > -DEVIATION_THRESHOLD) {
                prob = xexp(prob);  // Convert log prob to regular prob
                if (prob > 1.001) {
                    fprintf(stderr,
                            "[LinearPartition Warning] BPP value too high, something is wrong! bpp(%d, %d): %.5f\n", i,
                            j, prob);
                }
                prob = std::min(prob, 1.0);  // Clamp the probability to [0, 1]
                this->bpp[j][i] = prob;      // Set the bpp value
            }
        }
    }
}

double Partition::get_bpp(int i, int j) const {
    if (bpp == nullptr) {
        throw std::runtime_error(
            "[Error (LinearPartition)] BPP matrix not computed yet! You must run compute_bpp_matrix() first.");
    }
    // if (i > j) {
    //     std::swap(i, j);
    // }
    const auto &item = bpp[j].find(i);
    if (item == bpp[j].end()) {
        return 0.0;
    }
    return item->second;
}

double Partition::get_ensemble_energy() { return -kT * (bestC[seq->size() - 1].alpha) / 100.0; }  // -kT log(Q(x))

std::string Partition::get_mfe_structure() {
    std::string structure(seq->size(), '.');
    mfe_backtrack(0, seq->size() - 1, C, structure);
    return structure;
}

void Partition::print_alpha_beta() {
    std::cerr << "\nAlpha(C, n - 1): " << bestC[seq->size() - 1].alpha << std::endl;
    std::cerr << "Beta(C, -1): " << bestC[-1].beta << std::endl << std::endl;
}

std::string &Partition::get_threshknot_structure(float threshknot_threshold, int min_helix_size) const {
    std::vector<double> best_prob(seq->size(), 0.0);
    Structure structure(seq->size());
    std::set<int> visited;

    for (int j = 0; j < seq->size(); j++) {
        for (const auto &item : bpp[j]) {
            const int i = item.first;
            if (j - i < min_helix_size + 1) {
                continue;
            }
            const double prob = item.second;
            best_prob[i] = std::max(best_prob[i], prob);
            best_prob[j] = std::max(best_prob[j], prob);
        }
    }

    for (int j = 0; j < seq->size(); j++) {
        for (const auto &item : bpp[j]) {
            const int i = item.first;
            if (j - i < min_helix_size + 1) {
                continue;
            }

            const double prob = item.second;
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
    std::string &dotBracket = structure.getDotBracket();

    return dotBracket;
}

void Partition::dump_bpp(const std::string &filepath) const {
    if (bpp == nullptr) {
        throw std::runtime_error(
            "[LinearPartition: Error] BPP matrix not computed yet! You must run compute_bpp_matrix() first.");
    }

    // open the file for writing
    std::ofstream file(filepath);
    if (!file) {
        std::cout
            << "[Hint] The directory for the output file may not exist. Please create it before running the method."
            << std::endl;
        throw std::runtime_error("[LinearPartition: Error] Unable to open the file " + filepath +
                                 " for writing BPP matrix.");
    }

    // dump the BPP matrix to the file
    for (int j = 0; j < seq->size(); ++j) {
        for (const auto &item : bpp[j]) {
            const int i = item.first;
            const double prob = item.second;

            // output i, j, and the probability to the file
            file << i << " " << j << " " << std::fixed << std::setprecision(8) << prob << std::endl;
        }
    }

    // file automatically closes when it goes out of scope
    std::cout << "[LinearPartition] BPP matrix dumped to " << filepath << std::endl;
}
