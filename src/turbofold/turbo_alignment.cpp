#include <linearx/turbofold/linear_turbofold.hpp>

using namespace linearx::utils;
using namespace linearx::constants::math;
using namespace std;

TurboAlignment::TurboAlignment(const LinearTurbofold &turbofold, Sequence &seq1, Sequence &seq2, const float alpha1,
                               const float alpha2, const float alpha3)
    : LinearAlignment(seq1, seq2, alpha1, alpha2, alpha3), turbofold(turbofold) {}

void TurboAlignment::reset_saved_beams() {
    bestALN.clear();
    bestINS1.clear();
    bestINS2.clear();
    bestALN.reserve(seq_len_sum + 3);
    bestINS1.reserve(seq_len_sum + 1);
    bestINS2.reserve(seq_len_sum + 1);
}

void TurboAlignment::save_partition_function(const bool move) {
    // if move is true, the original beam data will be cleared
    reset_saved_beams();
    auto save = [this, move](auto &src, auto &dest) {
        if (move) {
            dest = std::move(src);
        } else {
            dest = src;
        }
    };
    save(bestALN, saved_bestALN);
    save(bestINS1, saved_bestINS1);
    save(bestINS2, saved_bestINS2);
}
