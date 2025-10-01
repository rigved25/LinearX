#include <linearx/turbofold/linear_turbofold.hpp>

using namespace linearx::utils;
using namespace linearx::constants::math;
using namespace std;

TurboAlignment::TurboAlignment(const LinearTurbofold& turbofold, Sequence& seq1, Sequence& seq2, const float alpha1,
                               const float alpha2, const float alpha3)
    : LinearAlignmentInterface<TurboAlignment>(seq1, seq2, alpha1, alpha2, alpha3), turbofold(turbofold) {
    reset_beams(turbofold.alignment_beam_size);
}

// void TurboAlignment::reset_saved_beams(const unsigned beam_size) {
//     reset_beam_vector(saved_bestALN, seq_len_sum + 3, beam_size);
//     reset_beam_vector(saved_bestINS1, seq_len_sum + 1, beam_size);
//     reset_beam_vector(saved_bestINS2, seq_len_sum + 1, beam_size);
// }

void TurboAlignment::save_partition_function(const bool move) {
    // if move is true, the original beam data will be cleared
    // reset_saved_beams(run_beam_size_);
    auto save = [this, move](auto& src, auto& dest) {
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
