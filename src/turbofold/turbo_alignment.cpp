#include <linearx/turbofold/linear_turbofold.hpp>
#include <linearx/constants.hpp>

using namespace linearx::utils;
using namespace linearx::constants::math;
using namespace linearx::constants::alignment;
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
    if (seq_len_sum > linearx::constants::alignment::SAVE_SWAP_THRESHOLD) {
        // Long sequences: O(1) swap to avoid large copy/dealloc cost. Caller must call reset_beams() after.
        std::swap(bestALN, saved_bestALN);
        std::swap(bestINS1, saved_bestINS1);
        std::swap(bestINS2, saved_bestINS2);
    } else {
        // Medium/short: copy or move to preserve cache locality for next pair.
        if (move) {
            saved_bestALN = std::move(bestALN);
            saved_bestINS1 = std::move(bestINS1);
            saved_bestINS2 = std::move(bestINS2);
        } else {
            saved_bestALN = bestALN;
            saved_bestINS1 = bestINS1;
            saved_bestINS2 = bestINS2;
        }
    }
}

void TurboAlignment::save_partition_function_old(const bool move) {
    // if move is true, the original beam data will be cleared
    // reset_saved_beams(run_beam_size_);
    auto save = [this, move](auto& src, auto& dest) {
        if (move) {
            dest = std::move(src);  // O(1) move + destructor
        } else {
            dest = src;  // deep copy - O(data size)
        }
    };
    save(bestALN, saved_bestALN);
    save(bestINS1, saved_bestINS1);
    save(bestINS2, saved_bestINS2);
}
