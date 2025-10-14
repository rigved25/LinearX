// src/turbofold/turbo_partition.cpp
#include <chrono>
#include <linearx/partition/linear_partition.hpp>
#include <linearx/turbofold/linear_turbofold.hpp>

using namespace linearx::math;
using namespace linearx::utils;

TurboPartition::TurboPartition(LinearTurbofold &turbofold, const Sequence &seq, const EnergyModel &energy_model,
                               const bool allow_sharp_turn)
    : LinearPartitionInterface<TurboPartition>(seq, energy_model, allow_sharp_turn), turbofold(turbofold) {
    prob_accm.upstrm.resize(seq.length());
    prob_accm.dwnstrm.resize(seq.length());
    reset_beams(turbofold.folding_beam_size);
}

// void TurboPartition::reset_saved_beams(const unsigned beam_size) {
//     reset_beam_vector(saved_bestH, seq_length, beam_size);
//     reset_beam_vector(saved_bestP, seq_length, beam_size);
//     reset_beam_vector(saved_bestM, seq_length, beam_size);
//     reset_beam_vector(saved_bestM2, seq_length, beam_size);
//     reset_beam_vector(saved_bestMulti, seq_length, beam_size);
// }

void TurboPartition::save_partition_function(const bool move, const unsigned beam_size) {
    // if move is true, the original beam data will be cleared
    // reset_saved_beams(beam_size);
    auto save = [this, move](auto &src, auto &dest) {
        if (move) {
            dest = std::move(src);
        } else {
            dest = src;
        }
    };
    save(bestH, saved_bestH);
    save(bestP, saved_bestP);
    save(bestM, saved_bestM);
    save(bestM2, saved_bestM2);
    save(bestMulti, saved_bestMulti);
    total_inside = bestC[seq_length - 1].alpha;
}

void TurboPartition::calc_prob_accm() {
    fill(prob_accm.upstrm.begin(), prob_accm.upstrm.end(), 0.0);
    fill(prob_accm.dwnstrm.begin(), prob_accm.dwnstrm.end(), 0.0);
    for (int j = 0; j < seq.length(); ++j) {
        for (const auto &item : this->bpp[j]) {
            const int i = item.first;
            const value_type prob = item.second;

            prob_accm.upstrm[j] += prob;
            prob_accm.dwnstrm[i] += prob;
        }
    }
    for (int j = 0; j < seq.length(); ++j) {
        if (prob_accm.upstrm[j] > 1.0) {
            prob_accm.upstrm[j] = 1.0;
            prob_accm.dwnstrm[j] = 0.0;
        }
        if (prob_accm.dwnstrm[j] > 1.0) {
            prob_accm.upstrm[j] = 0.0;
            prob_accm.dwnstrm[j] = 1.0;
        }
    }
}
