// src/turbofold/turbo_partition.cpp
#include <linearx/turbofold/linear_turbofold.hpp>

using namespace linearx::math;

TurboPartition::TurboPartition(const LinearTurbofold &turbofold, const Sequence &seq, const EnergyModel &energy_model,
                               const bool allow_sharp_turn = false)
    : LinearPartition(seq, energy_model, allow_sharp_turn), turbofold(turbofold) {}

inline __attribute__((always_inline)) bool TurboPartition::check_state(StateType type, int i, int j) {
    // Add extrinsic information to State P
    if (type == StateType::P) {
        auto &mp = bestP[j];
        auto item = mp.find(i);
        if (item == mp.end()) {
            return false;
        }
        State &state = item->second;
        const value_type ext_info = turbofold.get_extrinsic_info(seq, i, j);

        if (ext_info <= LOG_ZERO) {
            return false;
        } else {
            state.alpha = LOG_MUL(state.alpha, ext_info * turbofold.lambda);
        }
    }

    return true;
}

void TurboPartition::reset_saved_beams() {
    auto reset = [this](auto &vec) {
        vec.clear();
        vec.reserve(seq_length);
    };
    reset(saved_bestH);
    reset(saved_bestP);
    reset(saved_bestM);
    reset(saved_bestM2);
    reset(saved_bestMulti);
}

void TurboPartition::save_partition_function(const bool move) {
    // if move is true, the original beam data will be cleared
    reset_saved_beams();
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
}

void TurboPartition::calc_prob_accm() {
    prob_accm.upstrm.assign(seq.length(), 0);
    prob_accm.dwnstrm.assign(seq.length(), 0);

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
