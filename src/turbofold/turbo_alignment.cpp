#include <linearx/turbofold/linear_turbofold.hpp>

using namespace linearx::utils;
using namespace linearx::constants::math;
using namespace std;

TurboAlignment::TurboAlignment(const LinearTurbofold &turbofold, Sequence &seq1, Sequence &seq2, const float alpha1,
                               const float alpha2, const float alpha3)
    : LinearAlignment(
          seq1, seq2, [this](int i, int j) { return this->get_match_score(i, j); }, alpha1, alpha2, alpha3),
      turbofold(turbofold) {}

void TurboAlignment::reset_saved_beams(const unsigned beam_size) {
    reset_beam_vector(saved_bestALN, seq_len_sum + 3, beam_size);
    reset_beam_vector(saved_bestINS1, seq_len_sum + 1, beam_size);
    reset_beam_vector(saved_bestINS2, seq_len_sum + 1, beam_size);
}

void TurboAlignment::set_prob_accm(ProbAccm &prob_accm1, ProbAccm &prob_accm2) {
    pm1 = &prob_accm1;
    pm2 = &prob_accm2;
}

value_type TurboAlignment::get_match_score(const int i, const int j) {
    if (i >= seq1.length() || j >= seq2.length()) {
        return 0;
    }
    const value_type t1 = sqrt(pm1->upstrm[i] * pm2->upstrm[j]);
    const value_type t2 = sqrt(pm1->dwnstrm[i] * pm2->dwnstrm[j]);
    const value_type t3 =
        sqrt(max(1 - pm1->upstrm[i] - pm1->dwnstrm[i], 0.0) * max(1 - pm2->upstrm[j] - pm2->dwnstrm[j], 0.0));
    value_type output = ((t1 + t2) * alpha1) + (t3 * alpha2) + (alpha3);
    output = LOG(output);
    return output;
}

void TurboAlignment::save_partition_function(const bool move) {
    // if move is true, the original beam data will be cleared
    reset_saved_beams(run_beam_size_);
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
