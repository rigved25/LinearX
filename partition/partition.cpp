#include "partition.hpp"

float Partition::get_bpp(const int i, const int j) {
    if (this->bpp)
        return this->bpp[j][i];
    if (bestP[j].find(i) == bestP[j].end())
        return 0;

    State &state = bestP[j][i];
    float bpp = state.alpha + state.beta - viterbi->alpha;

    // if(bpp <= -DEVIATION_THRESHOLD) {
    //     return 0;
    // }

    bpp = Fast_Exp(bpp);

    // if (bpp > 1.00001) {
    //     printf("[WARNING] BPP value too high, something is wrong! bpp(%d, %d): %.5f\n", i, j, bpp);
    // }

    return std::max(bpp, 1.0f);
}

void Partition::compute_bpp_matrix() {
    delete[] this->bpp;                                          // delete any existing memory
    this->bpp = new std::unordered_map<int, float>[seq->size()]; // reallocate memory
    for (int j = 0; j < seq->size(); ++j) {
        for (auto &item : bestP[j]) {
            const int i = item.first;
            State &state = item.second;
            const float prob = Fast_Exp(state.alpha + state.beta - viterbi->alpha);
            this->bpp[j][i] = prob;
            if (prob > 1) {
                printf("[WARNING] BPP value too high, something is wrong! i: %d, j: %d\n", i, j);
            }
        }
    }
}

void Partition::calc_prob_accm() {
    prob_accm.upstrm.resize(seq->size(), 0);
    prob_accm.dwnstrm.resize(seq->size(), 0);

    for (int j = 0; j < seq->size(); ++j) {
        for (const auto &item : bestP[j]) {
            const int i = item.first;
            const float score = get_bpp(i, j);
            prob_accm.upstrm[j] += score;
            prob_accm.dwnstrm[i] += score;
        }
    }

    // print scores
    // for (int i = 0; i < seq.size(); ++i) {
    //     printf("i: %d, prob_accm_upstrm: %.5f, prob_accm_dwnstrm: %.5f\n", i, prob_accm.upstrm[i],
    //            prob_accm.dwnstrm[i]);
    // }
}