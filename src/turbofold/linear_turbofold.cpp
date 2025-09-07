#include <linearx/turbofold/linear_turbofold.hpp>

LinearTurbofold::LinearTurbofold(MultiSeq &msa, EnergyModel &energy_model, float lambda,
                                 float alignment_pruning_threshold, float folding_pruning_threshold,
                                 float threshknot_threshold, float min_helix_size, bool allow_sharp_turn, float alpha1,
                                 float alpha2, float alpha3)
    : msa(msa),
      energy_model(energy_model),
      lambda(lambda),
      alignment_pruning_threshold(alignment_pruning_threshold),
      folding_pruning_threshold(folding_pruning_threshold),
      threshknot_threshold(threshknot_threshold),
      min_helix_size(min_helix_size) {
    unsigned num_pairs = msa.size() * (msa.size() - 1) / 2;
    // reserve space for sequence pairs and sequence identities
    alns.reserve(num_pairs);
    pfs.reserve(msa.size());
    seq_identities.reserve(num_pairs);

    for (int i = 0; i < msa.size(); i++) {
        // create partition function object for each sequence
        pfs.emplace_back(*this, msa[i], energy_model, allow_sharp_turn);

        // enumerate all possible k^2 sequence pairs and create LinearAlign objects
        for (int j = i + 1; j < msa.size(); j++) {
            alns.emplace_back(*this, msa[i], msa[j], alpha1, alpha2, alpha3);
        }
    }
}

int LinearTurbofold::get_seq_pair_index(const int k1, const int k2) {
    // ensure the indices are within bounds
    if (k1 >= msa.size() || k2 >= msa.size()) {
        throw std::out_of_range("Seq index k out of range");
    }

    // ensure k1 != k2 and k1 < k2 by swapping if necessary
    int a = k1, b = k2;
    if (k1 == k2) {
        throw std::invalid_argument("k1 and k2 must be different");
    }
    if (k1 > k2) {
        std::swap(a, b);
    }

    // number of sequences
    int n = msa.size();

    // compute the index using the triangular indexing formula
    int index = a * (2 * n - a - 1) / 2 + (b - a - 1);

    return index;
}

value_type LinearTurbofold::get_extrinsic_info(const Sequence &x, int i, int j) const {
    value_type ext_info = 1.0;
    return ext_info;
}