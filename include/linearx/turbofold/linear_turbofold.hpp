// linearx/turbofold/linear_turbofold.hpp
#pragma once
#include <linearx/alignment/linear_align.hpp>
#include <linearx/partition/linear_partition.hpp>
#include <linearx/utility.hpp>

class LinearTurbofold {
   private:
    MultiSeq &msa;
    EnergyModel &energy_model;

    std::vector<float> seq_identities;  // sequence identities for all k*(k-1)/2 pairs
    std::vector<TurboPartition> pfs;    // contains partition function objects for k sequences
    std::vector<TurboAlignment> alns;   // contains alignment objects for all k*(k-1)/2 sequence pairs

   public:
    const float lambda;  // extrinsic information weight (contribution relative to intrinsic information)
    const float alignment_pruning_threshold;
    const float folding_pruning_threshold;
    const float threshknot_threshold;
    const float min_helix_size;

    LinearTurbofold(MultiSeq &msa, EnergyModel &energy_model, float lambda = 0.3,
                    float alignment_pruning_threshold = linearx::constants::limits::DEVIATION_THRESHOLD,
                    float folding_pruning_threshold = linearx::constants::limits::DEVIATION_THRESHOLD,
                    float threshknot_threshold = 0.3, float min_helix_size = 3, bool allow_sharp_turn = false,
                    float alpha1 = 1.0, float alpha2 = 0.8, float alpha3 = 0.5);

    int get_seq_pair_index(const int k1, const int k2);
    value_type get_extrinsic_info(const Sequence &x, int i, int j) const;
    // void reset_extinf_cache();
    void run();
};

class TurboPartition final : public LinearPartition {
   private:
    const LinearTurbofold &turbofold;
    std::vector<std::unordered_map<int, State>> saved_bestH;
    std::vector<std::unordered_map<int, State>> saved_bestP;
    std::vector<std::unordered_map<int, State>> saved_bestM;
    std::vector<std::unordered_map<int, State>> saved_bestM2;
    std::vector<std::unordered_map<int, State>> saved_bestMulti;

   public:
    friend class LinearTurbofold;
    linearx::utils::ProbAccm prob_accm;

    TurboPartition(const LinearTurbofold &turbofold, const Sequence &seq, const EnergyModel &energy_model,
                   const bool allow_sharp_turn = false);

    inline bool check_state(StateType type, int i, int j) override;
    void calc_prob_accm();

    void reset_saved_beams();
    void save_partition_function(const bool move);
};

class TurboAlignment final : public LinearAlignment {
   private:
    const LinearTurbofold &turbofold;
    std::vector<std::unordered_map<std::pair<int, int>, HState, linearx::utils::PairHash>> saved_bestALN;
    std::vector<std::unordered_map<std::pair<int, int>, HState, linearx::utils::PairHash>> saved_bestINS1;
    std::vector<std::unordered_map<std::pair<int, int>, HState, linearx::utils::PairHash>> saved_bestINS2;

   public:
    friend class LinearTurbofold;

    TurboAlignment(const LinearTurbofold &turbofold, Sequence &seq1, Sequence &seq2, const float alpha1 = 1.0,
                   const float alpha2 = 0.8, const float alpha3 = 0.5);

    void reset_saved_beams();
    void save_partition_function(const bool move);
};
