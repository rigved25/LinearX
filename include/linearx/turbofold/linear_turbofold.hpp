// linearx/turbofold/linear_turbofold.hpp
#pragma once
#include <linearx/alignment/linear_align.hpp>
#include <linearx/partition/linear_partition.hpp>
#include <linearx/turbofold/utility.hpp>
#include <linearx/utility.hpp>

// forward declarations
class TurboPartition;
class TurboAlignment;

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
    int curr_itr = 0;

    LinearTurbofold(MultiSeq &msa, EnergyModel &energy_model, float lambda = 0.3,
                    float alignment_pruning_threshold = linearx::constants::limits::DEVIATION_THRESHOLD,
                    float folding_pruning_threshold = linearx::constants::limits::DEVIATION_THRESHOLD,
                    float threshknot_threshold = 0.3, float min_helix_size = 3, bool allow_sharp_turn = false,
                    float alpha1 = 1.0, float alpha2 = 0.8, float alpha3 = 0.5);

    int get_seq_pair_index(const int k1, const int k2) const;
    value_type get_extrinsic_info(const Sequence &x, const int i, const int j);
    void run(const int num_itr = 3, const bool use_lazy_outside = true, const bool use_prev_itr_beta = false,
             const bool restrict_search = false, const bool verbose_output = false, const bool save_logs = false,
             const bool save_probs = false, const std::string out_dir = "./ltf_output/");
};

class TurboPartition final : public LinearPartition {
   private:
    LinearTurbofold &turbofold;
    std::vector<std::unordered_map<int, State>> saved_bestH;
    std::vector<std::unordered_map<int, State>> saved_bestP;
    std::vector<std::unordered_map<int, State>> saved_bestM;
    std::vector<std::unordered_map<int, State>> saved_bestM2;
    std::vector<std::unordered_map<int, State>> saved_bestMulti;

   public:
    friend class LinearTurbofold;
    linearx::utils::ProbAccm prob_accm;

    TurboPartition(LinearTurbofold &turbofold, const Sequence &seq, const EnergyModel &energy_model,
                   const bool allow_sharp_turn = false);

    inline unsigned long beam_prune(const int j, const StateType type,
                                    std::unordered_map<int, State> &beamstep) override final {
        if (turbofold.curr_itr > 0 && type == StateType::P) {
            // Add extrinsic information to State P
            auto it = bestP[j].begin();
            while (it != bestP[j].end()) {
                const int i = it->first;
                State &state = it->second;
                const double ext_info = turbofold.get_extrinsic_info(seq, i, j);
                if (ext_info <= linearx::math::LOG_ZERO) {
                    it = bestP[j].erase(it);  // erase the element and update the iterator
                } else {
                    state.alpha = LOG_MUL(state.alpha, ext_info * turbofold.lambda);  // adjust the weight as needed
                    ++it;                                                             // only increment if not erased
                }
            }
        }
        return LinearPartition::beam_prune(j, type, beamstep);
    }
    inline void update_score(const StateType type, const int i, const int j, State &state, int new_score,
                             value_type prev_score = 0) override final {
        state.alpha = LOG_SUM(state.alpha, prev_score + (new_score * linearx::constants::energy::INV_KT));
    }
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
