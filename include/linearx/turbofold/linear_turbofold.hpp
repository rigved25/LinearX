// linearx/turbofold/linear_turbofold.hpp
#pragma once
#include <linearx/alignment/linear_align.hpp>
#include <linearx/partition/linear_partition.hpp>
#include <linearx/turbofold/utility.hpp>
#include <linearx/utility.hpp>
#include <linearx/turbofold/utils/ProbabilisticModel.hpp>
#include <iostream>

// forward declarations
class TurboPartition;
class TurboAlignment;

class LinearTurbofold {
   private:
    MultiSeq& multi_seq;
    MultiSeq* multi_alignment = nullptr;
    EnergyModel& energy_model;

    std::vector<float> seq_identities;  // sequence identities for all k*(k-1)/2 pairs
    std::vector<TurboPartition> pfs;    // contains partition function objects for k sequences
    std::vector<TurboAlignment> alns;   // contains alignment objects for all k*(k-1)/2 sequence pairs

   public:
    const unsigned alignment_beam_size;
    const unsigned folding_beam_size;
    const float alignment_pruning_threshold;
    const float folding_pruning_threshold;
    const float lambda;  // extrinsic information weight (contribution relative to intrinsic information)
    const float threshknot_threshold;
    const float min_helix_size;
    int curr_itr = 0;

    // algorithm parameters
    bool use_lazy_outside_ = false;
    bool restrict_search_ = false;
    bool astar_viterbi_ = false;
    bool max_marginal_ = false;
    
    // run parameters
    int num_itr_ = 3;
    bool verbose_output_ = false;
    bool save_logs_ = false;
    bool save_probs_ = false;
    std::string out_dir_ = "";

    std::vector<std::vector<std::unordered_map<int, value_type>*>> posterior;

    LinearTurbofold(MultiSeq& multi_seq, EnergyModel& energy_model, const unsigned alignment_beam_size = 100,
                    const unsigned folding_beam_size = 100,
                    const float alignment_pruning_threshold = linearx::constants::limits::DEVIATION_THRESHOLD,
                    const float folding_pruning_threshold = 2 * linearx::constants::limits::DEVIATION_THRESHOLD,
                    const float lambda = 0.3, const float threshknot_threshold = 0.3, const float min_helix_size = 3,
                    const bool allow_sharp_turn = false, const float alpha1 = 1.0, const float alpha2 = 0.8,
                    const float alpha3 = 0.5, const int num_itr = 3, const bool verbose_output = false,
                    const bool save_logs = false, const bool save_probs = false, const std::string out_dir = "");
    ~LinearTurbofold();
    int get_seq_pair_index(int k1, int k2) const;
    
    value_type get_extrinsic_info(const Sequence& x, const int i, const int j) const;
    void dump_coinc_probs(const std::string &filepath, std::unordered_map<int, value_type>* posterior, int seqlen, int k1, int k2);
    // void dump_coinc_aln_probs(const std::string &filepath, const float threshold, std::unordered_map<int, AlnProbs>* coinc_prob, int seqlen);
    // void dump_extrinsic_info(const std::string &filepath, int iteration);
    // void dump_match_scores(const std::string &filepath, int iteration);
    // void dump_pairwise_match_scores(const std::string &filepath, int iteration, int sample_rate = 10);
    
    void run_phmm_alignment(TurboFoldLog& log);
    void multiple_sequence_alignment(TurboFoldLog& log, unsigned int beam_size = 100);
    void run(const bool use_lazy_outside = false, const bool restrict_search = false, const bool astar_viterbi = false, const bool max_marginal = false);
         
};

class TurboPartition final : public LinearPartitionInterface<TurboPartition> {
   private:
    LinearTurbofold& turbofold;
    std::vector<std::unordered_map<int, value_type>> extinfo_cache;
    std::vector<std::unordered_map<int, State>> saved_bestH;
    std::vector<std::unordered_map<int, State>> saved_bestP;
    std::vector<std::unordered_map<int, State>> saved_bestM;
    std::vector<std::unordered_map<int, State>> saved_bestM2;
    std::vector<std::unordered_map<int, State>> saved_bestMulti;
    value_type total_inside;

   public:
    friend class LinearTurbofold;
    linearx::utils::ProbAccm prob_accm;

    TurboPartition(LinearTurbofold& turbofold, const Sequence& seq, const EnergyModel& energy_model,
                   const bool allow_sharp_turn = false);

    State* get_saved_state(const StateType type, const int i, const int j) noexcept {
        std::unordered_map<int, State>* beam = nullptr;
        switch (type) {
            case H:
                beam = &saved_bestH[j];
                break;
            case Multi:
                beam = &saved_bestMulti[j];
                break;
            case P:
                beam = &saved_bestP[j];
                break;
            case M2:
                beam = &saved_bestM2[j];
                break;
            case M:
                beam = &saved_bestM[j];
                break;
            case C:
                return nullptr;  // control will never reach here
            default:
                return nullptr;
        }
        const auto it = beam->find(i);
        return it == beam->end() ? nullptr : &it->second;
    }

    unsigned long beam_prune(std::unordered_map<int, State>& beamstep, const unsigned beam_size,
                                    const StateType type, const int j) override final {
        if (turbofold.curr_itr > 0 && type == StateType::P) {
            // Add extrinsic information to State P
            auto it = bestP[j].begin();
            while (it != bestP[j].end()) {
                const int i = it->first;
                State& state = it->second;
                double ext_info = turbofold.get_extrinsic_info(seq, i, j);
                if (ext_info <= linearx::math::LOG_ZERO) {
                    it = bestP[j].erase(it);  // erase the element and update the iterator
                } else {
                    ext_info *= turbofold.lambda;
                    extinfo_cache[j][i] = ext_info;
                    state.alpha = LOG_MUL(state.alpha, ext_info);  // adjust the weight as needed
                    ++it;                                          // only increment if not erased
                }
            }
        }
        return LinearPartitionInterface<TurboPartition>::beam_prune(beamstep, beam_size, type, j);
    }

    State* check_state(const StateType type, const int i, const int j) {
        if (turbofold.restrict_search_ && turbofold.curr_itr > 0 && type != StateType::C) {
            State* state = TurboPartition::get_saved_state(type, i, j);
            if (!state) {
                return nullptr;
            }

            if (type != StateType::H) {  // no pruning for H states
                if (turbofold.use_lazy_outside_) {
                    if (state->beta <= linearx::math::LOG_ZERO) {
                        return nullptr;
                    }
                } else if (LOG_DIV(LOG_MUL(state->alpha, state->beta), total_inside) <=
                           -turbofold.folding_pruning_threshold) {
                    return nullptr;
                }
            }
        }

        return LinearPartitionInterface<TurboPartition>::get_state<true>(type, i, j);
    }

    template <Mode mode, typename F>
    void update_state_alpha(F&& get_score, const State* left, const State* right, const StateType type,
                                   const unsigned i, const unsigned j) {
        if constexpr (mode == Mode::PARTITION_INSIDE) {
            State* next_state = TurboPartition::check_state(type, i, j);
            if (next_state) {
                const value_type weight = get_score() * linearx::constants::energy::INV_KT;
                next_state->alpha =
                    LOG_SUM(next_state->alpha, (left ? left->alpha : 0) + (right ? right->alpha : 0) + weight);
            }
        }
    }

    template <typename F>
    void update_state_beta(F&& get_score, State* left, State* right, const StateType type, const unsigned i,
                                  const unsigned j) {
        const State* next_state = LinearPartitionInterface<TurboPartition>::get_state<false>(type, i, j);
        if (next_state) {
            value_type weight = get_score() * linearx::constants::energy::INV_KT;
            if (turbofold.curr_itr > 0 && type == StateType::P) {
                weight = LOG_MUL(weight, extinfo_cache[j][i]);  // adjust weight with extrinsic info
            }
            if (!right) {
                left->beta = LOG_SUM(left->beta, weight + next_state->beta);
            } else {
                left->beta = LOG_SUM(left->beta, right->alpha + weight + next_state->beta);
                right->beta = LOG_SUM(right->beta, left->alpha + weight + next_state->beta);
            }
        }
    }

    void calc_prob_accm();
    // void reset_saved_beams(const unsigned beam_size);
    void save_partition_function(const bool move, const unsigned beam_size);
};

class TurboAlignment final : public LinearAlignmentInterface<TurboAlignment> {
   private:
    const LinearTurbofold& turbofold;
    std::vector<std::unordered_map<std::pair<int, int>, HState, linearx::utils::PairHash>> saved_bestALN;
    std::vector<std::unordered_map<std::pair<int, int>, HState, linearx::utils::PairHash>> saved_bestINS1;
    std::vector<std::unordered_map<std::pair<int, int>, HState, linearx::utils::PairHash>> saved_bestINS2;
    value_type total_inside_alignment = linearx::math::LOG_ZERO;

   public:
    friend class LinearTurbofold;

    TurboAlignment(const LinearTurbofold& turbofold, Sequence& seq1, Sequence& seq2, const float alpha1 = 1.0,
                   const float alpha2 = 0.8, const float alpha3 = 0.5);

    // void reset_saved_beams(const unsigned beam_size);
    void save_partition_function(const bool move);

    void save_partition_function_old(const bool move);
    
    std::vector<std::unordered_map<std::pair<int, int>, HState, linearx::utils::PairHash>>& get_saved_beams(
        HStateType type) {
        switch (type) {
            case ALN:
                return saved_bestALN;
            case INS1:
                return saved_bestINS1;
            case INS2:
                return saved_bestINS2;
            default:
                return saved_bestALN;  // unreachable
        }
    }

    HState* get_saved_state(const HStateType type, const int i, const int j) noexcept {
        auto& beams = get_saved_beams(type);
        const int s = i + j;
        if (s >= beams.size()) {
            return nullptr;  // Out of bounds - no saved states on first iteration
        }
        auto& beam = beams[s];
        const std::pair<int, int> key = {i, j};
        const auto it = beam.find(key);
        if (it != beam.end()) {
            return &it->second;
        }
        return nullptr;
    }

    HState* check_state(const HStateType type, const int i, const int j) {
        if (turbofold.max_marginal_ && turbofold.curr_itr == 1) {
            // Filter partition mode using BEST outside results.
            // Since LHS_best < LHS_partition, we CANNOT use LHS_best > -threshold to filter.
            // Instead, just check if the state was visited during BEST outside (has valid beta).
            // The BEST outside threshold already does the pruning.
            const HState* state = TurboAlignment::get_saved_state(type, i, j);
            if (!state || state->beta <= linearx::math::LOG_ZERO) {
                return nullptr;
            }
        }

        if (turbofold.restrict_search_ && turbofold.curr_itr > 1) {
            HState* state = TurboAlignment::get_saved_state(type, i, j);
            if (!state) {
                return nullptr;
            }
            
            if (turbofold.use_lazy_outside_ && state->beta <= linearx::math::LOG_ZERO) {
                return nullptr;
            } 
            else if (LOG_DIV(LOG_MUL(state->alpha, state->beta),
                    total_inside_alignment) <= -turbofold.alignment_pruning_threshold) {
                return nullptr;
            }
            
        }

        return LinearAlignmentInterface<TurboAlignment>::get_state(type, i, j, true);
    }

    HState* check_state_AStar(const HStateType type, const int i, const int j) {
        if(turbofold.curr_itr > 1){
            if (turbofold.use_lazy_outside_ || turbofold.restrict_search_) {
                HState* state = TurboAlignment::get_saved_state(type, i, j);
                if (!state) {
                    return nullptr;
                }

                // With AStar: make lazy outside work irrespective of posterior pruning
                if (turbofold.use_lazy_outside_ && state->beta <= linearx::math::LOG_ZERO) {
                    return nullptr;
                }
                if (turbofold.restrict_search_ && LOG_DIV(LOG_MUL(state->alpha, state->beta),
                        total_inside_alignment) <= -turbofold.alignment_pruning_threshold) {
                    return nullptr;
                }
            }            
        }
        return LinearAlignmentInterface<TurboAlignment>::get_state(type, i, j, true);
    }

    // HState* check_state_best(const HStateType type, const int i, const int j) {


    //     // Only filter in PARTITION mode for iterations > 1.
    //     // BEST mode (A*) in iteration 2+ should NOT be filtered because:
    //     // - Probability weights change between iterations
    //     // - States pruned by lazy outside in iter N may be needed in iter N+1
    //     // - A* needs freedom to find the optimal path with new weights
    //     if (turbofold.restrict_search_ && turbofold.curr_itr > 1) {
    //         static unsigned long missing_end_reports = 0;
    //         static unsigned long missing_state_reports = 0;
    //         HState* end_state = TurboAlignment::get_saved_state(HStateType::ALN, seq1.length() + 1, seq2.length() + 1);
    //         if (!end_state) {
    //             return LinearAlignmentInterface<TurboAlignment>::get_state(type, i, j, true);
    //         }
    //         HState* state = TurboAlignment::get_saved_state(type, i, j);
    //         if (!state) {
    //             return LinearAlignmentInterface<TurboAlignment>::get_state(type, i, j, true);
    //         }
    //         if (turbofold.use_lazy_outside_) {
    //             if (state->beta <= linearx::math::LOG_ZERO) {
    //                 return nullptr;
    //             }
    //         } else if (LOG_DIV(LOG_MUL(state->alpha, state->beta), end_state->alpha) <=
    //                    -turbofold.alignment_pruning_threshold) {
    //             return nullptr;
    //         }
    //     }
    //     return LinearAlignmentInterface<TurboAlignment>::get_state(type, i, j, true);
    // }
};
