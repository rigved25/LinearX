// linearx/partition/linear_partition.hpp
#pragma once
#include <linearx/energy/energy_model.hpp>
#include <linearx/partition/utility.hpp>
#include <linearx/sequence/sequence.hpp>
#include <linearx/sequence/structure.hpp>
#include <linearx/utility.hpp>
#include <unordered_set>

class LinearPartition {
    inline const static StateType state_types[6] = {H, Multi, P, M2, M, C};

   private:
    std::vector<int> if_tetraloops;
    std::vector<int> if_hexaloops;
    std::vector<int> if_triloops;
    Mode run_mode_;
    int run_beam_size_;
    std::vector<std::pair<value_type, int>> beam_scores;
    std::vector<HEdge> incoming_hedges;
    std::vector<HEdge *> saved_hedges;
    TraceInfo best_trace;
    HEdge best_hedge;
    value_type _last_inside_exec_time = 1.0;

    // methods declared in file backward.cpp
    void update_best_trace(const HEdge &new_hedge, const TraceInfo &new_trace);

    void mfe_backtrack(const int i, const int j, const StateType type, Structure &structure);

   protected:
    std::vector<std::unordered_map<int, State>> bestH;
    std::vector<std::unordered_map<int, State>> bestP;
    std::vector<std::unordered_map<int, State>> bestM;
    std::vector<std::unordered_map<int, State>> bestM2;
    std::vector<std::unordered_map<int, State>> bestMulti;
    VectorWithNegOneIndex<State> bestC;

   public:
    friend struct PartitionFunctionBeam;

    std::vector<std::unordered_map<int, value_type>> bpp;

    std::array<std::vector<int>, 5> next_pair;
    std::array<std::vector<int>, 5> prev_pair;

    const Sequence &seq;
    const int seq_length;
    const EnergyModel &energy_model;
    const bool allow_sharp_turn;

    LinearPartition(const Sequence &seq, const EnergyModel &energy_model, const bool allow_sharp_turn = false);
    State &get_viterbi();

    void reset_beams();
    PartitionInsideLog compute_inside(const Mode mode, const unsigned beam_size = 100,
                                      const bool verbose_output = true);
    PartitionOutsideLog compute_outside(
        const bool use_lazy_outside = true,
        const value_type deviation_threshold = linearx::constants::limits::DEVIATION_THRESHOLD,
        const bool verbose_output = true);

    void compute_bpp_matrix();
    value_type get_ensemble_energy() const;
    void print_alpha_beta() const;
    Structure get_mfe_structure();

    inline value_type get_bpp(const int i, const int j) const {
        const auto &bpp_j = bpp[j];
        const auto item = bpp_j.find(i);
        return item == bpp_j.end() ? 0.0 : item->second;
    }
    inline virtual unsigned long beam_prune(const int j, const StateType type,
                                            std::unordered_map<int, State> &beamstep) {
        if (run_beam_size_ == 0 || beamstep.size() <= run_beam_size_) {
            return 0;
        }
        unsigned long num_pruned = 0;
        beam_scores.clear();
        beam_scores.reserve(beamstep.size());
        for (auto &item : beamstep) {
            const int i = item.first;
            const State &cand = item.second;
            const int k = i - 1;
            const value_type newalpha = ((k >= 0 ? bestC[k].alpha : 0) + cand.alpha);
            beam_scores.emplace_back(newalpha, i);
        }
        const value_type threshold =
            linearx::utils::quickselect(beam_scores, 0, beam_scores.size() - 1, beam_scores.size() - run_beam_size_);
        for (auto &p : beam_scores) {
            if (p.first < threshold) {
                beamstep.erase(p.second);
                num_pruned++;
            }
        }
        return num_pruned;
    }
    inline virtual void update_score(const StateType type, const int i, const int j, State &state, int new_score,
                                     value_type prev_score = 0) {
        if (run_mode_ == BEST) {
            state.alpha = std::max(state.alpha, prev_score + new_score);
        } else {
            state.alpha = LOG_SUM(state.alpha, prev_score + (new_score * linearx::constants::energy::INV_KT));
        }
    }

    // methods declared in file forward.cpp
    unsigned long beamstep_H(const int j);
    unsigned long beamstep_Multi(const int j);
    unsigned long beamstep_P(const int j);
    unsigned long beamstep_M2(const int j);
    unsigned long beamstep_M(const int j);
    void beamstep_C(const int j);

    // methods declared in file backward.cpp
    std::pair<unsigned long, unsigned long> backward_update(const int i, const int j, State &state,
                                                            const StateType type, const value_type edge_threshold);

    void get_incoming_edges_state(const int i, const int j, const StateType type);

    template <Mode mode>
    void get_incoming_hedges_C(const int j);
    template <Mode mode>
    void get_incoming_hedges_M(const int i, const int j);
    template <Mode mode>
    void get_incoming_hedges_M2(const int i, const int j);
    template <Mode mode>
    void get_incoming_hedges_P(const int i, const int j);
    template <Mode mode>
    void get_incoming_hedges_Multi(int i, int j);

    std::string get_threshknot_structure(float threshknot_threshold = 0.3f,
                                         int min_helix_size = linearx::constants::energy::MIN_HELIX_SIZE) const;

    void dump_bpp(const std::string &out_dir) const;

    inline void debug_states() const {
        // printf("\n");
        // for (int j = 0; j < seq.length(); ++j) {
        //     for (const auto &item : bestP[j]) {
        //         int i = item.first;
        //         const State &state = item.second;
        //         printf("P[%d][%d]: %.5f\t%.5f\n", i, j, state.alpha, state.beta);
        //     }
        // }
        // printf("\n");

        // print all p for last j
        for (const auto &item : bestP[seq_length - 1]) {
            int i = item.first;
            const State &state = item.second;
            printf("P[%d][%d]: %.5f\t%.5f\n", i, seq_length - 1, state.alpha, state.beta);
        }
        // for (int j = -1; j < (int)seq.length(); ++j) {
        //     printf("C[%d]: %.5f\t%.5f\n", j, bestC[j].alpha, bestC[j].beta);
        // }
        // printf("\n");
        // for (int j = 0; j < seq.length(); ++j) {
        //     for (const auto &item : bestM[j]) {
        //         int i = item.first;
        //         const State &state = item.second;
        //         printf("M[%d][%d]: %.5f\n", i, j, state.alpha);
        //     }
        // }
    }
};
