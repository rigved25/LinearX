// linearx/partition/linear_partition.hpp
#pragma once
#include <linearx/energy/energy_model.hpp>
#include <linearx/partition/utils.hpp>
#include <linearx/sequence/sequence.hpp>
#include <linearx/sequence/structure.hpp>
#include <linearx/utility/utils.hpp>

class LinearPartition {
    inline const static StateType state_types[6] = {H, Multi, P, M2, M, C};

   private:
    std::vector<int> if_tetraloops;
    std::vector<int> if_hexaloops;
    std::vector<int> if_triloops;
    std::vector<std::pair<double, int>> beam_scores;
    std::vector<HEdge> incoming_hedges;
    std::vector<HEdge *> saved_hedges;
    TraceInfo best_trace;
    HEdge best_hedge;
    double _last_inside_exec_time = 1.0;

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

    virtual bool check_state(const StateType type, const int i, const int j) const;

   public:
    friend struct PartitionFunctionBeam;

    std::vector<std::unordered_map<int, double>> bpp;

    std::array<std::vector<int>, 5> next_pair;
    std::array<std::vector<int>, 5> prev_pair;

    const Sequence &seq;
    const int seq_length;
    const EnergyModel &energy_model;
    const bool allow_sharp_turn;

    // std::vector<double> inside_exec_times;
    // std::vector<double> outside_exec_times;

    LinearPartition(const Sequence &seq, const EnergyModel &energy_model, const bool allow_sharp_turn = false);
    State &get_viterbi();

    void reset_beams(bool freeMemory = true);
    virtual PartitionInsideLog compute_inside(const Mode mode, const unsigned beam_size = 100,
                                              const bool verbose_output = true);
    PartitionOutsideLog compute_outside(
        const double deviation_threshold = linearx::constants::limits::DEVIATION_THRESHOLD,
        const bool verbose_output = true);

    void compute_bpp_matrix();
    double get_ensemble_energy() const;
    void print_alpha_beta() const;
    Structure get_mfe_structure();
    inline double get_bpp(const int i, const int j) const {
        const auto &bpp_j = bpp[j];
        const auto &item = bpp_j.find(i);
        if (item == bpp_j.end()) {
            return 0.0;
        }
        return item->second;
    }

    // methods declared in file forward.cpp
    unsigned beam_prune(std::unordered_map<int, State> &beamstep, const int beam_size);
    void update_score(const Mode mode, State &state, const int new_score, const double prev_score = 0);
    void beamstep_H(const int j, const Mode mode);
    void beamstep_Multi(const int j, const Mode mode);
    void beamstep_P(const int j, const Mode mode);
    void beamstep_M2(const int j, const Mode mode);
    void beamstep_M(const int j, const Mode mode);
    void beamstep_C(const int j, const Mode mode);

    // methods declared in file backward.cpp
    std::pair<unsigned long, unsigned long> backward_update(const int i, const int j, State &state,
                                                            const StateType type, const double edge_threshold);

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

    void dump_bpp(const std::string &filename) const;

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
