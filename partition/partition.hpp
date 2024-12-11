#ifndef PARTITION_HPP
#define PARTITION_HPP

#include <cctype>    // for std::toupper and std::tolower
#include <fstream>   // for file I/O
#include <iostream>  // [DEBUG] for debugging, remove later
#include <set>
#include <unordered_map>
#include <vector>

#include "./../energy_model.hpp"
#include "./../sequence/seq.hpp"
#include "./../shared.hpp"
#include "./../structure/structure.hpp"
#include "./../utility/log_math.hpp"
#include "./../utility/utility.hpp"
#include "./utility.hpp"

class Partition;

struct PartitionFunctionBeam {
    int length;
    bool has_data;
    std::vector<std::unordered_map<int, State>> bestH;
    std::vector<std::unordered_map<int, State>> bestP;
    std::vector<std::unordered_map<int, State>> bestM;
    std::vector<std::unordered_map<int, State>> bestM2;
    std::vector<std::unordered_map<int, State>> bestMulti;
    VectorWithNegOneIndex<State> bestC;

    PartitionFunctionBeam(int length) : length(length), bestC(length) { reset(true); }

    void reset(bool force = false);
    void save(std::unordered_map<int, State> *bestH, std::unordered_map<int, State> *bestP,
              std::unordered_map<int, State> *bestM, std::unordered_map<int, State> *bestM2,
              std::unordered_map<int, State> *bestMulti, VectorWithNegOneIndex<State> &bestC);
    void save(Partition &pf);
};

class Partition {
    inline const static StateType state_types[6] = {H, Multi, P, M2, M, C};

   private:
    std::vector<int> if_tetraloops;
    std::vector<int> if_hexaloops;
    std::vector<int> if_triloops;

   protected:
    std::unordered_map<int, State> *bestH = nullptr;
    std::unordered_map<int, State> *bestP = nullptr;
    std::unordered_map<int, State> *bestM = nullptr;
    std::unordered_map<int, State> *bestM2 = nullptr;
    std::unordered_map<int, State> *bestMulti = nullptr;
    VectorWithNegOneIndex<State> bestC;

   public:
    friend struct PartitionFunctionBeam;

    std::unordered_map<int, double> *bpp = nullptr;
    State *viterbi = nullptr;

    std::vector<int> next_pair[5];
    std::vector<int> prev_pair[5];

    const Seq *sequence;          // pointer to actual Seq object
    const std::vector<int> *seq;  // pointer to encoded sequence in Seq object

    EnergyModel *energy_model;
    InsideMode mode;

    bool allow_sharp_turn;
    bool verbose_output;

    double inside_execution_time = 0.0;
    double outside_execution_time = 0.0;

    Partition(const Seq *sequence, EnergyModel &energy_model, const InsideMode mode, bool allow_sharp_turn = false,
              bool verbose_output = true)
        : bestC(sequence->length() + 1),
          sequence(sequence),
          seq(&sequence->enc_seq),
          energy_model(&energy_model),
          mode(mode),
          allow_sharp_turn(allow_sharp_turn),
          verbose_output(verbose_output) {
        for (int nuc = 0; nuc < 5; ++nuc) {
            prev_pair[nuc].resize(seq->size(), -1);
            next_pair[nuc].resize(seq->size(), seq->size());

            if (nuc == 0) {
                continue; // skip N
            }

            int prev = -1;
            int next = seq->size();

            for (int j = 0; j < seq->size(); ++j) {
                prev_pair[nuc][j] = prev;
                if (Utility::check_valid_pair(nuc, seq->at(j))) {
                    prev = j;
                }
            }
            for (int j = seq->size() - 1; j >= 0; --j) {
                next_pair[nuc][j] = next;
                if (Utility::check_valid_pair(nuc, seq->at(j))) {
                    next = j;
                }
            }
        }

        reset_beams();

        std::string tmp = sequence->sequence;
        std::transform(tmp.begin(), tmp.end(), tmp.begin(), ::toupper);

        energy_model.init_tetra_hex_tri(tmp, seq->size(), if_tetraloops, if_hexaloops, if_triloops);
    }

    void reset_beams();

    virtual void compute_inside(int beam_size = 100);
    virtual void compute_outside(bool use_lazy_outside = true);
    virtual void compute_bpp_matrix();
    double get_bpp(int i, int j) const;
    double get_ensemble_energy();
    std::string get_mfe_structure();
    void print_alpha_beta();

    // methods declared in file forward.cpp
    double beam_prune(std::unordered_map<int, State> &beamstep, int beam_size);
    void update_score(State &state, const int new_score, const double prev_score);
    virtual bool check_state(const StateType type, const int i, const int j);
    void beamstep_H(const int j, const std::vector<int> *next_pair);
    void beamstep_Multi(const int j, const std::vector<int> *next_pair);
    void beamstep_P(const int j, const std::vector<int> *next_pair);
    void beamstep_M2(const int j, const std::vector<int> *next_pair);
    void beamstep_M(const int j);
    void beamstep_C(const int j);

    // methods declared in file backward.cpp
    void hedges_trace_update_helper(std::vector<HEdge> *incoming_hedges, HEdge &best_hedge, const HEdge &new_hedge,
                                    TraceInfo &best_trace, const TraceInfo &new_trace);
    void mfe_backtrack(int i, int j, StateType type, std::string &structure);
    virtual std::pair<int, int> backward_update(const int i, const int j, State &state, const StateType type,
                                                const double edge_threshold);

    virtual void get_incoming_edges_state(const int i, const int j, const StateType type,
                                          std::vector<HEdge> &incoming_hedges);
    TraceInfo get_incoming_hedges_C(const int j, std::vector<HEdge> *incoming_hedges);
    TraceInfo get_incoming_hedges_M(const int i, const int j, std::vector<HEdge> *incoming_hedges);
    TraceInfo get_incoming_hedges_M2(const int i, const int j, std::vector<HEdge> *incoming_hedges);
    TraceInfo get_incoming_hedges_P(const int i, const int j, std::vector<HEdge> *incoming_hedges);
    TraceInfo get_incoming_hedges_Multi(const int i, const int j, std::vector<HEdge> *incoming_hedges);

    std::string &get_threshknot_structure(float threshknot_threshold = 0.3f, int min_helix_size = MIN_HELIX_SIZE) const;

    void dump_bpp(const std::string &filename) const;

    void debug_states() {
        for (int j = 0; j < seq->size(); ++j) {
            for (auto item : bestH[j]) {
                int i = item.first;
                State &state = item.second;
                printf("H[%d][%d]: %.5f\t%.5f\n", i, j, state.alpha, state.beta);
            }
        }
        printf("\n");
        for (int j = 0; j < seq->size(); ++j) {
            for (auto item : bestP[j]) {
                int i = item.first;
                State &state = item.second;
                printf("P[%d][%d]: %.5f\t%.5f\n", i, j, state.alpha, state.beta);
            }
        }
        printf("\n");
        for (int j = -1; j < (int)seq->size(); ++j) {
            printf("C[%d]: %.5f\t%.5f\n", j, bestC[j].alpha, bestC[j].beta);
        }
        printf("\n");
        for (int j = 0; j < seq->size(); ++j) {
            for (auto item : bestM[j]) {
                int i = item.first;
                State &state = item.second;
                printf("M[%d][%d]: %.5f\n", i, j, state.alpha);
            }
        }
    }
};

#endif  // PARTITION_HPP
