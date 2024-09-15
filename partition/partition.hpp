#ifndef PARTITION_HPP
#define PARTITION_HPP

#include "./../energy_model.hpp"
#include "./../sequence/seq.hpp"
#include "./../utility/log_math.hpp"
#include "./../utility/utility.hpp"
#include "./../shared.hpp"
#include <iostream> // [DEBUG] for debugging, remove later
#include <unordered_map>
#include <vector>

struct State {
    double alpha;
    double beta;
    State() : alpha(VALUE_MIN), beta(VALUE_MIN) {};
};

enum StateType {
    H,
    Multi,
    P,
    M2,
    M,
    C,
};

enum InsideMode {
    MFE,
    PARTITION,
};

struct TraceInfo {
    int i;
    int j;
    int t; // split point
    StateType type_left;
    StateType type_right;

    TraceInfo() : i(-1), j(-1), t(-1), type_left(H), type_right(H) {};

    void set(const int i, const int j, const int t, const StateType type_left, const StateType type_right) {
        this->i = i;
        this->j = j;
        this->t = t;
        this->type_left = type_left;
        this->type_right = type_right;
    }
};

struct HEdge {
    float weight;
    State *left;
    State *right; // right=null <=> Edge

    HEdge() : weight(VALUE_MIN), left(nullptr), right(nullptr) {};
    HEdge(float weight, State *left, State *right) : weight(weight), left(left), right(right) {};

    void set(float weight, State *left, State *right) {
        this->weight = weight;
        this->left = left;
        this->right = right;
    }
};

class Partition {

    inline const static StateType state_types[6] = {H, Multi, P, M2, M, C};

  public:
    std::unordered_map<int, State> *bestH;
    std::unordered_map<int, State> *bestP;
    std::unordered_map<int, State> *bestM;
    std::unordered_map<int, State> *bestM2;
    std::unordered_map<int, State> *bestMulti;
    std::unordered_map<int, State> bestC;

    std::unordered_map<int, float> *bpp;
    State *viterbi = nullptr;

    std::vector<int> next_pair[5];
    std::vector<int> prev_pair[5];

    const Seq *sequence; // pointer to actual Seq object
    const std::vector<int> *seq; // pointer to encoded sequence in Seq object
    
    EnergyModel *energy_model;
    InsideMode mode;

    const int beam_size;
    bool allow_sharp_turn;
    bool verbose_output;

    float inside_execution_time = 0.0;
    float outside_execution_time = 0.0;

    ProbAccm prob_accm;

    Partition(const Seq *sequence, EnergyModel &energy_model, const InsideMode mode, const int beam_size = 100,
              bool allow_sharp_turn = false, bool verbose_output = true)
        : sequence(sequence), seq(&sequence->enc_seq), energy_model(&energy_model), mode(mode),
          beam_size(beam_size), allow_sharp_turn(allow_sharp_turn), verbose_output(verbose_output) {

        bestH = new std::unordered_map<int, State>[seq->size()];
        bestP = new std::unordered_map<int, State>[seq->size()];
        bestM = new std::unordered_map<int, State>[seq->size()];
        bestM2 = new std::unordered_map<int, State>[seq->size()];
        bestMulti = new std::unordered_map<int, State>[seq->size()];
        // bestC = new State[seq.size()];

        bpp = nullptr; // initialize later if required

        for (int nuc = 1; nuc < 5; ++nuc) {
            prev_pair[nuc].resize(seq->size(), -1);
            next_pair[nuc].resize(seq->size(), seq->size());

            int prev = -1;
            int next = seq->size();

            for (int j = 0; j < seq->size(); ++j) {
                prev_pair[nuc][j] = prev;
                if (check_valid_pair(nuc, seq->at(j))) {
                    prev = j;
                }
            }
            for (int j = seq->size() - 1; j >= 0; --j) {
                next_pair[nuc][j] = next;
                if (check_valid_pair(nuc, seq->at(j))) {
                    next = j;
                }
            }
        }

        bestC[-1].alpha = 0;
        bestC[seq->size() - 1].beta = 0;
        if (seq->size() > 0)
            bestC[0].alpha = 0.0;
        if (seq->size() > 1)
            bestC[1].alpha = 0.0;
    }

    virtual void compute_inside();
    virtual void compute_outside();
    void compute_bpp_matrix();

    std::string get_mfe_structure() {
        std::string structure(seq->size(), '.');
        mfe_backtrack(0, seq->size() - 1, C, structure);
        return structure;
    }

    double get_ensemble_energy() { return -kT * (bestC[seq->size() - 1].alpha) / 100.0; } // -kT log(Q(x))

    void print_alpha_beta() {
        std::cout << "\nAlpha(C, n - 1): " << bestC[seq->size() - 1].alpha << std::endl;
        std::cout << "Beta(C, -1): " << bestC[-1].beta << std::endl;
    }

    void calc_prob_accm();

    void debug_states() {

        // for (int j = 0; j < seq->size(); ++j) {
        //     for (auto item : bestH[j]) {
        //         int i = item.first;
        //         State &state = item.second;
        //         printf("H[%d][%d]: %.5f\t%.5f\n", i, j, state.alpha, state.beta);
        //     }
        // }
        // printf("\n");
        // for (int j = 0; j < seq->size(); ++j) {
        //     for (auto item : bestP[j]) {
        //         int i = item.first;
        //         State &state = item.second;
        //         printf("P[%d][%d]: %.5f\t%.5f\n", i, j, state.alpha, state.beta);
        //     }
        // }
        printf("\n");
        for (int j = -1; j < (int)seq->size(); ++j) {
            printf("C[%d]: %.5f\t%.5f\n", j, bestC[j].alpha, bestC[j].beta);
        }
        // printf("\n");
        // for (int j = 0; j < seq->size(); ++j) {
        //     for (auto item : bestM[j]) {
        //         int i = item.first;
        //         State &state = item.second;
        //         printf("M[%d][%d]: %.5f\n", i, j, state.alpha);
        //     }
        // }
    }

    float get_bpp(const int i, const int j);

    // methods declared in file forward.cpp
    float beam_prune(std::unordered_map<int, State> &beamstep);
    void update_score(State &state, const int new_score, const double prev_score);
    void beamstep_H(const int j, const std::vector<int> *next_pair);
    void beamstep_Multi(const int j, const std::vector<int> *next_pair);
    void beamstep_P(const int j, const std::vector<int> *next_pair);
    void beamstep_M2(const int j, const std::vector<int> *next_pair);
    void beamstep_M(const int j);
    void beamstep_C(const int j);

    // methods declared in file backward.cpp
    void hedges_trace_helper(std::vector<HEdge> *incoming_hedges, HEdge &best_hedge, HEdge &new_hedge,
                             TraceInfo &best_trace, TraceInfo &new_trace);
    void mfe_backtrack(int i, int j, StateType type, std::string &structure);
    std::pair<int, int> backward_update(const int i, const int j, State &state, const StateType type,
                                        const float edge_threshold);
    TraceInfo get_incoming_hedges_C(const int j, std::vector<HEdge> *incoming_hedges);
    TraceInfo get_incoming_hedges_M(const int i, const int j, std::vector<HEdge> *incoming_hedges);
    TraceInfo get_incoming_hedges_M2(const int i, const int j, std::vector<HEdge> *incoming_hedges);
    TraceInfo get_incoming_hedges_P(const int i, const int j, std::vector<HEdge> *incoming_hedges);
    TraceInfo get_incoming_hedges_Multi(const int i, const int j, std::vector<HEdge> *incoming_hedges);
};

#endif // PARTITION_HPP