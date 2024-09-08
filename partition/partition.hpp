#ifndef PARTITION_HPP
#define PARTITION_HPP

#include "./../energy_model.hpp"
#include "./../sequence/seq.hpp"
#include "./../utility/log_math.hpp"
#include "./../utility/utility.hpp"
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
    // State *bestC;
    std::unordered_map<int, State> bestC;

    std::vector<int> next_pair[5];
    std::vector<int> prev_pair[5];

    const std::vector<int> seq;
    EnergyModel *energy_model;
    InsideMode mode;
    const int beam_size;
    bool allow_sharp_turn;

    int num_alpha_updates = 0;
    int num_beta_updates = 0;

    Partition(const Seq &sequence, EnergyModel &energy_model, const InsideMode mode, const int beam_size = 100,
              bool allow_sharp_turn = false)
        : seq(Seq::encode(sequence, nuc_encoding_scheme)), energy_model(&energy_model), mode(mode),
          beam_size(beam_size), allow_sharp_turn(allow_sharp_turn) {

        bestH = new std::unordered_map<int, State>[seq.size()];
        bestP = new std::unordered_map<int, State>[seq.size()];
        bestM = new std::unordered_map<int, State>[seq.size()];
        bestM2 = new std::unordered_map<int, State>[seq.size()];
        bestMulti = new std::unordered_map<int, State>[seq.size()];
        // bestC = new State[seq.size()];

        for (int nuc = 1; nuc < 5; ++nuc) {
            prev_pair[nuc].resize(seq.size(), -1);
            next_pair[nuc].resize(seq.size(), seq.size());

            int prev = -1;
            int next = seq.size();

            for (int j = 0; j < seq.size(); ++j) {
                prev_pair[nuc][j] = prev;
                if (check_valid_pair(nuc, seq[j])) {
                    prev = j;
                }
            }
            for (int j = seq.size() - 1; j >= 0; --j) {
                next_pair[nuc][j] = next;
                if (check_valid_pair(nuc, seq[j])) {
                    next = j;
                }
            }
        }

        if (seq.size() > 0)
            bestC[0].alpha = 0.0;
        if (seq.size() > 1)
            bestC[1].alpha = 0.0;

        // compute inside scores
        for (int j = 0; j < seq.size(); j++) {
            std::cout << j << std::endl;

            // beam of H
            beam_prune(bestH[j]);
            beamstep_H(j, next_pair);
            if (j == 0)
                continue;

            // beam of Multi
            beam_prune(bestMulti[j]);
            beamstep_Multi(j, next_pair);

            // beam of P
            beam_prune(bestP[j]);
            beamstep_P(j, next_pair);

            // beam of M2
            beam_prune(bestM2[j]);
            beamstep_M2(j, next_pair);

            // beam of M
            beam_prune(bestM[j]);
            beamstep_M(j);

            beamstep_C(j); // beam of C
        }

        // debug_states();
    }

    void lazy_outside();

    void print_mfe_structure() {
        std::string structure(seq.size(), '.');
        mfe_backtrack(0, seq.size() - 1, C, structure);
        std::cout << structure << std::endl;
        std::cout << "MFE: " << bestC[seq.size() - 1].alpha / -100.0 << std::endl;
    }

    double get_ensemble_energy() { return -kT * (bestC[seq.size() - 1].alpha) / 100.0; } // -kT log(Q(x))

    void print_alpha_beta() {
        std::cout << "Alpha(n - 1): " << bestC[seq.size() - 1].alpha << std::endl;
        std::cout << "Beta(0): " << bestC[-1].beta << std::endl;
        std::cout << "Beta(-1): " << bestC[-1].beta << std::endl;
        std::cout << "Num Alpha Updates: " << num_alpha_updates << std::endl;
        std::cout << "Num Beta Updates: " << num_beta_updates << std::endl;
    }

    void debug_states() {

        // for (int j = 0; j < seq.size(); ++j) {
        //     for (auto item : bestH[j]) {
        //         int i = item.first;
        //         State &state = item.second;
        //         printf("H[%d][%d]: %.5f\t%.5f\n", i, j, state.alpha, state.beta);
        //     }
        // }
        printf("\n");
        for (int j = 0; j < seq.size(); ++j) {
            for (auto item : bestP[j]) {
                int i = item.first;
                State &state = item.second;
                printf("P[%d][%d]: %.5f\t%.5f\n", i, j, state.alpha, state.beta);
            }
        }
        printf("\n");
        for (int j = -1; j < (int)seq.size(); ++j) {
            printf("C[%d]: %.5f\t%.5f\n", j, bestC[j].alpha, bestC[j].beta);
        }
        printf("\n");
        for (int j = 0; j < seq.size(); ++j) {
            for (auto item : bestM[j]) {
                int i = item.first;
                State &state = item.second;
                printf("M[%d][%d]: %.5f\n", i, j, state.alpha);
            }
        }
    }

  private:
    double beam_prune(std::unordered_map<int, State> &beamstep);
    void update_score(State &state, const int new_score, const double prev_score);

    void beamstep_H(const int j, const std::vector<int> *next_pair);
    void beamstep_Multi(const int j, const std::vector<int> *next_pair);
    void beamstep_P(const int j, const std::vector<int> *next_pair);
    void beamstep_M2(const int j, const std::vector<int> *next_pair);
    void beamstep_M(const int j);
    void beamstep_C(const int j);

    std::pair<int, int> backward_update(const int i, const int j, State &state, const StateType type,
                                        const float edge_threshold);
    TraceInfo get_incoming_hedges_C(const int j, std::vector<HEdge> *incoming_hedges);
    TraceInfo get_incoming_hedges_M(const int i, const int j, std::vector<HEdge> *incoming_hedges);
    TraceInfo get_incoming_hedges_M2(const int i, const int j, std::vector<HEdge> *incoming_hedges);
    TraceInfo get_incoming_hedges_P(const int i, const int j, std::vector<HEdge> *incoming_hedges);
    TraceInfo get_incoming_hedges_Multi(const int i, const int j, std::vector<HEdge> *incoming_hedges);
    // void get_incoming_hedges_H(const int j, std::vector<HEdge> &incoming_hedges);

    void hedges_trace_helper(std::vector<HEdge> *incoming_hedges, HEdge &best_hedge, HEdge &new_hedge,
                             TraceInfo &best_trace, TraceInfo &new_trace) {
        if (incoming_hedges) {
            incoming_hedges->push_back(new_hedge);
        }
        const float best_value = best_hedge.weight + (best_hedge.left ? best_hedge.left->alpha : 0) +
                                 (best_hedge.right ? best_hedge.right->alpha : 0);
        const float new_value =
            new_hedge.weight + new_hedge.left->alpha + (new_hedge.right ? new_hedge.right->alpha : 0);
        if (new_value >= best_value) {
            best_hedge = new_hedge;
            best_trace = new_trace;
        }
    }

    void mfe_backtrack(int i, int j, StateType type, std::string &structure) {
        if (i >= j)
            return;
        
        // std::cout << i << " " << j << " " << type << std::endl;

        TraceInfo trace;
        switch (type) {
        case H:
            return;
        case Multi:
            trace = get_incoming_hedges_Multi(i, j, nullptr);
            break;
        case P:
            structure[i] = '(';
            structure[j] = ')';
            trace = get_incoming_hedges_P(i, j, nullptr);
            break;
        case M2:
            trace = get_incoming_hedges_M2(i, j, nullptr);
            break;
        case M:
            trace = get_incoming_hedges_M(i, j, nullptr);
            break;
        case C:
            trace = get_incoming_hedges_C(j, nullptr);
            break;
        }

        mfe_backtrack(trace.i, trace.t, trace.type_left, structure);
        if (trace.type_left != trace.type_right) {
            mfe_backtrack(trace.t + 1, trace.j, trace.type_right, structure);
        }
    }
};

#endif // PARTITION_HPP