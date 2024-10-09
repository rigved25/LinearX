#ifndef PHMM_HPP
#define PHMM_HPP

#include <iostream>
#include <unistd.h>

#include "./../utility/utility.hpp"
#include "./utility.hpp"

#define N_STATES (3)
#define N_OUTPUTS (27)
#define N_BINZ (10)

class Phmm {
  public:
    inline static double EMIT_PROBS[27][3] = {
        {0.000000, 0.000000, 0.134009}, // AA
        {0.000000, 0.000000, 0.027164}, // AC
        {0.000000, 0.000000, 0.049659}, // AG
        {0.000000, 0.000000, 0.028825}, // AU
        {0.211509, 0.000000, 0.000000}, // A.
        {0.000000, 0.000000, 0.027164}, // CA
        {0.000000, 0.000000, 0.140242}, // CC
        {0.000000, 0.000000, 0.037862}, // CG
        {0.000000, 0.000000, 0.047735}, // CU
        {0.257349, 0.000000, 0.000000}, // C.
        {0.000000, 0.000000, 0.049659}, // GA
        {0.000000, 0.000000, 0.037862}, // GC
        {0.000000, 0.000000, 0.178863}, // GG
        {0.000000, 0.000000, 0.032351}, // GU
        {0.271398, 0.000000, 0.000000}, // G.
        {0.000000, 0.000000, 0.028825}, // UA
        {0.000000, 0.000000, 0.047735}, // UC
        {0.000000, 0.000000, 0.032351}, // UG
        {0.000000, 0.000000, 0.099694}, // UU
        {0.259744, 0.000000, 0.000000}, // U.
        {0.000000, 0.211509, 0.000000}, // .A
        {0.000000, 0.257349, 0.000000}, // .C
        {0.000000, 0.271398, 0.000000}, // .G
        {0.000000, 0.259744, 0.000000}, // .U
        {0.000000, 0.000000, 0.000000}, // ..
        {0.000000, 0.000000, 1.000000}, // START
        {0.000000, 0.000000, 1.000000}  // END
    };

    inline static double TRANS_PROBS[3][3] = {
        {0.666439, 0.041319, 0.292242}, // INS1
        {0.041319, 0.666439, 0.292242}, // INS2
        {0.022666, 0.022666, 0.954668}  // ALIGN
    };

  private:
    double **emission_probs;
    double **trans_probs;

    double *fam_hmm_pars;
    double *fam_thresholds;
  public:
    // Replace the emission and transition probabilities.
    Phmm(double new_emit_probs[N_OUTPUTS][N_STATES], double new_trans_probs[N_STATES][N_STATES]);
    Phmm(const char *phmm_pars_file);
    ~Phmm();

    float similarity = -1.0f;

    void set_parameters_by_sim(double similarity);
    int get_bin_index(double similarity, int n_bins);
    double get_fam_threshold();

    double get_emit_prob(int sym_index, int state);
    double get_trans_prob(int prev, int next);

    void alloc_init_params();
    void free_params();

    void print_parameters();
};

#endif // PHMM_HPP
