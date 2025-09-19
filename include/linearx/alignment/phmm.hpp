// linearx/linear_align/phmm.hpp
#pragma once

#include <unistd.h>

#include <iostream>
#include <linearx/config.hpp>

#define N_STATES (3)
#define N_OUTPUTS (27)
#define N_BINZ (10)

class Phmm {
   public:
    inline static value_type EMIT_PROBS[27][3] = {
        {0.000000, 0.000000, 0.134009},  // AA
        {0.000000, 0.000000, 0.027164},  // AC
        {0.000000, 0.000000, 0.049659},  // AG
        {0.000000, 0.000000, 0.028825},  // AU
        {0.211509, 0.000000, 0.000000},  // A.
        {0.000000, 0.000000, 0.027164},  // CA
        {0.000000, 0.000000, 0.140242},  // CC
        {0.000000, 0.000000, 0.037862},  // CG
        {0.000000, 0.000000, 0.047735},  // CU
        {0.257349, 0.000000, 0.000000},  // C.
        {0.000000, 0.000000, 0.049659},  // GA
        {0.000000, 0.000000, 0.037862},  // GC
        {0.000000, 0.000000, 0.178863},  // GG
        {0.000000, 0.000000, 0.032351},  // GU
        {0.271398, 0.000000, 0.000000},  // G.
        {0.000000, 0.000000, 0.028825},  // UA
        {0.000000, 0.000000, 0.047735},  // UC
        {0.000000, 0.000000, 0.032351},  // UG
        {0.000000, 0.000000, 0.099694},  // UU
        {0.259744, 0.000000, 0.000000},  // U.
        {0.000000, 0.211509, 0.000000},  // .A
        {0.000000, 0.257349, 0.000000},  // .C
        {0.000000, 0.271398, 0.000000},  // .G
        {0.000000, 0.259744, 0.000000},  // .U
        {0.000000, 0.000000, 0.000000},  // ..
        {0.000000, 0.000000, 1.000000},  // START
        {0.000000, 0.000000, 1.000000}   // END
    };

    inline static value_type TRANS_PROBS[3][3] = {
        {0.666439, 0.041319, 0.292242},  // INS1
        {0.041319, 0.666439, 0.292242},  // INS2
        {0.022666, 0.022666, 0.954668}   // ALIGN
    };

   private:
    value_type** emission_probs;
    value_type** trans_probs;

    value_type* fam_hmm_pars;
    value_type* fam_thresholds;

   public:
    // Replace the emission and transition probabilities.
    Phmm(value_type new_emit_probs[N_OUTPUTS][N_STATES], value_type new_trans_probs[N_STATES][N_STATES]);
    Phmm(const char* phmm_pars_file);
    ~Phmm();

    float similarity = -1.0f;

    void set_parameters_by_sim(float similarity);
    int get_bin_index(float similarity, int n_bins);
    value_type get_fam_threshold();

    value_type get_emit_prob(int sym_index, int state);
    value_type get_trans_prob(int prev, int next);

    void alloc_init_params();
    void free_params();

    void print_parameters();
};
