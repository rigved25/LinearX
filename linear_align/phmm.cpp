#include "./phmm.hpp"

Phmm::Phmm(double new_emission_probs[N_OUTPUTS][N_STATES], double new_trans_probs[N_STATES][N_STATES]) {
    this->alloc_init_params();

    // copy transition matrix
    for (int cnt1 = 0; cnt1 < N_STATES; cnt1++) {
        for (int cnt2 = 0; cnt2 < N_STATES; cnt2++) {
            this->trans_probs[cnt1][cnt2] = LOG(new_trans_probs[cnt1][cnt2]);
        } // cnt2 loop
    } // cnt1 loop

    // copy emission probabilities
    for (int cnt1 = 0; cnt1 < N_OUTPUTS; cnt1++) {
        for (int cnt2 = 0; cnt2 < N_STATES; cnt2++) {
            this->emission_probs[cnt1][cnt2] = LOG(new_emission_probs[cnt1][cnt2]);
        } // cnt2 loop
    } // cnt1 loop
}

Phmm::Phmm(const char *phmm_pars_file) {
    this->alloc_init_params();

    // read the parameters file
    FILE *fam_par_file = Utility::open_f(phmm_pars_file, "r");

    if (fam_par_file == NULL) {
        printf("Cannot find phmm parameters file, exiting @ %s(%d).\n", __FILE__, __LINE__);
        exit(0);
    }

    // load all parameters from file
    for (int cnt = 0; cnt < N_BINZ * (N_STATES + N_OUTPUTS) * N_STATES; cnt++) {
        fscanf(fam_par_file, "%lf", &fam_hmm_pars[cnt]);
    }

    for (int cnt = 0; cnt < N_BINZ; cnt++) {
        fscanf(fam_par_file, "%lf", &fam_thresholds[cnt]);
    }

    fclose(fam_par_file);
}

void Phmm::alloc_init_params() {
    // Copy transition matrix.
    this->trans_probs = (double **)malloc(sizeof(double *) * (N_STATES + 2));
    for (int cnt1 = 0; cnt1 < N_STATES; cnt1++) {
        this->trans_probs[cnt1] = (double *)malloc(sizeof(double) * (N_STATES + 2));
        for (int cnt2 = 0; cnt2 < N_STATES; cnt2++) {
            trans_probs[cnt1][cnt2] = LOG(0.0f);
        } // cnt2 loop
    } // cnt1 loop

    // Copy emission probabilities.
    this->emission_probs = (double **)malloc(sizeof(double *) * (N_OUTPUTS + 2));
    for (int cnt1 = 0; cnt1 < N_OUTPUTS; cnt1++) {
        this->emission_probs[cnt1] = (double *)malloc(sizeof(double) * (N_STATES + 2));
        for (int cnt2 = 0; cnt2 < N_STATES; cnt2++) {
            emission_probs[cnt1][cnt2] = LOG(0.0f);
        } // cnt2 loop
    } // cnt1 loop

    this->fam_hmm_pars = (double *)malloc(sizeof(double) * (N_BINZ * (N_STATES + N_OUTPUTS) * N_STATES + 2));
    this->fam_thresholds = (double *)malloc(sizeof(double) * (N_BINZ + 2));
}

void Phmm::free_params() {
    // free transition matrix
    for (int cnt1 = 0; cnt1 < N_STATES; cnt1++) {
        free(this->trans_probs[cnt1]);
    } // cnt1 loop
    free(this->trans_probs);

    // free emission probabilities
    for (int cnt1 = 0; cnt1 < N_OUTPUTS; cnt1++) {
        free(this->emission_probs[cnt1]);
    } // cnt1 loop
    free(this->emission_probs);

    free(this->fam_hmm_pars);

    free(this->fam_thresholds);
}

Phmm::~Phmm() { this->free_params(); }

void Phmm::print_parameters() { 
    // dump emission probabilities
    for (int cnt1 = 0; cnt1 < N_OUTPUTS; cnt1++) {
        for (int cnt2 = 0; cnt2 < N_STATES; cnt2++) {
            printf("%.3f ", EXP(emission_probs[cnt1][cnt2]));
        }
        printf("\n");
    }

    // dump transition probabilities
    printf("\n");
    for (int cnt1 = 0; cnt1 < N_STATES; cnt1++) {
        for (int cnt2 = 0; cnt2 < N_STATES; cnt2++) {
            printf("%.3f ", EXP(trans_probs[cnt1][cnt2]));
        }
        printf("\n");
    }
}

void Phmm::set_parameters_by_sim(double similarity) {
    this->similarity = similarity;
    int fam_par_set_index = get_bin_index(similarity, N_BINZ);

    // load emission probabilities
    int start_linear_index = (N_STATES + N_OUTPUTS) * N_STATES * get_bin_index(similarity, N_BINZ);
    double *par_ptr = fam_hmm_pars + start_linear_index;

    for (int cnt1 = 0; cnt1 < N_OUTPUTS; cnt1++) {
        for (int cnt2 = 0; cnt2 < N_STATES; cnt2++) {            
            emission_probs[cnt1][cnt2] = LOG(par_ptr[cnt1 * N_STATES + cnt2]);
        }
    }

    start_linear_index = (N_STATES + N_OUTPUTS) * N_STATES * fam_par_set_index + N_STATES * N_OUTPUTS;
    par_ptr = fam_hmm_pars + start_linear_index;

    // load transition probabilities
    for (int cnt1 = 0; cnt1 < N_STATES; cnt1++) {
        for (int cnt2 = 0; cnt2 < N_STATES; cnt2++) {
            trans_probs[cnt1][cnt2] = LOG(par_ptr[cnt1 * N_STATES + cnt2]);
        }
    }
}

// get index of bin of parameters for a sequence alignment
int Phmm::get_bin_index(double similarity, int n_bins) {
    if (similarity == 1.0) {
        return (n_bins - 1);
    } else {
        return ((int)(n_bins * similarity));
    }
}

double Phmm::get_fam_threshold() {
    int bin_index = get_bin_index(this->similarity, N_BINZ);
    return (fam_thresholds[bin_index]);
}

double Phmm::get_trans_prob(int prev, int next) { return (this->trans_probs[prev][next]); }

double Phmm::get_emit_prob(int sym_index, int state) { return (this->emission_probs[sym_index][state]); }
