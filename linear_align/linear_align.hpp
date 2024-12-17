#include <fstream>   // for file I/O
#include <iostream>  // NOTE: for debugging, remove later
#include <unordered_map>
#include <vector>

#include "./../sequence/multi_seq.hpp"
#include "./../sequence/seq.hpp"
#include "./../shared.hpp"
#include "./phmm.hpp"
#include "./utility.hpp"

class LinearAlign;

struct AlignBeam {
    int seq1_size;
    int seq2_size;

    double total_alpha;

    std::unordered_map<std::pair<int, int>, HState, PairHash> *bestALN = nullptr;
    std::unordered_map<std::pair<int, int>, HState, PairHash> *bestINS1 = nullptr;
    std::unordered_map<std::pair<int, int>, HState, PairHash> *bestINS2 = nullptr;

    AlignBeam(int seq1_size, int seq2_size) : seq1_size(seq1_size), seq2_size(seq2_size) {}
    ~AlignBeam();

    void free();

    void save(std::unordered_map<std::pair<int, int>, HState, PairHash> *bestALN,
              std::unordered_map<std::pair<int, int>, HState, PairHash> *bestINS1,
              std::unordered_map<std::pair<int, int>, HState, PairHash> *bestINS2);

    void save(LinearAlign &la);
};

class LinearAlign {
    inline const static std::vector<HStateType> hstate_types = {INS1, INS2, ALN};

   private:
    ProbAccm *pm1 = nullptr, *pm2 = nullptr;

    std::unordered_map<std::pair<int, int>, HState, PairHash> *bestALN = nullptr;
    std::unordered_map<std::pair<int, int>, HState, PairHash> *bestINS1 = nullptr;
    std::unordered_map<std::pair<int, int>, HState, PairHash> *bestINS2 = nullptr;

    inline void update_state_alpha(HState &state, const double new_score, const HStateType h, const HStateType pre,
                                   const bool best_only);
    inline void update_state_beta(HState &state, const double new_score);

    inline void edges_trace_update_helper(std::vector<AlnEdge> *incoming_edges, const AlnEdge &new_edge,
                                          AlnEdge &best_edge, const HStateType &new_trace, HStateType &best_trace);

    std::pair<int, int> backward_update(const int i, const int j, const HState &state, const HStateType type,
                                        const double edge_threshold);
    HStateType get_incoming_edges(const int i, const int j, const HStateType type, std::vector<AlnEdge> *incoming_edges,
                                  bool use_match_score);

    void run_normal_outside(bool verbose_output);

   protected:
    virtual double beam_prune(std::unordered_map<std::pair<int, int>, HState, PairHash> &beamstep, HStateType h,
                              int beam_size);

    virtual inline bool check_state(const int i, const int j, const HStateType h) { return true; }

   public:
    friend struct AlignBeam;

    Phmm *phmm = nullptr;
    Seq *sequence1 = nullptr;  // pointer to actual Seq object
    Seq *sequence2 = nullptr;
    std::vector<int> *seq1 = nullptr;  // pointer to encoded sequence in Seq object
    std::vector<int> *seq2 = nullptr;
    int seq_len_sum;

    std::unordered_map<int, double> *coinc_prob = nullptr;
    std::vector<int> *prob_rev_idx = nullptr;

    double alpha1, alpha2, alpha3;
    double inside_execution_time = 0.0;
    double outside_execution_time = 0.0;

    LinearAlign(Seq *sequence1, Seq *sequence2, bool verbose = false, double alpha1 = 1.0, double alpha2 = 0.8,
                double alpha3 = 0.5)
        : sequence1(sequence1),
          sequence2(sequence2),
          seq1(&sequence1->enc_seq),
          seq2(&sequence2->enc_seq),
          seq_len_sum(seq1->size() + seq2->size()),
          alpha1(alpha1),
          alpha2(alpha2),
          alpha3(alpha3) {
        reset_beams(false);

        for (auto &x : *seq1) {
            if (x == 0) {
                x = (rand() % 4) + 1;
            }
        }

        for (auto &x : *seq2) {
            if (x == 0) {
                x = (rand() % 4) + 1;
            }
        }

        if (verbose) {
            print_seqs();
        }
    }

    ~LinearAlign() {
        // std::cerr << "[INFO] Free Alignment Memory" << std::endl;
        delete[] bestALN;
        delete[] bestINS1;
        delete[] bestINS2;
        delete[] coinc_prob;
        delete[] prob_rev_idx;
        // delete pm1;
        // delete pm2;
        // delete phmm;
        // delete seq1;
        // delete seq2;
        // delete sequence1;
        // delete sequence2;
    }

    void prob_set1() {
        delete phmm;
        phmm = new Phmm(Phmm::EMIT_PROBS, Phmm::TRANS_PROBS);
    }

    void prob_set2(float similarity) {
        delete phmm;
        std::string phmm_pars_fp =
            "/Users/malikap/ltf-project/LinearAlifold/linearx/linear_align/parameters/fam_hmm_pars.dat";
        phmm = new Phmm(phmm_pars_fp.c_str());
        phmm->set_parameters_by_sim(similarity);
    }

    void reset_beams(bool freeMemory = true);

    void compute_inside(bool best_only = false, int beam_size = 100, bool verbose_output = true);
    void compute_outside(bool use_lazy_outside = true, bool verbose_output = true);
    void compute_coincidence_probabilities(bool verbose_output = true);

    MultiSeq get_alignment();
    double get_trans_emit_prob(const int i, const int j, const HStateType h, const HStateType h_prev);
    double get_match_score(const int i, const int j);
    void set_prob_accm(ProbAccm &prob_accm1, ProbAccm &prob_accm2);
    double get_bpp(const int i, const int j) const;

    void dump_coinc_probs(const std::string &filepath, const float threshold = 0.001f) const;

    inline std::unordered_map<std::pair<int, int>, HState, PairHash> *get_beam(HStateType type) {
        switch (type) {
            case HStateType::INS1:
                return bestINS1;
            case HStateType::INS2:
                return bestINS2;
            case HStateType::ALN:
                return bestALN;
        }
    }

    inline void print_alpha_beta() {
        std::cerr << "Alpha(ALN, n1 + 1, n2 + 1): "
                  << bestALN[seq_len_sum + 2][{seq1->size() + 1, seq2->size() + 1}].alpha << std::endl;
        std::cerr << "Beta(ALN, 0, 0): " << bestALN[0][{0, 0}].beta << std::endl << std::endl;
    }

    inline void print_seqs() {
        printf("%s (%d):\n%s", sequence1->id.c_str(), sequence1->k_id, sequence1->sequence.c_str());
        for (auto &nuc : *seq1) {
            std::cout << nuc;
        }
        std::cout << std::endl;

        printf("%s (%d):\n%s", sequence2->id.c_str(), sequence2->k_id, sequence2->sequence.c_str());
        for (auto &nuc : *seq2) {
            std::cout << nuc;
        }
        std::cout << std::endl;
    }

    inline void print_beam() {
        for (int s = 0; s < seq_len_sum + 3; ++s) {
            for (auto &item : bestALN[s]) {
                std::cout << "ALN: (" << item.first.first << ", " << item.first.second << ") : " << item.second.alpha
                          << " " << item.second.beta << std::endl;
            }
        }
        std::cout << "--------------------------------" << std::endl;
        for (int s = 0; s < seq_len_sum + 1; ++s) {
            for (auto &item : bestINS1[s]) {
                std::cout << "INS1: (" << item.first.first << ", " << item.first.second << ") : " << item.second.alpha
                          << " " << item.second.beta << std::endl;
            }
        }
        std::cout << "--------------------------------" << std::endl;
        for (int s = 0; s < seq_len_sum + 1; ++s) {
            for (auto &item : bestINS2[s]) {
                std::cout << "INS2: (" << item.first.first << ", " << item.first.second << ") : " << item.second.alpha
                          << " " << item.second.beta << std::endl;
            }
        }
    }

    // legacy methods -----------------------------------------
    // void run_backward_phase(bool verbose_output = false);
    // MultiSeq old_traceback();
};
