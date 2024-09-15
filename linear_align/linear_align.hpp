#include <unordered_map>
#include <vector>

#include "./../sequence/seq.hpp"
#include "./../utility/constants.hpp"
#include "./../utility/log_math.hpp"
#include "./../utility/utility.hpp"
#include <cmath>    // not really needed if working with log scale
#include <iostream> // NOTE: for debugging, remove later

#include "./../shared.hpp"

enum HStateType {
    ALN,  // 0
    INS1, // 1
    INS2, // 2
};

struct HState {
    float alpha;
    float beta;
    HStateType pre;

    HState() : alpha(VALUE_MIN), beta(VALUE_MIN) {};
};

struct AlnEdge {
    float weight;
    HState *prev;

    void set(HState *prev, const float weight = 0) {
        this->prev = prev;
        this->weight = weight;
    }
};

class LinearAlign {
    inline const static std::vector<HStateType> hstate_types = {ALN, INS1, INS2};

  private:
    void edges_helper(std::vector<AlnEdge> *incoming_edges, const AlnEdge &new_edge);

    std::pair<int, int> backward_update(const int i, const int j, const HState &state, const HStateType type,
                                        const float edge_threshold);
    void get_incoming_edges(const int i, const int j, const HState &state, const HStateType type,
                            std::vector<AlnEdge> *incoming_edges);

  public:
    const Seq *sequence1; // pointer to actual Seq object
    const Seq *sequence2;
    const std::vector<int> *seq1; // pointer to encoded sequence in Seq object
    const std::vector<int> *seq2;
    int seq_len_sum;
    const double (*trans_probs_)[3][3];
    const double (*emit_probs_)[27][3];

    std::unordered_map<std::pair<int, int>, HState, PairHash> *bestALN;
    std::unordered_map<std::pair<int, int>, HState, PairHash> *bestINS1;
    std::unordered_map<std::pair<int, int>, HState, PairHash> *bestINS2;
    std::unordered_map<std::pair<int, int>, float, PairHash> coinc_prob;

    ProbAccm *pm1, *pm2;
    float alpha1 = 1.0, alpha2 = 0.8, alpha3 = 0.5;

    int beam_size;
    bool verbose_output;

    float inside_execution_time = 0.0;
    float outside_execution_time = 0.0;

    LinearAlign(const Seq *sequence1, const Seq *sequence2, const double (*transition_probs)[3][3] = &trans_probs,
                const double (*emission_probs)[27][3] = &emit_probs, const int beam_size = 100,
                const bool verbose_output = true)
        : sequence1(sequence1), sequence2(sequence2), seq1(&sequence1->enc_seq), seq2(&sequence2->enc_seq),
          seq_len_sum(seq1->size() + seq2->size()), trans_probs_(transition_probs), emit_probs_(emission_probs),
          beam_size(beam_size), verbose_output(verbose_output) {

        bestALN = new std::unordered_map<std::pair<int, int>, HState, PairHash>[seq_len_sum + 3];
        bestINS1 = new std::unordered_map<std::pair<int, int>, HState, PairHash>[seq_len_sum + 1];
        bestINS2 = new std::unordered_map<std::pair<int, int>, HState, PairHash>[seq_len_sum + 1];

        bestALN[0][{0, 0}].alpha = 1.0;
        bestALN[seq_len_sum + 2][{seq1->size() + 1, seq2->size() + 1}].beta = 1.0;
    }

    void compute_inside();
    void compute_outside();

    void update_state_alpha(HState &state, const double new_score, const HStateType h, const HStateType pre,
                            const bool best_only);
    void update_state_beta(HState &state, const double new_score);
    void update_coinc_prob(const std::pair<int, int> &pos, const double prob);

    double get_trans_emit_prob(const int i, const int j, const HStateType h, const HStateType h_prev);
    float get_match_score(const int i, const int j);
    void set_prob_accm(ProbAccm &prob_accm1, ProbAccm &prob_accm2);
    float get_bpp(int i, int j);

    // legacy methods
    void run_backward_phase(bool compute_coinc_probs);
    void traceback();

    inline std::unordered_map<std::pair<int, int>, HState, PairHash> *get_beam(HStateType type) {
        switch (type) {
        case HStateType::ALN:
            return bestALN;
        case HStateType::INS1:
            return bestINS1;
        case HStateType::INS2:
            return bestINS2;
        }
    }

    inline void print_alpha_beta() {
        std::cout << "\n(ALN, n1 + 1, n2 + 1): " << bestALN[seq_len_sum + 2][{seq1->size() + 1, seq2->size() + 1}].alpha << std::endl;
        std::cout << "(ALN, 0, 0): " << bestALN[0][{0, 0}].beta << std::endl;
    }

    inline void print_seqs() {
        std::cout << "Seq1: ";
        for (auto &nuc : *seq1) {
            std::cout << nuc << " ";
        }
        std::cout << std::endl;

        std::cout << "Seq2: ";
        for (auto &nuc : *seq2) {
            std::cout << nuc << " ";
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
};
