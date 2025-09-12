// linearx/alignment/linear_align.hpp
#pragma once
#include <linearx/alignment/phmm.hpp>
#include <linearx/alignment/utility.hpp>
#include <linearx/sequence/multi_sequence.hpp>
#include <linearx/utility.hpp>
#include <unordered_set>

class LinearAlignment {
    inline const static std::vector<HStateType> hstate_types = {INS1, INS2, ALN};

   private:
    std::vector<std::pair<value_type, std::pair<int, int>>> beam_scores;
    std::vector<AlnEdge> incoming_edges;
    std::vector<AlnEdge *> saved_edges;
    AlnEdge best_edge;
    HStateType best_trace;
    value_type _last_inside_exec_time = 1.0;

    void update_best_trace(const AlnEdge &new_edge, const HStateType &new_trace);

    std::pair<unsigned long, unsigned long> backward_update(const int i, const int j, const HState &state,
                                                            const HStateType type, const value_type edge_threshold);

    template <Mode mode>
    void get_incoming_edges(const int i, const int j, const HStateType type);

    AlignmentLog run_normal_outside(const bool verbose_output);

   protected:
    std::vector<std::unordered_map<std::pair<int, int>, HState, linearx::utils::PairHash>> bestALN;
    std::vector<std::unordered_map<std::pair<int, int>, HState, linearx::utils::PairHash>> bestINS1;
    std::vector<std::unordered_map<std::pair<int, int>, HState, linearx::utils::PairHash>> bestINS2;
    unsigned run_beam_size_;

    inline virtual HState *check_state(const HStateType type, const int i, const int j) {
        return get_state(type, i, j, true);
    }

    inline virtual unsigned beam_prune(
        std::unordered_map<std::pair<int, int>, HState, linearx::utils::PairHash> &beamstep) {
        if (run_beam_size_ == 0 || beamstep.size() <= run_beam_size_) {
            return 0;
        }
        unsigned num_pruned = 0;
        beam_scores.clear();
        beam_scores.reserve(beamstep.size());
        for (auto &item : beamstep) {
            const std::pair<int, int> aij = item.first;
            HState &cand = item.second;
            beam_scores.emplace_back(cand.alpha, aij);
        }
        value_type threshold =
            linearx::utils::quickselect(beam_scores, 0, beam_scores.size() - 1, beam_scores.size() - run_beam_size_);
        for (auto &p : beam_scores) {
            if (p.first < threshold) {
                beamstep.erase(p.second);
                num_pruned++;
            }
        }
        return num_pruned;
    }

   public:
    friend struct AlignBeam;

    Phmm *phmm = nullptr;
    Sequence &seq1;
    Sequence &seq2;
    int seq_len_sum;
    std::function<value_type(const unsigned, const unsigned)> get_match_score;

    std::vector<std::unordered_map<int, value_type>> coinc_prob;
    std::vector<std::vector<int>> prob_rev_idx;

    float alpha1, alpha2, alpha3;

    LinearAlignment(
        Sequence &seq1, Sequence &seq2,
        std::function<value_type(int, int)> get_match_score = [](const unsigned i,
                                                                 const unsigned j) { return linearx::math::LOG_ONE; },
        float alpha1 = 1.0, float alpha2 = 0.8, float alpha3 = 0.5);

    void use_prob_set1();
    void use_prob_set2(const float similarity);

    void reset_beams(const unsigned beam_size);
    inline std::vector<std::unordered_map<std::pair<int, int>, HState, linearx::utils::PairHash>> &get_beams(
        HStateType type) noexcept {
        switch (type) {
            case ALN:
                return bestALN;
            case INS1:
                return bestINS1;
            case INS2:
                return bestINS2;
        }
        throw std::invalid_argument("Invalid HStateType");
    };
    inline HState *get_state(const HStateType type, const int i, const int j, const bool create = false) noexcept {
        auto &beam = get_beams(type)[i + j];
        const std::pair<int, int> key = {i, j};
        const auto it = beam.find(key);
        if (it != beam.end()) {
            return &it->second;
        } else if (create) {
            return &beam[key];
        } else {
            return nullptr;
        }
    }

    value_type get_trans_emit_prob(const int i, const int j, const HStateType h, const HStateType h_prev) const;
    MultiSeq get_alignment();  // [TODO] Can make it more efficient, by avoiding copy, or using move semantics

    template <Mode mode>
    AlignmentLog compute_inside(const unsigned beam_size = 100, bool verbose_output = true);
    AlignmentLog compute_outside(const bool use_lazy_outside,
                                 const value_type deviation_threshold = linearx::constants::limits::DEVIATION_THRESHOLD,
                                 const bool verbose_output = true);

    void compute_coincidence_probabilities(const bool verbose_output = true);
    void dump_coinc_probs(const std::string &out_dir) const;
    inline value_type get_bpp(const int i, const int j) const {
        const auto &coinc_prob_i = coinc_prob[i];
        const auto it = coinc_prob_i.find(j);
        return it == coinc_prob_i.end() ? 0.0 : it->second;
    }

    void print_alpha_beta() const;
    void print_seqs() const;
    void print_beams() const;

    // legacy methods -----------------------------------------
    // void run_backward_phase(bool verbose_output = false);
    // MultiSeq old_traceback();
};
