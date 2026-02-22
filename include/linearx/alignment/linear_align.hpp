// linearx/alignment/linear_align.hpp
#pragma once
#include <linearx/alignment/phmm.hpp>
#include <linearx/alignment/utility.hpp>
#include <linearx/sequence/multi_sequence.hpp>
#include <linearx/utility.hpp>
#include <unordered_set>

template <typename T>
class LinearAlignmentInterface {
    inline const static std::vector<HStateType> hstate_types = {INS1, INS2, ALN};

   private:
    linearx::utils::ProbAccm *pm1 = nullptr, *pm2 = nullptr;
    std::vector<std::pair<value_type, std::pair<int, int>>> beam_scores;
    std::vector<AlnEdge> incoming_edges;
    std::vector<AlnEdge*> saved_edges;
    AlnEdge best_edge;
    HStateType best_trace;

    void update_best_trace(const AlnEdge& new_edge, const HStateType& new_trace);

    std::pair<unsigned long, unsigned long> backward_update(const int i, const int j, const HState& state,
                                                            const HStateType type, const value_type edge_threshold);
    std::pair<unsigned long, unsigned long> backward_update_BEST(const int i, const int j, const HState& state,
                                                                const HStateType type, const value_type edge_threshold);

    template <Mode mode>
    void get_incoming_edges(const int i, const int j, const HStateType type);

    void run_normal_outside(const bool verbose_output);
    void run_normal_outside_best(const bool verbose_output);

   protected:
    std::vector<std::unordered_map<std::pair<int, int>, HState, linearx::utils::PairHash>> bestALN;
    std::vector<std::unordered_map<std::pair<int, int>, HState, linearx::utils::PairHash>> bestINS1;
    std::vector<std::unordered_map<std::pair<int, int>, HState, linearx::utils::PairHash>> bestINS2;
    unsigned run_beam_size_;

    void set_prob_accm(linearx::utils::ProbAccm& prob_accm1, linearx::utils::ProbAccm& prob_accm2);

    value_type get_match_score(const unsigned i, const unsigned j) const {
        if (!pm1 || !pm2 || i >= seq1.length() || j >= seq2.length()) {
            return linearx::math::LOG_ONE;
        }
        const value_type t1 = sqrt(pm1->upstrm[i] * pm2->upstrm[j]);
        const value_type t2 = sqrt(pm1->dwnstrm[i] * pm2->dwnstrm[j]);
        const value_type t3 = sqrt(std::max(1 - pm1->upstrm[i] - pm1->dwnstrm[i], (value_type)0.0) *
                                   std::max(1 - pm2->upstrm[j] - pm2->dwnstrm[j], (value_type)0.0));
        value_type output = ((t1 + t2) * alpha1) + (t3 * alpha2) + (alpha3);
        output = LOG(output);
        return output;
    }

    inline HState* get_saved_state(const HStateType type, const int i, const int j) {
        return static_cast<T*>(this)->get_saved_state(type, i, j);
    }

    HState* check_state(const HStateType type, const int i, const int j) {
        return static_cast<T*>(this)->check_state(type, i, j);
    }

    HState* check_state_best(const HStateType type, const int i, const int j) {
        return static_cast<T*>(this)->check_state_best(type, i, j);
    }

    HState* check_state_AStar(const HStateType type, const int i, const int j) {
        return static_cast<T*>(this)->check_state_AStar(type, i, j);
    }

    unsigned beam_prune(std::unordered_map<std::pair<int, int>, HState, linearx::utils::PairHash>& beamstep) {
        if (run_beam_size_ == 0 || beamstep.size() <= run_beam_size_) {
            return 0;
        }
        unsigned num_pruned = 0;
        beam_scores.clear();
        beam_scores.reserve(beamstep.size());
        for (auto& item : beamstep) {
            const std::pair<int, int> aij = item.first;
            HState& cand = item.second;
            beam_scores.emplace_back(cand.alpha, aij);
        }
        value_type threshold =
            linearx::utils::quickselect(beam_scores, 0, beam_scores.size() - 1, beam_scores.size() - run_beam_size_);
        for (auto& p : beam_scores) {
            if (p.first < threshold) {
                beamstep.erase(p.second);
                num_pruned++;
            }
        }
        return num_pruned;
    }

   public:
    Phmm* phmm = nullptr;
    Sequence& seq1;
    Sequence& seq2;
    int seq_len_sum;
    std::vector<std::unordered_map<int, value_type>> cp;  // coincidence probabilities
    std::vector<std::vector<int>> prob_rev_idx;
    float alpha1, alpha2, alpha3;
    AlignmentLog log;

    LinearAlignmentInterface(Sequence& seq1, Sequence& seq2, float alpha1 = 1.0, float alpha2 = 0.8,
                             float alpha3 = 0.5);

    void use_prob_set1();
    void use_prob_set2(const float similarity);

    void reset_beams(const unsigned beam_size);

    std::vector<std::unordered_map<std::pair<int, int>, HState, linearx::utils::PairHash>>& get_beams(
        HStateType type) noexcept {
        switch (type) {
            case ALN:
                return bestALN;
            case INS1:
                return bestINS1;
            case INS2:
                return bestINS2;
            default:
                return bestALN;  // unreachable
        }
    }

    HState* get_state(const HStateType type, const int i, const int j, const bool create = false) noexcept {
        auto& beam = get_beams(type)[i + j];
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
    void compute_inside(const unsigned beam_size = 100, bool verbose_output = true);
    void compute_inside_Astar(const bool use_lazy_outside, const unsigned beam_size = 0, bool verbose_output = true);
    void compute_inside_Astar_lazy(const bool use_lazy_outside, const unsigned beam_size = 0, bool verbose_output = true);
    void compute_outside(const bool use_lazy_outside,
                         const value_type deviation_threshold = linearx::constants::limits::DEVIATION_THRESHOLD,
                         const bool verbose_output = true);
    void compute_outside_BEST(const bool use_lazy_outside,
                             const value_type deviation_threshold = linearx::constants::limits::DEVIATION_THRESHOLD,
                             const bool verbose_output = true);

    void compute_coincidence_probabilities(const bool verbose_output = true);
    void compute_posterior(std::vector<std::vector<std::unordered_map<int, value_type>*>>& posterior, const bool verbose_output = true);
    void dump_coinc_probs(const std::string& out_dir) const;
    value_type get_cp(const int i, const int j) const {
        const auto& cp_i = cp[i];
        const auto it = cp_i.find(j);
        return it == cp_i.end() ? 0.0 : it->second;
    }

    void print_alpha_beta() const;
    void print_seqs() const;
    void print_beams() const;
    void dump_beams(const std::string& out_dir, const std::string& prefix = "") const;
};

class LinearAlignment final : public LinearAlignmentInterface<LinearAlignment> {
   public:
    LinearAlignment(Sequence& seq1, Sequence& seq2, float alpha1 = 1.0, float alpha2 = 0.8, float alpha3 = 0.5)
        : LinearAlignmentInterface<LinearAlignment>(seq1, seq2, alpha1, alpha2, alpha3) {}

    HState* check_state(const HStateType type, const int i, const int j) {
        return LinearAlignmentInterface<LinearAlignment>::get_state(type, i, j, true);
    }

    HState* check_state_best(const HStateType type, const int i, const int j) {
        return LinearAlignmentInterface<LinearAlignment>::get_state(type, i, j, true);
    }

    HState* check_state_AStar(const HStateType type, const int i, const int j) {
        return LinearAlignmentInterface<LinearAlignment>::get_state(type, i, j, true);
    }

    inline HState* get_saved_state(const HStateType type, const int i, const int j) {
        (void)type; (void)i; (void)j;
        return nullptr;  // LinearAlignment has no saved states
    }
};
