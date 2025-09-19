// linearx/partition/linear_partition.hpp
#pragma once
#include <linearx/energy/energy_model.hpp>
#include <linearx/partition/utility.hpp>
#include <linearx/sequence/sequence.hpp>
#include <linearx/sequence/structure.hpp>
#include <linearx/utility.hpp>
#include <unordered_set>

template <typename T>
class LinearPartitionInterface {
    inline const static StateType state_types[6] = {H, Multi, P, M2, M, C};

   private:
    std::vector<int> if_tetraloops;
    std::vector<int> if_hexaloops;
    std::vector<int> if_triloops;
    std::vector<std::pair<value_type, int>> beam_scores;
    std::vector<HEdge> incoming_hedges;
    std::vector<HEdge*> saved_hedges;
    TraceInfo best_trace;
    HEdge best_hedge;
    value_type _last_inside_exec_time = 1.0;

    // methods declared in file backward.cpp
    void update_best_trace(const HEdge& new_hedge, const TraceInfo& new_trace);
    void mfe_backtrack(const int i, const int j, const StateType type, Structure& structure);
    PartitionLog run_regular_outside(const bool verbose_output);

   protected:
    std::vector<std::unordered_map<int, State>> bestH;
    std::vector<std::unordered_map<int, State>> bestP;
    std::vector<std::unordered_map<int, State>> bestM;
    std::vector<std::unordered_map<int, State>> bestM2;
    std::vector<std::unordered_map<int, State>> bestMulti;
    VectorWithNegOneIndex<State> bestC;

    inline State* check_state(const StateType type, const int i, const int j) {
        return static_cast<T*>(this)->check_state(type, i, j);
    }

    inline void update_state_beta(HEdge hedge, State* state, const StateType type, const int i, const int j) {
        return static_cast<T*>(this)->update_state_beta(hedge, state, type, i, j);
    }

    inline virtual unsigned long beam_prune(std::unordered_map<int, State>& beamstep, const unsigned beam_size,
                                            const StateType type, const int j) {
        if (beam_size == 0 || beamstep.size() <= beam_size) {
            return 0;
        }
        unsigned long num_pruned = 0;
        beam_scores.clear();
        beam_scores.reserve(beamstep.size());
        for (auto& item : beamstep) {
            const int i = item.first;
            const State& cand = item.second;
            const int k = i - 1;
            const value_type newalpha = ((k >= 0 ? bestC[k].alpha : 0) + cand.alpha);
            beam_scores.emplace_back(newalpha, i);
        }
        const value_type threshold =
            linearx::utils::quickselect(beam_scores, 0, beam_scores.size() - 1, beam_scores.size() - beam_size);
        for (auto& p : beam_scores) {
            if (p.first < threshold) {
                beamstep.erase(p.second);
                num_pruned++;
            }
        }
        return num_pruned;
    }

   public:
    friend struct PartitionFunctionBeam;

    std::vector<std::unordered_map<int, value_type>> bpp;
    std::array<std::vector<int>, 5> next_pair;
    std::array<std::vector<int>, 5> prev_pair;

    const Sequence& seq;
    const int seq_length;
    const EnergyModel& energy_model;
    const bool allow_sharp_turn;

    LinearPartitionInterface(const Sequence& seq, const EnergyModel& energy_model, const bool allow_sharp_turn = false);
    State& get_viterbi();

    void reset_beams(const unsigned beam_size);
    PartitionLog compute_inside(const Mode mode, const unsigned beam_size = 100, const bool verbose_output = true);
    PartitionLog compute_outside(const bool use_lazy_outside = true,
                                 const value_type deviation_threshold = linearx::constants::limits::DEVIATION_THRESHOLD,
                                 const bool verbose_output = true);

    void compute_bpp_matrix(const unsigned beam_size);
    value_type get_ensemble_energy() const;
    void print_alpha_beta() const;
    Structure get_mfe_structure();

    inline State* get_state(const StateType type, const int i, const int j, const bool create = false) noexcept {
        std::unordered_map<int, State>* beam = nullptr;
        switch (type) {
            case H:
                beam = &bestH[j];
                break;
            case Multi:
                beam = &bestMulti[j];
                break;
            case P:
                beam = &bestP[j];
                break;
            case M2:
                beam = &bestM2[j];
                break;
            case M:
                beam = &bestM[j];
                break;
            case C:
                return &bestC[j];  // special case
            default:
                return nullptr;
        }
        const auto it = beam->find(i);
        if (it != beam->end()) {
            return &it->second;
        } else if (create) {
            return &(*beam)[i];
        } else {
            return nullptr;
        }
    }

    inline value_type get_bpp(const int i, const int j) const {
        const auto& bpp_j = bpp[j];
        const auto item = bpp_j.find(i);
        return item == bpp_j.end() ? 0.0 : item->second;
    }

    // methods declared in file forward.cpp
    template <Mode mode>
    unsigned long beamstep_H(const int j);
    template <Mode mode>
    unsigned long beamstep_Multi(const int j);
    template <Mode mode>
    unsigned long beamstep_P(const int j);
    template <Mode mode>
    unsigned long beamstep_M2(const int j);
    template <Mode mode>
    unsigned long beamstep_M(const int j);
    template <Mode mode>
    void beamstep_C(const int j);

    // methods declared in file backward.cpp
    std::pair<unsigned long, unsigned long> backward_update(const int i, const int j, State& state,
                                                            const StateType type, const value_type edge_threshold);

    void get_incoming_edges_state(const int i, const int j, const StateType type);

    template <Mode mode>
    void get_incoming_hedges_C(const int j);
    template <Mode mode>
    void get_incoming_hedges_M(const int i, const int j);
    template <Mode mode>
    void get_incoming_hedges_M2(const int i, const int j);
    template <Mode mode>
    void get_incoming_hedges_P(const int i, const int j);
    template <Mode mode>
    void get_incoming_hedges_Multi(int i, int j);

    std::string get_threshknot_structure(float threshknot_threshold = 0.3f,
                                         int min_helix_size = linearx::constants::energy::MIN_HELIX_SIZE) const;

    void dump_bpp(const std::string& out_dir) const;
};

class LinearPartition final : public LinearPartitionInterface<LinearPartition> {
   public:
    LinearPartition(const Sequence& seq, const EnergyModel& energy_model, const bool allow_sharp_turn = false)
        : LinearPartitionInterface<LinearPartition>(seq, energy_model, allow_sharp_turn) {}

    inline State* check_state(const StateType type, const int i, const int j) {
        return LinearPartitionInterface<LinearPartition>::get_state(type, i, j, true);
    }

    inline void update_state_beta(HEdge hedge, State* state, const StateType type, const int i, const int j) {
        hedge.weight *= linearx::constants::energy::INV_KT;
        hedge.update_state_beta(*state);
    }
};
