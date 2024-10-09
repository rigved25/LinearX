#include <vector>

#include "./../linear_align/linear_align.hpp"
#include "./../partition/partition.hpp"
#include "./../sequence/multi_seq.hpp"
#include "./utility.hpp"

class TurboPartition;

class LinearTurboFold {

  private:
    MultiSeq *multi_seq;
    EnergyModel energy_model;

    std::vector<float> seq_identities; // contains pairwise sequence identities for all k^2 sequence pairs
    std::vector<TurboPartition> pfs;   // contains partition function objects for k sequences
    std::vector<LinearAlign> alns;     // contains alignment objects for all k^2 sequence pairs

    Cache3D ext_info_cache; // cache for extrinsic information

    int itr; // current iteration

  public:
    const int num_itr; // number of iterations
    const VerboseState verbose_state;

    LinearTurboFold(MultiSeq *multi_seq, const int num_itr, const VerboseState verbose_state)
        : multi_seq(multi_seq), energy_model(EnergyParamsType::VIENNA), num_itr(num_itr), verbose_state(verbose_state) {

        size_t num_pairs = (multi_seq->size() * (multi_seq->size() - 1)) / 2;

        // reserve space for sequence pairs and sequence identities
        alns.reserve(num_pairs);
        seq_identities.reserve(num_pairs);

        for (int i = 0; i < multi_seq->size(); i++) {
            multi_seq->at(i).k_id = i;                                 // set k_id for each sequence
            pfs.emplace_back(this, &(multi_seq->at(i)), energy_model); // better than pfs.push_back(),
                                                                       // creates Partition object directly
                                                                       // inside the container

            // enumerate all possible k^2 sequence pairs and create LinearAlign objects
            for (int j = i + 1; j < multi_seq->size(); j++) {
                alns.emplace_back(&(multi_seq->at(i)), &(multi_seq->at(j)));
                seq_identities.push_back(0.0f);
            }
        }
    }

    int get_seq_pair_index(const int k1, const int k2);
    float get_extrinsic_info(const Seq &x, int i, int j);
    void run();
};

class TurboPartition final : public Partition {
  private:
    LinearTurboFold *turbofold;

  public:
    TurboPartition(LinearTurboFold *turbofold, const Seq *sequence, EnergyModel &energy_model,
                   const int beam_size = 100, bool allow_sharp_turn = false, bool verbose_output = false)
        : Partition(sequence, energy_model, InsideMode::PARTITION, beam_size, allow_sharp_turn, verbose_output),
          turbofold(turbofold) {};

    void get_incoming_edges_state(const int i, const int j, const StateType type,
                                  std::vector<HEdge> &incoming_hedges) override;
    void compute_inside() override;
};
