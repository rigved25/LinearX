#include "./../linear_align/linear_align.hpp"
#include "./../partition/partition.hpp"
#include "./../sequence/multi_seq.hpp"
#include <vector>

class TurboPartition;

class LinearTurboFold {

  private:
    MultiSeq *multi_seq;
    EnergyModel energy_model;

    std::vector<TurboPartition> pfs; // contains partition function objects for k sequences
    std::vector<LinearAlign> alns;   // contains alignment objects for all k^2 sequence pairs

  public:
    const int num_itr; // number of iterations

    LinearTurboFold(MultiSeq *multi_seq, const int num_itr);

    LinearAlign* get_alignment_obj(int k1, int k2);
    float get_extrinsic_info(const Seq *x, int i, int j);
    void run();
};

class TurboPartition final : public Partition {
  private:
    LinearTurboFold *turbofold;

  public:
    TurboPartition(LinearTurboFold *turbofold, const Seq *sequence, EnergyModel &energy_model,
                   const int beam_size = 100, bool allow_sharp_turn = false, bool verbose_output = true)
        : Partition(sequence, energy_model, InsideMode::PARTITION, beam_size, allow_sharp_turn, verbose_output),
          turbofold(turbofold) {};

    float turbo_beam_prune(std::unordered_map<int, State> &beamstep, const int j);
    void compute_inside() override;
};
