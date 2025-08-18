#ifndef TURBOFOLD
#define TURBOFOLD

#include <vector>

#include "./../linear_align/linear_align.hpp"
#include "./../partition/partition.hpp"
#include "./utility.hpp"

#include "./utils/GuideTree.h"
#include "./utils/ProbabilisticModel.h"

class TurboAlign;
class TurboPartition;

class LinearTurboFold {
    const static int min_beam_size = 50;

   private:
    MultiSeq *multi_seq;
    MultiSeq *multi_alignment;
    EnergyModel energy_model;

    std::vector<float> seq_identities;  // contains pairwise sequence identities for all k^2 sequence pairs
    std::vector<TurboPartition> pfs;    // contains partition function objects for k sequences
    std::vector<TurboAlign> alns;       // contains alignment objects for all k^2 sequence pairs

   public:
    // const EnergyParamsType energy_params;
    const int num_itr;  // number of iterations
    const bool use_lazy_outside;
    const bool use_prev_outside_score;
    const bool shrink_beam;
    const float lambda;  // extrinsic information weight (contribution relative to intrinsic information)

    const float alignment_pruning_threshold;
    const float folding_pruning_threshold;
    const float threshknot_threshold;
    const float min_helix_size;

    const VerboseState verbose_state;

    int itr;                                                                 // current iteration
    int beam_size = 100;                                                     // current beam size for beam pruning
    std::vector<std::vector<std::unordered_map<int, double>>> extinf_cache;  // cache for extrinsic information
    vector<vector<unordered_map<int, double>*>> consistency_transform;

    LinearTurboFold(MultiSeq *multi_seq, const EnergyParamsType energy_params, const int num_itr,
                    const bool use_lazy_outside, const bool use_prev_outside_score, const bool shrink_beam,
                    const float lambda = 0.3, const float alignment_pruning_threshold = -DEVIATION_THRESHOLD,
                    const float folding_pruning_threshold = -DEVIATION_THRESHOLD, const float threshknot_threshold = 0.3,
                    const float min_helix_size = 3, VerboseState verbose_state = VerboseState::DEBUG)
        : multi_seq(multi_seq),
          energy_model(energy_params),
          num_itr(num_itr),
          use_lazy_outside(use_lazy_outside),
          use_prev_outside_score(use_prev_outside_score),
          shrink_beam(shrink_beam),
          lambda(lambda),
          alignment_pruning_threshold(alignment_pruning_threshold),
          folding_pruning_threshold(folding_pruning_threshold),
          threshknot_threshold(threshknot_threshold),
          min_helix_size(min_helix_size),
          verbose_state(verbose_state) {
        size_t num_pairs = (multi_seq->size() * (multi_seq->size() - 1)) / 2;

        // reserve space for sequence pairs and sequence identities
        alns.reserve(num_pairs);
        pfs.reserve(multi_seq->size());
        seq_identities.reserve(num_pairs);
        consistency_transform.resize(multi_seq->size());

        // reserve space extrinsic info cache
        extinf_cache.resize(multi_seq->size());

        for (int i = 0; i < multi_seq->size(); i++) {
            multi_seq->at(i).k_id = i;  // set k_id for each sequence
            pfs.emplace_back(this, &(multi_seq->at(i)), energy_model,
                             use_prev_outside_score);  // better than pfs.push_back(),
                                                       // creates Partition object directly
                                                       // inside the container
            
            consistency_transform[i].resize(multi_seq->size());

            // enumerate all possible k^2 sequence pairs and create LinearAlign objects
            for (int j = i + 1; j < multi_seq->size(); j++) {
                alns.emplace_back(this, &(multi_seq->at(i)), &(multi_seq->at(j)), use_prev_outside_score);
                seq_identities.push_back(0.0f);
            }

            // initialize extrinsic info cache for each sequence
            extinf_cache[i].resize(multi_seq->at(i).length());
        }
    }

    int get_seq_pair_index(const int k1, const int k2);
    double get_extrinsic_info(const Seq &x, int i, int j);
    void reset_extinf_cache();
    int multiple_sequence_alignment();
    void dump_coinc_probs2(const std::string &filepath, const float threshold, std::unordered_map<int, double>* coinc_prob, int seqlen);

    void run_phmm_alignment();
    void run();

    // ~LinearTurboFold() {
    //     std::cerr << "~LinearTurboFold: start\n";

    //     std::cerr << " clearing alns...\n";
    //     alns.clear();            // runs ~TurboAlign / ~LinearAlign now
    //     std::cerr << " cleared alns\n";

    //     std::cerr << " clearing pfs...\n";
    //     pfs.clear();             // runs ~TurboPartition / ~Partition now
    //     std::cerr << " cleared pfs\n";

    //     std::cerr << " deleting extinf_cache\n";

    //     // extinf first too
    //     extinf_cache.clear();

    //     std::cerr << " deleting seq_identities\n";

    //     seq_identities.clear();

    //     std::cerr << " deleting consistency_transform\n";

    //     // // only after alns/pfs are gone, free CT
    //     // int n = multi_seq ? multi_seq->size() : 0;
    //     // for (int i = 0; i < n; ++i) {
    //     //     for (int j = 0; j < n; ++j) {
    //     //         std::cerr << " deleting {i, j} " << i << " " << j << "\n";
    //     //         delete[] consistency_transform[i][j];   // ok if nullptr
    //     //         consistency_transform[i][j] = nullptr;
    //     //     }
    //     // }

    //     // std::cerr << " consistency_transform.clear()\n";

    //     // consistency_transform.clear();

    //     std::cerr << " deleting multi_alignment\n";
    //     delete multi_alignment; 
    //     multi_alignment = nullptr;

    //     std::cerr << "~LinearTurboFold: end\n";
    // }

};

class TurboPartition final : public Partition {
   private:
    LinearTurboFold *turbofold;
    PartitionFunctionBeam pfb;
    bool use_prev_outside_score;

   public:
    friend class LinearTurboFold;

    ProbAccm prob_accm;

    TurboPartition(LinearTurboFold *turbofold, const Seq *sequence, EnergyModel &energy_model,
                   bool use_prev_outside_score, bool allow_sharp_turn = false, bool verbose_output = true)
        : Partition(sequence, energy_model, InsideMode::PARTITION, allow_sharp_turn, verbose_output),
          turbofold(turbofold),
          pfb(use_prev_outside_score ? sequence->length() : 0),
          use_prev_outside_score(use_prev_outside_score) {}

    ~TurboPartition();

    double beam_prune(StateType type, int j, int beam_size);
    inline bool check_state(StateType type, int i, int j) const override;
    void compute_inside(int beam_size = 100) override;
    void calc_prob_accm();
};

class TurboAlign final : public LinearAlign {
   private:
    LinearTurboFold *turbofold;
    AlignBeam ab;
    bool use_prev_outside_score;

   public:
    friend class LinearTurboFold;

    TurboAlign(LinearTurboFold *turbofold, Seq *sequence1, Seq *sequence2, bool use_prev_outside_score,
               bool verbose = false, double alpha1 = 1.0, double alpha2 = 0.8, double alpha3 = 0.5)
        : LinearAlign(sequence1, sequence2, verbose, alpha1, alpha2, alpha3),
          turbofold(turbofold),
          ab(use_prev_outside_score ? sequence1->length() : 0, use_prev_outside_score ? sequence2->length() : 0),
          use_prev_outside_score(use_prev_outside_score) {}

    ~TurboAlign();

    double beam_prune(std::unordered_map<std::pair<int, int>, HState, PairHash> &beamstep, HStateType h,
                      int beam_size) override;

    inline bool check_state(const int i, const int j, const HStateType h) override;
};

#endif  // TURBOFOLD
