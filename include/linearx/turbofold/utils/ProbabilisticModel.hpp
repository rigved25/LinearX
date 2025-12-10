#include <list>
#include <cmath>
#include <cstdio>
#include <unordered_map>
#include <vector>
#include <fstream>
#include <string>
#include <utility>
#include <limits>

#include "GuideTree.h"
#include "linearx/utility.hpp"
#include "linearx/math/log_math.hpp"
#include "linearx/sequence/multi_sequence.hpp"
#include "linearx/turbofold/utility.hpp"

using namespace linearx::utils;

struct MEAState {
    value_type score;      // MEA reward (standard DP)
    char traceback;        // 'D'=diagonal, 'U'=up, 'L'=left

    MEAState()
        : score(linearx::math::LOG_ZERO),
          traceback('\0') {}

    MEAState(value_type s, char t) : score(s), traceback(t) {}
};

struct DijkstraLazyState {
    value_type score;      // reward accumulated along path
    value_type cost;       // priority value (1 - prob accum)
    char traceback;

    DijkstraLazyState()
        : score(linearx::math::LOG_ZERO),
          cost(std::numeric_limits<value_type>::infinity()),
          traceback('\0') {}

    DijkstraLazyState(value_type s, value_type c, char t) : score(s), cost(c), traceback(t) {}
};

struct DijkstraState {
    value_type cost;
    char traceback;

    DijkstraState()
        : cost(std::numeric_limits<value_type>::infinity()),
          traceback('\0') {}

    DijkstraState(value_type c, char t) : cost(c), traceback(t) {}
};

class ProbabilisticModel {
   public:
    std::string out_dir_;
    unsigned int skip_beam_prune_ = 6;
    unsigned int num_iterative_refinement_reps_ = 100;
    unsigned int num_consistency_reps_ = 2;

    ProbabilisticModel(const std::string& out_dir) : out_dir_(out_dir + "/"),
                                                    num_iterative_refinement_reps_(100),
                                                    skip_beam_prune_(6),
                                                    num_consistency_reps_(2)
    {}

    inline unsigned beam_prune(std::unordered_map<std::pair<int, int>, MEAState, linearx::utils::PairHash>& beamstep, int run_beam_size) {
        if (run_beam_size == 0 || beamstep.size() <= run_beam_size) {
            return 0;
        }

        std::vector<std::pair<value_type, std::pair<int, int>>> beam_scores;
        unsigned num_pruned = 0;
        beam_scores.clear();
        beam_scores.reserve(beamstep.size());
        for (auto& item : beamstep) {
            const std::pair<int, int> aij = item.first;
            MEAState& cand = item.second;
            beam_scores.emplace_back(cand.score, aij);
        }
        value_type threshold =
            linearx::utils::quickselect(beam_scores, 0, beam_scores.size() - 1, beam_scores.size() - run_beam_size);
        for (auto& p : beam_scores) {
            if (p.first < threshold) {
                beamstep.erase(p.second);
                num_pruned++;
            }
        }
        return num_pruned;
    }

    pair<string *, value_type> LinearComputeAlignment(int hmmBeam, int seq1Length, int seq2Length, const unordered_map<int, value_type>* posterior);
    pair<string *, value_type> LinearComputeAlignmentDijkstraScore(int hmmBeam, int seq1Length, int seq2Length, const unordered_map<int, value_type>* posterior);
    pair<string *, value_type> LinearComputeAlignmentDijkstraCost(int hmmBeam, int seq1Length, int seq2Length, const unordered_map<int, value_type>* posterior);
    pair<string *, value_type> LinearComputeAlignmentDijkstra(int hmmBeam, int seq1Length, int seq2Length, const unordered_map<int, value_type>* posterior);
    pair<string *, value_type> LinearComputeAlignmentDijkstraLogs(int hmmBeam, int seq1Length, int seq2Length, const unordered_map<int, value_type>* posterior);

    // pair<string *, value_type> LinearComputeAlignmentDijkstraEfficient(int hmmBeam, int seq1Length, int seq2Length, const unordered_map<int, value_type>* posterior);

    vector<vector<unordered_map<int, value_type>*>> LinearMultiConsistencyTransform(MultiSeq &sequences, vector<vector<unordered_map<int, value_type>*>> &posterior);

    void LinearConsistencyTransform(int seq1Length, unordered_map<int, value_type>* &xz_posterior, unordered_map<int, value_type>* &zy_posterior, unordered_map<int, value_type>* &new_xy_posterior);

    unordered_map<int, value_type>* LinearMultiAlnResults(MultiSeq* align1, MultiSeq* align2, const vector<vector<unordered_map<int, value_type>*>> &posterior, value_type cutoff) const;

    pair<MultiSeq*, value_type> LinearAlignAlignments(MultiSeq* align1, MultiSeq* align2, const vector<vector<unordered_map<int, value_type>*>> &posterior, int hmmBeam);

    pair<MultiSeq*, value_type> LinearProcessTree(const TreeNode *tree, MultiSeq* sequences, const vector<vector<unordered_map<int, value_type>*>> &posterior, int hmmBeam);

    pair<MultiSeq*, value_type> LinearDoIterativeRefinement (const vector<vector<unordered_map<int, value_type>*>> &posterior, MultiSeq* alignment, int i, int hmmBeam);

    MultiSeq* LinearComputeFinalAlignment(const TreeNode *tree, MultiSeq* sequences, const vector<vector<unordered_map<int, value_type>*>> &posterior, int hmmBeam, TurboFoldLog* log = nullptr);

};