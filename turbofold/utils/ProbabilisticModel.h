#ifndef _TURBOFOLD_PROBABILISTICMODEL_
#define _TURBOFOLD_PROBABILISTICMODEL_

#include <list>
#include <cmath>
#include <cstdio>
#include <unordered_map>
#include <vector>

#include "../../sequence/multi_seq.hpp"
#include "./../../utility/log_math.hpp"
#include "GuideTree.h"

#ifdef TURBOHOMOLOGY
	#define ChooseBestOfThree ChooseBestOfThreeturbohomology
	#define ProbabilisticModel ProbabilisticModelturbohomology
#endif


using namespace std;

//! ProbabilisticModel Class
/*!
    The ProbabilisticModel Class stores the parameters of a probabilistic model.

*/

#define NumMatrixTypes 3   // One match state and two insert states.

//! Store the largest of three values x1, x2, and x3 in *x.  
//! If xi is the largest value, then store bi in *b.
void ChooseBestOfThree(float x1, float x2, float x3, char b1, char b2, char b3, float *x, char *b);

class ProbabilisticModel{

public:

    //! Computes an alignment based on given posterior matrix.
    //! This is done by finding the maximum summing path (or maximum weight trace) through the posterior matrix.
    //! The final alignment is returned as a pair consisting of:
    //! (1) a string (e.g., XXXBBXXXBBBBBBYYYYBBB) where X's and denote insertions in one of the two sequences and B's denote that both sequences are present (i.e. matches).
    //! (2) a float indicating the sum achieved.
    vector<float> *LinearBuildPosterior(MultiSeq *align1, MultiSeq *align2,
                      const vector<vector<unordered_map<int, double>*>> &mul_aln_rets, float cutoff = 0.0f) const;

    unordered_map<int, double> *LinearMultiAlnResults(MultiSeq *align1, MultiSeq *align2,
                   const vector<vector<unordered_map<int, double>*>> &mul_aln_rets, float cutoff = 0.0f) const;

    pair<vector<char> *, float> LinearComputeAlignment(int hmmBeam, int seq1Length, int seq2Length, const unordered_map<int, double>* posterior) const;

    //! This function computes the consistency transformation for sequence Z by taking two posterior probabilities matrices of alignments between X-Z and Z-Y.
    //! For the case that sequence Z's index is larger than that of X.
    //! The transformed matrix is added to \param posterior.
    void LinearConsistencyTransform(int lengthX, unordered_map<int, double>* &xz_aln_ret, unordered_map<int, double>* &zy_aln_ret, unordered_map<int, double>* &new_xy_ret);

    //! This function takes multiple sequences and posterior probability matices to perform three-way probabilistic consistency transformation.
    //! Returns new re-estimated alignment score matrices. 
    //! The formula is: P'(x[i]-y[j])=(1/|S|)*sum_z_in_S{ sum_k{ P(x[i]-z[k]) * P(z[k]-y[j]) } }
    vector<vector<unordered_map<int, double>*>> LinearMultiConsistencyTransform(MultiSeq *sequences, vector<vector<unordered_map<int, double>*>> &consistency_transform);

    //! This function takes two multiple sequence alignments as input.
    //! Returns the alignment of the two MultiSequence objects.
    MultiSeq *LinearAlignAlignments (MultiSeq *align1, MultiSeq *align2,
                                    const vector<vector<unordered_map<int, double>*>> &consistency_transform,
                                    const ProbabilisticModel &model, int hmmBeam);

    //! This function takes guide tree (computed by distance) as input.
    //! Returns the aligned sequences corresponding to a node or leaf of a guide tree.
    MultiSeq *LinearProcessTree (const TreeNode *tree, MultiSeq *sequences,
                                const vector<vector<unordered_map<int, double>*>> &consistency_transform,
                                const ProbabilisticModel &model, int hmmBeam);

    //! This function computes the final alignment by calling ProcessTree() and performing iterative refinement.
    MultiSeq *LinearComputeFinalAlignment (const TreeNode *tree, MultiSeq *sequences,
                                        const vector<vector<unordered_map<int, double>*>> &consistency_transform,
                                        const ProbabilisticModel &model, int hmmBeam);

    //! This function performs randomized partitioning iterative refinement. 
    //! Taking posterior probability matrices, parameters of probabilistic model, and multiple sequence alignments.
    //! Returns a new multiple sequence alignment.
    void LinearDoIterativeRefinement (const vector<vector<unordered_map<int, double>*>> &consistency_transform,
                                const ProbabilisticModel &model, MultiSeq* &alignment, int i, int hmmBeam);
                                
};

struct AlignState {
	unsigned i;
	unsigned k;
    double alpha;
    double beta;
    int manner;
	int pre;
	unsigned step;
	bool close; 

    AlignState(): alpha(xlog(0)), beta(xlog(0)), manner(0), pre(0), close(false), step(0), i(0), k(0) {};

    void set(double score_, int pre_manner_, int manner_, unsigned step_, unsigned i_, unsigned k_) {
        alpha = score_; pre = pre_manner_; manner = manner_; step = step_ ; i = i_; k = k_;
    }

    void set(double score_, int manner_) {
        alpha = score_; manner = manner_;
    }

	void set_beta(double score_, int pre_manner_, int manner_, unsigned step_, unsigned i_, unsigned k_) {
        beta = score_; pre = pre_manner_; manner = manner_; step = step_ ; i = i_; k = k_;
    }
};

#endif