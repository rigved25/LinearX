#include "ProbabilisticModel.h"
#include <list>
#include <cmath>
#include <cstdio>
#include <string>
#include <sstream>
#include <iomanip>
#include <iostream>
#include <set>
#include <algorithm>
#include <cstdlib>
#include <cerrno>
#include <iomanip>
#include <unordered_map>
#include "./../../linear_align/utility.hpp"

using namespace std;

//! ProbabilisticModel Class
/*!
    The ProbabilisticModel Class stores the parameters of a probabilistic model.

*/


//! Store the largest of three values x1, x2, and x3 in *x.  
//! If xi is the largest value, then store bi in *b.
// void ChooseBestOfThree(float x1, float x2, float x3, char b1, char b2, char b3, float *x, char *b){
//     if (x1 >= x2){
//         if (x1 >= x3){
//             *x = x1;
//             *b = b1;
//             return;
//         }
//         *x = x3;
//         *b = b3;
//         return;
//     }
//     if (x2 >= x3){
//         *x = x2;
//         *b = b2;
//         return;
//     }
//     *x = x3;
//     *b = b3;
// }

unsigned long tmp_quickselect_partition(vector<pair<double, int>>& scores, unsigned long lower, unsigned long upper) {
    double pivot = scores[upper].first;
    while (lower < upper) {
        while (scores[lower].first < pivot) ++lower;
        while (scores[upper].first > pivot) --upper;
        if (scores[lower].first == scores[upper].first) ++lower;
        else if (lower < upper) swap(scores[lower], scores[upper]);
    }
    return upper;
}

double tmp_quickselect(vector<pair<double, int>>& scores, unsigned long lower, unsigned long upper, unsigned long k) {
    if ( lower == upper ) return scores[lower].first;
    unsigned long split = tmp_quickselect_partition(scores, lower, upper);
    unsigned long length = split - lower + 1;
    if (length == k) return scores[split].first;
    else if (k  < length) return tmp_quickselect(scores, lower, split-1, k);
    else return tmp_quickselect(scores, split+1, upper, k - length);
}

double tmp_beam_prune(std::unordered_map<int, AlignState> &beamstep, int beamsize){
    vector<pair<double, int>> scores;
    for (auto &item : beamstep) {
        int ik = item.first;
        AlignState &cand = item.second;
        scores.push_back(make_pair(cand.alpha, ik));
        // scores.push_back(make_pair(cand.alpha, ik));
    }
    if (scores.size() <= beamsize) return VALUE_MIN;
    double threshold = tmp_quickselect(scores, 0, scores.size() - 1, scores.size() - beamsize);
    for (auto &p : scores) {
        if (p.first < threshold) beamstep.erase(p.second);
    }
    // cout << "threshold: " << threshold << endl;
    return threshold;
}

pair<vector<char> *, float> ProbabilisticModel::LinearComputeAlignment(int hmmBeam, int seq1Length, int seq2Length, const unordered_map<int, double>* posterior) const {
    unordered_map<int, AlignState>* max_exp_acc = new unordered_map<int, AlignState>[seq1Length + seq2Length + 3];
    unsigned max_len = seq1Length > seq2Length ? seq1Length : seq2Length;
    // step 0
    max_exp_acc[0][0].alpha = 0;
    max_exp_acc[0][0].manner = 3; // ALIGN_ALN;
    max_exp_acc[0][0].step = 0;
    max_exp_acc[0][0].i = 0;
    max_exp_acc[0][0].k = 0;

    for(int s = 0; s < seq1Length + seq2Length + 1; ++s){
        int beamsize = hmmBeam;
        if (beamsize > 0 && max_exp_acc[s].size() > beamsize) tmp_beam_prune(max_exp_acc[s], beamsize);
        for (auto &item : max_exp_acc[s]) {
            AlignState &state = item.second;
            unsigned i = state.i;
            unsigned k = state.k;
            unsigned step = state.step;
            unsigned next_key;
            int manner = state.manner;
            int next_manner;
            unsigned next_i, next_k, next_step;
            for (int m = 3; m >= 1; m--){
                next_step = step;
                switch (m)
                {
                case 3: // ALIGN_ALN:
                    next_manner = m; // ALIGN_ALN;
                    next_i = i + 1;
                    next_k = k + 1;
                    next_step += 2;
                    next_key = next_i * max_len + next_k;

                    // if ((next_i <= seq1Length) && (next_k <= seq2Length)) {
                    //     if (posterior[next_i].find(next_k) == posterior[next_i].end()) continue;
                    //     max_exp_acc[next_step][next_key].alpha = state.alpha + posterior[next_i].at(next_k);
                    //     max_exp_acc[next_step][next_key].i = next_i;
                    //     max_exp_acc[next_step][next_key].k = next_k;
                    //     max_exp_acc[next_step][next_key].step = next_step;
                    //     max_exp_acc[next_step][next_key].manner = m;
                    //     // cout << i << " " << k << " " << next_i << " " << next_k << " " << tmp_aln_rets[next_step][next_key].alpha << endl;
                    // }

                    if (next_i > 0 && next_i <= seq1Length && next_k > 0 && next_k <= seq2Length) {
                        // 0-based indexed in posterior 
                        if (posterior[next_i-1].find(next_k-1) == posterior[next_i-1].end()) continue;
                        max_exp_acc[next_step][next_key].alpha = state.alpha + posterior[next_i-1].at(next_k-1);
                        max_exp_acc[next_step][next_key].i = next_i;
                        max_exp_acc[next_step][next_key].k = next_k;
                        max_exp_acc[next_step][next_key].step = next_step;
                        max_exp_acc[next_step][next_key].manner = m;
                        // cout << i << " " << k << " " << next_i << " " << next_k << " " << tmp_aln_rets[next_step][next_key].alpha << endl;
                    }

                    // if (next_i > 0 && next_i <= seq1Length && next_k > 0 && next_k <= seq2Length) {
                    //     auto &rowMap = posterior[next_i - 1];    // safe array access
                    //     auto it = rowMap.find(next_k - 1);       // safe map lookup
                    //     if (it == rowMap.end()) 
                    //         continue;
                    //     double p = it->second;

                    //     auto &cell = max_exp_acc[next_step][ next_i * max_len + next_k ];
                    //     cell.alpha   = state.alpha + p;
                    //     cell.i       = next_i;
                    //     cell.k       = next_k;
                    //     cell.step    = next_step;
                    //     cell.manner  = m;
                    // }
                    break;

                case 1: // ALIGN_INS1:
                    next_manner = m; // ALIGN_INS1;
                    next_i = i + 1;
                    next_k = k;
                    next_step += 1;
                    next_key = next_i * max_len + next_k;
                    if ((next_i <= seq1Length) && (next_k <= seq2Length)){
                        if (state.alpha > max_exp_acc[next_step][next_key].alpha){
                            max_exp_acc[next_step][next_key].alpha = state.alpha;
                            max_exp_acc[next_step][next_key].i = next_i;
                            max_exp_acc[next_step][next_key].k = next_k;
                            max_exp_acc[next_step][next_key].step = next_step;
                            max_exp_acc[next_step][next_key].manner = m;
                            // cout << next_i << " " << next_k << " " << max_exp_acc[next_step][next_key].alpha << endl;
                        }
                    }
                    break;

                case 2: // ALIGN_INS2:
                    next_manner = m; // ALIGN_INS2;
                    next_i = i;
                    next_k = k + 1;
                    next_step += 1;
                    next_key = next_i * max_len + next_k;
                    if ((next_i <= seq1Length) && (next_k <= seq2Length)) {
                        if (state.alpha > max_exp_acc[next_step][next_key].alpha){
                            max_exp_acc[next_step][next_key].alpha = state.alpha;
                            max_exp_acc[next_step][next_key].i = next_i;
                            max_exp_acc[next_step][next_key].k = next_k;
                            max_exp_acc[next_step][next_key].step = next_step;
                            max_exp_acc[next_step][next_key].manner = m;
                            // cout << next_i << " " << next_k << " " << max_exp_acc[next_step][next_key].alpha << endl;
                        }
                    }
                    break;
                default:
                    break;
                }
            }
        }
    }

    // compute traceback
    vector<char> *alignment = new vector<char>; 
    assert (alignment);

    int seq1_pos = seq1Length;
    int seq2_pos = seq2Length;
    int step = seq1_pos + seq2_pos;
    int key = seq1_pos * max_len + seq2_pos;
    int pre_manner = max_exp_acc[seq1Length + seq2Length][key].manner;
    double total = max_exp_acc[seq1Length + seq2Length][key].alpha;
    // cout << "max_exp_acc[seq1Length + seq2Length][key].manner" << pre_manner << endl;
    while (true) {
        if ((seq1_pos == 0) && (seq2_pos == 0)) break;
        // cout << step << " " << key << " " << seq1_pos << " " << seq2_pos << " " << pre_manner  << endl;
        switch (pre_manner)
        {
        case 3:
            alignment->push_back ('B');
            // cout << "B";
            seq1_pos--;
            seq2_pos--;
            step -= 2;
            key = seq1_pos * max_len + seq2_pos;
            pre_manner = max_exp_acc[step][key].manner;
            break;
        case 2: // INS2
            alignment->push_back ('Y');
            // cout << "Y";
            seq2_pos--;
            step--;
            key = seq1_pos * max_len + seq2_pos;
            pre_manner = max_exp_acc[step][key].manner;
            break;
        case 1: // INS1
            alignment->push_back ('X');
            // cout << "X";
            seq1_pos--;
            step--;
            key = seq1_pos * max_len + seq2_pos;
            pre_manner = max_exp_acc[step][key].manner;
            break;
        case 0:
            cout << "ComputeAlignment traceback error !!" << endl;
            cout << seq1_pos << " " << seq2_pos << endl;
            exit(0);
        default:
            break;
        }
    }
    // cout << endl;
    // alignment->push_back('X');
    reverse(alignment->begin(), alignment->end());
    // for(auto i = alignment->begin(); i!=alignment->end(); i++){
    //     cout << *i;
    // }
    // cout << endl;

    // cout << "total: " << total - posterior[0][0].value << endl;
    delete[] max_exp_acc;
    return make_pair(alignment, total);
}

/// Equivalent to BuildPosterior from the Probcons 
unordered_map<int, double> * ProbabilisticModel::LinearMultiAlnResults(MultiSeq *align1, MultiSeq *align2, const vector<vector<unordered_map<int, double>*>> &consistency_transform, float cutoff) const {
    const int seq1Length = align1->at(0).length();
    const int seq2Length = align2->at(0).length();
    
    unordered_map<int, double>* sum_aln_ret = new unordered_map<int, double>[seq1Length + 1];
    // cout << seq1Length << " " << seq2Length << endl;
    for (int i = 0; i < align1->size(); i++){
        int first = align1->at(i).k_id;
        vector<int> *mapping1 = align1->at(i).get_mapping();

        // Loops through align2
        for (int j = 0; j < align2->size(); j++){
            int second = align2->at(j).k_id;
            vector<int> *mapping2 = align2->at(j).get_mapping();
            // cout << "seqs: " << i << " " << j << " " << align1->GetSequence(i)->GetLength() << " " << align2->GetSequence(j)->GetLength() << endl;

            if(first < second){
                unordered_map<int, double>* aln_ret = consistency_transform[first][second];

                int seq1len = mapping1->size() - 1;
                for (int ii = 0; ii < seq1len; ii++){
                    int ibase = (*mapping1)[ii];

                    for (auto &item : aln_ret[ii]) {
                        // cout << "ii: " << ii << " " << item.first << endl;
                        int jbase = (*mapping2)[item.first];
                        if (item.second < 0.01) continue;

                        sum_aln_ret[ibase][jbase] += item.second;
                        // cout << i <<  " "  << j <<  " "  << first << " "  << second << " " << ii << " "  << item.first << " " << ibase  << " " << jbase  << " " <<  item.second.value << " " << sum_aln_ret[ibase][jbase].value << endl;
                    }
                }
            } else {
                unordered_map<int, double>* aln_ret = consistency_transform[second][first];

                int seq2len = mapping2->size() - 1;
                for (int jj = 0; jj <= seq2len; jj++){
                    int jbase = (*mapping2)[jj];
                    for (auto &item : aln_ret[jj]) {
                        int ibase = (*mapping1)[item.first];
                        if (item.second < 0.01) continue;
                        sum_aln_ret[ibase][jbase] += item.second;
                        // cout << i <<  " "  << j <<  " "  <<  first << " "  << second << " " << ibase  << " " << jbase <<  item.second.value << " " << sum_aln_ret[ibase][jbase].value << endl;
                    }
                }
            }
            delete mapping2;
        }
        delete mapping1;
    }

    
    return sum_aln_ret;
}

void ProbabilisticModel::LinearConsistencyTransform(int lengthX, unordered_map<int, double>* &xz_consistency_transform, unordered_map<int, double>* &zy_consistency_transform, unordered_map<int, double>* &new_xy_consistency_transform){
    // int lengthX = xz_consistency_transform->size();
    // cout << "lengthX: " << lengthX << endl;
    // int lengthZ = zy_consistency_transform->size();

    for(int i = 0; i < lengthX; i++){

        //cerr << " i " << i << " xz_CT[i].size() " << xz_consistency_transform[i].size()  << endl;  
        for(auto &xz_cand : xz_consistency_transform[i]){
            int k = xz_cand.first;

            //cerr << " k " << k << "zy_CT[k].size() " << zy_consistency_transform[k].size()  << endl;  
            for(auto &zy_cand : zy_consistency_transform[k]){
                int j = zy_cand.first;
                
                //cerr << " j " << j << endl;
                new_xy_consistency_transform[i][j] += xz_cand.second * zy_cand.second;
                // cerr << xz_cand.second * zy_cand.second << " ";
            }
        }
        // cerr << endl;
    }
}

// linearTurboFold
vector<vector<unordered_map<int, double>*>> ProbabilisticModel::LinearMultiConsistencyTransform(MultiSeq *sequences, vector<vector<unordered_map<int, double>*>> &consistency_transform){
    const int numSeqs = sequences->size();

    // For every pair of sequences
    for (int i = 0; i < numSeqs; i++){
        for (int j = i+1; j < numSeqs; j++){
            Seq seq1 = sequences->at(i);
            Seq seq2 = sequences->at(j);

            const int seq1Length = seq1.length();
            const int seq2Length = seq2.length();

            // allocate space for temporary results
            unordered_map<int, double>* temp_pair_CT = new unordered_map<int, double>[seq1Length];

            // Get the original alignment result
            unordered_map<int, double>* &pair_CT = consistency_transform[i][j];

            // Contribution from the summation where z = x and z = y
            // cout << "seq1length: " << seq1Length << endl;
            for (int k = 0; k < seq1Length; k++){
                for(auto &item : pair_CT[k]){
                    int l = item.first;
                    temp_pair_CT[k][l] = 2 * item.second;
                    // cout << i << " " << j << " " << k << " " << l << " " << new_aln_ret[k][l].value << endl;
                }
            }

            // Contribution from all other sequences
            for (int k = 0; k < numSeqs; k++) {
                if (k == i || k == j) continue;

                // unordered_map<int, double> *xzCT, *zyCT;
                // if (k < i) {
                //     xzCT = consistency_transform[k][i];
                //     zyCT = consistency_transform[k][j];
                // }
                // else if (k > j) {
                //     xzCT = consistency_transform[i][k];
                //     zyCT = consistency_transform[j][k];
                // }
                // else{
                //     xzCT = consistency_transform[i][k];
                //     zyCT = consistency_transform[k][j];
                // }
                // cerr << "seqs: " << i << " " << j << " " << k << endl;
                // LinearConsistencyTransform(seq1Length, xzCT, zyCT, temp_pair_CT);

                // cout << "seqs: " << i << " " << j << " " << k << endl;
                LinearConsistencyTransform(seq1Length, consistency_transform[i][k], consistency_transform[k][j], temp_pair_CT);
            }

            // Renormalization
            for (int k = 0; k < seq1Length; k++){
                for(auto &item : temp_pair_CT[k]){
                    int l = item.first;
                    temp_pair_CT[k][l] /= numSeqs;
                }
            }

            // Mask out positions not originally in the posterior matrix
            for (int k = 0; k < seq1Length; k++){
                for(auto &item : pair_CT[k]){
                    int l = item.first;
                    if (temp_pair_CT[k].find(l) == temp_pair_CT[k].end()) continue; // N.B.
                    if (temp_pair_CT[k][l] >= 0.01){
                        consistency_transform[i][j][k][l] = temp_pair_CT[k][l];
                        consistency_transform[j][i][l][k] = temp_pair_CT[k][l];
                        // cout << i << " " << j << " " << k << " " << l << " " << temp_pair_CT[k][l].value << endl;
                    }
                }
            }
            delete[] temp_pair_CT;
        }
    }

    return consistency_transform;
}

MultiSeq* ProbabilisticModel::LinearAlignAlignments (MultiSeq *align1, MultiSeq *align2,
                                const vector<vector<unordered_map<int, double>*>> &consistency_transform,
                                const ProbabilisticModel &model, int hmmBeam){

    // Print some info about the alignment
    // vector<float> *posterior = model.LinearBuildPosterior (align1, align2, consistency_transform);
    // model.ComputeAlignment (align1->GetSequence(0)->GetLength(), align2->GetSequence(0)->GetLength(), *posterior);
    // delete posterior;

    // Choose the alignment routine depending on the "cosmetic" gap penalties used
    const unordered_map<int, double> *sum_aln_ret = model.LinearMultiAlnResults(align1, align2, consistency_transform);

    pair<vector<char> *, float> alignment = LinearComputeAlignment(hmmBeam, align1->at(0).length(), align2->at(0).length(), sum_aln_ret);
    delete[] sum_aln_ret;

    // Build final alignment
    MultiSeq *result = new MultiSeq();
    for (int i = 0; i < align1->size(); i++)
        result->add_sequence (  *(align1->at(i).add_gaps(alignment.first, 'X')) );
    for (int i = 0; i < align2->size(); i++)
        result->add_sequence (  *(align2->at(i).add_gaps(alignment.first, 'Y')) );
    // result->SortByLabel();  // naukarkr

    // Free temporary alignment
    delete alignment.first;

    return result;
}

MultiSeq* ProbabilisticModel::LinearProcessTree (const TreeNode *tree, MultiSeq *sequences,
                            const vector<vector<unordered_map<int, double>*>> &consistency_transform,
                            const ProbabilisticModel &model, int hmmBeam){
    MultiSeq *result;
    
    // Check if this is an internal node of the alignment tree
    if (tree->GetSequenceLabel() == -1){
        MultiSeq *alignLeft = LinearProcessTree (tree->GetLeftChild(), sequences, consistency_transform, model, hmmBeam);
        MultiSeq *alignRight = LinearProcessTree (tree->GetRightChild(), sequences, consistency_transform, model, hmmBeam);

        assert (alignLeft);
        assert (alignRight);

        result = LinearAlignAlignments (alignLeft, alignRight, consistency_transform, model, hmmBeam);
        assert (result);

        delete alignLeft;
        delete alignRight;
    }

    // Otherwise, this is a leaf of the alignment tree
    else {
        result = new MultiSeq(); assert (result);
        result->add_sequence( *(sequences->at(tree->GetSequenceLabel()).clone()) );
    }

    return result;
}

void ProbabilisticModel::LinearDoIterativeRefinement (const vector<vector<unordered_map<int, double>*>> &consistency_transform,
                            const ProbabilisticModel &model, MultiSeq* &alignment, int i, int hmmBeam){
    set<int> groupOne, groupTwo;

    // Create two separate groups
    for (int i = 0; i < alignment->size(); i++){
        if (rand() % 2)
          groupOne.insert (i);
        else
          groupTwo.insert (i);
    }

    if (groupOne.empty() || groupTwo.empty()) return;

    // Project into the two groups
    MultiSeq *groupOneSeqs = alignment->Project (groupOne); assert (groupOneSeqs);
    MultiSeq *groupTwoSeqs = alignment->Project (groupTwo); assert (groupTwoSeqs);
    delete alignment;

    // Realign
    alignment = LinearAlignAlignments (groupOneSeqs, groupTwoSeqs, consistency_transform, model, hmmBeam);    

    delete groupOneSeqs;
    delete groupTwoSeqs;
}

MultiSeq* ProbabilisticModel::LinearComputeFinalAlignment (const TreeNode *tree, MultiSeq *sequences,
                                      const vector<vector<unordered_map<int, double>*>> &consistency_transform,
                                      const ProbabilisticModel &model, int hmmBeam){
    
    unsigned int num_iterative_refinement_reps = 100;

    cerr << endl << "[Multi Seq Align] Processing the Guide Tree for initial alignment " << endl;

    MultiSeq *alignment = LinearProcessTree (tree, sequences, consistency_transform, model, hmmBeam);

    cout << endl << "[Multi Seq Align] Alignment from the processing of Guide Tree " << endl;
    alignment->print_sequences();

    // Iterative refinement
    for (int i = 0; i < num_iterative_refinement_reps; i++)
        LinearDoIterativeRefinement (consistency_transform, model, alignment, i, hmmBeam);

    return alignment;
}