#include "linearx/turbofold/utils/ProbabilisticModel.hpp"
#include <linearx/turbofold/utils/random.h>
#include <algorithm>
#include <vector>


/////////////////////////////////////////////////////////////////
// ProbabilisticModel::LinearComputeAlignment()
//
// Computes an alignment based on the given posterior matrix using
// a push-based DP approach (like Viterbi). This finds the maximum
// summing path through the posterior matrix. The final alignment
// is returned as a pair consisting of:
//    (1) a string (e.g., XXXBBXXXBBBBBBYYYYBBB) where X's and
//        denote insertions in one of the two sequences and
//        B's denote that both sequences are present (i.e.
//        matches).
//    (2) a value_type indicating the sum achieved
/////////////////////////////////////////////////////////////////

pair<string *, value_type> ProbabilisticModel::LinearComputeAlignment(int hmmBeam, int seq1Length, int seq2Length, const unordered_map<int, value_type>* posterior) {
    vector<unordered_map<pair<int, int>, MEAState, PairHash>> beam(seq1Length + seq2Length + 3);

    // Initialize
    beam[0][{0, 0}] = MEAState(0.0, '\0');

    for (int s = 0; s <= seq1Length + seq2Length; ++s) {
        // Beam pruning
        if (s % skip_beam_prune_ == 0 || s % skip_beam_prune_ == 1) beam_prune(beam[s], hmmBeam);
        // if (beam[s].size() > hmmBeam) beam_prune(beam[s], hmmBeam);


        // Process each state in current step
        for (auto& item : beam[s]) {
            int i = item.first.first;
            int j = item.first.second;
            value_type score = item.second.score;

            // ALN
            if (i + 1 <= seq1Length && j + 1 <= seq2Length) {
                int next_i = i + 1;
                int next_j = j + 1;
                auto row_it = posterior[i].find(j);
                if (row_it != posterior[i].end()) {
                    value_type new_score = score + row_it->second;
                    MEAState& next_state = beam[s + 2][{next_i, next_j}];
                    // Update always: it is always the best
                    next_state.score = new_score;
                    next_state.traceback = 'D'; // Diagonal (match)
                }
            }

            // INS1
            if (i + 1 <= seq1Length) {
                MEAState& next_state = beam[s + 1][{i + 1, j}];

                // Update if better
                if (score > next_state.score) {
                    next_state.score = score;
                    next_state.traceback = 'U'; // Up (ins1)
                }
            }

            // INS2
            if (j + 1 <= seq2Length) {
                MEAState& next_state = beam[s + 1][{i, j + 1}];

                // Update if better
                if (score > next_state.score) {
                    next_state.score = score;
                    next_state.traceback = 'L'; // Left (ins2)
                }
            }
        }
    }

    // Find the best score at the final position
    value_type total = 0.0;
    unsigned final_step = seq1Length + seq2Length;

    if (beam[final_step].find({seq1Length, seq2Length}) != beam[final_step].end()) {
        total = beam[final_step][{seq1Length, seq2Length}].score;
    }

    // Traceback to reconstruct alignment
    string* alignment = new string();
    unsigned current_i = seq1Length;
    unsigned current_j = seq2Length;
    unsigned current_step = final_step;

    while (current_i > 0 || current_j > 0) {
        auto state_it = beam[current_step].find({current_i, current_j});
        if (state_it == beam[current_step].end()) break;

        char ch = state_it->second.traceback;
        switch (ch) {
        case 'D': // Diagonal (match)
            alignment->push_back('B');
            current_i--;
            current_j--;
            current_step -= 2;
            break;
        case 'U': // Up (ins1)
            alignment->push_back('X');
            current_i--;
            current_step -= 1;
            break;
        case 'L': // Left (ins2)
            alignment->push_back('Y');
            current_j--;
            current_step -= 1;
            break;
        default:
            current_step = 0; // safety
            current_i = current_j = 0;
            break;
        }
    }

    // Reverse the alignment string
    reverse(alignment->begin(), alignment->end());

    return make_pair(alignment, total);
}

void ProbabilisticModel::LinearConsistencyTransform(int seq1Length, unordered_map<int, value_type>* &xz_posterior, unordered_map<int, value_type>* &zy_posterior, unordered_map<int, value_type>* &new_xy_posterior){
    for(int i = 0; i < seq1Length; i++){
        for(auto &xz_cand : xz_posterior[i]){
            int k = xz_cand.first;

            for(auto &zy_cand : zy_posterior[k]){
                int j = zy_cand.first;
                new_xy_posterior[i][j] += xz_cand.second * zy_cand.second;
            }
        }
    }
}

vector<vector<unordered_map<int, value_type>*>> ProbabilisticModel::LinearMultiConsistencyTransform(MultiSeq &sequences, vector<vector<unordered_map<int, value_type>*>> &posterior){
    const int numSeqs = sequences.size();
    
    // allocate space for the new posterior matrix
    vector<vector<unordered_map<int, value_type>*>> new_posterior(numSeqs);
    for (int i = 0; i < numSeqs; i++){
        new_posterior[i].resize(numSeqs);
        for (int j = 0; j < numSeqs; j++){
            if (i == j) continue;
            new_posterior[i][j] = new unordered_map<int, value_type>[sequences.at(i).length()];
        }
    }

    // For every pair of sequences
    for (int x = 0; x < numSeqs; x++){
        for (int y = x+1; y < numSeqs; y++){
            Sequence seq1 = sequences.at(x);
            Sequence seq2 = sequences.at(y);

            const int seq1Length = seq1.length();
            const int seq2Length = seq2.length();

            // allocate space for temporary results
            unordered_map<int, value_type>* transformation = new unordered_map<int, value_type>[seq1Length];

            // Get the original alignment result
            unordered_map<int, value_type>* &original = posterior[x][y];

            // Contribution from the summation where z = x and z = y
            for (int i = 0; i < seq1Length; i++){
                for(auto &item : original[i]){
                    int j = item.first;
                    transformation[i][j] = 2 * item.second;
                }
            }

            // Contribution from all other sequences
            for (int z = 0; z < numSeqs; z++) {
                if (z == x || z == y) continue;
                LinearConsistencyTransform(seq1Length, posterior[x][z], posterior[z][y], transformation);
            }

            // Renormalization
            for (int i = 0; i < seq1Length; i++){
                for(auto &item : transformation[i]){
                    int j = item.first;
                    transformation[i][j] /= numSeqs;
                }
            }

            // // Renormalization and mask out smaller values
            // for (int i = 0; i < seq1Length; ++i) {
            //     auto& cons_trans_i = transformation[i];

            //     for (auto it = cons_trans_i.begin(); it != cons_trans_i.end(); ) {
            //         int j = it->first;
            //         value_type prob = it->second;
            //         prob /= numSeqs;

            //         if (prob < 0.01) {
            //             it = cons_trans_i.erase(it); // returns iterator to next
            //         } else {
            //             posterior[x][y][i][j] = prob;
            //             posterior[y][x][j][i] = prob;
            //             ++it;
            //         }
            //     }
            // }

            // Mask out positions not originally in the posterior matrix
            for (int i = 0; i < seq1Length; i++){
                for(auto &item : original[i]){
                    int j = item.first;
                    if (transformation[i].find(j) == transformation[i].end()) continue; // N.B.
                    if (transformation[i][j] >= 0.01){
                        new_posterior[x][y][i][j] = transformation[i][j];
                        new_posterior[y][x][j][i] = transformation[i][j];
                        // cout << i << " " << j << " " << k << " " << l << " " << temp_pair_CT[k][l].value << endl;
                    }
                }
            }
            delete[] transformation;
        }
    }

    return new_posterior;
    // return posterior;
}

/// Equivalent to BuildPosterior from the Probcons 
/////////////////////////////////////////////////////////////////
// ProbabilisticModel::BuildPosterior()
//
// Builds a posterior probability matrix needed to align a pair
// of alignments.  Mathematically, the returned matrix M is
// defined as follows:
//    M[i,j] =     sum          sum      f(s,t,i,j)
//             s in align1  t in align2
// where
//                  [  P(s[i'] <--> t[j'])
//                  [       if s[i'] is a letter in the ith column of align1 and
//                  [          t[j'] it a letter in the jth column of align2
//    f(s,t,i,j) =  [
//                  [  0    otherwise
//
/////////////////////////////////////////////////////////////////
unordered_map<int, value_type>* ProbabilisticModel::LinearMultiAlnResults(MultiSeq* align1, MultiSeq* align2, const vector<vector<unordered_map<int, value_type>*>> &posterior, value_type cutoff) const {
    const int seq1Length = align1->at(0).length();
    const int seq2Length = align2->at(0).length();
    
    unordered_map<int, value_type>* transformation = new unordered_map<int, value_type>[seq1Length];
    for (int i = 0; i < align1->size(); i++){
        int first = align1->at(i).id;
        vector<int>* mapping1 = align1->at(i).get_mapping();

        // Loops through align2
        for (int j = 0; j < align2->size(); j++){
            int second = align2->at(j).id;
            vector<int>* mapping2 = align2->at(j).get_mapping();

            if(first < second){
                unordered_map<int, value_type>* original = posterior[first][second];

                int seq1len = mapping1->size();
                for (int ii = 0; ii < seq1len; ii++){
                    int ibase = (*mapping1)[ii];

                    for (auto &item : original[ii]) {
                        int jbase = (*mapping2)[item.first];
                        if (item.second < cutoff) continue;

                        transformation[ibase][jbase] += item.second;
                    }
                }
            } else {
                unordered_map<int, value_type>* original = posterior[second][first];

                int seq2len = mapping2->size();
                for (int jj = 0; jj < seq2len; jj++){
                    int jbase = (*mapping2)[jj];

                    for (auto &item : original[jj]) {
                        int ibase = (*mapping1)[item.first];
                        if (item.second < cutoff) continue;

                        transformation[ibase][jbase] += item.second;
                    }
                }
            }

            delete mapping2;
        }

        delete mapping1;
    }
    return transformation;
}

MultiSeq* ProbabilisticModel::LinearAlignAlignments (MultiSeq* align1, MultiSeq* align2, const vector<vector<unordered_map<int, value_type>*>> &posterior, int hmmBeam){

    // Choose the alignment routine depending on the "cosmetic" gap penalties used
    unordered_map<int, value_type>* posterior_matrix = LinearMultiAlnResults(align1, align2, posterior, 0.01f);

    pair<string*, value_type> alignment = LinearComputeAlignment(hmmBeam, align1->at(0).length(), align2->at(0).length(), posterior_matrix);
    delete[] posterior_matrix;

    // Build final alignment
    MultiSeq *result = new MultiSeq();
    for (int i = 0; i < align1->size(); i++)
        result->add_sequence (  *(align1->at(i).add_gaps(alignment.first, 'X')) );
    for (int i = 0; i < align2->size(); i++)
        result->add_sequence (  *(align2->at(i).add_gaps(alignment.first, 'Y')) );

    delete alignment.first;

    return result;
}

MultiSeq* ProbabilisticModel::LinearProcessTree (const TreeNode *tree, MultiSeq* sequences, const vector<vector<unordered_map<int, value_type>*>> &posterior, int hmmBeam){
    MultiSeq *result;

    // Check if this is an internal node of the alignment tree
    if (tree->GetSequenceLabel() == -1){
        MultiSeq *alignLeft = LinearProcessTree (tree->GetLeftChild(), sequences, posterior, hmmBeam);
        MultiSeq *alignRight = LinearProcessTree (tree->GetRightChild(), sequences, posterior, hmmBeam);

        assert (alignLeft);
        assert (alignRight);
        
        result = LinearAlignAlignments (alignLeft, alignRight, posterior, hmmBeam);
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

void ProbabilisticModel::LinearDoIterativeRefinement (const vector<vector<unordered_map<int, value_type>*>> &posterior, MultiSeq* alignment, int i, int hmmBeam){
    set<int> groupOne, groupTwo;

    randomnumber rn;
    rn.seed(1234+i);

    // Create two separate groups
    for (int i = 0; i < alignment->size(); i++){
        int x = rn.roll_int(1,10);

        if (x % 2)
            groupOne.insert (i);
        else
            groupTwo.insert (i);
    }

    if (groupOne.empty() || groupTwo.empty()) return;

    // Project into the two groups
    MultiSeq *groupOneSeqs = alignment->Project (groupOne); assert (groupOneSeqs);
    MultiSeq *groupTwoSeqs = alignment->Project (groupTwo); assert (groupTwoSeqs);

    // Realign
    alignment = LinearAlignAlignments (groupOneSeqs, groupTwoSeqs, posterior, hmmBeam); 

    delete groupOneSeqs;
    delete groupTwoSeqs;
}

MultiSeq* ProbabilisticModel::LinearComputeFinalAlignment (const TreeNode *tree, MultiSeq* sequences, const vector<vector<unordered_map<int, value_type>*>> &posterior, int hmmBeam){
    MultiSeq* alignment = LinearProcessTree (tree, sequences, posterior, hmmBeam);

    // alignment->print();

    // Iterative refinement
    for (int i = 0; i < num_iterative_refinement_reps_; i++)
        LinearDoIterativeRefinement (posterior, alignment, i, hmmBeam);

    return alignment;
}

