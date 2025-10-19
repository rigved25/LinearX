#include "linearx/turbofold/utils/GuideTree.h"
#include <string>
#include <list>
#include <stdio.h>
#include <iostream>
#include <vector>
#include <utility>
#include <algorithm>

//! TreeNode Class
/*!
    The TreeNode Class provides unit for representing an alignment tree. 
*/

//! Constructor for a tree node. 
//! Note that sequenceLabel = -1 implies that the current node is not a leaf in the tree.
TreeNode::TreeNode(int sequenceLabel) : sequenceLabel(sequenceLabel), left(NULL), right(NULL), parent(NULL) {
    assert(sequenceLabel >= -1);
}

//! Destructor for a tree node.  Recursively deletes all children.
TreeNode::~TreeNode() {
    if(left){
        delete left; 
        left = NULL;
    }
    if(right){
        delete right;
        right = NULL;
    }
    parent = NULL;
}

// Getters
int TreeNode::GetSequenceLabel() const { return sequenceLabel; }
TreeNode * TreeNode::GetLeftChild() const { return left; }
TreeNode * TreeNode::GetRightChild() const { return right; }
TreeNode * TreeNode::GetParent() const { return parent; }

// Setters
void TreeNode::SetSequenceLabel(int sequenceLabel){ 
    this->sequenceLabel = sequenceLabel; 
    assert (sequenceLabel >= -1); 
}
void TreeNode::SetLeftChild(TreeNode *left){ this->left = left; }
void TreeNode::SetRightChild(TreeNode *right){ this->right = right; }
void TreeNode::SetParent(TreeNode *parent){ this->parent = parent; }

//! Computes a guide tree based on a given distance matrix. 
TreeNode * TreeNode::ComputeTree(const vector<vector<value_type> > &distMatrix){

    // Number of sequences
    int numSeqs = distMatrix.size();
    // A copy of distance matrix.
    vector<vector<value_type> > distances(numSeqs, vector<value_type> (numSeqs));
    vector<TreeNode *> nodes (numSeqs, NULL);
    // valid[i] tells whether or not the ith nodes in the distances and nodes array are valid.
    vector<int> valid (numSeqs, 1);

    // initialization: make a copy of the distance matrix
    for (int i = 0; i < numSeqs; i++)
        for (int j = 0; j < numSeqs; j++)
            distances[i][j] = distMatrix[i][j];

    // initialization: create all the leaf nodes
    for (int i = 0; i < numSeqs; i++){
        nodes[i] = new TreeNode (i);
        assert (nodes[i]);
    }

    // repeat until only a single node left
    for (int numNodesLeft = numSeqs; numNodesLeft > 1; numNodesLeft--){
        value_type bestProb = -1;
        pair<int,int> bestPair;

        // find the closest pair
        for (int i = 0; i < numSeqs; i++) if (valid[i]){
            for (int j = i+1; j < numSeqs; j++) if (valid[j]){
                if (distances[i][j] > bestProb){
                    bestProb = distances[i][j];
                    bestPair = make_pair(i, j);
                }
            }
        }

        // merge the closest pair
        cerr << "[GuideTree] merging: " << bestPair.first << " and " << bestPair.second
             << " (distance=" << bestProb << ")" << endl;
        TreeNode *newParent = new TreeNode (-1);
        newParent->SetLeftChild (nodes[bestPair.first]);
        newParent->SetRightChild (nodes[bestPair.second]);
        nodes[bestPair.first]->SetParent (newParent);
        nodes[bestPair.second]->SetParent (newParent);
        nodes[bestPair.first] = newParent;
        nodes[bestPair.second] = NULL;

        // now update the distance matrix
        for (int i = 0; i < numSeqs; i++) if (valid[i]){
            distances[bestPair.first][i] = distances[i][bestPair.first]
                = (distances[i][bestPair.first] + distances[i][bestPair.second]) * bestProb / 2;
        }

        // mark the second node entry as no longer valid
        valid[bestPair.second] = 0;
    }

    assert (nodes[0]);
    return nodes[0];
}

//! Computes a guide tree based on a given distance matrix using a max Fibonacci heap.
TreeNode * TreeNode::FastComputeTree(const vector<vector<value_type> > &distMatrix){
    using std::pair;
    using std::make_pair;
    using std::cerr;
    using std::endl;

    struct PairEntry {
        value_type distance;
        int a;
        int b;
    };
    struct PairCompare {
        bool operator()(const PairEntry &lhs, const PairEntry &rhs) const {
            // Max-heap by distance
            return lhs.distance < rhs.distance;
        }
    };
    typedef boost::heap::fibonacci_heap<PairEntry, boost::heap::compare<PairCompare> > PairHeap;
    typedef PairHeap::handle_type PairHandle;

    const int numSeqs = static_cast<int>(distMatrix.size());
    if (numSeqs == 0) return NULL;

    // Flattened NxN distance matrix to simplify indexing and in-place updates.
    std::vector<value_type> distances1D(static_cast<size_t>(numSeqs) * static_cast<size_t>(numSeqs), static_cast<value_type>(0));
    auto idx = [numSeqs](int i, int j) -> size_t {
        return static_cast<size_t>(i) * static_cast<size_t>(numSeqs) + static_cast<size_t>(j);
    };
    for (int i = 0; i < numSeqs; ++i) {
        for (int j = 0; j < numSeqs; ++j) {
            distances1D[idx(i, j)] = distMatrix[i][j];
        }
    }

    // Tree nodes and active flags.
    std::vector<TreeNode*> nodes(numSeqs, NULL);
    std::vector<int> isActive(numSeqs, 1);
    for (int i = 0; i < numSeqs; ++i) {
        nodes[i] = new TreeNode(i);
        assert(nodes[i]);
    }

    // Heap of all pairs with handles to support key updates.
    PairHeap heap;
    std::vector<char> handleExists(static_cast<size_t>(numSeqs) * static_cast<size_t>(numSeqs), 0);
    std::vector<PairHandle> handleOf(static_cast<size_t>(numSeqs) * static_cast<size_t>(numSeqs));
    auto canon = [numSeqs](int a, int b) -> size_t {
        if (a > b) std::swap(a, b);
        return static_cast<size_t>(a) * static_cast<size_t>(numSeqs) + static_cast<size_t>(b);
    };

    for (int i = 0; i < numSeqs; ++i) {
        for (int j = i + 1; j < numSeqs; ++j) {
            PairEntry entry{ distances1D[idx(i, j)], i, j };
            PairHandle h = heap.push(entry);
            handleOf[canon(i, j)] = h;
            handleExists[canon(i, j)] = 1;
        }
    }

    for (int numNodesLeft = numSeqs; numNodesLeft > 1; --numNodesLeft) {
        // Pop until we find a valid active pair
        PairEntry top;
        bool found = false;
        while (!heap.empty()) {
            top = heap.top();
            heap.pop();
            if (isActive[top.a] && isActive[top.b]) {
                found = true;
                break;
            }
        }
        assert(found && "Heap underflow while building guide tree");
        const int i = top.a;
        const int j = top.b;
        const value_type bestProb = distances1D[idx(i, j)];

        // Merge i and j under a new parent; keep merged cluster id as i, deactivate j
        cerr << "[GuideTree] merging: " << i << " and " << j
             << " (distance=" << bestProb << ")" << endl;
        TreeNode *newParent = new TreeNode(-1);
        newParent->SetLeftChild(nodes[i]);
        newParent->SetRightChild(nodes[j]);
        nodes[i]->SetParent(newParent);
        nodes[j]->SetParent(newParent);
        nodes[i] = newParent;
        nodes[j] = NULL;

        // For all active k, update distance(i,k) using the original formula
        for (int k = 0; k < numSeqs; ++k) {
            if (!isActive[k] || k == i || k == j) continue;
            const value_type dik = distances1D[idx(i, k)];
            const value_type djk = distances1D[idx(j, k)];
            const value_type newDist = (dik + djk) * bestProb / static_cast<value_type>(2);

            // Update flattened matrix symmetrically
            distances1D[idx(i, k)] = newDist;
            distances1D[idx(k, i)] = newDist;

            // Update heap key for pair (min(i,k), max(i,k))
            size_t keyIK = canon(i, k);
            if (handleExists[keyIK]) {
                PairHandle &h = handleOf[keyIK];
                (*h).distance = newDist;
                heap.update(h);
            } else {
                // Should not happen, but be robust
                int a = std::min(i, k);
                int b = std::max(i, k);
                PairEntry updated{ newDist, a, b };
                PairHandle h = heap.push(updated);
                handleOf[keyIK] = h;
                handleExists[keyIK] = 1;
            }

            // Erase pairs that involve the deactivated node j
            size_t keyJK = canon(j, k);
            if (handleExists[keyJK]) {
                heap.erase(handleOf[keyJK]);
                handleExists[keyJK] = 0;
            }
        }

        // Remove handle for (i,j) if it still exists (it was popped)
        handleExists[canon(i, j)] = 0;
        isActive[j] = 0;
    }

    assert (nodes[0]);
    return nodes[0];
}

