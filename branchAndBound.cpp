#include "common.h"
#include "branchAndBound.h"

// FUNCTION DEFINITION
std::map<int,int> branchAndBoundLossMinimizationBlockPermutation(BlockNormsViewType& blockNorms, bool verbose){
    // Create the object
    branchAndBoundPermutationSearch bb = branchAndBoundPermutationSearch(blockNorms);

    // Run reordering
    bb.solveExhuastive();

    // Print results
    if (verbose){
        bb.minLossOrder.print();
    }


    // Convert to map. Keys are ints from 0, values are 
    std::map<int,int> orderMap;
    int blockIdx=0;
    for (const auto& block : bb.minLossOrder.list) {
        for (const int& member : block){
            orderMap[member]=blockIdx;
        }
        blockIdx++;
    }
    return orderMap;
}

// CLASS DEFINITION
branchAndBoundPermutationSearch::branchAndBoundPermutationSearch(BlockNormsViewType& blockNorms)
    : blockNorms(blockNorms)
{
    n = blockNorms.extent(0);
}

void branchAndBoundPermutationSearch::solveExhuastive(){
    // Loop through all permutations

    // Start with a vector of counting numbers
    std::vector<int> orderVector(n);
    std::iota(orderVector.begin(), orderVector.end(), 0);

    do {
        // Turn the vector into a partial order in the list of lists format
        partialOrder newOrder;
        for (const int& member : orderVector) {
            std::list<int> block; // Create a new sublist
            block.push_back(member);    // Add the current number to the sublist
            newOrder.list.push_back(block);
        }

        // Score it
        updateLoss(newOrder);
        // Compare it
        compareCandidateOrder(newOrder);

    } while (std::next_permutation(orderVector.begin(), orderVector.end()));
    return;
}

void branchAndBoundPermutationSearch::updateLoss(partialOrder& order){
    order.loss=0; // Starting from scratch

    // Loop through each block (as if they were getting inserted sort of)
    int blockIdx=0;
    for (const auto& block : order.list){ // Each block
        // Add loss associated with merging

        // Add the cost associated with each other member
        for (const int& member : block){ // Each member
            int otherBlockIdx=0;
            for (const auto& otherBlock : order.list){ // For each other block
                for (const int& otherMember : otherBlock){ // Each other member
                    double v=0.0;
                    if (otherBlockIdx < blockIdx){ // In lower triangle?
                        v = blockNorms(member, otherMember); // Look to the left of (blockIdx,blockIdx)
                        order.loss += v * v;
                    }
                }
                otherBlockIdx++;
            }
        }
        blockIdx++;
    }
}

partialOrder branchAndBoundPermutationSearch::insertIntoOrder(const partialOrder& oldOrder, const int newMember, const int intoBlock, const bool merge){
    // Make insertion
    partialOrder newOrder(oldOrder.loss);
    int blockIdx=0;
    for (const auto& blockList : oldOrder.list) {
        std::list<int> innerListCopy = blockList;
        if (blockIdx == intoBlock){
            if (merge){
                innerListCopy.push_back(newMember); // Add an extra member to the existing block
            }
            else{
                newOrder.list.push_back(std::list<int>(1,newMember)); // Add an extra block
            }
        }
        newOrder.list.push_back(innerListCopy); 
        blockIdx++;
    }

    // Add to score based on insertion
    int otherBlockIdx=0;
    for (const auto& otherBlockList : newOrder.list){ // For each other block
        for (const int& otherMember : otherBlockList){ // Each other member
            int v=0;
            if (otherBlockIdx < intoBlock){
                v = blockNorms(newMember, otherMember); // To the left of (blockIdx, blockIdx)
            }
            else if (otherBlockIdx > intoBlock){
                v = blockNorms(otherMember, newMember); // To the below (blockIdx, blockIdx)
            }
            newOrder.loss += v * v;
        }
        otherBlockIdx++;
    }
    return newOrder;
}

bool branchAndBoundPermutationSearch::compareCandidateOrder(partialOrder& newOrder){
    if (newOrder.loss < minLossOrder.loss){
        minLossOrder=newOrder;
        return true;
    }
    return false;
}