#include "common.h"
#include "branchAndBound.h"
#include <algorithm>
#include <cassert>


// FUNCTION DEFINITION
std::map<int,int> branchAndBoundLossMinimizationBlockPermutation(BlockNormsViewType& blockNorms, bool verbose){
    // Create the object
    BranchAndBoundPermutationSearch bb = BranchAndBoundPermutationSearch(blockNorms);
    bb.allowBranchCutting=true;

    // Run reordering
    bb.solve();

    // Print results
    if (verbose){
        bb.minLossOrder.print();
    }

    return bb.minLossOrder.getMap();
}

// CLASS DEFINITION
BranchAndBoundPermutationSearch::BranchAndBoundPermutationSearch(BlockNormsViewType& blockNorms)
    : blockNorms(blockNorms)
{
    n = blockNorms.extent(0);
}

void BranchAndBoundPermutationSearch::solve(){
    // Precompute the lower bounds
    // Find the presort order to add members
    auto memberImportance=std::vector<double>(n,0);
    for (int member=0; member<n; member++){
        for (int otherMember=0; otherMember<n; otherMember++){
            memberImportance[member]+=std::abs(blockNorms(member,otherMember)-blockNorms(otherMember,member));
        }
    }
    // Argsort
    memberPresort=std::vector<int>(n);
    std::iota(memberPresort.begin(), memberPresort.end(), 0);
    std::sort(memberPresort.begin(),memberPresort.end(), [&memberImportance](const int memberA, const int memberB){return memberImportance[memberA] > memberImportance[memberB];}); // More important members get a true and go first

    // Lower bound based on add order
    remainingLossLowerBound=std::vector<double>(n);
    double totalLossLowerBound=0;
    for (int memberIdx=n-1; memberIdx>=0; memberIdx--){ // Loop backwards for cumsum
        // What is the min loss relative to all of the members before it
        for (int otherMemberIdx=0; otherMemberIdx<memberIdx; otherMemberIdx++){ // Min cost due to interaction w all previous members on list
            // Find the min loss of adding it
            double additionalLoss = std::min(blockNorms(memberPresort[memberIdx],memberPresort[otherMemberIdx]),blockNorms(memberPresort[otherMemberIdx],memberPresort[memberIdx]));
            if (allowMerge){
                // Test option to merge
                // TESTING ONLY TODO
                //assert(false);
                double mergeLoss=0;
                additionalLoss = std::min(additionalLoss,mergeLoss);
            }
            totalLossLowerBound+=additionalLoss;
        }
        // Save running total
        remainingLossLowerBound[memberIdx]=totalLossLowerBound;
    }

    // Make recursive call to add to order
    PartialOrder base = PartialOrder(n);
    addToOrder(base, true);
    if (earlyTerminationWarning){
        std::cout << "WARNING: search terminated early after checking " << terminateSearchAfterNumLeafNodes << " orders." << std::endl;
    }
}

void BranchAndBoundPermutationSearch::addToOrder(const PartialOrder& order, bool prevBest){
    numInternalNodes++;
    // Allocate for branches
    int nBranches=order.nBlocks+1; // Insert between existing blocks and before and after
    bool preventFrontBackBranches = allowBranchCutting && !prevBest;
    nBranches+=preventFrontBackBranches ? -2 : 0; // No need for branches for new member inserted before or after all previous if it wasn't the best performer previously
    nBranches+=allowMerge ? order.nBlocks : 0; // Additional branches to allow for merging
    std::vector<PartialOrder> branchOrders = std::vector<PartialOrder>(nBranches,PartialOrder(order)); // Fill with copies. Preallocated for speed

    // Insert into branches
    int branchIdx=0;
    for (int addAsBlock=(preventFrontBackBranches ? 1 : 0); addAsBlock < order.nBlocks+(preventFrontBackBranches ? 0 : 1); addAsBlock++){
        insertMember(branchOrders[branchIdx],memberPresort[order.nMembers],addAsBlock,false);
        branchIdx++;
    }
    if (allowMerge){
        for (int addIntoBlock=0; addIntoBlock < order.nBlocks; addIntoBlock++){
            insertMember(branchOrders[branchIdx],memberPresort[order.nMembers],addIntoBlock,true);
            branchIdx++;
        }
    }
    assert(branchIdx==nBranches); // Make sure we got them all
    // Sort branches
    std::sort(branchOrders.begin(), branchOrders.end(),[](const PartialOrder& branchA, const PartialOrder& branchB){return branchA.loss < branchB.loss;});

    // Are we at the bottom of the recursion?
    if (branchOrders[0].nMembers == n){
        compareCandidateOrder(branchOrders[0]);
        numLeafNodes+=nBranches; // All of these branches are terminal
        return;
    }

    // Loop through the sorted branches and make recursive call
    bool best=true;
    for (auto branchOrder : branchOrders){
        // Make the recursive call
        if (!allowBranchCutting || branchOrder.loss + remainingLossLowerBound[branchOrder.nMembers] < minLossOrder.loss){ // Test if branch should be cut
            if (numLeafNodes>terminateSearchAfterNumLeafNodes){ // TODO add time limit
                earlyTerminationWarning=true;
                return;
            }
            addToOrder(branchOrder,best);
        }
        best=false;
    }
}

void BranchAndBoundPermutationSearch::solveExhuastive(){
    // Loop through all permutations

    // Start with a vector of counting numbers
    PartialOrder order(n,0,n);
    std::iota(order.blocks.begin(), order.blocks.end(), 0);

    do {
        // Score
        updateLoss(order);
        // Compare
        compareCandidateOrder(order);
        numLeafNodes++;
        // Permute
    } while (std::next_permutation(order.blocks.begin(), order.blocks.end()));
    return;
}

void BranchAndBoundPermutationSearch::updateLoss(PartialOrder& order){
    order.loss=0; // Starting from scratch
    // Loop through potential vector
    for (int i=0; i<n; i++){ // Loop through rows
        for (int j=0; j<n; j++){ // Loop through cols
            if (order.blocks[i] && order.blocks[j]){
                if (order.blocks[i] > order.blocks[j]){
                    // Physics j with lower block number cannot impact physics i with higher block number, so add it to cost
                    double v = blockNorms(i, j); // Look to the left of ()
                    order.loss += v * v;
                }
            }
        }
    }
}

void BranchAndBoundPermutationSearch::insertMember(PartialOrder& order, const int newMember, const int intoBlock, const bool merge){
    // MAKE INSERTION
    // Update with new member
    assert(!order.blocks[newMember]); // Make sure this member has not already been assigned
    if (!merge){ // Bump higher blocks forward
        order.nBlocks++;
        for (int member=0; member<n; member++){
            if (!order.blocks[member]) continue;
            if (order.blocks[member] >= intoBlock){
                order.blocks[member].value()++; // Bumped up by new insertion
            }
        }
    }
    order.blocks[newMember]=intoBlock;
    order.nMembers++;

    // Update score due to new member
    if (merge){
        order.blocks[newMember]=intoBlock;
        // TODO update score due to merge
        //assert(false);
        double mergeLoss=0; // TESTING ONLY
        order.loss += mergeLoss;
    }
    for (int member=0; member<n; member++){ // Loop through members
        if (!order.blocks[member]) continue; // Has this member been assigned yet?
        double v=0; // Default, used within the same block
        if (order.blocks[newMember] > order.blocks[member]){
            v = blockNorms(newMember, member); // To the left
        }
        else if (order.blocks[newMember] < order.blocks[member]){
            v = blockNorms(member, newMember); // Below
        }
        order.loss += v * v;
    }
}

bool BranchAndBoundPermutationSearch::compareCandidateOrder(PartialOrder& newOrder){
    // Check that the order is complete
    assert(newOrder.nBlocks=n);
    for (auto block : newOrder.blocks){assert(block);}
    // Compare
    if (newOrder.loss < minLossOrder.loss){
        minLossOrder=newOrder;
        return true;
    }
    return false;
}