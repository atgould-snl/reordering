#pragma once
#include "common.h"
#include "branchAndBound.h"
#include <iostream>
#include <optional>
#include <vector>


using BlockNormsViewType = Kokkos::View<double**>; //GaussSeidelBlockReorderPreconditionerFactory::BlockNormsViewType;
// using PermutationScoreType = std::pair<PermutationType, double>;
// using PermutationType = std::map<int, int>;

// Function prototypes
//std::map<int,int>  doNothing(int a);

std::map<int,int> branchAndBoundLossMinimizationBlockPermutation(BlockNormsViewType& blockNorms, bool verbose=false); // Function interface, this makes the reordering object, runs it, and extracts the map

struct PartialOrder
{
    // FIELDS
    std::vector<std::optional<int>> blocks; // List of lists of doubles
    double loss;                        // Loss value
    int nBlocks;
    int nMembers;
    
    // CONSTRUCTION
    PartialOrder(int initialn=0, double initialLoss=0.0, int initialnBlock=0, int initialnMember=0) : blocks(std::vector<std::optional<int>>(initialn)), loss(initialLoss) , nBlocks(initialnBlock), nMembers(initialnMember){} 

    // PRINTING
    void print() const {
        std::cout << "Order Details:" << std::endl;
        std::cout << "Loss: " << loss << std::endl;
        std::cout << "Blocks IDs by member:" << std::endl;
        if (blocks.empty()) {
            std::cout << "Empty" << std::endl;
        } else {
            int member=0;
            for (const auto& block : blocks) {
                std::cout << member << ": " << block.value();
                if (member < blocks.size() - 1) {
                    std::cout << std::endl; // Print a comma for all but the last element
                }
                member++;
            }
            std::cout << std::endl;
        }
    }

    // TO MAP
    std::map<int,int> getMap(){ // Convert to map. Keys are ints from 0, values are 
        std::map<int,int> orderMap;
        for (int member=0; member<blocks.size(); member++){
            assert(blocks[member]); // All members should map to a block
            orderMap[member]=blocks[member].value();
        }
        return orderMap;
    }
};

class BranchAndBoundPermutationSearch
{
public:
    BranchAndBoundPermutationSearch(BlockNormsViewType& blockNorms);
    BlockNormsViewType& blockNorms;
    int n=0;
    bool allowMerge=false;
    bool allowBranchCutting=false;

    // Results tracking
    PartialOrder minLossOrder = PartialOrder(n, std::numeric_limits<double>::max());
    int numInternalNodes = 0;
    int numLeafNodes = 0;
    int terminateSearchAfterNumLeafNodes = 1E6;
    bool earlyTerminationWarning=false;

    // User interface methods
    void solve(); // After problem is set up, calls the private recursive function
    void solveExhuastive();

    private:
        // Methods
        void addToOrder(const PartialOrder& order, bool prevBest=true);
        void insertMember(PartialOrder& order, const int newMember, const int intoBlock, const bool merge);
        void updateLoss(PartialOrder& order);
        bool compareCandidateOrder(PartialOrder& newOrder); // Compares orders and replaces if the new one is lower cost

        // Internal vars
        std::vector<int> memberPresort;
        std::vector<double> remainingLossLowerBound;
};