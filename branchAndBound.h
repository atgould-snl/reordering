#pragma once
#include "common.h"
#include "branchAndBound.h"
#include <iostream>
#include <optional>
#include <vector>


using BlockNormsViewType = Kokkos::View<double**>; //GaussSeidelBlockReorderPreconditionerFactory::BlockNormsViewType;

// Function prototypes
//std::map<int,int>  doNothing(int a);

std::map<int,int> branchAndBoundLossMinimizationBlockPermutation(BlockNormsViewType& blockNorms, bool verbose=false); // Function interface, this makes the reordering object, runs it, and extracts the map

struct partialOrder
{
    // FIELDS
    std::vector<std::optional<int>> blocks; // List of lists of doubles
    double loss;                        // Loss value
    double nBlocks;
    double nMembers;
    
    // CONSTRUCTION
    partialOrder(int initialn=0, double initialLoss=0.0, int initialnBlock=0, int initialnMember=0) : blocks(std::vector<std::optional<int>>(initialn)), loss(initialLoss) , nBlocks(initialnBlock), nMembers(initialnMember){} 
    // As copy
    partialOrder(const partialOrder& other) 
        : blocks(other.blocks), loss(other.loss), nBlocks(other.nBlocks), nMembers(other.nMembers) {}

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

class branchAndBoundPermutationSearch
{
private:
    //void add_to_order( , bool prevMin=true);
    //void 
    void addToOrder(const partialOrder& order, bool prevBest=true);
    void insertMember(partialOrder& order, const int newMember, const int intoBlock, const bool merge);
    void updateLoss(partialOrder& order);
    bool compareCandidateOrder(partialOrder& newOrder); // Compares orders and replaces if the new one is lower cost
    //double lowerTriangularLoss(BlockNormsViewType& blockNorms, partialOrder);
public:
    branchAndBoundPermutationSearch(BlockNormsViewType& blockNorms);
    ~branchAndBoundPermutationSearch()=default; // Destructor
    BlockNormsViewType& blockNorms;
    int n=0;
    bool allowMerge=false;
    bool allowBranchCutting=false;
    std::vector<int> memberPresort;
    std::vector<double> remainingLossLowerBound;

    // Results tracking
    partialOrder minLossOrder = partialOrder(n, std::numeric_limits<double>::max());
    //double minLoss = std::numeric_limits<double>::max();
    int numInternalNodes = 0;
    int numLeafNodes = 0;

    void solve(); // After problem is set up, calls the private recursive function
    void solveExhuastive();
    
};


//void update_score_vector_mat(std::vector<double>& T, int n, order& order_now);
//Kokkos::View<double**> create_random_T(int n);
//void print_matrix(const Kokkos::View<double**>& T);
//std::vector<double> kokkos_to_vector(const Kokkos::View<double**>& T);
//order find_optimal_order_exhaustive(Kokkos::View<double**>& T);
