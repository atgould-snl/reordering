#pragma once
#include "common.h"
#include "branchAndBound.h"


using BlockNormsViewType = Kokkos::View<double**>; //GaussSeidelBlockReorderPreconditionerFactory::BlockNormsViewType;

// Function prototypes
//std::map<int,int>  doNothing(int a);

std::map<int,int> branchAndBoundLossMinimizationBlockPermutation(BlockNormsViewType& blockNorms, bool verbose=false); // Function interface, this makes the reordering object, runs it, and extracts the map

struct partialOrder
{
    // FIELDS
    std::list<std::list<int>> list; // List of lists of doubles
    double loss;                        // Loss value

    // CONSTRUCTION
    partialOrder(double initialLoss=0.0) : loss(initialLoss) {} 
    // Copy constructor
    partialOrder(const partialOrder& other) 
        : loss(other.loss) { // Copy the loss value
        // Perform a deep copy of the list of lists
        for (const auto& block : other.list) {
            std::list<int> blockCopy = block;
            list.push_back(blockCopy); // Copy each inner list
        }
    }

    // PRINTING
    void print() const {
        std::cout << "Order Details:" << std::endl;
        std::cout << "Loss: " << loss << std::endl;
        std::cout << "Blocks: ";
        if (list.empty()) {
            std::cout << "Empty" << std::endl;
        } else {
            int blockIdx=0;
            for (const auto& block : list) {
                // Print members of the block with + between
                int memberIdx=0;
                for (const int& member : block){
                    std::cout << member;
                    if (memberIdx>0){
                        std::cout << '+';
                    }
                    memberIdx++;
                }
                if (blockIdx < list.size() - 1) {
                    std::cout << ", "; // Print a comma for all but the last element
                }
            }
            std::cout << std::endl;
        }
    }
};

class branchAndBoundPermutationSearch
{
private:
    //void add_to_order( , bool prevMin=true);
    //void 
    partialOrder insertIntoOrder(const partialOrder& oldOrder, const int newMember, const int intoBlock, const bool merge);
    void updateLoss(partialOrder& order);
    bool compareCandidateOrder(partialOrder& newOrder); // Compares orders and replaces if the new one is lower cost
    //double lowerTriangularLoss(BlockNormsViewType& blockNorms, partialOrder);
public:
    branchAndBoundPermutationSearch(BlockNormsViewType& blockNorms);
    ~branchAndBoundPermutationSearch()=default; // Destructor
    BlockNormsViewType& blockNorms;
    int n;

    // Results tracking
    partialOrder minLossOrder = partialOrder(std::numeric_limits<double>::max());
    //double minLoss = std::numeric_limits<double>::max();
    int numInternalNodes = 0;
    int numLeafNodes = 0;

    //void solve(); // After problem is set up, calls the private recursive function
    void solveExhuastive();
    
};


//void update_score_vector_mat(std::vector<double>& T, int n, order& order_now);
//Kokkos::View<double**> create_random_T(int n);
//void print_matrix(const Kokkos::View<double**>& T);
//std::vector<double> kokkos_to_vector(const Kokkos::View<double**>& T);
//order find_optimal_order_exhaustive(Kokkos::View<double**>& T);
