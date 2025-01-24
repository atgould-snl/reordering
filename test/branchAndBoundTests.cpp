#include <gtest/gtest.h>
#include "../common.h"
#include "../exhaustive.h"
#include "../branchAndBound.h"
#include "Kokkos_Core.hpp"

// Function to convert std::vector<std::vector<double>> to Kokkos::View<double**>
Kokkos::View<double**> vectorToKokkosView(const std::vector<std::vector<double>>& vec) {
    // Get the dimensions of the input vector
    size_t numRows = vec.size();
    size_t numCols = (numRows > 0) ? vec[0].size() : 0;

    // Create a Kokkos View with the same dimensions
    Kokkos::View<double**> kokkosView("kokkosView", numRows, numCols);

    // Copy data from std::vector to Kokkos View
    for (size_t i = 0; i < numRows; ++i) {
        for (size_t j = 0; j < numCols; ++j) {
            kokkosView(i, j) = vec[i][j];
        }
    }

    return kokkosView;
}

TEST(branchAndBoundTests, all){
    // Define a 4x4 matrix using a vector of vectors
    std::vector<std::vector<double>> T_vec_of_vec = {
        {2, 10.1, 1, 10.2},
        {1, 2, 1, 10.3},
        {10.4, 10.5, 2, 10.6},
        {1, 1, 1, 2}
    };
    // Expected result is 1 2 0 3
    //partialOrder expectedOrder(4,6,4);
    //std::vector<std::optional<int>> solnBlocks = {1,2,0,3};
        std::map<int, int> solnMap = {
        {0, 1},   // Eq 0 maps to block 1
        {1, 2},
        {2, 0},
        {3, 3}
    };
    double solnLoss=6;
    //expectedOrder.blocks=solnBlocks;

    // Load into kokkos
    Kokkos::initialize();
    Kokkos::View<double**> T = vectorToKokkosView(T_vec_of_vec);
    //Kokkos::View<double**> T = create_random_T(4);
    std::cout << "Test matrix:" << std::endl;
    print_matrix(T);
    //std::map<int,int> foundMap = branchAndBoundLossMinimizationBlockPermutation(T, true);

    branchAndBoundPermutationSearch bb = branchAndBoundPermutationSearch(T);
    bb.solveExhuastive();
    bb.minLossOrder.print();
    std::map<int,int> foundMap=bb.minLossOrder.getMap();
    EXPECT_TRUE(foundMap==solnMap);

    branchAndBoundPermutationSearch bb2 = branchAndBoundPermutationSearch(T);
    bb2.solve();
}