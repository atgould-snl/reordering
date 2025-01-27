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

Kokkos::View<double**> getRandomT(int n, double maxExp = 3.) {
    // Create a Kokkos View with the same dimensions
    Kokkos::View<double**> kokkosView("kokkosView", n, n);

    // Seed the random number generator
    std::random_device rd;  // Obtain a random number from hardware
    std::mt19937 gen(rd());  // Seed the generator
    std::uniform_real_distribution<> dis(0.0, maxExp);  // Define the range

    // Fill the Kokkos View with random values
    for (size_t i = 0; i < n; ++i) {
        for (size_t j = 0; j < n; ++j) {
            // Generate a random exponent between 1 and maxExp
            double exp=dis(gen);
            // Generate a random number in the range [10^1, 10^maxExp]
            kokkosView(i, j) = std::pow(10, exp);
        }
    }

    // Ensure the matrix has highest value (1) on diag
    for (size_t i = 0; i < n; ++i) {
        // Find the maximum value in the row
        double maxVal = kokkosView(i, 0);
        size_t maxIndex = 0;

        for (size_t j = 1; j < n; ++j) {
            if (kokkosView(i, j) > maxVal) {
                maxVal = kokkosView(i, j);
                maxIndex = j;
            }
        }

        // Swap the maximum value with the diagonal element
        if (maxIndex != i) {
            std::swap(kokkosView(i, i), kokkosView(i, maxIndex));
        }
    }

    // Normalize the matrix based on the diagonal
    for (size_t i = 0; i < n; ++i) {
        double diagonalValue = kokkosView(i, i);
        if (diagonalValue != 0) { // Avoid division by zero
            for (size_t j = 0; j < n; ++j) {
                kokkosView(i, j) /= diagonalValue;
            }
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

    // EXHAUSTIVE SOLN
    branchAndBoundPermutationSearch exhaustiveObj = branchAndBoundPermutationSearch(T);
    exhaustiveObj.solveExhuastive();
    exhaustiveObj.minLossOrder.print();
    std::map<int,int> exhuastiveMap=exhaustiveObj.minLossOrder.getMap();
    EXPECT_EQ(exhuastiveMap,solnMap);

    // BB WITH NO BRANCH CUTTING NO MERGING
    branchAndBoundPermutationSearch bbObj = branchAndBoundPermutationSearch(T);
    bbObj.solve();
    exhaustiveObj.minLossOrder.print();
    std::map<int,int> bbMap=exhaustiveObj.minLossOrder.getMap();
    EXPECT_EQ(bbMap,solnMap);
    EXPECT_NEAR(exhaustiveObj.minLossOrder.loss,bbObj.minLossOrder.loss,1E-8);
    // Check total scanned
    EXPECT_EQ(bbObj.n,4);
    EXPECT_EQ(bbObj.numLeafNodes,24);

    // BB WITH BRANCH CUTTING
    branchAndBoundPermutationSearch bbObj_cutting = branchAndBoundPermutationSearch(T);
    bbObj_cutting.allowBranchCutting=true;
    bbObj_cutting.solve();
    exhaustiveObj.minLossOrder.print();
    std::map<int,int> bbMap_cutting=exhaustiveObj.minLossOrder.getMap();
    EXPECT_EQ(bbMap_cutting,solnMap);
    EXPECT_NEAR(exhaustiveObj.minLossOrder.loss,bbObj_cutting.minLossOrder.loss,1E-8);
    // Check total scanned, should have gone down
    EXPECT_EQ(bbObj_cutting.n,4);
    EXPECT_LE(bbObj_cutting.numLeafNodes, 23);

    // Test in loop with random T
    std::cout << "Testing T random test: " << std::endl;
    T = getRandomT(6);
    print_matrix(T);
    EXPECT_EQ(T(3,3), 1.);
    EXPECT_GE(T(3,3), T(3,0)); // Diag dom
    EXPECT_GE(T(3,3), T(3,1)); // Diag dom
    EXPECT_GE(T(3,3), T(3,2)); // Diag dom

    // Test a bunch of T to make sure exhaustive and branch cutting give the same result
    for (int i=0; i<10; i++){
        T = getRandomT(7);
        branchAndBoundPermutationSearch exhaustiveObj = branchAndBoundPermutationSearch(T);
        exhaustiveObj.solveExhuastive();
        exhaustiveObj.minLossOrder.print();
        std::map<int,int> exhuastiveMap=exhaustiveObj.minLossOrder.getMap();

        branchAndBoundPermutationSearch bbObj = branchAndBoundPermutationSearch(T);
        bbObj.solve();
        exhaustiveObj.minLossOrder.print();
        std::map<int,int> bbMap=exhaustiveObj.minLossOrder.getMap();

        EXPECT_EQ(exhuastiveMap,bbMap);
        EXPECT_NEAR(exhaustiveObj.minLossOrder.loss,bbObj.minLossOrder.loss,1E-8);
    }
    // Check with merging. Want to total to be as estimated from formula (see python).
    // 3 -> 13   7 -> 47293
    T = getRandomT(7);
    branchAndBoundPermutationSearch bbObj_merging = branchAndBoundPermutationSearch(T);
    bbObj_merging.allowBranchCutting=false;
    bbObj_merging.allowMerge=true;
    bbObj_merging.solve();
    // Check total scanned
    EXPECT_EQ(bbObj_merging.numLeafNodes, 47293);

    for (int i=0; i<10; i++){
        T = getRandomT(20);
        easy_timer et1 = easy_timer();
        branchAndBoundPermutationSearch speedObj = branchAndBoundPermutationSearch(T);
        speedObj.allowBranchCutting=true;
        speedObj.allowMerge=false;
        speedObj.solve();
        //speedObj.minLossOrder.print();
        et1.print_time();
    }
}