//#include "exhaustive.h"
#include "common.h"
#include <vector>
// Ordering 3124 example:
// a is down, b is across
//  3                "what was in 3 is now in spot 1"
//  x  1
//  x  x  2
//  x  x  x  4
//   ->
// i is down j is across. The xs are the entries than need to be summed on the origional unsorted matrix
// 1     x
// x  2  x
//       3
// x  x  x  4
// Note: for entry in Aria the input deck, the argsort must be taken of the ordering. This would be input as 2 3 1 4
void update_score(Kokkos::View<double**>& T, order& order_now){
    order_now.loss=0; // Default is inf so must set to zero to start
    for (int b=0; b<order_now.vect.size(); b++){ // Old mat is i,j new mat is a,b, but new mat is never formed. //vect maps from a,b to i,j. We loop through a,b lower triangle 
        int j=order_now.vect[b]; // Column is determined in order of the list
        for (int a=b+1; a<order_now.vect.size(); a++){ // Now we loop through and sum and indicies AFTER the one were on 
            int i=order_now.vect[a];
            order_now.loss+=T(i, j);
        }
    }
}

void update_score_vector_mat(std::vector<double>& T_vector, int n, order& order_now){ // This version uses a vector instead. It is faster.
    order_now.loss=0; // Default is inf so must set to zero to start
    for (int b=0; b<order_now.vect.size(); b++){ // Old mat is i,j new mat is a,b, but new mat is never formed. //vect maps from a,b to i,j. We loop through a,b lower triangle 
        int j=order_now.vect[b]; // Column is determined in order of the list
        for (int a=b+1; a<order_now.vect.size(); a++){ // Now we loop through and sum and indicies AFTER the one were on 
            int i=order_now.vect[a];
            order_now.loss+=T_vector[(i + j * n)];
        }
    }
}

Kokkos::View<double**> create_random_T(int n) {
    // Create a Kokkos view
    Kokkos::View<double**> randomMatrix("RandomMatrix", n, n);

    // Fill the Kokkos view with random numbers
    Kokkos::parallel_for("FillRandomMatrix", Kokkos::RangePolicy<>(0, n), KOKKOS_LAMBDA(int i) {
        // Create a random number generator and distribution for each thread
        std::default_random_engine generator(i+10); // Seed with the row index for reproducibility
        std::uniform_real_distribution<double> distribution(0.0, 1.0);
        for (int j = 0; j < n; ++j) {
            // Generate a random number and assign it to the matrix
            randomMatrix(i, j) = distribution(generator);
        }
    });

    // Wait for the parallel operations to finish
    Kokkos::fence();

    return randomMatrix;
}

void print_matrix(const Kokkos::View<double**>& T) {
    const int numRows = T.extent(0);
    const int numCols = T.extent(1);

    // Print the matrix
    for (int i = 0; i < numRows; ++i) {
        for (int j = 0; j < numCols; ++j) {
            std::cout << T(i, j) << " "; // Accessing T(i, j)
        }
        std::cout << std::endl;
    }
}

std::vector<double> kokkos_to_vector(const Kokkos::View<double**>& T) {
    const int numRows = T.extent(0);
    const int numCols = T.extent(1);
    std::vector<double> T_vector(numRows*numCols);
    // Print the matrix
    for (int i = 0; i < numRows; ++i) {
        for (int j = 0; j < numCols; ++j) {
            T_vector[i + j * numCols]=T(i, j);
        }
    }
    return T_vector;
}

order find_optimal_order_exhaustive(Kokkos::View<double**>& T) {
    // Find n from T. Assume T is square
    int n = T.extent(0);
    assert(n == T.extent(0) && "Error: T must be square");

    // Create a vector with numbers from 1 to n
    std::vector<int> numbers(n);
    for (int i = 0; i < n; ++i) { numbers[i] = i + 1; } // Fill the vector with 1, 2, ..., n

    // Assign it as the vect in a temporary ordering
    order temp_order(numbers);
    order best_order;

    // Go through all orderings
    std::vector<double> T_vector=kokkos_to_vector(T); // Convert to vector for speed
    do {
        // Test new ordering against old
        //update_score(T,temp_order); // This uses the kokkos matrix dirrectly. It is slower by a factor of 1.6 or so from limited testing
        update_score_vector_mat(T_vector,n,temp_order);
        //temp_order.print();
        if (temp_order.loss<best_order.loss){
            best_order=temp_order; // We found a new best order so replace
        }
    } while (std::next_permutation(temp_order.vect.begin(), temp_order.vect.end()));
    return best_order;
}