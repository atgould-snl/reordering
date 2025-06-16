#include "branchAndBound.h"
#include "branchAndBound.cpp" // TODO: THIS SHOULD NOT BE INCLUDED HERE! SHOULD GO IN THE CMAKE
#include <iostream>
#include <cmath>
#include <math.h>
#include <ostream>
#include <string>
#include <vector>
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

Kokkos::View<double**> getRandomT(int n, bool expMode = true, double randMax = 3., double sparsity = 0.5) { // Exponential mode simulates variability over orders of magnitude
    // Create a Kokkos View with the same dimensions
    Kokkos::View<double**> kokkosView("kokkosView", n, n);

    // Seed the random number generator
    std::random_device rd;  // Obtain a random number from hardware
    //std::mt19937 gen(5);  // Seed the generator // rd();
    std::mt19937 gen(rd());
    std::uniform_real_distribution<> dis(0.0, randMax);  // Define the range
    std::uniform_real_distribution<> dis_sparse(0.0, 1.0);  // Define the range

    // Fill the Kokkos View with random values
    for (size_t i = 0; i < n; ++i) {
        for (size_t j = 0; j < n; ++j) {
            // Generate a random exponent between 1 and maxExp
            double rand=dis(gen);
            // Generate a random number in the range [10^1, 10^maxExp]
            if (dis_sparse(gen) > sparsity){
                kokkosView(i, j) = expMode ? std::pow(10, rand) : rand;
            }
            else{
                kokkosView(i, j) = 0.0;
            }
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


Kokkos::View<double**> getUniformMergeCostsLike(Kokkos::View<double**> T) { // Exponential mode simulates variability over orders of magnitude
    // Create a Kokkos View with the same dimensions
    int n = T.extent(0);
    Kokkos::View<double**> kokkosView("kokkosView", n, n);
    // Fill the Kokkos View with 0.5 values
    for (size_t i = 0; i < n; ++i) {
        for (size_t j = 0; j < n; ++j) {
            kokkosView(i, j) = 0.5;
        }
    }
    return kokkosView;
}

void print_matrix(const Kokkos::View<double**>& T) {
    const int numRows = T.extent(0);
    const int numCols = T.extent(1);

    std::cout << "\n";
    std::cout << "T = {\n";
    for (int i = 0; i < numRows; i++){
    std::cout << "    {";
    for (int j = 0; j < numCols; j++){ 
        // Print each element, followed by a comma if it's not the last element
        std::cout << std::fixed << std::setprecision(12) << T(i,j);
        if (j < numCols - 1) {
            std::cout << ", ";
        }
    }
    std::cout << "},\n";
    }
    std::cout << "};\n";
}

void print_map(const std::map<int, int>& myMap) {
    std::cout << "Map contents:\n";
    for (const auto& pair : myMap) {
        std::cout << "Key: " << pair.first << ", Value: " << pair.second << '\n';
    }
}

void print_order(const std::vector<int>& values) {
    std::cout << "Order contents: ";
    for (const int& value : values) {
        std::cout << value << " ";
    }
    std::cout << std::endl;
}

#include <iostream>
#include <fstream>
#include <vector>
#include <string>
void writeVectorToCSV(const std::vector<std::vector<double>>& data, const std::string& filename = "cost_over_its.csv") {
    // Open a file in write mode
    std::ofstream file(filename);

    // Check if the file is open
    if (!file.is_open()) {
        std::cerr << "Error opening file: " << filename << std::endl;
        return;
    }

    // Find the maximum length of the inner vectors
    size_t max_length = 0;
    for (const auto& vec : data) {
        if (vec.size() > max_length) {
            max_length = vec.size();
        }
    }

    // Write the data to the CSV file
    for (size_t i = 0; i < max_length; ++i) {
        for (size_t j = 0; j < data.size(); ++j) {
            if (i < data[j].size()) {
                file << data[j][i]; // Write the element
            }
            // Add a comma if it's not the last column
            if (j < data.size() - 1) {
                file << ",";
            }
        }
        file << "\n"; // New line after each row
    }

    // Close the file
    file.close();
    std::cout << "Data written to " << filename << " successfully." << std::endl;
}

// Function to calculate Stirling numbers of the second kind
int stirling_second_kind(int n, int k) { // https://brilliant.org/wiki/distinct-objects-into-identical-bins/
    // Base cases
    if (n == k || k == 1) {
        return 1;
    }
    if (k == 0 || n == 0) {
        return 0;
    }
    // Recursive calculation
    return k * stirling_second_kind(n - 1, k) + stirling_second_kind(n - 1, k - 1);
}

// Function to calculate total grouping configurations
int total_grouping_configs(int n) {
    int total = 0;
    for (int r = 1; r <= n; ++r) {
        total += stirling_second_kind(n, r);
    }
    return total;
}

////////////////////////////////
//////////// TESTS /////////////
////////////////////////////////


std::map<int,int> run_test(Kokkos::View<double**> T, double costTarget = 0, bool exhuastive = false, std::string test_name = "default_name"){
    BucketOrderingSolver soln = BucketOrderingSolver(T,getUniformMergeCostsLike(T),costTarget,1,exhuastive);
    soln.solve();
    std::cout << "Ran test: " << test_name << std::endl;
    std::cout << "Cost: " << soln.get_best().get_cost() << std::endl;
    std::cout << "Score (sqrt cost): " << std::sqrt(soln.get_best().get_cost()) << std::endl;
    std::cout << "Merge cost: " << soln.get_best().get_mergeCost() << std::endl;
    std::cout << "Leaf nodes : " << soln.get_leafNodes() << std::endl;
    std::cout << "Inner nodes: " << soln.get_internalNodes() << std::endl;
    auto map = soln.get_best().get_map();
    print_map(map);
    std::cout << std::endl;
    return map;
}

std::map<int,int> run_test(std::vector<std::vector<double>> T_vec_of_vec, double costTarget = 0, bool exhuastive = false, std::string test_name = "default_name"){
    auto T = vectorToKokkosView(T_vec_of_vec);
    return run_test(T, costTarget, exhuastive, test_name);
}

double frobeniusNorm(const std::vector<std::vector<double>>& matrix) {
    double sum = 0.0;

    // Iterate through each element of the matrix
    for (const auto& row : matrix) {
        for (double value : row) {
            sum += value * value; // Sum of squares
        }
    }

    // Return the square root of the sum of squares
    return std::sqrt(sum);
}

BlockNormsViewType get_tournament(const Kokkos::View<double**>& T) {
  const int N = T.extent(0);
  BlockNormsViewType tournament("normalForm", N, N);
  // Tournament based on the buckets
  for (int a = 0; a < N; a++)
  {
    for (int b = 0; b < N; b++)
    {
      if (a > b) continue; // Block is below diag, filter out
      if (a == b)
      {
        tournament(a, a) = 0;
        continue;
      }
      // Sum block above diag
      double Total = 0;
      Total += T(a, b)*T(a, b);

      // Take away sum of block below diag
      Total -= T(b, a)*T(b, a);
      tournament(a, b) = Total;
      tournament(b, a) = -Total;
    }
  }
  return tournament;
}

std::tuple<BlockNormsViewType, double, double> sortByOrder(std::vector<int> order, BlockNormsViewType H)// H mats are not squares, T mats are
{
    int N=order.size();
    double norm_A=0;
    // Make squared copy
    BlockNormsViewType T("newSortedMat", N, N);
    for (int i=0; i<N; i++){
      for (int j=0; j<N; j++){
        T(i,j)=H(i,j)*H(i,j);
        norm_A += T(i,j);
      }
    }
    norm_A = std::sqrt(norm_A);

    // Get vector from order along with new size

    int N_new = *std::max_element(order.begin(), order.end())+1;

    // Make sorted merged matrix explicitly
    BlockNormsViewType T_sorted("newSortedMat", N_new, N_new);
    for (int i=0; i<N; i++){
      for (int j=0; j<N; j++){
        T_sorted(order[i],order[j]) += T(i,j);
      }
    }

    // Take sqaure roots
    double norm_A_sorted = 0;
    BlockNormsViewType H_sorted("newSortedMat", N_new, N_new);
    for (int i=0; i<N_new; i++){
        for (int j=0; j<N_new; j++){
            norm_A_sorted += T_sorted(i,j);
            H_sorted(i,j) = std::sqrt(T_sorted(i,j));
        }
    }
    norm_A_sorted = std::sqrt(norm_A_sorted);

    // Take LD sum
    double total=0;
    for (int i=0; i<N_new; i++){
        for (int j=0; j<N_new; j++){
            if (i>j){
                total += T_sorted(i,j);
            }
        }
    }
    double finalScore = sqrt(total);
    double finalScoreTimesNormA = finalScore * norm_A;

  return {H_sorted, finalScore, finalScoreTimesNormA};
}

double get_LowerDiagSumByOrder(BlockNormsViewType matrix, std::vector<int> order){ // Find the score based on the order
  double cost=0;
  int N = order.size();
  for (int a=0; a < int(N); a++){
    for (int b=0; b < int(N); b++){
      cost += order[a]>order[b] ? matrix(a,b) : 0;
    }
  }
  return cost;
}

// Or use map overload
std::tuple<BlockNormsViewType, double, double> sortByMap(std::map<int, int> map, BlockNormsViewType H)// H mats are not squares, T mats are
{
    // Get vector from map along with new size
    std::vector<int> values;
    for (const auto& pair : map) {
      values.push_back(pair.second);
    }
    return sortByOrder(values, H);
}


//////// LOP ISSUE DEMO ////////////
int main_2(int argc, char* argv[]) {
    // Load into kokkos
    Kokkos::initialize();
    std::cout << "Running test serries..." << std::endl;
    std::vector<std::vector<double>> T_vec_of_vec;
    std::vector<std::vector<double>> Tournament;

    T_vec_of_vec = {
    {0.315781, 0.027645, 0.018343, 0.007621, 0.024263, 0.015416, 0.019881, 0.012559},
    {0.025879, 0.507779, 0.030138, 0.012945, 0.037173, 0.019956, 0.031640, 0.023114},
    {0.022455, 0.031327, 0.377431, 0.010560, 0.028490, 0.016374, 0.019379, 0.012608},
    {0.007493, 0.009718, 0.011001, 0.119354, 0.011526, 0.005281, 0.006693, 0.006526},
    {0.023247, 0.038834, 0.028760, 0.011672, 0.454487, 0.019286, 0.027340, 0.017193},
    {0.013233, 0.023079, 0.015198, 0.007801, 0.020683, 0.238708, 0.016659, 0.010489},
    {0.018528, 0.034186, 0.025299, 0.011024, 0.028358, 0.017391, 0.377431, 0.015701},
    {0.013561, 0.022519, 0.014828, 0.005540, 0.016741, 0.006081, 0.017621, 0.238708},
    };

    Tournament = {
    {0.000000000000, 0.000038265851, 0.000168820955, 0.000159902185, 0.000070610834, 0.000067523839, 0.000006216579, -0.000109547992},
    {-0.000038265851, 0.000000000000, 0.000061041767, -0.000049882820, 0.000024115617, -0.000007155300, -0.000051749280, 0.000032089360},
    {-0.000168820955, -0.000061041767, 0.000000000000, -0.000015048915, 0.000007640325, -0.000077430560, -0.000054874368, -0.000036466776},
    {-0.000159902185, 0.000049882820, 0.000015048915, 0.000000000000, 0.000041383280, -0.000009953568, -0.000027760617, 0.000031531200},
    {-0.000070610834, -0.000024115617, -0.000007640325, -0.000041383280, 0.000000000000, -0.000001934592, -0.000032966640, 0.000011897076},
    {-0.000067523839, 0.000007155300, 0.000077430560, 0.000009953568, 0.000001934592, 0.000000000000, 0.000062540767, -0.000026172240},
    {-0.000006216579, 0.000051749280, 0.000054874368, 0.000027760617, 0.000032966640, -0.000062540767, 0.000000000000, 0.000073040560},
    {0.000109547992, -0.000032089360, 0.000036466776, -0.000031531200, -0.000011897076, 0.000026172240, -0.000073040560, 0.000000000000},
    };

    std::cout << "\nWe experiment on tournament: \n";
    print_matrix(vectorToKokkosView(Tournament));
    LinearOrderingSolver solver{vectorToKokkosView(Tournament), 1}; // Is communicator pulled in from the namespace??
    solver.solve();
    std::vector<int> order_of_rows = solver.order();
    std::vector<int> order = std::vector<int>(order_of_rows.size());
    std::iota( order.begin(), order.end(), 0);
    std::sort(order.begin(), order.end(), [order_of_rows](int a, int b) { return order_of_rows[a] < order_of_rows[b];});

    std::cout << "Running the LOP code on this order/map we get a score of: " << get_LowerDiagSumByOrder(vectorToKokkosView(Tournament), order) << " with \n";
    print_order(order);

    std::cout << "However, we can hardcode this order/map \n";
    std::vector<int> better_order = {0, 4, 7, 2, 5, 3, 1, 6};
    //better_order = {0, 1, 2, 3, 4, 5, 6, 7};
    auto [mergedMatrix_2, scoreUnscaledRecomputed_2, scoreRecomputedTimesNormA_2] =sortByOrder(better_order, vectorToKokkosView(T_vec_of_vec));
    print_order(better_order);
    std::cout << "The issue is that this hard code order gets us a slightly better score of: " << get_LowerDiagSumByOrder(vectorToKokkosView(Tournament), better_order) << "\n";

    return 0;
}

int main(int argc, char* argv[]) {
    // Load into kokkos
    Kokkos::initialize();
    std::cout << "Running test serries..." << std::endl;
    std::vector<std::vector<double>> T_vec_of_vec;

    //// T1 TESTING ////
    T_vec_of_vec = {
        {1000}
    };
    //run_test(T_vec_of_vec,100);

    ///// T2 TESTING ///


    T_vec_of_vec = {
        {0.329500, 0.020269, 0.022056, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000},
        {0.019202, 0.409623, 0.028737, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000},
        {0.019759, 0.025354, 0.444298, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000},
        {0.000000, 0.000000, 0.000000, 0.140499, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000},
        {0.000000, 0.000000, 0.000000, 0.000000, 0.140499, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000},
        {0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.140499, 0.000000, 0.000000, 0.000000, 0.000000},
        {0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.140499, 0.000000, 0.000000, 0.000000},
        {0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.476457, 0.000000, 0.000000},
        {0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.280999, 0.000000},
        {0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.371727},
    };

    T_vec_of_vec = {
        {0.329500, 0.020269, 0.000000, 0.000000},
        {0.019202, 0.409623, 0.000000 ,0.000000},
        {0.000000, 0.000000, 0.444298, 0.000000},
        {0.000000, 0.000000, 0.000000, 0.140499}
    };

    T_vec_of_vec = {
        {0.197756, 0.014164, 0.010598, 0.010881, 0.011217, 0.012621, 0.014232, 0.005635, 0.007046, 0.005192},
        {0.014733, 0.364644, 0.015606, 0.021709, 0.020916, 0.024087, 0.027712, 0.007428, 0.011036, 0.006702},
        {0.010944, 0.017466, 0.293320, 0.018680, 0.018676, 0.021544, 0.020389, 0.007171, 0.003552, 0.008424},
        {0.012608, 0.022486, 0.018679, 0.364644, 0.022650, 0.026436, 0.029524, 0.007473, 0.008716, 0.009953},
        {0.015170, 0.023547, 0.014865, 0.020995, 0.395512, 0.026741, 0.029068, 0.009233, 0.009797, 0.008176},
        {0.015016, 0.026666, 0.019795, 0.023168, 0.025560, 0.424140, 0.033202, 0.011390, 0.009841, 0.011845},
        {0.014653, 0.027151, 0.019600, 0.023683, 0.027258, 0.031372, 0.450953, 0.010544, 0.008380, 0.009057},
        {0.004841, 0.006627, 0.007130, 0.009934, 0.011210, 0.006135, 0.012443, 0.125072, 0.004660, 0.005308},
        {0.004758, 0.008725, 0.004402, 0.011756, 0.008029, 0.010360, 0.010426, 0.003622, 0.125072, 0.001013},
        {0.006988, 0.007196, 0.006739, 0.008619, 0.009940, 0.011953, 0.011980, 0.003016, 0.000602, 0.125072},
    };

    T_vec_of_vec = {
        {2, 10.1, 1, 10.2},
        {1, 2, 1, 10.3},
        {10.4, 10.5, 2, 10.6},
        {1, 1, 1, 2}
    };

    T_vec_of_vec = {
        {0, 10, 0.1, 0},
        {10, 1, 0.1, 0},
        {1, 1 ,  2, 10},
        {10, 10 , 0, 3}
    };

    T_vec_of_vec = {
      {1, 1, 0.01 , 0},
      {1, 1, 0, 0},
      {1, 1, 1, 1},
      {1, 1, 0, 1}
    };

    T_vec_of_vec = {
    {0.315781, 0.015670, 0.015786, 0.016417, 0.018343, 0.007621, 0.024263, 0.015416, 0.019881, 0.012559},
    {0.015440, 0.279910, 0.015204, 0.013552, 0.017976, 0.008121, 0.019864, 0.009482, 0.016608, 0.014276},
    {0.013106, 0.013043, 0.279910, 0.017161, 0.015132, 0.005997, 0.022091, 0.011521, 0.017257, 0.011905},
    {0.016111, 0.015282, 0.017594, 0.315781, 0.018873, 0.008103, 0.022343, 0.013252, 0.020675, 0.013738},
    {0.022455, 0.016729, 0.018186, 0.019256, 0.377431, 0.010560, 0.028490, 0.016374, 0.019379, 0.012608},
    {0.007493, 0.006468, 0.005322, 0.004927, 0.011001, 0.119354, 0.011526, 0.005281, 0.006693, 0.006526},
    {0.023247, 0.021185, 0.022992, 0.023036, 0.028760, 0.011672, 0.454487, 0.019286, 0.027340, 0.017193},
    {0.013233, 0.011902, 0.013697, 0.014261, 0.015198, 0.007801, 0.020683, 0.238708, 0.016659, 0.010489},
    {0.018528, 0.017411, 0.018007, 0.023266, 0.025299, 0.011024, 0.028358, 0.017391, 0.377431, 0.015701},
    {0.013561, 0.013104, 0.013349, 0.012538, 0.014828, 0.005540, 0.016741, 0.006081, 0.017621, 0.238708},
    };

    //run_test(T_vec_of_vec,100);
    run_test(T_vec_of_vec,0.5);

    return 0;

    // Find map
    auto map=run_test(T_vec_of_vec,0.1*0.1,false, "2x2");
    int N=T_vec_of_vec.size();

    // Get vector from map along with new size
    std::vector<int> values;
    for (const auto& pair : map) {
        values.push_back(pair.second);
    }

    int N_new = *std::max_element(values.begin(), values.end())+1;

    // Make sorted merged matrix explicitly
    std::vector<std::vector<double>> T_vec_of_vec_sorted = std::vector<std::vector<double>>(N_new, std::vector<double>(N_new,0));
    for (int i=0; i<N; i++){
        for (int j=0; j<N; j++){
                double v_old = T_vec_of_vec_sorted[values[i]][values[j]];
                double v_add = T_vec_of_vec[i][j];
                double v_new = std::sqrt(v_old*v_old + v_add*v_add);
                T_vec_of_vec_sorted[values[i]][values[j]] = v_new;
        }
    }

    // Check that origional has approximantly the same 
    double f_old = frobeniusNorm(T_vec_of_vec);
    double f_new = frobeniusNorm(T_vec_of_vec_sorted);

    STK_ThrowRequire( (f_old*0.99 < f_new) && (f_new*0.99 < f_old));

    // Take LD sum
    double total=0;
    for (int i=0; i<N_new; i++){
        for (int j=0; j<N_new; j++){
            if (i>j){
                double v = T_vec_of_vec_sorted[i][j];
                total += v*v;
            }
        }
    }
    print_matrix(vectorToKokkosView(T_vec_of_vec_sorted));
    double finalScore = sqrt(total);
    std::cout << "Total cost recalced: " << finalScore << std::endl;

    return 0;
    

    ///// T4 TESTING ///
    T_vec_of_vec = {
        {2, 10.1, 1, 10.2},
        {1, 2, 1, 10.3},
        {10.4, 10.5, 2, 10.6},
        {1, 1, 1, 2}
    };
    // run_test(T_vec_of_vec,1000,true);
    //run_test(T_vec_of_vec,1000,false);


    //// CORRECTNESS CHECK ////
    int n_trials = 10;
    int n = 6;
    for (int i=0; i<n_trials; i++){
        auto T = getRandomT(n);
        auto m1 = run_test(T,0.01,false,"BB to 0.01");
        auto m1_ex = run_test(T,0.01,true,"Exhuastive to 0.01");
        STK_ThrowRequire(m1==m1_ex);
        auto m2 = run_test(T,0.1,false, "BB to 0.1");
        auto m2_ex = run_test(T,0.1,true,"Exhuastive to 0.1");
        STK_ThrowRequire(m2==m2_ex);
    }

    //// COUNT CHECK ///
    int n_max = 7;
    int n_min = 3;
    for (int i=n_min; i<=n_max; i++){
        auto T = getRandomT(i);
        BucketOrderingSolver soln = BucketOrderingSolver(T,getUniformMergeCostsLike(T),1E6,1,true);
        soln.solve();
        STK_ThrowRequire(soln.get_leafNodes() == total_grouping_configs(i));
    }

    return 0;

    auto T = getRandomT(8);
    print_matrix(T);
    global_timer_all.reset();
    global_timer_LOP.reset();
    global_LOP_call_counter=0;

    run_test(T, 0.3,false, "Big mat test");
    std::cout << "LOP Time: " << global_timer_LOP.time() << std::endl;
    std::cout << "LOP Calls: " << global_LOP_call_counter << std::endl;
    std::cout << "Total Time: " << global_timer_all.time() << std::endl;
    //run_test(T, 0.1,true);

    //// CONVERGNCE TEST ////

    n_trials = 10;
    n = 8;
    for (int i=0; i<n_trials; i++){
        global_best_cost_record.push_back(std::vector<double>());
        auto T = getRandomT(n);
        auto m = run_test(T,0.1,false);
    }
    writeVectorToCSV(global_best_cost_record);
    return 0;
}



    /*

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


    // EXHAUSTIVE SOLN
    BranchAndBoundPermutationSearch exhaustiveObj = BranchAndBoundPermutationSearch(T);
    exhaustiveObj.solveExhuastive();
    exhaustiveObj.minLossOrder.print();
    std::map<int,int> exhuastiveMap=exhaustiveObj.minLossOrder.getMap();
    EXPECT_EQ(exhuastiveMap,solnMap);

    // BB WITH NO BRANCH CUTTING NO MERGING
    BranchAndBoundPermutationSearch bbObj = BranchAndBoundPermutationSearch(T);
    bbObj.solve();
    exhaustiveObj.minLossOrder.print();
    std::map<int,int> bbMap=exhaustiveObj.minLossOrder.getMap();
    EXPECT_EQ(bbMap,solnMap);
    EXPECT_NEAR(exhaustiveObj.minLossOrder.loss,bbObj.minLossOrder.loss,1E-8);
    // Check total scanned
    EXPECT_EQ(bbObj.n,4);
    EXPECT_EQ(bbObj.numLeafNodes,24);

    // BB WITH BRANCH CUTTING
    BranchAndBoundPermutationSearch bbObj_cutting = BranchAndBoundPermutationSearch(T);
    bbObj_cutting.allowBranchCutting=true;
    bbObj_cutting.solve();
    bbObj_cutting.minLossOrder.print();
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
    EXPECT_GE(T(3,3), T(3,0)); // Diag largest
    EXPECT_GE(T(3,3), T(3,1)); // Diag largest
    EXPECT_GE(T(3,3), T(3,2)); // Diag largest

    // Test a bunch of T to make sure exhaustive and branch cutting give the same result
    for (int i=0; i<100; i++){
        T = getRandomT(7);
        BranchAndBoundPermutationSearch exhaustiveObj = BranchAndBoundPermutationSearch(T);
        exhaustiveObj.solveExhuastive();
        exhaustiveObj.minLossOrder.print();
        std::map<int,int> exhuastiveMap=exhaustiveObj.minLossOrder.getMap();

        BranchAndBoundPermutationSearch bbObj = BranchAndBoundPermutationSearch(T);
        bbObj.allowBranchCutting=true;
        bbObj.solve();
        bbObj.minLossOrder.print();
        std::map<int,int> bbMap=exhaustiveObj.minLossOrder.getMap();

        EXPECT_EQ(exhuastiveMap,bbMap);
        EXPECT_NEAR(exhaustiveObj.minLossOrder.loss,bbObj.minLossOrder.loss,1E-8);
    }
    
    // Check with merging. Want to total to be as estimated from formula (see python).
    // 3 -> 13   7 -> 47293
    T = getRandomT(7);
    BranchAndBoundPermutationSearch bbObj_merging = BranchAndBoundPermutationSearch(T);
    bbObj_merging.allowBranchCutting=false;
    bbObj_merging.allowMerge=true;
    bbObj_merging.solve();
    // Check total scanned
    EXPECT_EQ(bbObj_merging.numLeafNodes, 47293);

    easy_timer et1 = easy_timer();
    int trials=10;
    int n=13;
    int total_nodes = 0;
    for (int i=0; i<n; i++){

        T = getRandomT(16,true);
        easy_timer et2 = easy_timer();
        BranchAndBoundPermutationSearch speedObj = BranchAndBoundPermutationSearch(T);
        speedObj.allowBranchCutting=true;
        speedObj.allowMerge=false;
        speedObj.solve();

        // Report results
        std::cout << "Solved problem with:" << std::endl;
        et2.print_time();
        std::cout << "Leaves:" << std::endl;
        std::cout << speedObj.numLeafNodes << std::endl;
        std::cout << "Forks:" << std::endl;
        std::cout << speedObj.numInternalNodes;
        std::cout << std::endl;
        total_nodes+=speedObj.numInternalNodes+speedObj.numLeafNodes;
        et2.restart();
        //print_matrix(T);

        //speedObj.minLossOrder.print();
    }
    std::cout << std::endl << "Ave time for " << trials << " trials at "<< n <<"x"<< n <<": " << et1.time() /trials << std::endl;
    std::cout << "Ave nodes : " << total_nodes/trials << std::endl;
    EXPECT_LE(et1.time() , 10.); // Expect less than a second
}

*/



/*
Failed example: 

T_vec_of_vec = {
    {0.226279, 0.016248, 0.006368, 0.010865, 0.008612, 0.016572, 0.014919, 0.016260, 0.008149, 0.016284},
    {0.016589, 0.383674, 0.015264, 0.015110, 0.009293, 0.029911, 0.022622, 0.028333, 0.015424, 0.027723},
    {0.006234, 0.010748, 0.178889, 0.007710, 0.006107, 0.013070, 0.010551, 0.013701, 0.005816, 0.015335},
    {0.009474, 0.017704, 0.004697, 0.226279, 0.004024, 0.017521, 0.014766, 0.015778, 0.010741, 0.016080},
    {0.007326, 0.009766, 0.005830, 0.003364, 0.113139, 0.008140, 0.007094, 0.009506, 0.004939, 0.011047},
    {0.018681, 0.027312, 0.015172, 0.014085, 0.009174, 0.407930, 0.023032, 0.027849, 0.019686, 0.026296},
    {0.013098, 0.023481, 0.011708, 0.014949, 0.007953, 0.026772, 0.329855, 0.024353, 0.018153, 0.022958},
    {0.016841, 0.026471, 0.013829, 0.017085, 0.012334, 0.029668, 0.022964, 0.407930, 0.019651, 0.031745},
    {0.011978, 0.017768, 0.007117, 0.008510, 0.007758, 0.018227, 0.016067, 0.021345, 0.265335, 0.018717},
    {0.014141, 0.027953, 0.012369, 0.017727, 0.011182, 0.033509, 0.026915, 0.028921, 0.019580, 0.407930},
};
Greedy score: 0.091026
Greedy perm: 
Map contents:
Key: 0, Value: 3
Key: 1, Value: 0
Key: 2, Value: 6
Key: 3, Value: 1
Key: 4, Value: 5
Key: 5, Value: 2
Key: 6, Value: 2
Key: 7, Value: 0
Key: 8, Value: 4
Key: 9, Value: 0

BB score: 0.099435
BB map: 
Map contents:
Key: 0, Value: 3
Key: 1, Value: 0
Key: 2, Value: 6
Key: 3, Value: 1
Key: 4, Value: 5
Key: 5, Value: 2
Key: 6, Value: 2
Key: 7, Value: 0
Key: 8, Value: 4
Key: 9, Value: 0
*** Failure on proc 0 in /ceeGPFS/fgs/atgould/code/tftk/tftk_linsolv/unit_tests/src/UnitTestBlockReorder.C:697
Expected: (scoreBranchAndBound*0.999) <= (scoreGreedy), actual: 0.0993358 vs 0.0910256
[  FAILED  ] on 1 of 1 procs 
BlockReorderingFixture.BranchAndBoundOutperformsGreedyMerging (71 ms)
[  PASSED  ] 0 tests.
[  FAILED  ] 1 tests:
[  FAILED  ] --gtest_filter=BlockReorderingFixture.BranchAndBoundOutperformsGreedyMerging
*/