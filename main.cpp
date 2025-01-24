#include "exhaustive.h"
#include "branchAndBound.h"



int main(int argc, char* argv[]) {
    int n=9;
    Kokkos::initialize(argc, argv);
    Kokkos::View<double**> T=create_random_T(n);
    print_matrix(T);
    easy_timer main_time=easy_timer();
    //doNothing(1);
    std::map<int,int> orderMap = branchAndBoundLossMinimizationBlockPermutation(T, true);
    main_time.print_time();
    return 0;
}


int main_(int argc, char* argv[]) {
    int n=9;
    Kokkos::initialize(argc, argv);
    Kokkos::View<double**> T=create_random_T(n);
    print_matrix(T);
    easy_timer main_time=easy_timer();
    order best_order=find_optimal_order_exhaustive(T);
    main_time.print_time();
    best_order.print();
    return 0;
}


int main_old2(int argc, char* argv[]) {
    Kokkos::initialize(argc, argv);

    {
        Kokkos::parallel_for("HelloWorld", 1, KOKKOS_LAMBDA(const int i) {
            printf("Hello, World from Kokkos!\n");
        });
    }

    Kokkos::finalize();

    return 0;
}