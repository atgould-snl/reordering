#include "Kokkos_Core.hpp"
#include <iostream>

int main(int argc, char* argv[]) {
    Kokkos::initialize(argc, argv);

    {
        Kokkos::parallel_for("HelloWorld", 1, KOKKOS_LAMBDA(const int i) {
            printf("Hello, World from Kokkos!\n");
        });
    }
    Kokkos::finalize();

    return 0;
}