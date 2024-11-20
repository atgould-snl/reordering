#pragma once
#include "common.h"
#include "exhaustive.cpp"

// Function prototypes
void update_score(Kokkos::View<double**>& T, order& order_now);
void update_score_vector_mat(std::vector<double>& T, int n, order& order_now);
Kokkos::View<double**> create_random_T(int n);
void print_matrix(const Kokkos::View<double**>& T);
std::vector<double> kokkos_to_vector(const Kokkos::View<double**>& T);
order find_optimal_order_exhaustive(Kokkos::View<double**>& T);