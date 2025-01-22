#pragma once
#include <iostream>
#include <vector>
#include <list>
#include <limits>
#include <random>

#include <chrono> // easy timer

#include "Kokkos_Core.hpp"

struct order{
    double loss=std::numeric_limits<double>::infinity();
    std::vector<int> vect;

    // Default constructor
    order() : loss(std::numeric_limits<double>::infinity()), vect() {}

    // Order specified constructor
    order(const std::vector<int>& vals): vect(vals){}

    // Print method
    void print() const {
        std::cout << "Order Details:" << std::endl;
        std::cout << "Loss: " << loss << std::endl;
        std::cout << "Vector: ";
        if (vect.empty()) {
            std::cout << "Empty" << std::endl;
        } else {
            for (size_t i = 0; i < vect.size(); ++i) {
                std::cout << vect[i];
                if (i < vect.size() - 1) {
                    std::cout << ", "; // Print a comma for all but the last element
                }
            }
            std::cout << std::endl;
        }
    }

};



////////////////// EASY TIMER ///////////////////
class easy_timer{ // Use as  easy_timer main_time=easy_timer();   then     main_time.print_time();    e.g.
public:
    double t;
    std::chrono::_V2::system_clock::time_point start, end;
    easy_timer(){
        restart();
    }
    void restart(){
        start  = std::chrono::high_resolution_clock::now();
    }
    double time(){
        end  = std::chrono::high_resolution_clock::now();
        t = std::chrono::duration_cast<std::chrono::duration<double>>(end-start).count();
        return t;
    }
    void print_time(){
        t=time();
        std::cout << t << " s elapsed" << std::endl;
    }
};