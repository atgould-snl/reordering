#pragma once
#include <iostream>
#include <vector>
#include <list>
#include <limits>
#include <random>
#include <algorithm>
#include <cassert>
#include <cmath>
#include <utility>
#include <chrono> // easy timer

#include "Kokkos_Core.hpp"


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