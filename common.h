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
#include <set>
#include <iomanip>

#include "Kokkos_Core.hpp"

using BlockNormsViewType = Kokkos::View<double**>; //GaussSeidelBlockReorderPreconditionerFactory::BlockNormsViewType;

////////////////// EASY TIMER ///////////////////
class easy_timer{ // Use as  easy_timer main_time=easy_timer();   then     main_time.print_time();    e.g.
public:
    double t;
    bool running = false;
    std::chrono::_V2::system_clock::time_point start_time, end_time;
    easy_timer(){
        reset();
    }
    void reset(){
        t = 0;
        start();
    }
    void start(){
        start_time  = std::chrono::high_resolution_clock::now();
        running = true;
    }
    void pause(){
        if (!running){return;}
        end_time  = std::chrono::high_resolution_clock::now();
        t+= std::chrono::duration_cast<std::chrono::duration<double>>(end_time-start_time).count();
        running = false;
    }
    double time(){
        if (running){
            end_time  = std::chrono::high_resolution_clock::now();
            return t+std::chrono::duration_cast<std::chrono::duration<double>>(end_time-start_time).count();
        }
        return t;
    }
    void print_time(){
        double t_now=time();
        std::cout << t_now << " s elapsed" << std::endl;
    }
};

easy_timer global_timer_LOP;
easy_timer global_timer_all;
int global_LOP_call_counter;
std::vector<std::vector<double>> global_best_cost_record = std::vector<std::vector<double>>(1,std::vector<double>());
std::vector<double> global_best_cost_record_time;