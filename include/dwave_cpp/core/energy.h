//
// Created by Humberto Munoz Bauza on 9/15/19.
//

#include "dwave_cpp_core.h"

#ifndef DW_CPP_UTIL_ENERGY_H
#define DW_CPP_UTIL_ENERGY_H

namespace dwave_cpp{
    struct EnergyEval{
        const vector<ProblemEntry>& problem;

        double operator()(const vector<int8_t>& ising_spins);
    };
}

#endif //DW_CPP_UTIL_ENERGY_H
