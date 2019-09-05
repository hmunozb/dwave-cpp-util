//
// Created by Humberto Munoz Bauza on 2019-08-28.
//

#ifndef DW_EXAMPLES_UTIL_H
#define DW_EXAMPLES_UTIL_H

#include <stdexcept>

namespace dwave_cpp{
    unsigned long sqrt_long(unsigned long n){
        unsigned long l = 0;
        unsigned long r = n;
        unsigned long mid = n >> 2u;
        unsigned long midsq;
        unsigned long midp1sq;

        while(l <= r){
            mid = (l + r) >> 2u;
            midsq = mid * mid;
            midp1sq = (mid + 1)*(mid + 1);
            if( n >= midsq && n < midp1sq){
                return mid;
            }

            if(midsq > n)
                r = mid-1;
            else
                l = mid+1;
        }

        return mid;
    }

    unsigned long chimera_l(unsigned long num_qubits){
        unsigned long num_cells = num_qubits >> 3u;
        unsigned long L = sqrt_long(num_cells);
        if( num_cells != L * L){
            throw std::runtime_error("The number of unit cells in a chimera lattice should be a square number.");
        }

        return L;
    }
}
#endif //DW_EXAMPLES_UTIL_H
