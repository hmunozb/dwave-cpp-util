//
// Created by Humberto Munoz Bauza on 2019-05-17.
//

#ifndef DW_EXAMPLES_CHIMERA_H
#define DW_EXAMPLES_CHIMERA_H

#include <cstdint>
#include <array>
namespace dwave_cpp{
    // Chimera cell specification by X - Y coordinates of cell, bipartition (0, 1) of the cell, and relative
    // qubit (0-3) in the bipartition
     struct chimera_cell{
         const uint16_t L;
         uint16_t x;
         uint16_t y;
         uint8_t lo;
         uint8_t i;

         inline uint32_t n(){
             return 8*c() + 4*lo + i;
         }
         inline uint32_t c(){
             return x + L*y;
         }
         inline uint32_t b(){
             return 2*c() + lo;
         }

     };

     // Specifies the physical qubits and the penalty qubit for a QAC encoding
     struct qac_spec{
         array<int32_t, 3> phys_qs;
         int32_t pen_q;
         inline int32_t& operator[](size_t n);
     };

    inline int32_t& qac_spec::operator[](size_t n) {
        return phys_qs[n];
    }

    // Get the QAC qubits for the corresponding chimera bipartition (ignores i in the chimera_cell struct)
    qac_spec qac_qubits(chimera_cell cell){
        cell.i = 0;
        int32_t n0 = cell.n();
         if(cell.lo == 0){ // First bipartition
              int32_t pen = n0 + 7;
              return qac_spec{{n0, n0 + 1, n0 + 2}, pen};
         } else {
            int32_t pen = n0 - 1;
            return qac_spec{{n0, n0 + 1, n0 + 2}, pen};
         }

    }
    class Chimera{

    };
}


#endif //DW_EXAMPLES_CHIMERA_H
