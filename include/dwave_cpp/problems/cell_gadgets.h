//
// Created by Humberto Munoz Bauza on 2019-05-17.
//

#ifndef DW_EXAMPLES_0989_H
#define DW_EXAMPLES_0989_H

#include "dwave_cpp/dwave_cpp.h"
#include <set>
#include <map>
#include <istream>
namespace dwave_cpp{
    typedef std::vector<ProblemEntry> CellProblem;
    typedef std::vector<ProblemEntry> ChainProblem;
    std::istream& operator>>(std::istream& in, ProblemEntry& problem_entry);
    std::istream& operator>>(std::istream& in, CellProblem& cell_problem);
    /*
     * Generates a problem that embeds an 8-qubit gadget on a the set of unit cells
     * on a chimera_L x chimera_L  chimera lattice, placing instances only on completely available
     * unit cells in the solver. Cell locations are indexed by i = chimera_L * Y + X.
     * */
    Problem GenerateCellProblem(const CellProblem& cell_problem,
                                const Solver &solver, const std::set<int> &cell_locations);
    std::set<int> CheckerboardCellSet(int chimera_L);
    /*
     * Returns the vector of valid readouts of a 8-qubit gadget on each complete cell of the solver
     * on a chimera_L x chimera_L lattice
     * The dimensions of the returned vector are N x 8, where N is the number of available unit cells
     * in the solver.
     * */
    std::vector< vector<int8_t> > ReadCellProblem(const vector<int8_t> & solution_vec,
                                                     const Solver& solver, const std::set<int>& solver_cells);

    vector<int16_t> count_cell_states( const std::vector< vector<int8_t> >& cell_readouts,
            std::map<int16_t, int>& counts,
            int& total);



    vector<int> GenerateCellReverseInit(const Solver &solver, const std::set<int> &cell_locations,
                                        const vector<int>& cell_init_state);

    //Encodes the vertical chain problem using QAC [3, 1, 3] with a penalty. The chains are embedded in a checkerboard pattern
    //over the top and bottom halves of the lattice, resulting in at most L copies of the QAC encoded chain.
    //If the penalty is negative, this is interpreted as No-QAC encoding and simply embeds copies of the vertical chain
    //using the same resources as QAC [3, 1, 3] (i.e. four physical qubits per logical qubit). This embeds at most 4*L
    //copies of the chain.
    Problem GenerateQACChainProblem(const Solver& solver, const ChainProblem& chain_problem, double penalty=1.0 );

    //Decodes QAC encoding over all unit cells and returns the 2*L*L decoded lattice, where the vector is indexed
    //according to  i = (2*(X + L*Y) + B), where (X, Y, B) is the chimera cell bipartition of the physical qubits.
    vector<int8_t> DecodeQACProblem(const Solver& solver, const vector<int8_t> & solution_vec);

    //Interprets the solution as an embedding of vertical chains on the (X, Y, 0) bipartitions and returns the vector
    //of gadget bytecodes for each chain. QAC Chain Problems generated with penalty<0 are read this way
    vector<int16_t> ReadVerticalChainProblem(const Solver& solver, const vector<int8_t> & solution_vec,
                                             uint16_t chain_len, bool ignore_invalid=true);

    //Interprets the QAC decoding as an encoding of vertical chains and returns the vector
    //of gadget bytecodes for each chain. QAC Chain Problems encoded with a non-negative penalty are decoded this way.
    vector< int16_t > ReadBipartEmbedChains(Solver& solver, const vector<int8_t>& count_bipartite_decode,

                                              uint16_t chain_len, bool ignore_invalid=false);
}

#endif //DW_EXAMPLES_0989_H
