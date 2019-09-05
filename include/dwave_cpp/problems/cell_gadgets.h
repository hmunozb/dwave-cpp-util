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

    Problem GenerateQACChainProblem(const ChainProblem& chain_problem,
                                    const Solver& solver);
}

#endif //DW_EXAMPLES_0989_H
