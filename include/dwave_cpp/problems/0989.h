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

    /*
     * Generates a problem that embeds Instance 0989 on a checkerboard pattern of unit cells
     * on a chimera_L x chimera_L  chimera lattice, placing instances only on completely available
     * unit cells in the solver
     * */
    typedef std::vector<ProblemEntry> CellProblem;
    typedef std::vector<ProblemEntry> ChainProblem;
    //extern const ProblemEntry[] array_0989_V1;
    //extern const ProblemEntry[] array_0989_V6;
    std::istream& operator>>(std::istream& in, ProblemEntry& problem_entry);
    std::istream& operator>>(std::istream& in, CellProblem& cell_problem);
    CellProblem C0989_V1();
    CellProblem C0989_V6();
    Problem GenerateCellProblem(const CellProblem& cell_problem,
                                const Solver &solver, const std::set<int> &cell_locations);
    vector<int> GenerateCellReverseInit(const Solver &solver, const std::set<int> &cell_locations,
                                        const vector<int>& cell_init_state);
    std::vector< vector<int8_t> > ReadCellProblem(const vector<int8_t> & solution_vec,
                                                     const Solver& solver, const std::set<int>& solver_cells);

    vector<int16_t> count_cell_states( const std::vector< vector<int8_t> >& cell_readouts,
            std::map<int16_t, int>& counts,
            int& total);

    std::set<int> CheckerboardCellSet(int chimera_L);


    Problem GenerateQACChainProblem(const ChainProblem& chain_problem,
                                    const Solver& solver);

    Problem generate_0989_v1_array(const Solver& solver, int chimera_L);
    Problem generate_0989_array(const Solver& solver, int chimera_L);

    Problem generate_0989_cells(const Solver& solver, const std::set<int>& cells);

    /*
     * Returns the vector of valid readouts of a 0989 on each complete cell of the solver
     * on a chimera_L x chimera_L lattice
     * The dimensions of the returned vector are N x 8, where N is the number of available unit cells
     * in the solver.
     * */
    std::vector< vector<short int> > read_0989(const vector<vector<short int> >& results_vec,
                                               const Solver& solver, int chimera_L);
    std::vector< vector<short int> > read_0989(const  vector<short int>& solution_vec,
                                               const Solver& solver, int chimera_L);
}

#endif //DW_EXAMPLES_0989_H
