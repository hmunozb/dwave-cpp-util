//
// Created by Humberto Munoz Bauza on 2019-05-16.
//

#ifndef DW_EXAMPLES_SOLVE_H
#define DW_EXAMPLES_SOLVE_H

#include "dwave_cpp_core.h"
#include "parameters.h"
#include "properties.h"

namespace dwave_cpp{
    struct sapi_SubmittedProblem;
    struct sapi_IsingResult;
    typedef vector<vector<vector<short> > > ResultsVec;

    class ProblemSubmission{
        //todo: add ability to submit a batch of problems
    public:
        ~ProblemSubmission();
        void asyncSolveIsing(Solver& solver, vector<Problem>& problem, QuantumSolverParameters& params);
        void await(double timeout);
        void fetch_done();
        vector<sapi_IsingResult*>& get_results();
        ResultsVec results_vector();

    private:
        vector<sapi_SubmittedProblem*> submitted_problem_ptr;
        vector<sapi_IsingResult*> ising_results;
        ErrorMessage err_msg;

    };
    ProblemSubmission asyncSolveIsing(Solver& solver, Problem& problem,
            QuantumSolverParameters& params, ProblemSubmission& submission);
}
#endif //DW_EXAMPLES_SOLVE_H
