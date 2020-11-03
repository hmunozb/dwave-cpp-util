//
// Created by Humberto Munoz Bauza on 2019-05-16.
//

#ifndef DW_EXAMPLES_SOLVE_H
#define DW_EXAMPLES_SOLVE_H

#include "dwave_cpp/core/dwave_cpp_core.h"
#include "parameters.h"
#include "properties.h"
#include <vector>

using namespace std;

namespace dwave_cpp{
    struct sapi_SubmittedProblem;
    struct sapi_IsingResult;
    typedef vector<vector<vector<int8_t> > > ResultsVec;

    class ProblemSubmission{
        //todo: add ability to submit a batch of problems
    public:
        ~ProblemSubmission();
        void asyncSolveIsing(Solver& solver, vector<Problem>& problem, QuantumSolverParameters& params);
        int await(double timeout);
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

    dwave_cpp::ProblemSubmission run_problem_vector(dwave_cpp::Solver& solver, vector<dwave_cpp::Problem>& prob_vec,
                                             dwave_cpp::QuantumSolverParameters& params, double timeout);

}
#endif //DW_EXAMPLES_SOLVE_H
