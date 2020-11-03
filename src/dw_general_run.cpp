//
// Created by Humberto Munoz Bauza on 1/14/20.
//

#include <dwave_cpp/dwave_cpp.h>
#include <iostream>
#include <algorithm>
#include <boost/program_options.hpp>
#include <chrono>

#include "dw_prog.h"

using namespace std;
using namespace dwave_cpp;

namespace po = boost::program_options;
namespace fs = boost::filesystem;

struct dw_general_run : public advanced_schedule_program{
    string input_file;
    ProblemAdj  problem_adjacency;
    void run();
};

void dw_general_run::run(){
    if(verbose) cout << "Importing problem from" << input_file << endl;
    problem_adjacency = import_problem(input_file);
    generate_schedule();

    cout << sched << endl;
    dwave_cpp::Connection connection(url, token, prompt_retries);
    dwave_cpp::Solver solver = connection.get_solver(solver_str);
    dwave_cpp::QuantumSolverParameters params;
    set_parameters(params);


}