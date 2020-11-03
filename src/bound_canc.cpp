//
// Created by Humberto Munoz Bauza on 2019-05-17.
//

#include <dwave_cpp/dwave_cpp.h>
#include <dwave_cpp/problems/cell_gadgets.h>
#include <dwave_cpp/schedules/beta.h>
#include <iostream>
#include <algorithm>
#include <boost/program_options.hpp>
#include <chrono>

#include "dw_prog.h"
#include "dwave_cpp/dwave_cpp.h"

using namespace std;
namespace po = boost::program_options;
namespace fs = boost::filesystem;

struct rev_gadget_program: public gadget_program{
    std::set<int> SolverCells;
    vector<int> ReverseInitVec;

    po::options_description bound_canc_opts;

    rev_gadget_program();

    void set_schedule_parameters(dwave_cpp::Solver &solver, dwave_cpp::QuantumSolverParameters &params) override;
    dwave_cpp::Problem encode_problem( dwave_cpp::Solver& solver,
                                       dwave_cpp::ProblemAdj &cell_problem) override;
    vector<int16_t> decode_problem(
            dwave_cpp::Solver& solver, const vector<int8_t>& readout) override;
};

void rev_gadget_program::set_schedule_parameters(dwave_cpp::Solver &solver,
                                                 dwave_cpp::QuantumSolverParameters &params) {
    if(sched_type == rev){
        ReverseInitVec = dwave_cpp::GenerateCellReverseInit(solver, SolverCells, reverse_init_cell_vec);
//        int n = ReverseInitVec.size();
//        for(int i = 0; i < n; ++i){
//            int c = ReverseInitVec[i];
//            if (c!= -1 && c != 1) ReverseInitVec[i] = 1;
//
//        }
        params.set_reverse_anneal(ReverseInitVec);
    }
    gadget_program::set_schedule_parameters(solver, params);
}

dwave_cpp::Problem rev_gadget_program::encode_problem(dwave_cpp::Solver& solver,
                                                      dwave_cpp::ProblemAdj &cell_problem) {
    if( cell_locations_vec.empty() ){
        SolverCells = dwave_cpp::CheckerboardCellSet(16);
    }
    else {
        for(int cell : cell_locations_vec) SolverCells.insert(cell);
    }
    if(verbose){
        cout << "Solver Cells used:\n";
        for(int cell : SolverCells)
            cout << cell << ", ";
        cout << endl;
    }

    return dwave_cpp::GenerateCellProblem(cell_problem, solver, SolverCells);
}

vector<int16_t> rev_gadget_program::decode_problem(
        dwave_cpp::Solver& solver, const vector<int8_t>& readout) {
    int num_readouts = 0;
    auto sep_readout = dwave_cpp::ReadCellProblem(readout, solver, SolverCells);
    return dwave_cpp::count_cell_states(sep_readout, counts, num_readouts);
}

rev_gadget_program::rev_gadget_program() : gadget_program(){
    positional_options.add("sched", 1)
            .add("output", 1);

}

int main(int argc, const char* argv[]){
    rev_gadget_program Prog;
    if(Prog.parse_all_options(argc, argv))
        return 0;
    if(Prog.check_options())
        return 1;
    try{
        Prog.run();
    } catch( std::runtime_error& e){
        cerr << "An unexpected error occured:\n\t" << e.what() << endl;
        return 1;
    }

    return 0;
}