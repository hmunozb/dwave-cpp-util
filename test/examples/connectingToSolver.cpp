//
// Created by Humberto Munoz Bauza on 2019-05-15.
//

#include "dwave_cpp/dwave_cpp.h"
#include <iostream>
#include <string>
#include <vector>
#include <dwave_sapi.h>
#include <dwave_cpp/core/solve.h>

using namespace std;
int main(int argc, const char* argv[]) {
    string url;
    string token;
    string solver_name;
    if (! (argc == 3 || argc == 4)){
        cout << "must have 2 or 3 arguments: url token [solver_name]" << endl;
        return 1;
    }
    url = argv[1];
    token = argv[2];
    dwave_cpp::Connection connection(url, token);

    if(argc == 4){
        solver_name = argv[3];
    } else{
        vector<string> solver_list = connection.solver_list();
        cout << "List of available solvers: \n";
        for(string s : solver_list){
            cout << "\t" << s << "\n";
        }
        return 0;
    }

    dwave_cpp::Solver solver = connection.get_solver(solver_name);
    const dwave_cpp::sapi_SolverProperties* props = solver.get_solver_properties();
    auto anneal_schedule = props->anneal_schedule;
    auto anneal_offset = props->anneal_offset;
    if(anneal_schedule){
        cout << "Max number of points:\n\t" <<
            anneal_schedule->max_points << "\n"
            << "Anneal time range:\n\t"
            << "[" << anneal_schedule->min_annealing_time << ", " << anneal_schedule->max_annealing_time << "]\n"
            ;
    }
    cout << "Anneal offset properties:" << endl;
    if(anneal_offset){
        cout << "\tStep:" << anneal_offset->step << endl;
        cout <<"\tOffset steps: ";
        for(int i = 0; i < anneal_offset->ranges_len; ++i)
            cout << i << ": [" << anneal_offset->ranges[i].min << ", " << anneal_offset->ranges[i].max << "],  ";

        cout << endl;
    } else{
        cout << "\t Not found\n";
    }
    const auto& broken_cells = solver.get_broken_cells();
    cout << "List of Broken Cells:\n\t";
    for( const auto& cell : broken_cells){
        cout << cell << "  ";
    } 
    cout << endl;
    
    dwave_cpp::QuantumSolverParameters params;

    vector< dwave_cpp::AnnealSchedulePoint> schedule{{0.0, 0.0},
                                                     {10.0, 0.2},
                                                     {20.0, 0.4},
                                                     {30.0, 0.6},
                                                     {40.0, 0.8},
                                                     {50.0, 0.9},
                                                     {60.0, 0.95},
                                                     {100.0, 1.00}};
    params.set_anneal_schedule(schedule);

    //dwave_cpp::sapi_Problem problems[2] = {{NULL, 0}, {NULL, 0}};
    vector<dwave_cpp::Problem> problem_array(1);
    vector<dwave_cpp::ProblemEntry>& problem_entries = problem_array[0].entries();

    cout << "Constructing Problem..." << endl;
    for(int j = 0; j < props->quantum_solver->couplers_len; j++)
    {
        dwave_cpp::ProblemEntry p;
        p.i = props->quantum_solver->couplers[j].q1;
        p.j = props->quantum_solver->couplers[j].q2;
        p.value = 1;
        problem_entries.push_back(p);
    }
    params.get()->num_reads = 20;

    dwave_cpp::ProblemSubmission problem_submission;
    cout << "Submitting..." << endl;
    problem_submission.asyncSolveIsing(solver, problem_array, params);

    cout << "Waiting for results..." << endl;
    problem_submission.await(30.0);
    problem_submission.fetch_done();
    vector<dwave_cpp::sapi_IsingResult*>& ising_results = problem_submission.get_results();
    cout << "Done." << endl;

    for(auto res : ising_results){
        if(res != nullptr){
            for(int j = 0; j < res->num_solutions; ++j){
                cout << "Solution " << j << "\n";
                cout << "\tEnergy: "<< res->energies[j] << "\n";
                cout << "\tOccurences: " << res->num_occurrences[j] << "\n\n";
                for(int k =0; k< res->solution_len; ++k){
                    cout << res->solutions[j*res->solution_len + k];
                }
                cout << "\n";
            }
        }
    }

    return 0;

}