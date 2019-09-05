//
// Created by Humberto Munoz Bauza on 2019-05-17.
//
#include <dwave_cpp/dwave_cpp.h>
#include <iostream>

using namespace std;

int main(int argc, const char* argv[]){
    string url;
    string token;
    string solver_name;
    if (argc != 4){
        cout << "must have 3 arguments: url token solver_name" << endl;
        return 1;
    }
    url = argv[1];
    token = argv[2];
    solver_name = argv[3];

    dwave_cpp::Connection connection(url, token);
    dwave_cpp::Solver solver = connection.get_solver(solver_name);

    const std::set<int>& broken_qubits = solver.get_broken_qubits();
    const std::set<int>& broken_cells = solver.get_broken_cells();
    cout << "Solver " << solver_name << " Status:" <<endl;
    cout <<"\t Broken Qubits: ";
    for(int q : broken_qubits){
        cout << q << ", ";
    }
    cout << "\n";

    cout << "\t Broken Cells: ";
    for(int c : broken_cells){
        cout << c << ", ";
    }
    cout << "\n";

    const dwave_cpp::sapi_SolverProperties* props = solver.get_solver_properties();
    auto anneal_schedule = props->anneal_schedule;
    auto ising_ranges = props->ising_ranges;
    if(anneal_schedule){
        cout << "Max number of points:\n\t" <<
             anneal_schedule->max_points << "\n"
             << "Anneal time range:\n\t"
             << "[" << anneal_schedule->min_annealing_time << ", " << anneal_schedule->max_annealing_time << "]\n"
                ;
    }
    if(ising_ranges){
        cout << "h range\n\t [" << ising_ranges->h_min << ", " << ising_ranges->h_max << "]\n";
        cout << " J range\n\t [" << ising_ranges->j_min << ", " << ising_ranges->j_max << "]\n";
        cout << endl;
    }

    return 0;

}