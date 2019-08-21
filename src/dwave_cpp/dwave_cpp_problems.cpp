//
// Created by Humberto Munoz Bauza on 2019-05-17.
//

#include "dwave_cpp/problems/0989.h"

#include <vector>
#include <array>
#include <tuple>
#include <set>
#include <map>
#include <iostream>
#include <sstream>

namespace dwave_cpp{
#include <dwave_sapi.h>

    std::istream& operator>>(std::istream& in, ProblemEntry& problem_entry){
        in  >> problem_entry.i
            >> problem_entry.j
            >> problem_entry.value;
        return in;
    }
    std::istream& operator>>(std::istream& in, CellProblem& cell_problem){
        std::string line;
        cell_problem.clear();
        while( std::getline(in, line)){
            std::istringstream iss( line );
            ProblemEntry problem_entry;
            if(iss >> problem_entry)
                cell_problem.push_back(problem_entry);
        }
        in.eof();
        return in;
    }

    const int _num_cell_entries_0989 = 24;
    const ProblemEntry array_0989_V6[_num_cell_entries_0989] =
    {
            {0,	0, -3},
            {1,	1, -2},
            {2,	2,	2},
            {3,	3,	-3},
            {4,	4,	1},
            {5,	5,	3},
            {6,	6,	-3},
            {7,	7,	3},

            {0,	4,	3},
            {0,	5,	-3},
            {0,	6,	-3},
            {0,	7,	-3},
            {1,	4,	-3},
            {1,	5,	-3},
            {1,	6,	3},
            {1,	7,	-3},
            {2,	4,	-3},
            {2,	5,	-3},
            {2,	6,	-3},
            {2,	7,	-3},
            {3,	4,	-3},
            {3,	5,	-3},
            {3,	6,	-3},
            {3,	7,	-3}
    };

    const ProblemEntry array_0989_V1[_num_cell_entries_0989] =
    {       {0, 0, -2},
            {1, 1, -2},
            {2, 2, 2},
            {3, 3, -3},
            {4, 4, 1},
            {5, 5, 3},
            {6, 6, -3},
            {7, 7, 3},

            {0, 4, 1},
            {0, 5, -3},
            {0, 6, -1},
            {0, 7, -1},
            {1, 4, -1},
            {1, 5, -2},
            {1, 6, 2},
            {1, 7, -3},
            {2, 4, -2},
            {2, 5, -1},
            {2, 6, -3},
            {2, 7, -2},
            {3, 4, -3},
            {3, 5, -3},
            {3, 6, -1},
            {3, 7, -3}
    };
    CellProblem C0989_V1(){
        CellProblem cell_problem;
        cell_problem.resize(_num_cell_entries_0989);
        for(int i = 0; i < _num_cell_entries_0989; ++i){
            cell_problem[i] = array_0989_V1[i];
        }
        return cell_problem;
    }
    CellProblem C0989_V6(){
        CellProblem cell_problem;
        cell_problem.resize(_num_cell_entries_0989);
        for(int i = 0; i < _num_cell_entries_0989; ++i){
            cell_problem[i] = array_0989_V6[i];
        }
        return cell_problem;
    }



    template<typename _int>
    inline bool _cell_cond(_int x, _int y){
        return (x+y)%4 == 0 && y%2==0;
    }

    std::set<int> CheckerboardCellSet(int chimera_L){
        std::set<int> cell_locations;
        for(int x = 0; x < chimera_L; ++x){
            for(int y = 0; y < chimera_L; ++y){
                if(_cell_cond(x, y)){
                    int i = chimera_L * y + x;
                    cell_locations.insert(i);
                }
            }
        }
        return  cell_locations;
    }

    Problem GenerateCellProblem(const ProblemEntry *cell_problem, int size,
                                const Solver &solver, const std::set<int> &cell_locations){

        Problem problem;
        ProblemEntry P0;
        vector<ProblemEntry>& entries = problem.entries();
        const std::set<int>& broken_cells = solver.get_broken_cells();
        for(int i : cell_locations){
            if(broken_cells.find(i) == broken_cells.end()){
                int n0 = i*8;
                for( int k = 0; k < size; ++k){
                    P0 = cell_problem[k];
                    ProblemEntry entry = ProblemEntry{
                            n0 + P0.i,
                            n0 + P0.j,
                            P0.value};
                    entries.push_back(entry);
                }
            }
        }

        return problem;
    }



    Problem GenerateCellProblem(const CellProblem& cell_problem,
                                const Solver &solver, const std::set<int> &cell_locations){

        Problem problem;
        vector<ProblemEntry>& entries = problem.entries();

        const std::set<int>& broken_cells = solver.get_broken_cells();
        for(int i : cell_locations){
            if(broken_cells.find(i) == broken_cells.end()){
                int n0 = i*8;
                for( ProblemEntry P0 : cell_problem){
                    ProblemEntry entry = ProblemEntry{
                            n0 + P0.i,
                            n0 + P0.j,
                            P0.value};
                    entries.push_back(entry);
                }
            }
        }

        return problem;
    }

    vector<int> GenerateCellReverseInit(const Solver &solver, const std::set<int> &cell_locations,
                                const vector<int>& cell_init_state){
        if(cell_init_state.size() > 8){
            throw runtime_error("Cell initial state too large");
        }

        int num_qubits = solver.get_solver_properties()->quantum_solver->num_qubits;
        vector<int> init_state(num_qubits, 3); //initialize unspecified init states at 3

        const std::set<int>& broken_cells = solver.get_broken_cells();
        for(int i : cell_locations){
            if (broken_cells.find(i) == broken_cells.end()){
                int n0 = i*8;
                for( int j = 0; j < 8; ++j){
                    init_state[n0 + j] = cell_init_state[j];
                }
            }
        }

        return init_state;
    }

    std::vector< vector<short int> >_read_array(const vector<short int> & solution_vec,
                const Solver& solver, const std::set<int>& cell_locations){
        std::vector< vector<short int> > instance_results;

        for( int i : cell_locations){
            int n0 = i*8;
            instance_results.emplace_back(
                    &solution_vec[n0], &solution_vec[n0+8] );
        }

        return instance_results;

    }

    Problem _generate_array(const ProblemEntry cell_problem[], int size, const Solver& solver, int chimera_L){

        Problem problem;
        ProblemEntry P0;
        vector<ProblemEntry>& entries = problem.entries();
        const std::set<int>& broken_cells = solver.get_broken_cells();

        for(int x = 0; x < chimera_L; ++x){
            for(int y = 0; y < chimera_L; ++y){
                if(_cell_cond(x, y)){
                    int i = chimera_L * y + x;
                    if(broken_cells.find(i) == broken_cells.end()){
                        int n0 = i*8;

                        for( int k = 0; k < size; ++k){
                            P0 = cell_problem[k];
                            entries.push_back(ProblemEntry{
                                    n0 + P0.i,
                                    n0 + P0.j,
                                    P0.value});
                        }
                    }
                }
            }
        }
        return problem;
    }

    Problem generate_0989_array(const Solver& solver, int chimera_L){

//        Problem problem_0989;
//        vector<ProblemEntry>& entries = problem_0989.entries();
//        const std::set<int>& broken_cells = solver.get_broken_cells();

/*        for(int x = 0; x < chimera_L; ++x){
            for(int y = 0; y < chimera_L; ++y){
                if((x+y)%2 == 0){
                    int i = chimera_L * y + x;
                    if(broken_cells.find(i) == broken_cells.end()){
                        int n0 = i*8;
                        for( ProblemEntry P0 : array_0989){
                            entries.push_back(ProblemEntry{
                                n0 + P0.i,
                                n0 + P0.j,
                                P0.value});
                        }
                    }
                }
            }
        }*/

        Problem p = _generate_array(array_0989_V6, _num_cell_entries_0989, solver, chimera_L);
        return p;
    }

    Problem generate_0989_v1_array(const Solver& solver, int chimera_L){
        Problem p = _generate_array(array_0989_V1, _num_cell_entries_0989, solver, chimera_L);

        return p;
    }

    Problem generate_0989_cells(const Solver& solver, const std::set<int>& cells){
        return GenerateCellProblem(array_0989_V6, _num_cell_entries_0989, solver, cells);
    }


    std::vector< vector<short int> > ReadCellProblem(const vector<short int> & solution_vec,
            const Solver& solver, const std::set<int>& solver_cells){
        std::vector< vector<short int> > instance_results;
        const std::set<int>& broken_cells = solver.get_broken_cells();
        for(int i : solver_cells){
            if(broken_cells.find(i) == broken_cells.end()){
                int n0 = i*8;
                instance_results.emplace_back(
                        &solution_vec[n0], &solution_vec[n0+8] );
            }
        }
        return instance_results;
    }

    vector<short> count_cell_states(
            const std::vector< vector<short int> >& cell_readouts,
            std::map<short, int>& counts,
            int& total){
        unsigned short st = 0;
        short c;

        vector<short> arr(cell_readouts.size());

        int k = 0;
        for (const auto &cell : cell_readouts) {
            st = 0;
            for (short i = 0; i < 8; ++i) {
                c = cell[i];
                if (not(c == 1 or c == -1)) {
                    st = -1;
                    break;
                }
                unsigned q = (c == -1 ? 1 : 0);
                st |= (q << unsigned(7-i));
            }
            arr[k++] = st;
            counts[st] += 1;
            total += 1;
        }
        return arr;
    }

    std::vector< vector<short int> > read_0989(const vector<short int> & solution_vec,
                                               const Solver& solver, int chimera_L){
        std::vector< vector<short int> > instance_results;
        const std::set<int>& broken_cells = solver.get_broken_cells();

        for(int x = 0; x < chimera_L; ++x){
            for(int y = 0; y < chimera_L; ++y){
                if(_cell_cond(x, y)){
                    int i = chimera_L * y + x;
                    if(broken_cells.find(i) == broken_cells.end()){
                        int n0 = i*8;
                        instance_results.emplace_back(
                                &solution_vec[n0], &solution_vec[n0+8] );
                    }
                }
            }
        }


        return instance_results;

    }

    std::vector< vector<short int> > read_0989(const vector<vector<short int> >& results_vec,
            const Solver& solver, int chimera_L){
        std::vector< vector<short int> > instance_results;
        const std::set<int>& broken_cells = solver.get_broken_cells();

        for( const std::vector<short>& solution : results_vec){
            for(int x = 0; x < chimera_L; ++x){
                for(int y = 0; y < chimera_L; ++y){
                    if(_cell_cond(x, y)){
                        int i = chimera_L * y + x;
                        if(broken_cells.find(i) == broken_cells.end()){
                            int n0 = i*8;
                            instance_results.emplace_back(
                                    &solution[n0], &solution[n0+8] );
                        }
                    }
                }
            }
        }

        return instance_results;

    }
}