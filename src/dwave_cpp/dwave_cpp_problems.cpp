//
// Created by Humberto Munoz Bauza on 2019-05-17.
//

#include "dwave_cpp/problems/cell_gadgets.h"
#include "dwave_cpp/core/util.h"
#include "dwave_cpp/core/chimera.h"

#include <vector>
#include <array>
#include <tuple>
#include <set>
#include <map>
#include <iostream>
#include <sstream>

namespace dwave_cpp{

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

    void EmbedQACProblemEntry(double J, chimera_cell i_cell, chimera_cell j_cell,
                              vector<ProblemEntry>& physical_entries, double penalty=1.0){
        qac_spec qi = qac_qubits(i_cell);
        if( i_cell.b() == j_cell.b()){ // Just a field spec
            for(int32_t n : qi.phys_qs){
                physical_entries.push_back({n, n, J});
            }
            //Encode QAC for i_cell if penalty is non-negative
            if(penalty >= 0.0)
                for(int32_t n : qi.phys_qs){
                    physical_entries.push_back({n, qi.pen_q, penalty});
                }
            else{ //Otherwise, embed a 4th copy after the 3rd physical qubit
                physical_entries.push_back({qi[2]+1, qi[2]+1, J});
            }

        } else { // Embed two adjacent QAC qubits
            qac_spec qj = qac_qubits(j_cell);
            // Connect the corresponding physical qubits
            for(int n = 0; n < 3; ++n){
                physical_entries.push_back({qi[n], qj[n], J});
            }
            //Encode QAC for j_cell if penalty is non-negative
            if(penalty >= 0)
                for(int32_t n : qj.phys_qs){
                    physical_entries.push_back({n, qj.pen_q, penalty});
                }
            else //Otherwise, embed a 4th copy after the 3rd physical qubit
                physical_entries.push_back({qi[2]+1, qj[2]+1, J} );
        }
    }

    Problem GenerateQACChainProblem(const Solver& solver, const ChainProblem& chain_problem, double penalty){
        Problem problem;
        vector<ProblemEntry>& entries = problem.entries();
        const std::set<int>& broken_cells = solver.get_broken_cells();
        unsigned int num_qubits = solver.get_solver_properties()->quantum_solver->num_qubits;
        uint16_t L = chimera_l(num_qubits);
        uint16_t mid_L = L / 2;

        for( uint16_t x = 0; x < L; ++x){
            uint16_t y0 = (x % 2 ? mid_L : 0);
            // Check that the chain is uninterrutped
            bool chain_ok = true;
            for(uint16_t i = 0; i < mid_L; ++i){
                uint16_t yi = y0 + i;
                chimera_cell ci{L, x, yi, 0, 0};
                uint32_t cci = ci.c();
                if( broken_cells.find(cci) != broken_cells.end()){
                    cout << "Ignoring broken cell "<< cci << "\n";
                    cout << "The chain at x = " << x << ", "
                        << " y = [" << y0 << ", " << y0 + mid_L - 1 << "] will not be embedded\n";
                    chain_ok = false;
                    break;
                }
            }
            if( not chain_ok){
                continue;
            }
            for( const ProblemEntry& p : chain_problem){
                if( p.i >= 8 || p.j >= 8){
                    throw runtime_error("Cannot embed chains longer than 8 qubits with QAC");
                }
                uint16_t yi = y0 + uint16_t(p.i);
                uint16_t yj = y0 + uint16_t(p.j);

                chimera_cell ci{L, x, yi, 0, 0};
                chimera_cell cj{L, x, yj, 0, 0};
                uint32_t cci = ci.c();
                uint32_t ccj = cj.c();
                if( broken_cells.find(cci) != broken_cells.end()){
                    cout << "Ignoring broken cell "<< cci << "\n";
                    continue;
                }
                if( broken_cells.find(ccj) != broken_cells.end()){
                    cout << "Ignoring broken cell " << ccj << "\n";
                    continue;
                }
                EmbedQACProblemEntry(p.value, ci, cj, entries, penalty);
            }
        }

        return problem;
    }


    Problem GenerateCellProblem(const CellProblem& cell_problem,
                                const Solver &solver, const std::set<int> &cell_locations){

        Problem problem;
        vector<ProblemEntry>& entries = problem.entries();

        const std::set<int>& broken_cells = solver.get_broken_cells();
        vector<int> ignored_cells;
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
            } else {
                ignored_cells.push_back(i);

            }
        }
        if(! ignored_cells.empty()){
            cout << "Ignoring broken cells: ";
            for(int i : ignored_cells){
                cout << i << ", ";
            }
            cout << endl;
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


    std::vector< vector<int8_t> > ReadCellProblem(const vector<int8_t> & solution_vec,
            const Solver& solver, const std::set<int>& solver_cells){
        std::vector< vector<int8_t> > instance_results;
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

    // Decodes a solution vector in ising spins according to QAC decoding with 3 physical qubits
    // and 2 logical qubits per cell
    // Returns a vector of 2 * L * L decoded qubits by bipartition indexing (2*(X + L*Y) + B)
    vector<int8_t> DecodeQACProblem(const Solver& solver, const vector<int8_t> & solution_vec){
        unsigned int num_qubits = solver.get_solver_properties()->quantum_solver->num_qubits;
        uint16_t L = chimera_l(num_qubits);
        // Constructs the array with invalid qubits (q=3) by default
        vector<int8_t> instance_results(2 * L * L, 3);
        for( uint16_t x = 0; x < L; ++x){
            for( uint16_t y = 0; y < L; ++y){
                for(uint8_t b = 0; b < 2; ++b){
                    int8_t lq = 0;
                    chimera_cell c0{L, x, y, b, 0};
                    uint32_t n0 = c0.n();
                    bool all_qubits_valid = true;
                    for(uint8_t i = 0; i < 3; ++i){
                        short int qi = solution_vec[n0 + i];
                        if(qi == 3){
                            all_qubits_valid = false;
                            break;
                        }
                        lq += qi;
                    }
                    if( all_qubits_valid){ // All qubits in the bipartition were valid
                        instance_results[c0.b()] = (lq > 0 ? 1 : -1); //Write
                    }
                }
            }
        }
        return instance_results;
    }

    vector<int16_t> ReadVerticalChainProblem(const Solver& solver, const vector<int8_t> & solution_vec,
                                             uint16_t chain_len, bool ignore_invalid){
        if( chain_len > 8){
            chain_len = 8;
            cout << "Warning: chain_len set to 8 (the maximum allowed) for decoding";
        }

        unsigned int num_qubits = solver.get_solver_properties()->quantum_solver->num_qubits;
        uint16_t L = chimera_l(num_qubits);
        // Constructs the array with invalid qubits (q=3) by default
        vector<int16_t> instance_results;
        instance_results.reserve(2 * 4 * L);

        for(uint8_t i = 0; i < 3; ++i){
            for( uint16_t x = 0; x < L; ++x){
                for(uint8_t h = 0; h < 2; ++h){
                    int16_t st = 0;
                    uint16_t y0 = h * L / 2;
                    for( uint16_t dy = 0; dy < chain_len; ++dy){
                        chimera_cell c{L, x, uint16_t(y0 + dy), 0, i};
                        int8_t q = solution_vec[c.n()];
                        if( not(q == 1 or q == -1)) {
                            st = -1;
                            break;
                        }
                        uint b = ( q == -1 ? 1 : 0);
                        st |= (b << unsigned(7 - dy));
                    }
                    if( st == -1 and ignore_invalid)
                        continue;
                    instance_results.push_back(st);
                }
            }
        }
        return instance_results;
    }

    //Reads the results of vertically QAC embedded 8-qubit chains
    vector< int16_t > ReadBipartEmbedChains(Solver& solver, const vector<int8_t>& count_bipartite_decode,
                                            uint16_t chain_len, bool ignore_invalid){
        if( chain_len > 8){
            chain_len = 8;
            cout << "Warning: chain_len set to 8 (the maximum allowed) for decoding";
        }

        unsigned int num_qubits = solver.get_solver_properties()->quantum_solver->num_qubits;
        uint16_t L = chimera_l(num_qubits);
        uint16_t mid_L = L / 2;
        vector< int16_t > chain_counts;

        for( uint16_t x = 0; x < L; ++x){
            int16_t st = 0;
            uint16_t y0 = (x % 2 ? mid_L : 0);
            for( uint16_t dy = 0; dy < chain_len; ++dy){
                uint16_t y = y0 + dy;
                chimera_cell c{L, x, y, 0, 0};
                int8_t q = count_bipartite_decode[c.b()];
                if( not(q == 1 or q == -1)) {
                    st = -1;
                    break;
                }
                uint b = ( q == -1 ? 1 : 0);
                st |= (b << unsigned(7 - dy));
            }
            if( st == -1 and ignore_invalid)
                continue;
            chain_counts.push_back(st);
            //counts[st] += 1;
        }

        return chain_counts;
    }

    vector<int16_t> count_cell_states(
            const std::vector< vector<int8_t> >& cell_readouts,
            std::map<int16_t, int>& counts,
            int& total){
        int16_t st = 0;
        short c;

        vector<int16_t> arr(cell_readouts.size());

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
            //counts[st] += 1;
            total += 1;
        }
        return arr;
    }


}