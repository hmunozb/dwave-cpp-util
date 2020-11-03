//
// Created by Humberto Munoz Bauza on 2019-05-14.
//

#include "dwave_cpp/core/dwave_cpp_core.h"
#include "dwave_cpp/core/parameters.h"
#include "dwave_cpp/core/properties.h"
#include "dwave_cpp/core/solve.h"
#include "dwave_cpp/core/energy.h"
#include <iostream>
#include <sstream>
#include <fstream>
#include <string>
#include <tuple>
#include <cassert>

namespace dwave_cpp{

extern "C" {
    #include <dwave_sapi.h>
}


    ErrorMessage::ErrorMessage() {
        msg_buffer = new char[SAPI_ERROR_MESSAGE_MAX_SIZE];
    }

    char *ErrorMessage::get_str_buffer() {
        return msg_buffer;
    }

    const string &ErrorMessage::get_str() {
        msg = msg_buffer;
        return msg;
    }

    struct connection_result{
        sapi_Code  code;
        sapi_Connection* con_ptr;
    };

    connection_result create_connection(const char* url, const char* token, const char* proxy,
            ErrorMessage& err_msg){
        connection_result cr;
        cr.code =
                sapi_remoteConnection(url, token, proxy,
                                      &cr.con_ptr, err_msg.get_str_buffer() );
        return cr;
    }

    Connection::Connection(const string& url, const string& token, const string& proxy) {
        sapi_Code err_code;
        err_code = sapi_globalInit();
        ErrorMessage err_msg;
        if(err_code != SAPI_OK){
            //std::cout << "dwave_sapi failed to initialize" << std::endl;
            throw exception();
        }

        connection_result cr = create_connection(url.c_str(), token.c_str(), proxy.c_str(), err_msg);

        if(cr.code != SAPI_OK){
            throw exception();
        }
        connection = cr.con_ptr;

    }
    bool _handle_connection_failure(bool handle_interactive){
        if(not handle_interactive)
            return false;

        cout << "Failed to establish connection. Try again? (y/n) :";
        cout.flush();
        char c = 0;
        do{
            cin >> c;
            if(c == 'y' || c == 'Y')
                return true;
            else if (c == 'n' || c == 'N')
                return false;
            cin.ignore(INTMAX_MAX);
            cout << "(y/n)";
        } while ( not(c == 'y' || c == 'Y' || c == 'n' || c == 'N') );
        assert(false);
    }

    Connection::Connection(const string &url, const string &token, bool enable_interactive_retry) {
        sapi_Code err_code;
        err_code = sapi_globalInit();
        ErrorMessage err_msg;
        if(err_code != SAPI_OK){
            //std::cout << "dwave_sapi failed to initialize" << std::endl;
            throw std::runtime_error("Failed to initialize dwave_sapi");
        }

        do{
        //Attempt connection three times
            for(int i = 0; i < 3; ++i) {
                connection_result cr = create_connection(url.c_str(), token.c_str(), nullptr, err_msg);
                if(cr.code == SAPI_OK){
                    connection = cr.con_ptr;
/*
 * * * * * * * * * * * * * * * * * * * * *
 * NORMAL EXIT POINT FOR THIS CONSTRUCTOR*
 * * * * * * * * * * * * * * * * * * * * *
 */
                    return;
                }
                else
                    cout << "Connection failed. Retrying..." << endl;
            }
        } while ( _handle_connection_failure(enable_interactive_retry) );

        throw std::runtime_error("Failed to establish connection to " + url);
//        if (cr.code != SAPI_OK) {
//            sapi_globalCleanup(); //needs to be called since destructor is not called
//            throw std::runtime_error("Failed to establish connection to " + url);
//        }

    }
    Connection::~Connection() {
        sapi_freeConnection(connection);

        sapi_globalCleanup();
    }

    vector<string> Connection::solver_list(){
        const char ** lst = sapi_listSolvers(connection);
        vector<string> solver_vec;
        for (int i = 0; lst[i] != nullptr; ++i){
            solver_vec.emplace_back(lst[i]);
        }
        return solver_vec;
    }


    std::ostream& operator<<(std::ostream& os, const AnnealSchedule& sched){
        os  << '#'
            << "sched" << '\n';
        for( const auto& p : sched){
            os << os.fill() << ' ' << p.time;
        }
        os << '\n';
        for( const auto& p : sched){
            os << os.fill() << ' ' << p.relative_current;
        }

        return os;
    }

    Solver Connection::get_solver(const string &solver_name) {
        sapi_Solver* solver_ptr = sapi_getSolver(connection, solver_name.c_str());
        const sapi_SolverProperties* props_ptr = sapi_getSolverProperties(solver_ptr);

        Solver solver(solver_ptr, props_ptr);

        return solver;
    }



    Solver::Solver(sapi_Solver * sapi_solver_ptr, const dwave_cpp::sapi_SolverProperties * sapi_props_ptr)
        : solver_ptr(sapi_solver_ptr), props_ptr(sapi_props_ptr){
        calculate_broken_units();
    }
    const sapi_Solver* Solver::get() {return solver_ptr;}

    Solver::~Solver() {
        if(solver_ptr != nullptr)
            sapi_freeSolver(solver_ptr);
    }
    const sapi_SolverProperties* Solver::get_solver_properties() const{
        return props_ptr;
    }
    bool Solver::is_quantum() {
        return (bool)(props_ptr->quantum_solver);
    }


    const std::set<int>& Solver::get_broken_qubits() const{
        return broken_qubits;
    }
    const std::set<int>& Solver::get_broken_cells() const {
        return broken_cells;
    }

    void Solver::calculate_broken_units() {
        broken_qubits.clear();
        broken_cells.clear();

        const sapi_SolverProperties* props = props_ptr;
        const int* qubits = props->quantum_solver->qubits;
        int num_qubits = props->quantum_solver->qubits_len;

        //Gather the indices of broken qubits
        int expected_qubit = 0;

        for(int i = 0; i < num_qubits; ++i){
            int q = qubits[i];
            assert(q >= expected_qubit);
            if(q != expected_qubit){
                //Possibly more than one broken qubit in a row
                for( int bq = expected_qubit; bq < q; ++bq)
                    broken_qubits.insert(bq);
                expected_qubit = q + 1;
            } else{
                ++expected_qubit;
            }
        }

        //Calculate the indices of cells with  at least one broken qubit
        for(int q : broken_qubits){
            broken_cells.insert( q / 8 );
        }

    }




    Problem::Problem() : problem_ptr(new sapi_Problem) { }
    //Problem::Problem(unsigned int num_entries) : problem_entries(num_entries), problem_ptr(new sapi_Problem) {}

    vector<ProblemEntry>& Problem::entries() { return problem_entries;}
    sapi_Problem* Problem::prepare_problem() {
        problem_ptr->len = problem_entries.size();
        problem_ptr->elements = problem_entries.data();
        return problem_ptr.get();
    }

    Problem::Problem(ProblemAdj && problem_adj) :
        problem_entries(problem_adj), problem_ptr(new sapi_Problem){

    }

    Problem::Problem(const Problem& problem) : problem_entries(problem.problem_entries), problem_ptr(new sapi_Problem){

    }

    QuantumSolverParameters::QuantumSolverParameters()
            : q_sol_params_ptr(new sapi_QuantumSolverParameters(SAPI_QUANTUM_SOLVER_DEFAULT_PARAMETERS)){

    }

    sapi_QuantumSolverParameters* QuantumSolverParameters::get() {
        return q_sol_params_ptr.get();
    }
    sapi_QuantumSolverParameters* QuantumSolverParameters::operator->(){
        return q_sol_params_ptr.get();
    }
    void QuantumSolverParameters::set_anneal_schedule(const std::vector<AnnealSchedulePoint>& schedule){
        int num_pnts = schedule.size();
        if(num_pnts == 0){
            throw runtime_error("set_anneal_schedule(): Schedule cannot be empty.");
        }

        anneal_schedule_ptr = std::make_unique<sapi_AnnealSchedule>();
        anneal_schedule_pnts_ptr = std::make_unique<AnnealSchedulePoint[]>(num_pnts);
        for(int i = 0; i < num_pnts; ++i){
            anneal_schedule_pnts_ptr[i] = schedule[i];
        }
        anneal_schedule_ptr->len = num_pnts;
        anneal_schedule_ptr->elements = anneal_schedule_pnts_ptr.get();

        q_sol_params_ptr->anneal_schedule = anneal_schedule_ptr.get();
    }
    void QuantumSolverParameters::set_reverse_anneal(std::vector<int>& init_state, bool reinitialize) {
        reverse_anneal_ptr = std::make_unique<sapi_ReverseAnneal>();
        reverse_anneal_ptr->initial_state = init_state.data();
        reverse_anneal_ptr->initial_state_len = init_state.size();
        reverse_anneal_ptr->reinitialize_state = reinitialize ? 1 : 0;

        q_sol_params_ptr->reverse_anneal = reverse_anneal_ptr.get();
    }

// solve.h

    ProblemSubmission::~ProblemSubmission() {
        for(sapi_SubmittedProblem* ptr : submitted_problem_ptr){
            if(ptr != nullptr)
                sapi_freeSubmittedProblem(ptr);
        }
        for(sapi_IsingResult* ptr : ising_results){
            if(ptr != nullptr)
                sapi_freeIsingResult(ptr);
        }
    }

    void ProblemSubmission::asyncSolveIsing(Solver& solver, vector<Problem>& problem, QuantumSolverParameters& params){
        int num_problems = problem.size();
        submitted_problem_ptr.resize(num_problems);
        ising_results.resize(num_problems);

        const sapi_Solver* solver_ptr = solver.get();
        auto params_ptr = ( sapi_SolverParameters *) params.get();

        for(int i = 0; i < num_problems; ++i){
            sapi_Problem * problem_ptr = problem[i].prepare_problem();
            sapi_Code err_code;
            err_code = sapi_asyncSolveIsing(solver_ptr, problem_ptr, params_ptr,
                    & submitted_problem_ptr[i], err_msg.get_str_buffer());
            if(err_code != sapi_Code::SAPI_OK){
                cerr << "Problem Submission Failed: " << err_msg.get_str() << endl;
                throw runtime_error(err_msg.get_str());
            }
        }
    }

    int ProblemSubmission::await(double timeout) {
        return sapi_awaitCompletion((const sapi_SubmittedProblem **) submitted_problem_ptr.data(),
                submitted_problem_ptr.size(), submitted_problem_ptr.size(), timeout);
    }
    void ProblemSubmission::fetch_done() {
        int num_problems = submitted_problem_ptr.size();
        for(int k =0; k < num_problems; ++k){
            int done = sapi_asyncDone(submitted_problem_ptr[k]);
                sapi_Code code = sapi_asyncResult(submitted_problem_ptr[k], &ising_results[k], err_msg.get_str_buffer());
                //sapi_cancelSubmittedProblem(submitted_problem_ptr[k]);
                if(code != sapi_Code::SAPI_OK){
                    cerr << "Problem Fetch Failed: " << err_msg.get_str() << endl;
                    throw runtime_error(err_msg.get_str());
                }
        }
    }

    vector<sapi_IsingResult*>& ProblemSubmission::get_results() {
        return ising_results;
    }

    ResultsVec ProblemSubmission::results_vector() {
        ResultsVec results_vec;
        for(sapi_IsingResult* ir: ising_results){
            int num_sols = ir->num_solutions;
            int sol_len = ir->solution_len;
            vector<vector<int8_t> > sol_vec;
            for(int i = 0; i < num_sols; ++i){
                int* sol_begin = &ir->solutions[i*sol_len];
                int* sol_end = &ir->solutions[(i+1)*sol_len];

                sol_vec.emplace_back(sol_begin, sol_end);
            }
            results_vec.push_back(std::move(sol_vec));
        }
        return results_vec;
    }


    dwave_cpp::ProblemSubmission run_problem_vector(dwave_cpp::Solver& solver, vector<dwave_cpp::Problem>& prob_vec,
                                             dwave_cpp::QuantumSolverParameters& params, double timeout){
        dwave_cpp::ProblemSubmission problem_submission;
        problem_submission.asyncSolveIsing(solver, prob_vec, params);

        problem_submission.await(timeout);
        problem_submission.fetch_done();
        //auto results_vec = problem_submission.results_vector();

        return problem_submission;
    }

    std::istream& operator>>(std::istream& in, ProblemEntry& problem_entry){
        in  >> problem_entry.i
            >> problem_entry.j
            >> problem_entry.value;
        return in;
    }

    std::istream& operator>>(std::istream& in, ProblemAdj& cell_problem){
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

    void import_canary_problem(string file, ProblemAdj& canary_problem, AnnealSchedule& sched ){
        ifstream ifs(file);
        std::string line;

        sched.clear();
        sched.push_back({0.0, 0.0});

        if(std::getline(ifs, line)) {
            std::istringstream iss( line );
            double sq, tp;
            if(iss >> sq >> tp){
                if(sq > 0.0 && sq < 1.0 )
                {   sched.push_back({1.0, sq});
                    sched.push_back({1.0 + tp, sq});
                    sched.push_back({2.0+tp, 1.0});
                } else if (sq == 1.0){
                    sched.push_back({tp, 1.0});
                } else {
                    throw runtime_error("Invalid schedule for canary problem.");
                }

            } else {
                throw runtime_error("Failed to read canary problem.");
            }
        }
        ifs >> canary_problem;
    }

    double EnergyEval::operator()(const vector<int8_t> &ising_spins) {
        double e = 0.0;
        for(const ProblemEntry& p : problem){
            e +=  p.value * (p.i == p.j
                                ?   ising_spins.at(p.i)
                                :   ising_spins.at(p.i) * ising_spins.at(p.j) ) ;
        }

        return e;
    }

    ProblemAdj import_problem(string input_file, double J_scale){
        ProblemAdj problem_adj;
        ifstream ifs(input_file);
        ifs >> problem_adj;
        if(J_scale > 0){
            for(auto& p : problem_adj){
                p.value /= J_scale;
            }
        }
        return problem_adj;
    }
}