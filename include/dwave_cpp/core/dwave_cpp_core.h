//
// Created by Humberto Munoz Bauza on 7/13/18.
//

#ifndef DW_EXAMPLES_DWAVE_CPP_CORE_H
#define DW_EXAMPLES_DWAVE_CPP_CORE_H

#include <vector>
#include <string>
#include <memory>
#include <set>

using namespace std;

namespace dwave_cpp{

    struct sapi_Connection;
    struct sapi_Solver;
    struct sapi_SolverProperties;
    struct sapi_Problem;
    struct sapi_ProblemEntry;
    struct sapi_SupportedProblemTypeProperty;

    class Solver;
    class QuantumSolver;

    struct ErrorMessage{
        ErrorMessage();

        char * get_str_buffer();
        const string& get_str();
        char * msg_buffer;
        string msg;
        int err_code;
    };

    class Connection{
    public:
        Connection(const string& url, const string& token, const string& proxy );
        Connection(const string& url, const string& token, bool enable_interactive_retry=false);
        vector<string> solver_list();
        Solver get_solver(const string& solver_name);
        ~Connection();
        const sapi_Connection * get_ptr() const{
            return connection;
        }
    private:
        sapi_Connection* connection = nullptr;
    };

    typedef sapi_ProblemEntry ProblemEntry ;

    class Problem{
    public:
        Problem();
        Problem(Problem&&) = default;
        Problem(const Problem&) = delete;
        //Problem(unsigned int num_entries);
        sapi_Problem* prepare_problem();
        vector<ProblemEntry>& entries();
    private:
        vector<ProblemEntry> problem_entries;
        std::unique_ptr<sapi_Problem> problem_ptr;

    };



    class Solver{
    public:
        ~Solver();

        const sapi_SolverProperties* get_solver_properties() const;
        const sapi_Solver* get();
        bool is_quantum();

        const std::set<int>& get_broken_qubits() const;
        const std::set<int>& get_broken_cells() const;
        friend class Connection;

    protected:
        Solver(sapi_Solver * solver_ptr, const sapi_SolverProperties * props_ptr);

    private:
        void calculate_broken_units();

        sapi_Solver * solver_ptr;
        const sapi_SolverProperties * props_ptr;

        std::set<int> broken_qubits;
        std::set<int> broken_cells;
    };


}
#endif //DW_EXAMPLES_DWAVE_CPP_CORE_H
