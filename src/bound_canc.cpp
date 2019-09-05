//
// Created by Humberto Munoz Bauza on 2019-05-17.
//

#include <dwave_cpp/dwave_cpp.h>
#include <dwave_cpp/problems/0989.h>
#include <dwave_cpp/schedules/beta.h>
#include <iostream>
#include <fstream>
#include <map>
#include <algorithm>
#include <boost/program_options.hpp>
#include <boost/filesystem.hpp>
#include <boost/algorithm/string.hpp>
#include <boost/tokenizer.hpp>
#include <regex>
#include <chrono>

#include "dw_prog.h"

using namespace std;
namespace po = boost::program_options;
namespace fs = boost::filesystem;


vector<vector<int8_t>> run_single_problem(dwave_cpp::Solver& solver,dwave_cpp::Problem& problem,
        dwave_cpp::QuantumSolverParameters& params,double timeout){

    vector<dwave_cpp::Problem> prob_vec;

    auto results_vec = run_problem_vector(solver, prob_vec, params, timeout);

    return results_vec[0];
}

vector<vector<int8_t> > run_case_0989(dwave_cpp::Solver& solver,
        int n,
        const std::set<int>& SolverCells,
        const dwave_cpp::CellProblem& cell_problem,
        const vector<dwave_cpp::AnnealSchedulePoint>& sched,
        int num_spin_reversals=0,
        double timeout=60.0,
        bool autoscale=true){

    vector<dwave_cpp::Problem> prob_0989;
    prob_0989.push_back(dwave_cpp::GenerateCellProblem(cell_problem, solver, SolverCells));

    dwave_cpp::QuantumSolverParameters params;

    params->answer_mode = dwave_cpp::sapi_SolverParameterAnswerMode::SAPI_ANSWER_MODE_RAW;
    params->num_spin_reversal_transforms = num_spin_reversals;
    params->num_reads = n;
    params->auto_scale = (autoscale ? 1 : 0);
    //params->reduce_intersample_correlation = 0;
    params.set_anneal_schedule(sched);

    dwave_cpp::ProblemSubmission problem_submission;

    problem_submission.asyncSolveIsing(solver, prob_0989, params);

    problem_submission.await(timeout);
    problem_submission.fetch_done();
    auto results_vec = problem_submission.results_vector();

    return results_vec[0];
}

vector<vector<int8_t> > run_reverse_anneal_gadget(dwave_cpp::Solver& solver,
                                     int n,
                                     const std::set<int>& SolverCells,
                                     const dwave_cpp::CellProblem& cell_problem,
                                     const vector<dwave_cpp::AnnealSchedulePoint>& sched,
                                     const vector<int>& reverse_anneal_init_cell,
                                     int num_spin_reversals=0,
                                     double timeout=60.0,
                                     bool autoscale=true){

    vector<dwave_cpp::Problem> prob_gadget;
    vector<int> reverse_anneal_init_state;
    prob_gadget.push_back(dwave_cpp::GenerateCellProblem(cell_problem, solver, SolverCells));
    reverse_anneal_init_state = dwave_cpp::GenerateCellReverseInit(solver, SolverCells, reverse_anneal_init_cell);

    dwave_cpp::QuantumSolverParameters params;

    params->answer_mode = dwave_cpp::sapi_SolverParameterAnswerMode::SAPI_ANSWER_MODE_RAW;
    params->num_spin_reversal_transforms = num_spin_reversals;
    params->num_reads = n;
    params->auto_scale = (autoscale ? 1 : 0);
    //params->reduce_intersample_correlation = 0;
    params.set_anneal_schedule(sched);
    params.set_reverse_anneal(reverse_anneal_init_state);

    dwave_cpp::ProblemSubmission problem_submission;

    problem_submission.asyncSolveIsing(solver, prob_gadget, params);

    problem_submission.await(timeout);
    problem_submission.fetch_done();
    auto results_vec = problem_submission.results_vector();

    return results_vec[0];
}



vector<vector<int8_t> > run_qac_chain_gadget(dwave_cpp::Solver& solver,
        int n,
        const std::set<int>& SolverCells,
        const dwave_cpp::ChainProblem& chain_problem,
        const vector<dwave_cpp::AnnealSchedulePoint>& sched,
        const vector<int>& reverse_anneal_init_cell,
        int num_spin_reversals=0,
        double timeout=60.0,
        bool autoscale=true){

    vector<dwave_cpp::Problem> prob_vec;
    prob_vec.push_back(dwave_cpp::GenerateQACChainProblem(chain_problem, solver));

    dwave_cpp::QuantumSolverParameters params;
    params->answer_mode = dwave_cpp::sapi_SolverParameterAnswerMode::SAPI_ANSWER_MODE_RAW;
    params->num_spin_reversal_transforms = num_spin_reversals;
    params->num_reads = n;
    params->auto_scale = (autoscale ? 1 : 0);

    params.set_anneal_schedule(sched);

    dwave_cpp::ProblemSubmission problem_submission;

    problem_submission.asyncSolveIsing(solver, prob_vec, params);

    problem_submission.await(timeout);
    problem_submission.fetch_done();
    auto results_vec = problem_submission.results_vector();

    return results_vec[0];

}


struct bound_canc_0989_prog: public gadget_program{
//    string cell_file;
//    string cell_locations_file;
//    bool write_cells_states=false;
//    vector<int> cell_locations_vec;
//    vector<string> cell_tgts_strs;
//    //vector<int> cell_tgts_vecs;

    std::set<int> SolverCells;

    po::options_description desc;

    bound_canc_0989_prog();

    dwave_cpp::Problem encode_problem( dwave_cpp::Solver& solver,
                                       dwave_cpp::CellProblem &cell_problem);
    vector<int16_t> decode_problem(
            dwave_cpp::Solver& solver, const vector<int8_t>& readout);


    int main(int argc, const char* argv[]);
    void write_all(
            const vector<dwave_cpp::AnnealSchedulePoint>& sched,
            const std::map<short, int>& counts,
            const std::map<string, int>& tgt_counts,
            const std::vector< vector<short> >& read_arr,
            int num_readouts);
};

dwave_cpp::Problem bound_canc_0989_prog::encode_problem(dwave_cpp::Solver& solver,
        dwave_cpp::CellProblem &cell_problem) {
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

vector<int16_t> bound_canc_0989_prog::decode_problem(
        dwave_cpp::Solver& solver, const vector<int8_t>& readout) {
    int num_readouts = 0;
    auto sep_readout = dwave_cpp::ReadCellProblem(readout, solver, SolverCells);
    return dwave_cpp::count_cell_states(sep_readout, counts, num_readouts);
}

bound_canc_0989_prog::bound_canc_0989_prog() : gadget_program(){
    po::positional_options_description pd;

//    desc.add_options()
//            ("problem", po::value<string>(&cell_file), "Cell problem input file")
//            ("cellsfile,", po::value<string>(&cell_locations_file), "Cell locations input file")
//            ("cells,C", po::value<vector<int>>(&cell_locations_vec)->multitoken(), "List of cell locations")
//            ("tgts", po::value<vector<string>>(&cell_tgts_strs)->multitoken(), "List of cell locations")
//            ("write-cells", po::bool_switch(&write_cells_states), "Write individual cell states in output file "
//                                                                  "(May result in large output file)")
//            ;
//    //pd.add("sched", 1).add("output", 1);
//
//    program_options.add(desc);

    positional_options.add("sched", 1)
            .add("output", 1);

}


void bound_canc_0989_prog::write_all(
        const vector<dwave_cpp::AnnealSchedulePoint>& sched,
        const std::map<short, int>& counts,
        const std::map<string, int>& tgt_counts,
        const std::vector< vector<short> >& read_arr,
        int num_readouts) {
    {
        string s(1, output_separator);
        ofstream f(output_file_pref);
        //the fill character indicates to operator<<(sched) which separator to use
        f.fill(output_separator);
        f   << "#" << s
            << "n" << s
            << "tf" << s
            << "sq" << s
            << "sc" << s
            << "a" << s
            << "b" << '\n';
        f   << s
            << n << s
            << tf << s
            << sq << s
            << sc << s
            << a << s
            << b << '\n';

//        f   << '#'
//            << "sched" << '\n';
//        for( const auto& p : sched){
//            f << s << p.time;
//        }
//        f << '\n';
//        for( const auto& p : sched){
//            f << s << p.relative_current;
//        }
//        f << '\n';

        f << sched << '\n';

        f   << '#' << s
            << "tgt" << s
            << "count" << s
            << "frac" << '\n';
        for( const auto& p : tgt_counts)
            f << s << p.first
              << s  << p.second
              << s  << double(p.second)/num_readouts << "\n";

        f   << '#' << s
            << "state" << s
            << "count" << s
            << "frac" << '\n';
        for( const auto& p : counts){
            f << s << p.first
              << s  << p.second
              << s  << double(p.second)/num_readouts << "\n";
        }


        if(write_cells_states){
            f   << "# arrs" << '\n';

            for(const auto& result : read_arr) {
                for(short id: result){
                    f << s << id;
                }
                f << '\n';
            }
        }
    }

}

int bound_canc_0989_prog::main(int argc, const char* argv[]) {
    int parse_code = parse(argc, argv);
    if(parse_code ){
        return parse_code;
    }

//    if( cell_locations_vec.empty() ){
//        SolverCells = dwave_cpp::CheckerboardCellSet(16);
//    }
//    else {
//        for(int cell : cell_locations_vec) SolverCells.insert(cell);
//    }
//    if(verbose){
//        cout << "Solver Cells used:\n";
//        for(int cell : SolverCells)
//            cout << cell << ", ";
//        cout << endl;
//    }


    if(verbose) cout << "Import cell problem from " << cell_file << endl;
    dwave_cpp::CellProblem cell_problem = import_cell_problem();

    cout << "Generating schedule...\n";
    generate_schedule();

    cout << "\n";
    cout << sched << endl;

    dwave_cpp::Connection connection(url, token, prompt_retries);
    dwave_cpp::Solver solver = connection.get_solver(solver_str);
    dwave_cpp::QuantumSolverParameters params;
    set_parameters(params);
    set_schedule_parameters(params);

    vector<dwave_cpp::Problem> problem_vector;
    problem_vector.push_back(encode_problem(solver, cell_problem));

    int num_readouts = 0;
    int num_tgts = cell_tgts_strs.size();
    //vector< vector<int16_t> > read_arr;
    //std::map<int16_t, int> counts;
    //std::vector<int> tgt_count_vec(num_tgts, 0);
    //std::map<string, int> tgt_counts;

    for(int r = 0; r < reps; ++r) {
        cout << r+1 << "/" << reps << "..." << endl;
        vector<vector<int8_t> > results_vec;
        results_vec = run_problem_vector(solver, problem_vector, params, timeout)[0];

//        if(sched_type == rev)
//            results_vec = run_reverse_anneal_gadget(solver, n, SolverCells, cell_problem, sched, reverse_init_cell_vec);
//        else
//            results_vec = run_case_0989(solver, n, SolverCells, cell_problem, sched, rand_gauges ? 1 : 0, timeout,
//                    J_scale < 0.0);

        for (const auto &result : results_vec) {
            //auto readouts = decode_problem(solver, result);
            //auto readouts = dwave_cpp::ReadCellProblem(result, solver, SolverCells);

            vector<int16_t> arr = decode_problem(solver, result);
            num_readouts += arr.size();
            read_arr.push_back(std::move(arr));
        }
    }
    cout << endl;


    cout << "Counts: \n";
    for( const auto& p : counts){
        cout << "\t" << p.first
        << "\t\t"   << p.second
        << "\t\t"   << double(p.second)/num_readouts << "\n";
        for(int k = 0; k < num_tgts; ++k){
            auto& tgt = cell_tgts_strs[k];
            tgt_counts[tgt] += 0; //force the tgt_counts map to include zero-count targets
            if(cmp_tgts(p.first, tgt)){
                tgt_count_vec[k] += p.second;
                tgt_counts[tgt] += p.second;
            }
        }
    }
    cout << endl;
    cout << "Target Counts: \n" ;
    for(int k = 0; k < num_tgts; ++k){
        cout << "\t\t" << cell_tgts_strs[k]
            << ":\t\t" << tgt_count_vec[k]
            << "\t\t" << double(tgt_count_vec[k]) / num_readouts << "\n";
    }

    /*
     * OUTPUT TO FILE
     *
     * */
    write_all(sched, counts, tgt_counts, read_arr, num_readouts);
    return 0;
}

int main(int argc, const char* argv[]){
    bound_canc_0989_prog Prog;
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