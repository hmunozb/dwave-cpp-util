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

using namespace std;
namespace po = boost::program_options;
namespace fs = boost::filesystem;

const int CHIMERA_C16_L = 16;
auto CHIMERA_C16_SOLVER = "C16";

enum CellProb{
    none=0,
    c0989v1,
    c0989v6,
    fer
};
enum SchedType{
    def=0,
    lin,
    rev,
    beta,
    pl2b,
    pl3,
    ps3
};

bool cmp_tgts( short s, const string& tgt_str){
    bool eq = true;
    string rev = tgt_str;
    std::reverse(rev.begin(), rev.end());
    for(char c : rev){
        if( c == '0'){
            eq = eq && !(s & 1 );
        } else if ( c == '1'){
            eq = eq && (s & 1);
        } else if ( c == '3'){

        } else {
            throw runtime_error("Invalid tgt_str character in "+tgt_str);
        }

        if(!eq) break;
        s >>= 1;
    }

    return eq;
}

vector<vector<short> > run_case_0989(dwave_cpp::Solver& solver,
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

vector<vector<short> > run_reverse_anneal_gadget(dwave_cpp::Solver& solver,
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

struct bound_canc_0989_prog{
    string url;
    string token;
    bool verbose=false;
    string cell_file;
    string cell_locations_file;
    string reverse_init_file;
    string reverse_init_cell;
    bool write_cells_states=false;
    bool prompt_retries=false;
    vector<int> cell_locations_vec;
    vector<int> reverse_init_cell_vec;
    vector<string> cell_tgts_strs;
    vector<int> cell_tgts_vecs;
    string solver_str = CHIMERA_C16_SOLVER;
    string sched_str;
    SchedType sched_type=def;
    CellProb cell_prob = CellProb::c0989v1;
    double J_scale = -1.0;

    double timeout = 0.0;
    int a = 0;
    int b = 0;
    int n = 1;
    int reps = 1;
    double tf = 0.0; //Default anneal time
    double tw = 0.0;
    double t1 = 0.0;
    double s1 = 0.0;
    double sq = -1.0;
    double sc = 0.0;
    pair<double, double> pl1;
    pair<double, double> pl2;
    string output_file_pref;
    char output_separator = ' ';
    bool rand_gauges=false;

    int main(int argc, const char* argv[]);
    int parse(int argc, const char* argv[]);
    void write_all(
            const vector<dwave_cpp::AnnealSchedulePoint>& sched,
            const std::map<short, int>& counts,
            const std::map<string, int>& tgt_counts,
            const std::vector< vector<short> >& read_arr,
            int num_readouts);
};


int bound_canc_0989_prog::parse(int argc, const char **argv) {
    // Read and handle command line arguments

    po::options_description desc("Allowed options");
    po::positional_options_description pd;

    desc.add_options()
            ("help", "produce help message")
            ("prompt-retry,R", po::bool_switch(&prompt_retries), "Enable retry prompt if solver connection fails")
            ("verbose,v", po::bool_switch(&verbose), "Enable verbose printing to console")
            ("url", po::value<string>(&url), "Solver URL")
            ("token", po::value<string>(&token), "Solver Token")
            ("sched", po::value<string>(&sched_str), "Schedule Type")
            ("problem", po::value<string>(&cell_file), "Cell problem input file")
            //("reverse-init", po::value<string>(&reverse_init_file), "Initial state input file for reverse anneal")
            ("reverse-init-cell", po::value<string>(&reverse_init_cell), "Initial state input file for reverse anneal")
            ("cellsfile,", po::value<string>(&cell_locations_file), "Cell locations input file")
            ("cells,C", po::value<vector<int>>(&cell_locations_vec)->multitoken(), "List of cell locations")
            ("tgts", po::value<vector<string>>(&cell_tgts_strs)->multitoken(), "List of cell locations")
            ("rand-gauge", po::bool_switch(&rand_gauges), "Use random global gauge")
            ("write-cells", po::bool_switch(&write_cells_states), "Write individual cell states in output file "
                                                           "(May result in large output file)")
            //("a", po::value<int>(&a), "Beta Schedule a")
            //("b", po::value<int>(&b), "Beta Schedule b")
            ("tf", po::value<double>(&tf)->default_value(50.0), "Anneal time (Default 50.0 us)")
            //("sq", po::value<double>(&sq)->default_value(-1.0), "Quenching point (Default 1.0)")
            ("sc", po::value<double>(&sc)->default_value(0.9), "Beta schedule codomain fulcrum (Default 0.9)")
            ("n", po::value<int>(&n)->default_value(20),
             "Number of samples per repetition (Default 20)")
            ("max-j", po::value<double>(&J_scale)->default_value(-1.0), "Maximum energy scale for problem couplings."
                                                                        " Problem is autoscaled if not set.")
            ("reps,r", po::value<int>(&reps)->default_value(1),"Number of repetitions (Default 1)")
            ("output,o", po::value<string>(&output_file_pref), "Output file prefix")
            ("timeout", po::value<double>(&timeout)->default_value(60.0), "Solver Response Timeout (Default 1min)")
        //("compression", po::value<int>(), "set compression level")
            ;
    //pd.add("a", 1).add("b", 1).add("output", 1);
    pd.add("sched", 1).add("output", 1);
    po::variables_map vm;
    po::store(
            po::command_line_parser(argc, argv)
                    .options(desc).positional(pd).run(),
            vm);
    po::notify(vm);

    if (vm.count("help")) {
        cout << desc << "\n";
        return 1;
    }
    if (!vm.count("url")) {
        throw std::runtime_error("Must specify solver url with --url");
    }
    if (!vm.count("token")) {
        throw std::runtime_error("Must specify solver API token with --token");
    }
    if (!vm.count("output")) {
        throw std::runtime_error("Must specify an output file");
    }

    if(vm.count("tgts")){
        for(string& tgt : cell_tgts_strs){
            for(const char& c : tgt)
                if( !( (c == '0') || (c == '1') || (c == '3'))){
                    throw std::runtime_error("Invalid target string: "+tgt);
                }
            if(tgt.size() != 8){
                throw std::runtime_error("Invalid target string length: "+tgt);
            }
        }
    }

    fs::path output_path(output_file_pref);
    auto output_dir = output_path.parent_path();
    if (not(fs::exists(output_dir))) {
        throw std::runtime_error("Could not find output directory ");
    }
    auto output_ext = output_path.extension();
    if (output_ext.compare(".tsv") == 0 or output_ext.compare(".txt") == 0) {
        output_separator = '\t';
    } else if (output_ext.compare(".csv")) {
        output_separator = ',';
    } else {
        cout << "Note: file extension " << output_ext
             << "not recognized. Outputting tap-separated data to the destination." << endl;
        output_separator = '\t';
    }

    regex re(R"((\w+)(\{\s*([.[:s:][:alnum:]]*)\})?)");
    //regex args_re(R"(\{\})");
    smatch match;
    regex_match(sched_str,  match, re);
    if(!match.ready()){
        throw runtime_error("Failed to parse schedule specifier "+ sched_str);
    }
    string sched_type_str = match[1].str();
    auto sched_args_m = match[3];
    string sched_args_str;
    vector<string> sched_args_v;
    if(sched_args_m.matched){
        sched_args_str = sched_args_m.str();
        boost::char_separator<char> sep{" ", ","};
        boost::tokenizer<boost::char_separator<char> > tok(sched_args_str, sep);

        for( auto& arg_token : tok){
            sched_args_v.push_back(arg_token);
        }
    } else{
        throw runtime_error("Failed to parse schedule arguments"+ sched_str);
    }

    if(sched_type_str == "lin"){
        sched_type = lin;

    } else if(sched_type_str == "rev"){
        sched_type = rev;
        if(vm.count("reverse-init-cell")){
            long r = stol(reverse_init_cell, nullptr, 2);
            for(int i = 0; i < 8; ++i){
                reverse_init_cell_vec.push_back((-2)*(r & 1) + 1);
                r >>= 1;
            }
        } else {
            throw runtime_error("reverse-init-cell required for reverse anneal schedule");
        }
        int num_args = sched_args_v.size();
        if(num_args == 3){
            t1 = stod(sched_args_v[0]);
            s1 = stod(sched_args_v[1]);
            tw = stod(sched_args_v[2]);
        }else {
            throw runtime_error("Invalid arguments for reverse anneal schedule");
        }
    }
    else if(sched_type_str == "pl2b") {
        sched_type = pl2b;
        int num_args = sched_args_v.size();
        switch(num_args){
            case 8:
                sc = stod(sched_args_v[7]);
            case 7:
                sq = stod(sched_args_v[6]);
            case 6:
                a = stoi(sched_args_v[0]);
                b = stoi(sched_args_v[1]);
                pl1.first = stod(sched_args_v[2]);
                pl1.second = stod(sched_args_v[3]);
                pl2.first = stod(sched_args_v[4]);
                pl2.second = stod(sched_args_v[5]);
                break;
            default:
                throw runtime_error("Invalid arguments for beta schedule");
        }
    } else if (sched_type_str == "pl3" || sched_type_str=="ps3"){
        int num_args = sched_args_v.size();
        if (num_args == 5){
            if(sched_type_str == "pl3")
                sched_type = pl3;
            else
                sched_type = ps3;

            pl1.first = stod(sched_args_v[0]);
            pl1.second = stod(sched_args_v[1]);
            pl2.first = stod(sched_args_v[2]);
            pl2.second = stod(sched_args_v[3]);
            tf = stod(sched_args_v[4]);
        } else {
            throw runtime_error("pl3 requires 5 argumnts");
        }
    } else if (sched_type_str == "beta"){
        sched_type = beta;
        int num_args = sched_args_v.size();
        switch(num_args){
            case 4:
                sc = stod(sched_args_v[3]);
            case 3:
                sq = stod(sched_args_v[2]);
            case 2:
                a = stoi(sched_args_v[0]);
                b = stoi(sched_args_v[1]);
                break;
            default:
                throw runtime_error("Invalid arguments for beta schedule");
        }
    } else {
        throw runtime_error("Unrecognized schedule type "+sched_type_str);
    }


    return 0;
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
    if(parse_code != 0){
        return parse_code;
    }

    std::set<int> SolverCells;
    if( cell_locations_vec.empty() ){
        SolverCells = dwave_cpp::CheckerboardCellSet(CHIMERA_C16_L);
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


    //auto c0989_cell = dwave_cpp::C0989_V1();
    if(verbose) cout << "Import cell problem from " << cell_file << endl;
    dwave_cpp::CellProblem cell_problem;
    ifstream ifs(cell_file);
    ifs >> cell_problem;
    if(J_scale > 0){
        for(auto& p : cell_problem){
            p.value /= J_scale;
        }
    }

    cout << "Generating schedule...\n";
    vector<dwave_cpp::AnnealSchedulePoint> sched;
    switch(sched_type){
        case def:
            cout << "Schedule specification " << sched_str << " not recognized" << endl;
            return 1;
        case lin:
            cout << "\t***Linear***";
            sched = dwave_cpp::linear_schedule(tf);
            break;
        case rev:
            cout << "\t***Reverse Anneal**";
            cout << "\tInitial State: " << reverse_init_cell << "\n";
            sched = dwave_cpp::reverse_anneal_schedule(s1, t1, tw);
            break;
        case pl2b:
            cout << "\t***2 Piecewise Linear with Quenched Beta***";
            sched = dwave_cpp::pl2_beta_quenched_schedule(pl1, pl2, tf, a, b, sq, sc);
            break;
        case pl3:
            cout << "\t***3 Piecewise Linear***";
            sched = dwave_cpp::pl3_schedule(pl1, pl2, tf);
            break;
        case ps3:
            cout << "\t***3 Piecewise Slopes***\n";
            sched = dwave_cpp::ps3_schedule(pl1, pl2, tf);
            break;
        case beta:
            cout << "Schedule type: Beta " << a << ", " << b << endl;
            if(sq < 0.0){
                cout << "\tUnquenched, r=" << sc << endl;
                sched = dwave_cpp::beta_schedule(tf, sc, a, b);
            } else {
                cout << "\tQuenched at " << sq*100 << "% , r=" << sc <<endl;
                sched = dwave_cpp::beta_quenched_schedule(tf, sq, sc, a, b);
            }
            break;
        default:
            throw std::logic_error("Schedule type enum switch failed");
    }
    //auto sched = dwave_cpp::linear_schedule(Prog.tf);
    cout << "\n";
    cout << sched << endl;

    dwave_cpp::Connection connection(url, token, prompt_retries);
    dwave_cpp::Solver solver = connection.get_solver(CHIMERA_C16_SOLVER);



 /*   dwave_cpp::Connection connection(url, token);
    dwave_cpp::Solver solver = connection.get_solver(solver_str);
    //auto SolverCells = dwave_cpp::CheckerboardCellSet(CHIMERA_C16_L);
    std::set<int> SolverCells;
    SolverCells.insert(150);
    vector<dwave_cpp::Problem> prob_0989;
    //prob_0989.push_back(dwave_cpp::generate_0989_array(solver, CHIMERA_C16_L));
    prob_0989.push_back(dwave_cpp::GenerateCellProblem(dwave_cpp::C0989_V1(), solver, SolverCells));
    dwave_cpp::QuantumSolverParameters params;

    double anneal_time = tf;
    params.get()->answer_mode = dwave_cpp::sapi_SolverParameterAnswerMode::SAPI_ANSWER_MODE_RAW;
    params.get()->num_spin_reversal_transforms = 0;
    params.get()->num_reads = n;
   // auto sched = dwave_cpp::beta_schedule(anneal_time, 0.90, 3, 3);
    //auto sched = dwave_cpp::beta_quenched_schedule(anneal_time, sq, sc, a, b);
    auto sched = dwave_cpp::linear_schedule(anneal_time);
    cout << "Constructed schedule:\n";
    for( auto p : sched){
        cout << '\t' << p.time << '\t' << p.relative_current << '\n';
    }
    cout << endl;

    dwave_cpp::ProblemSubmission problem_submission;
    cout << "Submitting..." << endl;

    problem_submission.asyncSolveIsing(solver, prob_0989, params);

    cout << "Waiting for results..." << endl;
    problem_submission.await(timeout);
    problem_submission.fetch_done();
    //vector<dwave_cpp::sapi_IsingResult*>& ising_results = problem_submission.get_results();
    auto results_vec =problem_submission.results_vector();
    cout << "Done." << endl;
    cout << "Results: "<<endl;*/

    int num_readouts = 0;
    int num_tgts = cell_tgts_strs.size();
    vector< vector<short> > read_arr;
    std::map<short, int> counts;
    std::vector<int> tgt_count_vec(num_tgts, 0);
    std::map<string, int> tgt_counts;

    for(int r = 0; r < reps; ++r) {
        cout << r+1 << "/" << reps << "..." << endl;
        vector<vector<short> > results_vec;

        if(sched_type == rev)
            results_vec = run_reverse_anneal_gadget(solver, n, SolverCells, cell_problem, sched, reverse_init_cell_vec);
        else
            results_vec = run_case_0989(solver, n, SolverCells, cell_problem, sched, rand_gauges ? 1 : 0, timeout,
                    J_scale < 0.0);

        for (const auto &result : results_vec) {
            auto readouts = dwave_cpp::ReadCellProblem(result, solver, SolverCells);

            vector<short> arr = dwave_cpp::count_cell_states(readouts, counts, num_readouts);

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
    return Prog.main(argc, argv);
}