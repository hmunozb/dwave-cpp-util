//
// Created by Humberto Munoz Bauza on 2019-09-03.
//

#include "dw_prog.h"


const int CHIMERA_C16_L = 16;
auto CHIMERA_C16_SOLVER = "C16";

generic_dwave_program::generic_dwave_program() : program_base(), dwave_opts("D-Wave Solver Configuration"){
    dwave_opts.add_options()
            ("url", boost::program_options::value<string>(&url), "Solver URL")
            ("solver", boost::program_options::value<string>(&solver_str)->default_value(CHIMERA_C16_SOLVER), "Solver type (Default C16)")
            ("token", boost::program_options::value<string>(&token), "Solver Token")
            ("prompt-retry,R", boost::program_options::bool_switch(&prompt_retries),
             "Enable retry prompt if solver connection fails")
            ("timeout", boost::program_options::value<double>(&timeout)->default_value(60.0),
             "Solver Response Timeout (Default 1min)")
            ("n", boost::program_options::value<int>(&n)->default_value(20),
             "Number of samples per repetition (Default 20)")
            ("max-j", boost::program_options::value<double>(&J_scale)->default_value(-1.0), "Maximum energy scale for problem couplings."
                                                                                            " Problem is autoscaled if not set.")
            ("rand-gauge", boost::program_options::bool_switch(&rand_gauges), "Use random global gauge")
            ("reps,r", boost::program_options::value<int>(&reps)->default_value(1), "Number of repetitions (Default 1)")
            ;

    program_options.add(dwave_opts);
}

int advanced_schedule_program::parse_schedule(){
    regex re(R"((\w+)(\{\s*([.[:s:][:alnum:]]*)\})?)");
    //regex args_re(R"(\{\})");
    smatch match;
    regex_match(sched_str,  match, re);
    if(!match.ready()){
        cerr << "Failed to parse schedule specifier " << sched_str << endl;
        return 1;
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
        cerr << "Failed to parse schedule specifier " << sched_str << endl;
        return 1;
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
        }  else if(vm.count("reverse-init-file")){
            cerr << "Reverse anneal from file not implemented" << endl;
            return 1;
        }
        else {
            cerr << "reverse-init-cell required for reverse anneal schedule" << endl;
            return 1;
        }
        int num_args = sched_args_v.size();
        if(num_args == 3){
            t1 = stod(sched_args_v[0]);
            s1 = stod(sched_args_v[1]);
            tw = stod(sched_args_v[2]);
        }else {
            cerr << "Invalid arguments for reverse anneal schedule" << endl;
            return 1;
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
                cerr << "Invalid arguments for beta schedule" << endl;
                return 1;
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
            cerr << "pl3 requires 5 arguments" << endl;
            return 1;
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
                cerr << "Invalid arguments for beta schedule" << endl;
                return 1;
        }
    } else {
        cerr << "Unrecognized schedule type " << sched_type_str << endl;
        return 1;
    }

    return 0;
}

void advanced_schedule_program::generate_schedule(){

    switch(sched_type){
        case def:
//                cout << "Schedule specification " << sched_str << " not recognized" << endl;
            throw runtime_error("Schedule specification "+ sched_str + " not recognized");
        case lin:
            cout << "\t***Linear***" << endl;
            sched = dwave_cpp::linear_schedule(tf);
            break;
        case rev:
            cout << "\t***Reverse Anneal***";
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
            throw logic_error("Schedule type enum switch failed");
    }
}

gadget_program::gadget_program() : advanced_schedule_program(), gadget_opts("8-Qubit Gadget Options"){
    po::positional_options_description pd;

    gadget_opts.add_options()
            ("problem", po::value<string>(&cell_file), "Cell problem input file")
            ("cellsfile,", po::value<string>(&cell_locations_file), "Cell locations input file")
            ("cells,C", po::value<vector<int>>(&cell_locations_vec)->multitoken(), "List of cell locations")
            ("tgts", po::value<vector<string>>(&cell_tgts_strs)->multitoken(), "List of cell locations")
            ("write-cells", po::bool_switch(&write_cells_states), "Write individual cell states in output file "
                                                                  "(May result in large output file)")
            ;
    //pd.add("sched", 1).add("output", 1);

    program_options.add(gadget_opts);
    positional_options.add("sched", 1)
            .add("output", 1);

}

int gadget_program::check_options() {
    if(advanced_schedule_program::check_options())
        return 1;

    if(vm.count("tgts")){
        for(string& tgt : cell_tgts_strs){
            for(const char& c : tgt)
                if( !( (c == '0') || (c == '1') || (c == '3'))){
                    cerr << "Invalid target string: " << tgt << endl;
                    return 1;
                }
            if(tgt.size() != 8){
                cerr << "Invalid target string length: " <<  tgt << endl;
                return 1;
            }
        }
    }

    return 0;
}

dwave_cpp::CellProblem gadget_program::import_cell_problem() {
    dwave_cpp::CellProblem cell_problem;
    ifstream ifs(cell_file);
    ifs >> cell_problem;
    if(J_scale > 0){
        for(auto& p : cell_problem){
            p.value /= J_scale;
        }
    }
    return cell_problem;
}

void gadget_program::write_results(){
    string s(1, output_separator);
    ofstream f(output_file);
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
          << s  << double(p.second)/num_gadget_readouts << "\n";

    f   << '#' << s
        << "state" << s
        << "count" << s
        << "frac" << '\n';
    for( const auto& p : counts){
        f << s << p.first
          << s  << p.second
          << s  << double(p.second)/num_gadget_readouts << "\n";
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

void gadget_program::run(){
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

    num_gadget_readouts = 0;
    int num_tgts = cell_tgts_strs.size();
    tgt_count_vec.resize(num_tgts, 0);

    for(int r = 0; r < reps; ++r) {
        cout << r+1 << "/" << reps << "..." << endl;
        vector<vector<int8_t> > results_vec;
        results_vec = run_problem_vector(solver, problem_vector, params, timeout)[0];

        for (const auto &result : results_vec) {
            vector<int16_t> arr = decode_problem(solver, result);
            num_gadget_readouts += arr.size();
            for(int16_t st : arr){
                counts[st] += 1;
            }
            read_arr.push_back(std::move(arr));
        }
    }
    cout << endl;


    cout << "Counts: \n";
    for( const auto& p : counts){
        cout << "\t" << p.first
             << "\t\t"   << p.second
             << "\t\t"   << double(p.second)/num_gadget_readouts << "\n";
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
             << "\t\t" << double(tgt_count_vec[k]) / num_gadget_readouts << "\n";
    }

    /*
     * OUTPUT TO FILE
     *
     * */
    write_results();
}

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