//
// Created by Humberto Munoz Bauza on 2019-09-03.
//

#ifndef DW_EXAMPLES_DW_PROG_H
#define DW_EXAMPLES_DW_PROG_H

#include <dwave_cpp/dwave_cpp.h>
#include <dwave_cpp/problems/cell_gadgets.h>
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

namespace po = boost::program_options;

bool cmp_tgts( short s, const string& tgt_str);

enum SchedType{
    def=0,
    lin,
    rev,
    beta,
    pl2b,
    pl3,
    ps3
};

// Base class for all programs
// Only contains options for basic program functions (--help, --verbose, etc...)
struct program_base{
    bool verbose=false;
    string output_file;
    char output_separator = ' ';

    boost::program_options::options_description program_options;
    boost::program_options::positional_options_description positional_options;
    boost::program_options::variables_map vm;

    program_base(): program_options("Program Options"){
        program_options.add_options()
                ("help", "Print help message")
                ("verbose,v", boost::program_options::bool_switch(&verbose), "Enable verbose printing to console during runs")
                ("output,o", boost::program_options::value<string>(&output_file), "Output file [required]")
                ;
    }

    // Parses all command line options. Returns 1 if --help was passed
    // indicating that the derived program should terminate
    int parse_all_options(int argc, const char **argv){
        try {
            boost::program_options::store(
                    boost::program_options::command_line_parser(argc, argv)
                            .options(program_options)
                            .positional(positional_options).run(),
                    vm);
            boost::program_options::notify(vm);
        } catch(std::exception& e) {
            cerr << "Error: " << e.what() << endl;
            return -1;
        }

        if (vm.count("help")) {
            cout << program_options << "\n";
            return 1;
        }

        return 0;
    }

    // Checks that the passed options are valid
    // If this is overriden, each derived method should call its parent method before performing its own checks
    virtual int check_options(){
        if(!vm.count("output")){
            cout << "An output file (--output) is required." << endl;
            return 1;
        }
        boost::filesystem::path output_path(output_file);
        auto output_dir = output_path.parent_path();
        if (not(output_dir.empty()) and not(boost::filesystem::exists(output_dir))) {
            cout << "Could not find output directory " << output_dir;
            return 1;
        }
        auto output_ext = output_path.extension();
        if (output_ext.compare(".tsv") == 0 or output_ext.compare(".txt") == 0) {
            output_separator = '\t';
        } else if (output_ext.compare(".csv")) {
            output_separator = ',';
        } else if (output_ext.empty()){
            cout << "Note: no file extension. Outputting tap-separated text to the destination." << endl;
            output_separator = '\t';
        } else {
            cout << "Note: file extension " << output_ext
                 << "not recognized. Outputting tap-separated data to the destination." << endl;
            output_separator = '\t';
        }

        return 0;
    }
};

struct generic_dwave_program: public program_base{
    string url;
    string token;
    bool prompt_retries=false;
    string solver_str;
    double timeout = 0.0;

    bool rand_gauges=false;
    double J_scale = -1.0;
    int n=0;
    int reps=0;

    boost::program_options::options_description dwave_opts;

    generic_dwave_program() ;

    int check_options() override {
        if(program_base::check_options()){
            return 1;
        }

        if (!vm.count("url")) {
            cout << "Must specify solver url with --url" << endl;
            return 1;
        }
        if (!vm.count("token")) {
            cout << "Must specify solver API token with --token";
            return 1;
        }
        if (!vm.count("output")) {
            cout << "Must specify an output file";
            return 1;
        }

        return 0;
    }

    void set_parameters(dwave_cpp::QuantumSolverParameters& params){

        params->answer_mode = dwave_cpp::SAPI_ANSWER_MODE_RAW;
        params->num_spin_reversal_transforms = int(rand_gauges);
        params->num_reads = n;
        params->auto_scale = int(J_scale < 0);

    }


};

struct basic_schedule_program: public generic_dwave_program {
    double tf = 0.0; //Default anneal time
    boost::program_options::options_description sched_opts;

    basic_schedule_program(): generic_dwave_program(), sched_opts("Schedule Options"){
        sched_opts.add_options()
                ("tf", boost::program_options::value<double>(&tf)->default_value(50.0), "Anneal time (Default 50.0 us)")
                ;
        program_options.add(sched_opts);
    }

    int check_options() override{
        if(generic_dwave_program::check_options()){
            return 1;
        }


        return 0;
    }


};

struct advanced_schedule_program : public basic_schedule_program{
    boost::program_options::options_description adv_sched_opts;
    vector<dwave_cpp::AnnealSchedulePoint> sched;

    string sched_str;
    SchedType sched_type=def;
    //CellProb cell_prob = CellProb::c0989v1;
    string reverse_init_cell;
    vector<int> reverse_init_cell_vec;

    int a = 0;
    int b = 0;
    //double tf = 0.0; //Default anneal time
    double tw = 0.0;
    double t1 = 0.0;
    double s1 = 0.0;
    double sq = -1.0;
    double sc = 0.9;
    pair<double, double> pl1;
    pair<double, double> pl2;

    advanced_schedule_program(): basic_schedule_program(), adv_sched_opts("Advanced Schedule Options"){
        adv_sched_opts.add_options()
            ("sched", boost::program_options::value<string>(&sched_str), "Schedule Type")
            ;
        program_options.add(adv_sched_opts);
    }


    int parse_schedule();
    int check_options() override {
        if(basic_schedule_program::check_options()){
            return 1;
        }
        if(vm.count("tf") && !vm.count("sched")){
            sched_type = lin;
            return 0;
        }
        if(!vm.count("tf") && !vm.count("sched")){
            cout << "Error: Must specify an anneal time --tf or an advanced schedule --sched." << endl;
            return 1;
        }
        if(vm.count("tf") && vm.count("sched")){
            cout << "Note: option --sched has precedence over --tf for advanced schedules." << endl;
        }

        return parse_schedule();
    }

    void generate_schedule();

    void set_schedule_parameters(dwave_cpp::QuantumSolverParameters& params){
        if(sched_type == rev){
            params.set_reverse_anneal(reverse_init_cell_vec);
        }
        params.set_anneal_schedule(sched);
    }
};


struct gadget_program : public advanced_schedule_program{
// Parametric data
    string cell_file;
    string cell_locations_file;
    bool write_cells_states=false;
    vector<int> cell_locations_vec;
    vector<string> cell_tgts_strs;
    po::options_description gadget_opts;
// Result data
    vector< vector<int16_t> > read_arr;
    std::map<int16_t, int> counts;
    std::vector<int> tgt_count_vec;
    std::map<string, int> tgt_counts;
    int num_gadget_readouts=0;

    gadget_program();
    int check_options() override;
    dwave_cpp::CellProblem import_cell_problem();
    virtual dwave_cpp::Problem encode_problem(
            dwave_cpp::Solver& solver,
            dwave_cpp::CellProblem &cell_problem ) = 0;

    virtual vector<int16_t> decode_problem(
            dwave_cpp::Solver& solver,
            const vector<int8_t>& readout) = 0;

    void write_results();
    void run();
};

#endif //DW_EXAMPLES_DW_PROG_H
