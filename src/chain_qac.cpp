//
// Created by Humberto Munoz Bauza on 2019-09-03.
//
#include <dwave_cpp/dwave_cpp.h>
#include "dw_prog.h"
//#include "dwave_cpp_core.h"

namespace po = boost::program_options;

struct chain_qac_prog : public gadget_program {
    double qac_penalty=0.0;
    bool no_qac=false;
    string qac_strategy_str;
    enum qac_strategy{
        complete,
        classical,
        no_penalty,
        penalty_only,
        unprotected
    };
    qac_strategy strategy = complete;
    po::options_description chain_qac_opts;

    chain_qac_prog();
    int check_options() override;

    dwave_cpp::Problem encode_problem( dwave_cpp::Solver& solver,
                                       dwave_cpp::ProblemAdj &cell_problem) override;
    vector<int16_t> decode_problem(
            dwave_cpp::Solver& solver, const vector<int8_t>& readout) override;
};

chain_qac_prog::chain_qac_prog() : gadget_program(),chain_qac_opts("Chain QAC Options"){
    chain_qac_opts.add_options()("qac-penalty", po::value<double>(&qac_penalty)->default_value(1.0),
            "Penalty strength for QAC (default 1.0)")
            ("mode", po::value<string>(&qac_strategy_str),
                    R"(
Indicates an alternative QAC method is to be used. Available:
         unprotected   -  Embed uncorrelated chain over 4n qubits, using no penalty or decoding.
         classical     -  Embed uncorrelated chain over 4n qubits, using no penalty, and decode by lowest energy found
         no-penalty    -  Embed uncorrelated chain over 3n qubits, using no penalty, and decode by majority vote
         penalty-only  -  Embed the chains with QAC penalty, but read all chains directly (no decoding)

Performs ordinary QAC with the penalty parameter if not specified.
)")
            ("no-qac", po::bool_switch(&no_qac),
                    "Encode copies of the linear chain only and do not use QAC encoding [Deprecated]")
            ;

    program_options.add(chain_qac_opts);

}

int chain_qac_prog::check_options() {
    if(gadget_program::check_options())
        return 1;

    if(!vm.count("mode"))
        return 0;

    if(qac_strategy_str == "unprotected"){
        strategy = unprotected;
    } else if (qac_strategy_str == "classical"){
        strategy = classical;
    } else if (qac_strategy_str == "no-penalty"){
        strategy = no_penalty;
    } else if (qac_strategy_str == "penalty-only"){
        strategy = penalty_only;
    } else {
        cerr << "Invalid QAC mode \"" << qac_strategy_str << "\" " << endl;
        return 1;
    }

    cout << "QAC Mode Specified: " << qac_strategy_str << endl;

    return 0;
}


dwave_cpp::Problem chain_qac_prog::encode_problem(dwave_cpp::Solver &solver, dwave_cpp::ProblemAdj &cell_problem) {
    double penalty = 0.0;
    if(no_qac){
        penalty = -1.0;
    } else {
        switch(strategy){
            case complete:
            case penalty_only:
                penalty = qac_penalty;
                break;
            case no_penalty:
                penalty = 0.0;
                break;
            case classical:
            case unprotected:
                penalty = -1.0;
                break;
        }
    }
    return dwave_cpp::GenerateQACChainProblem(solver, cell_problem, penalty);
}

vector<int16_t> chain_qac_prog::decode_problem(dwave_cpp::Solver &solver, const vector<int8_t> &readout) {
    if(no_qac){
        return dwave_cpp::ReadVerticalChainProblem(solver, readout, 8);
    } else {
        switch(strategy){
            case complete:
            case no_penalty:
            {
                auto decode = dwave_cpp::DecodeQACProblem(solver, readout);
                return dwave_cpp::ReadBipartEmbedChains(solver, decode,  8, true);
            }
            case unprotected:
                return dwave_cpp::ReadVerticalChainProblem(solver, readout, 8);
            case classical:
            {
                auto vert_chains = dwave_cpp::ReadVerticalChains(solver, readout, 8, 2, 4);
                return dwave_cpp::ClassicalDecode(vert_chains, gadget_problem, true);
            }
            case penalty_only:
            {
                auto vert_chains = dwave_cpp::ReadVerticalChains(solver, readout, 8, 2, 3);
                return dwave_cpp::UnprotectedDecode(vert_chains, true);
            }
        }
//        auto decode = dwave_cpp::DecodeQACProblem(solver, readout);
//        return dwave_cpp::ReadBipartEmbedChains(solver, decode,  8, true);
    }
}

int main(int argc, const char** argv){
    chain_qac_prog Prog;
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