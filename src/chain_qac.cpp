//
// Created by Humberto Munoz Bauza on 2019-09-03.
//
#include <dwave_cpp/dwave_cpp.h>
#include "dw_prog.h"

namespace po = boost::program_options;

struct chain_qac_prog : public gadget_program {
    double qac_penalty=0.0;
    bool no_qac=false;
    po::options_description chain_qac_opts;

    chain_qac_prog();

    dwave_cpp::Problem encode_problem( dwave_cpp::Solver& solver,
                                       dwave_cpp::CellProblem &cell_problem) override;
    vector<int16_t> decode_problem(
            dwave_cpp::Solver& solver, const vector<int8_t>& readout) override;
};

chain_qac_prog::chain_qac_prog() : gadget_program(),chain_qac_opts("Chain QAC Options"){
    chain_qac_opts.add_options()("qac-penalty", po::value<double>(&qac_penalty)->default_value(1.0),
            "Penalty strength for QAC (default 1.0)")
            ("no-qac", po::bool_switch(&no_qac), 
                    "Encode copies of the linear chain only and do not use QAC encoding")
            ;
                                
    program_options.add(chain_qac_opts);

}
dwave_cpp::Problem chain_qac_prog::encode_problem(dwave_cpp::Solver &solver, dwave_cpp::CellProblem &cell_problem) {
    
    return dwave_cpp::GenerateQACChainProblem(solver, cell_problem, (no_qac ? -1.0 : qac_penalty));
}

vector<int16_t> chain_qac_prog::decode_problem(dwave_cpp::Solver &solver, const vector<int8_t> &readout) {
    if(no_qac){
        return dwave_cpp::ReadVerticalChainProblem(solver, readout, 8);
    } else{
        auto decode = dwave_cpp::DecodeQACProblem(solver, readout);
        return dwave_cpp::ReadBipartEmbedChains(solver, decode,  8, true);
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