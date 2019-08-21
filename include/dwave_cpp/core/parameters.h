//
// Created by Humberto Munoz Bauza on 2019-05-15.
//

#ifndef DW_EXAMPLES_PARAMETERS_H
#define DW_EXAMPLES_PARAMETERS_H
#include <memory>
#include <vector>
#include <tuple>
namespace  dwave_cpp{
    struct sapi_SolverParameters;
    struct sapi_QuantumSolverParameters;
    struct sapi_DoubleArray;
    struct sapi_AnnealSchedule;
    struct sapi_AnnealSchedulePoint;
    typedef sapi_AnnealSchedulePoint AnnealSchedulePoint ;
    typedef std::vector<AnnealSchedulePoint> AnnealSchedule;
    struct sapi_ReverseAnneal;
    struct sapi_DoubleArray;

    std::ostream& operator<<(std::ostream& ostream1, const AnnealSchedule& sched);

    class SolverParameters{
    public:
        SolverParameters() = default;
        //sapi_SolverParameters& get_params();
        //~SolverParameters();
    private:
        sapi_SolverParameters * sol_params_ptr;
    };

    class QuantumSolverParameters{
    public:
        QuantumSolverParameters();
        void set_anneal_schedule(const std::vector<sapi_AnnealSchedulePoint >& schedule);
        void set_reverse_anneal(std::vector<int>& init_state, bool reinitialize=true);
        sapi_QuantumSolverParameters* get();
        sapi_QuantumSolverParameters* operator->();
        //~QuantumSolverParameters();
    private:
        std::unique_ptr<sapi_QuantumSolverParameters> q_sol_params_ptr;
        std::unique_ptr<sapi_DoubleArray> anneal_offsets_ptr;
        std::unique_ptr<sapi_AnnealSchedule> anneal_schedule_ptr;
        std::unique_ptr<AnnealSchedulePoint[]> anneal_schedule_pnts_ptr;
        std::unique_ptr<sapi_ReverseAnneal> reverse_anneal_ptr;
        std::unique_ptr<sapi_DoubleArray> flux_biases_ptr;
    };
}
#endif //DW_EXAMPLES_PARAMETERS_H
