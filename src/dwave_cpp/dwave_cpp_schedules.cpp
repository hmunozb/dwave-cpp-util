//
// Created by Humberto Munoz Bauza on 2019-05-19.
//

#include <stdexcept>
#include "dwave_cpp/schedules/beta.h"
#include <boost/math/special_functions/beta.hpp>
#include <boost/format.hpp>
namespace dwave_cpp{

#include <dwave_sapi.h>

    std::vector<AnnealSchedulePoint> pl3_schedule(std::pair<double, double> pl1,
                                                  std::pair<double, double> pl2, double tf){
        std::vector<AnnealSchedulePoint> sched_pnts;
        sched_pnts.push_back({0.0, 0.0});
        sched_pnts.push_back({pl1.first, pl1.second});
        sched_pnts.push_back({pl2.first, pl2.second});
        sched_pnts.push_back({tf, 1.0});

        return sched_pnts;
    }

    std::vector<AnnealSchedulePoint> ps3_schedule(std::pair<double, double> pl1,
                                                  std::pair<double, double> pl2, double tf){
        std::vector<AnnealSchedulePoint> sched_pnts;
        sched_pnts.push_back({0.0, 0.0});
        double t1 = pl1.first * pl1.second;
        double t2 = t1 + pl2.first * (pl2.second - pl1.second);
        double t3 = t2 + tf*(1.0 - pl2.second);

        sched_pnts.push_back({t1, pl1.second});
        sched_pnts.push_back({t2, pl2.second});
        sched_pnts.push_back({t3, 1.0});

        return sched_pnts;
    }

    std::vector<AnnealSchedulePoint> _beta_schedule(double ctrl_s, int start_pnts, int end_pnts, unsigned a, unsigned b){
        double start_div = ctrl_s / (start_pnts + 1);
        int total_pnts = start_pnts + end_pnts + 3;
        std::vector<AnnealSchedulePoint> ctrl_pnts(total_pnts);

        ctrl_pnts[0].relative_current = 0.0;
        ctrl_pnts[0].time = 0.0;

        for(int i = 0; i < start_pnts ; ++i){
            ctrl_pnts[1 + i].relative_current = start_div * (i + 1);
        }

        ctrl_pnts[1 + start_pnts].relative_current = ctrl_s;

        double d = (1.0 - ctrl_s) / 2.0;
        for(int i = 0; i < end_pnts; ++i){
            ctrl_pnts[2 + start_pnts + i].relative_current = 1 - d;
            d = d / 2.0;
        }

        ctrl_pnts[total_pnts - 1].relative_current = 1.0;
        ctrl_pnts[total_pnts - 1].time = 1.0;

        //calculate the inverse ibeta for the intermediate control points
        for(int i = 1; i < total_pnts - 1; ++i){
            ctrl_pnts[i].time = boost::math::ibeta_inv(a, b, ctrl_pnts[i].relative_current);
        }

        //post condition: assert monotonicity
        for(int i = 0; i < total_pnts - 1; ++i){
            assert( (ctrl_pnts[i+1].time >= ctrl_pnts[i].time) );
            assert( (ctrl_pnts[i+1].relative_current >= ctrl_pnts[i].relative_current) );
        }

        return ctrl_pnts;
    }

    std::vector<AnnealSchedulePoint> linear_schedule(double anneal_time){
        std::vector<AnnealSchedulePoint> sched;
        sched.push_back({0,0});
        sched.push_back({anneal_time, 1});

        return sched;
    }

    std::vector<AnnealSchedulePoint> reverse_anneal_schedule(double reverse_s, double t1, double t_wait){
        std::vector<AnnealSchedulePoint> sched;
        sched.push_back({0.0, 1.0});
        sched.push_back({t1, reverse_s});
        sched.push_back({t1 + t_wait, reverse_s});
        sched.push_back({2.0*t1 + t_wait, 1.0});

        return sched;
    }

    std::vector<AnnealSchedulePoint> pl2_beta_quenched_schedule(
            std::pair<double, double> pl1,
            std::pair<double, double> pl2,
            double beta_segment_time,
            unsigned a, unsigned b,
            double quench_s,
            double ctrl_r)
    {
        int BETA_START_PNTS = 2;
        int BETA_END_PNTS = (quench_s < 0 ? 5 : 4);
        //10 points for the beta schedule
        std::vector<AnnealSchedulePoint> sched_pnts;
        sched_pnts.push_back({0.0, 0.0});
        sched_pnts.push_back({pl1.first, pl1.second});
        double beta_init_time = pl2.first;
        double beta_init_current = pl2.second;
        double rem_current = (quench_s < 0 ? 1.0 - beta_init_current : quench_s - beta_init_current);

        std::vector<AnnealSchedulePoint> ctrl_pnts = _beta_schedule(ctrl_r, BETA_START_PNTS, BETA_END_PNTS, a, b);

        for(AnnealSchedulePoint p : ctrl_pnts){
            sched_pnts.push_back({ beta_init_time + beta_segment_time * p.time,
                                   beta_init_current + rem_current * p.relative_current});
        }
        AnnealSchedulePoint last_p = sched_pnts.back();

        if(quench_s >= 0)  //add a quench if it is enabled
            sched_pnts.push_back({last_p.time + 1.0, 1.0});

        return sched_pnts;
    }

    std::vector<AnnealSchedulePoint> beta_quenched_schedule( double anneal_time,
            double quench_s, double ctrl_r, unsigned a, unsigned b) {
        if(!(ctrl_r <= 1 && ctrl_r >= 0)){
            throw std::invalid_argument(
                    boost::str(boost::format("Invalid beta schedule s_ctrl %1%")%ctrl_r) );
        }
        if( !(quench_s <= 1 && quench_s >= 0)){
            throw std::invalid_argument(
                    boost::str(boost::format("Invalid beta schedule quench s %1%")%quench_s) );
        }
        int START_PNTS = 4;
        int END_PNTS = 4;
        int TOTAL_PNTS = 12;

        //The first 11 points are the beta schedule points
        std::vector<AnnealSchedulePoint> ctrl_pnts = _beta_schedule(ctrl_r, START_PNTS, END_PNTS, a, b);

        for( AnnealSchedulePoint& p : ctrl_pnts){
            p.time *= anneal_time;
            p.relative_current *= quench_s;
        }
        //The final point is the quenche
        ctrl_pnts.push_back({anneal_time + 1.0, 1.0});

        //post condition: assert monotonicity
        for(int i = 0; i < TOTAL_PNTS - 1; ++i){
            assert( (ctrl_pnts[i+1].time >= ctrl_pnts[i].time) );
            assert( (ctrl_pnts[i+1].relative_current >= ctrl_pnts[i].relative_current) );
        }

        return ctrl_pnts;
    }

    std::vector<AnnealSchedulePoint> beta_schedule( double anneal_time,
            double ctrl_r, unsigned a, unsigned b){
        if(!(ctrl_r <= 1 && ctrl_r >= 0)){
            throw std::invalid_argument(
                    boost::str(boost::format("Invalid beta schedule s_ctrl %1%")%ctrl_r) );
        }

        int START_PNTS = 4;
        int END_PNTS = 5;

        std::vector<AnnealSchedulePoint> ctrl_pnts = _beta_schedule(ctrl_r, START_PNTS, END_PNTS, a, b);

        for( AnnealSchedulePoint& p : ctrl_pnts){
            p.time *= anneal_time;
        }

        //post condition: assert monotonicity
        for(int i = 0; i < 11; ++i){
            assert( (ctrl_pnts[i+1].time >= ctrl_pnts[i].time) );
            assert( (ctrl_pnts[i+1].relative_current >= ctrl_pnts[i].relative_current) );
        }

        return ctrl_pnts;

    }
}