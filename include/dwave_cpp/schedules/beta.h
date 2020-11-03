//
// Created by Humberto Munoz Bauza on 2019-05-17.
//

#ifndef DW_EXAMPLES_BETA_H
#define DW_EXAMPLES_BETA_H

#include <vector>
#include "../core/parameters.h"

namespace dwave_cpp{

    std::vector<AnnealSchedulePoint> linear_schedule(double anneal_time);
    std::vector<AnnealSchedulePoint> reverse_anneal_schedule(double reverse_s, double t1, double t_wait);
    std::vector<AnnealSchedulePoint> reverse_anneal_schedule(double reverse_s, double t1, double t2, double t_wait);
    std::vector<AnnealSchedulePoint> pl3_schedule(std::pair<double, double> pl1,
                                                  std::pair<double, double> pl2, double tf);
    std::vector<AnnealSchedulePoint> ps3_schedule(std::pair<double, double> pl1,
                                                  std::pair<double, double> pl2, double tf);

    std::vector<AnnealSchedulePoint> pl2_beta_quenched_schedule(
            std::pair<double, double> pl1,
            std::pair<double, double> pl2,
            double beta_segment_time,
            unsigned a, unsigned b,
            double quench_s,
            double ctrl_r);
    std::vector<AnnealSchedulePoint> beta_quenched_schedule( double anneal_time,
                                                             double quench_s, double ctrl_r, unsigned a, unsigned b);
    std::vector<AnnealSchedulePoint> beta_schedule(double anneal_time, double ctrl_r, unsigned a, unsigned b);
}

#endif //DW_EXAMPLES_BETA_H
