#ifndef PROBLEMS_TUMOR_VASCULATURE_CONFIG_HPP_
#define PROBLEMS_TUMOR_VASCULATURE_CONFIG_HPP_

namespace tumor {
namespace vasc {

struct config {
    double init_stability = 0.5;
    double degeneration = 0.05;

//    double t_ec_sprout = 0.5;
//    double t_ec_migr = 100;

    // double t_ec_sprout = 10;
    // double t_ec_sprout = 2;
    double t_ec_sprout = 0.3;


    double t_ec_migr = 2;


    // my values
    // double segment_length = 0.03;
    double segment_length = 0.03; // zmienic na 40 mikrometrow

    double vessel_diameter = 0; // 5um!!

    double t_ec_collapse = 10; // tak, zeby dt/t_ec_collapse wyszlo 0.01 dla zyl/tetnic, 0.004 dla kapilar
    double c_min = 0.003;
};


}
}

#endif /* PROBLEMS_TUMOR_VASCULATURE_CONFIG_HPP_ */
