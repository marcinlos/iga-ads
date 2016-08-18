#ifndef PROBLEMS_TUMOR_VASCULATURE_CONFIG_HPP_
#define PROBLEMS_TUMOR_VASCULATURE_CONFIG_HPP_

namespace tumor {
namespace vasc {

struct config {
    double init_stability = 0.5;
    double degeneration = 0.05;

//    double t_ec_sprout = 0.5;
//    double t_ec_migr = 100;

    double t_ec_sprout = 10;
    double t_ec_migr = 2;


    // my values
    double segment_length = 0.03;
    double t_ec_collapse = 10;
    double c_min = 0.2;
};


}
}

#endif /* PROBLEMS_TUMOR_VASCULATURE_CONFIG_HPP_ */
