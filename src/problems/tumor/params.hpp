#ifndef PROBLEMS_TUMOR_PARAMS_HPP_
#define PROBLEMS_TUMOR_PARAMS_HPP_

#include "problems/tumor/skin.hpp"

namespace ads {
namespace tumor {


struct params {
    double c_b_norm = 1; // normal concentration of tumor cells
    double c_b_max = 2; // maximal concentration of tumor cells
    double tau_b = 0.5; // ?

    double o_prol_TC = 0.1;   // oxygen proliferation threshold
    double o_death_TC = 0.01; // oxygen survival threshold
    double t_prol_TC = 10;
    double t_death_TC = 100;
    double P_b = 0.001; // stimulated mitosis rate
    double r_b = 1e-4;  // chemoattractant sensitivity, 1-3-5 x 10^-4

    double beta_m = 1;      // ???? some ECM decay coefficient
    double gamma_a = 0.5;   // production rate of attractants
    double chi_aA = 0.01;   // diffusion of degraded ECM
    double gamma_oA = 0.01; // decay of degraded ECM

    double rho_n = 0.28;  // haptotactic cell migration
    double D_n = 0.0003;  // diffusion of endothelial cells
    double chi_n = 0.38;  // chemotactic cell migration
    double delta_n = 0.6; // chemotactic constant

    double beta_f = 0.05; // production rate of fibronectin
    double gamma_f = 0.1; // degradation rate of fibronectin

    double alpha_m = 0.000001; // MDA generation rate
    double epsilon_m = 0.01;   // MDA diffusion coefficient
    double upsilon_m = 3;      // MDA degradation rage

    double diff_c = 0.01; // TAF diffusion rate
    double cons_c = 0.3;  // TAF consumption

    skin_model skin;

    double init_M = 0.015;
};


}
}


#endif /* PROBLEMS_TUMOR_PARAMS_HPP_ */
