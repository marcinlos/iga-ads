// SPDX-FileCopyrightText: 2015 - 2021 Marcin Łoś <marcin.los.91@gmail.com>
// SPDX-License-Identifier: MIT

#ifndef PROBLEMS_TUMOR_PARAMS_HPP_
#define PROBLEMS_TUMOR_PARAMS_HPP_

#include "skin.hpp"


namespace tumor {

struct params {
    double c_b_norm = 1; // normal concentration of tumor cells
    double c_b_max = 2; // maximal concentration of tumor cells
    double tau_b = 0.5; // ?

    double o_prol_TC = 10;   // oxygen proliferation threshold
    double o_death_TC = 2; // oxygen survival threshold
    double t_prol_TC = 10;
    double t_death_TC = 100;
    double P_b = 0.001; // stimulated mitosis rate
    double r_b = 0.3;  // chemoattractant sensitivity, 1-3-5 x 10^-4

    double beta_m = 1 / 16.0;      // ???? some ECM decay coefficient
    double gamma_a = 0.5 / 15.6;   // production rate of attractants
    double chi_aA = 0.01 / 15.6;   // diffusion of degraded ECM
    double gamma_oA = 0.01 / 15.6; // decay of degraded ECM

    double diff_c = 0.0000555 ; // TAF diffusion rate
    double cons_c = 0.01; // TAF consumption

    double o_max = 60.0;
    double alpha_0 = 0.0000555;
    double gamma_T = 0.01;
    double alpha_1 = 0.4;

    skin_model skin;

    double init_M = 0.015;
};


}

#endif /* PROBLEMS_TUMOR_PARAMS_HPP_ */
