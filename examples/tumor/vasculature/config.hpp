// SPDX-FileCopyrightText: 2015 - 2021 Marcin Łoś <marcin.los.91@gmail.com>
// SPDX-License-Identifier: MIT

#ifndef PROBLEMS_TUMOR_VASCULATURE_CONFIG_HPP_
#define PROBLEMS_TUMOR_VASCULATURE_CONFIG_HPP_


namespace tumor::vasc {

struct config {
    double init_stability = 0.5;
    double degeneration = 0.05;

    double t_ec_sprout = 0.3;

    double segment_length = 40.0 / 3000.0; // 40um

    double r_sprout = 5.0;
    double r_max = 25.0;

    double t_ec_switch = 24.0;
    double c_switch = 0.003; // ???????????
    double dilatation = 0.4; // 4um/h

    double t_ec_collapse = 10; // tak, zeby dt/t_ec_collapse wyszlo 0.01 dla zyl/tetnic, 0.004 dla kapilar
    double c_min = 0.003;

    // to remove
    double t_ec_migr = 2;
};


}

#endif /* PROBLEMS_TUMOR_VASCULATURE_CONFIG_HPP_ */
