// SPDX-FileCopyrightText: 2015 - 2023 Marcin Łoś <marcin.los.91@gmail.com>
// SPDX-License-Identifier: MIT

#ifndef TUMOR_VASCULATURE_CONFIG_HPP
#define TUMOR_VASCULATURE_CONFIG_HPP

namespace tumor::vasc {

struct config {
    double init_stability = 0.5;
    double degeneration = 0.05;

    double t_ec_sprout = 0.3;

    double segment_length = 40.0 / 3000.0;  // 40um

    double r_sprout = 5.0;
    double r_max = 25.0;

    double t_ec_switch = 24.0;
    double c_switch = 0.003;  // ???????????
    double dilatation = 0.4;  // 4um/h

    double t_ec_collapse = 10;  // dt/t_ec_collapse should be 0.01 for veins, 0.004 for capillaries
    double c_min = 0.003;

    // to remove
    double t_ec_migr = 2;
};

}  // namespace tumor::vasc

#endif  // TUMOR_VASCULATURE_CONFIG_HPP
