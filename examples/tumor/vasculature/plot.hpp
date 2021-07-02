// SPDX-FileCopyrightText: 2015 - 2021 Marcin Łoś <marcin.los.91@gmail.com>
// SPDX-License-Identifier: MIT

#ifndef TUMOR_VASCULATURE_PLOT_HPP
#define TUMOR_VASCULATURE_PLOT_HPP

#include <ostream>

#include "defs.hpp"


namespace tumor::vasc {

void plot(std::ostream& os, const val_array& v);

void plot(const std::string& file, const val_array& v);

}

#endif // TUMOR_VASCULATURE_PLOT_HPP
