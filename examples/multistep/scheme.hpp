// SPDX-FileCopyrightText: 2015 - 2021 Marcin Łoś <marcin.los.91@gmail.com>
// SPDX-License-Identifier: MIT

#ifndef MULTISTEP_SCHEME_HPP
#define MULTISTEP_SCHEME_HPP

#include <iostream>
#include <iterator>
#include <sstream>
#include <vector>

namespace ads {

struct scheme {
    int s;
    std::vector<double> as;
    std::vector<double> bs;
};

std::ostream& operator<<(std::ostream& out, const scheme& s);

ads::scheme get_scheme(const std::string& name);

scheme parse_scheme(const std::string& text);

}  // namespace ads

#endif  // MULTISTEP_SCHEME_HPP
