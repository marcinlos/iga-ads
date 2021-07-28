// SPDX-FileCopyrightText: 2015 - 2021 Marcin Łoś <marcin.los.91@gmail.com>
// SPDX-License-Identifier: MIT

#include "scheme.hpp"

#include <map>

namespace ads {

const std::map<std::string, ads::scheme> schemes{
    {"AM-0", {0, {-1}, {1}}},
    {"AM-1", {1, {-1}, {0.5, 0.5}}},
    {"AM-2", {2, {-1, 0}, {5. / 12, 2. / 3, -1. / 12}}},
    {"AM-3", {3, {-1, 0, 0}, {9. / 24, 19. / 24, -5. / 24, 1. / 24}}},
    {"AM-4", {4, {-1, 0, 0, 0}, {251. / 270, 646. / 720, -264. / 720, 106. / 720, -19. / 720}}},

    {"BDF-1", {1, {-1}, {1, 0}}},
    {"BDF-2", {2, {-4. / 3, 1. / 3}, {2. / 3, 0, 0}}},
    {"BDF-3", {3, {-18. / 11, 9. / 11, -2. / 11}, {6. / 11, 0, 0, 0}}},
    {"BDF-4", {4, {-48. / 25, 36. / 25, -16. / 25, 3. / 25}, {12. / 25, 0, 0, 0, 0}}},
    {"BDF-5",
     {5,
      {-300. / 137, 300. / 137, -200. / 137, 75. / 137, -12. / 137},
      {60. / 137, 0, 0, 0, 0, 0}}},
};

std::ostream& operator<<(std::ostream& out, const scheme& s) {
    out << "s = " << s.s << ", as = [ ";
    std::copy(begin(s.as), end(s.as), std::ostream_iterator<double>{out, " "});
    out << "], bs = [ ";
    std::copy(begin(s.bs), end(s.bs), std::ostream_iterator<double>{out, " "});
    out << "]";
    return out;
}

ads::scheme get_scheme(const std::string& name) {
    auto it = schemes.find(name);
    if (it != schemes.end()) {
        return it->second;
    } else {
        return parse_scheme(name);
    }
}

scheme parse_scheme(const std::string& text) {
    std::istringstream input{text};
    std::string buf;

    ads::scheme scm;

    input >> scm.s;
    if (!input)
        throw std::invalid_argument{"Invalid scheme (cannot read s value)"};

    auto eat_pipe = [&]() { input >> buf; };

    auto read_numbers = [&](auto out) {
        while (input) {
            double a;
            input >> a;
            if (input)
                *out++ = a;
        }
        input.clear();
    };

    eat_pipe();
    read_numbers(std::back_inserter(scm.as));
    eat_pipe();
    read_numbers(std::back_inserter(scm.bs));

    return scm;
}

}  // namespace ads
