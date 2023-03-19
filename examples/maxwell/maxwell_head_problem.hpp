// SPDX-FileCopyrightText: 2015 - 2023 Marcin Łoś <marcin.los.91@gmail.com>
// SPDX-License-Identifier: MIT

#ifndef MAXWELL_MAXWELL_HEAD_PROBLEM_HPP
#define MAXWELL_MAXWELL_HEAD_PROBLEM_HPP

#include <cmath>
#include <string_view>

#include "head_data.hpp"
#include "problems.hpp"

class maxwell_head_problem : public maxwell_problem<maxwell_head_problem> {
private:
    maxwell_manufactured1 init_state_{1, 1};
    head_data data_;

public:
    using point_type = head_data::point_type;

    explicit maxwell_head_problem(std::string_view data_path)
    : data_{read_density_data(data_path)} { }

    auto init_E1() const { return init_state_.init_E1(); }
    auto init_E2() const { return init_state_.init_E2(); }
    auto init_E3() const { return init_state_.init_E3(); }

    auto init_H1() const { return init_state_.init_H1(); }
    auto init_H2() const { return init_state_.init_H2(); }
    auto init_H3() const { return init_state_.init_H3(); }

    auto eps(point_type x) const -> double { return data_.eps(x); }

    auto mu(point_type x) const -> double { return data_.mu(x); }

    constexpr static value_type unknown{NAN, NAN, NAN, NAN};

    auto E1(point_type, double) const -> value_type { return unknown; }
    auto E2(point_type, double) const -> value_type { return unknown; }
    auto E3(point_type, double) const -> value_type { return unknown; }

    auto H1(point_type, double) const -> value_type { return unknown; }
    auto H2(point_type, double) const -> value_type { return unknown; }
    auto H3(point_type, double) const -> value_type { return unknown; }
};

#endif  // MAXWELL_MAXWELL_HEAD_PROBLEM_HPP
