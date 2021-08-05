// SPDX-FileCopyrightText: 2015 - 2021 Marcin Łoś <marcin.los.91@gmail.com>
// SPDX-License-Identifier: MIT

#ifndef MAXWELL_MAXWELL_HEAD_PROBLEM_HPP
#define MAXWELL_MAXWELL_HEAD_PROBLEM_HPP

#include <string>

#include "head_data.hpp"
#include "problems.hpp"

class maxwell_head_problem {
private:
    maxwell_manufactured1 init_state_{1, 1};
    head_data data_;

public:
    using point_type = head_data::point_type;

    explicit maxwell_head_problem(const std::string& data_path)
    : data_{read_density_data(data_path)} { }

    auto init_E1() const { return init_state_.E1_val_at(0); }
    auto init_E2() const { return init_state_.E2_val_at(0); }
    auto init_E3() const { return init_state_.E3_val_at(0); }

    auto init_H1() const { return init_state_.H1_val_at(0); }
    auto init_H2() const { return init_state_.H2_val_at(0); }
    auto init_H3() const { return init_state_.H3_val_at(0); }

    auto eps(point_type x) const -> double { return data_.eps(x); }

    auto mu(point_type x) const -> double { return data_.mu(x); }
};

#endif  // MAXWELL_MAXWELL_HEAD_PROBLEM_HPP
