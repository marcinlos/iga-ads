// SPDX-FileCopyrightText: 2015 - 2021 Marcin Łoś <marcin.los.91@gmail.com>
// SPDX-License-Identifier: MIT

#ifndef ADS_BASIS_DATA_HPP
#define ADS_BASIS_DATA_HPP

#include <utility>
#include <vector>

#include <boost/range/counting_range.hpp>

#include "ads/bspline/bspline.hpp"

namespace ads {

using element_id = int;
using dof_id = int;

// TODO: use proper multidimensional arrays
// NOLINTNEXTLINE(cppcoreguidelines-special-member-functions, hicpp-special-member-functions)
struct basis_data {
    using range_type = decltype(boost::counting_range(0, 0));

    std::vector<int> first_dofs;
    std::vector<std::pair<int, int>> element_ranges;

    int degree;
    int elements;
    int dofs;
    int derivatives;
    int quad_order;
    int elem_division;
    std::vector<double> points;
    bspline::basis basis;
    double**** b;
    double** x;
    const double* w;
    double* J;

    ~basis_data();

    basis_data(const basis_data& other);

    basis_data& operator=(const basis_data& other) {
        auto tmp = other;
        swap(*this, tmp);
        return *this;
    }

    basis_data(bspline::basis basis, int derivatives)
    : basis_data{std::move(basis), derivatives, basis.degree + 1, 1} { }

    basis_data(bspline::basis basis, int derivatives, int quad_order, int elem_division);

    int first_dof(element_id e) const { return first_dofs[e / elem_division]; }

    int last_dof(element_id e) const { return first_dof(e) + dofs_per_element() - 1; }

    int dofs_per_element() const { return degree + 1; }

    range_type dof_range(element_id e) const {
        return boost::counting_range(first_dof(e), last_dof(e) + 1);
    }

    range_type element_range(dof_id dof) const {
        return boost::counting_range(element_ranges[dof].first, element_ranges[dof].second + 1);
    }

    friend void swap(basis_data& a, basis_data& b) noexcept {
        using std::swap;
        swap(a.first_dofs, b.first_dofs);
        swap(a.element_ranges, b.element_ranges);
        swap(a.degree, b.degree);
        swap(a.elements, b.elements);
        swap(a.dofs, b.dofs);
        swap(a.derivatives, b.derivatives);
        swap(a.quad_order, b.quad_order);
        swap(a.elem_division, b.elem_division);
        swap(a.points, b.points);
        swap(a.basis, b.basis);
        swap(a.b, b.b);
        swap(a.x, b.x);
        swap(a.w, b.w);
        swap(a.J, b.J);
    }
};

}  // namespace ads

#endif  // ADS_BASIS_DATA_HPP
