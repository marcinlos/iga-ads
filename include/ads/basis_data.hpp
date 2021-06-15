#ifndef ADS_BASIS_DATA_HPP_
#define ADS_BASIS_DATA_HPP_

#include <utility>
#include <vector>

#include <boost/range/counting_range.hpp>

#include "ads/bspline/bspline.hpp"


namespace ads {

typedef int element_id;
typedef int dof_id;


struct basis_data {
    using range_type = decltype(boost::counting_range(0, 0));

    std::vector<int> first_dofs;
    std::vector<std::pair<int, int>> element_ranges;

    int degree;
    int elements;
    int dofs;
    int quad_order;
    int elem_division;
    std::vector<double> points;
    bspline::basis basis;
    double**** b;
    double** x;
    const double* w;
    double* J;

    basis_data(bspline::basis basis, int derivatives)
    : basis_data{ std::move(basis), derivatives, basis.degree + 1, 1 }
    { }

    basis_data(bspline::basis basis, int derivatives, int quad_order, int elem_division);

    int first_dof(element_id e) const {
        return first_dofs[e / elem_division];
    }

    int last_dof(element_id e) const {
        return first_dof(e) + dofs_per_element() - 1;
    }

    int dofs_per_element() const {
        return degree + 1;
    }

    range_type dof_range(element_id e) const {
        return boost::counting_range(first_dof(e), last_dof(e) + 1);
    }

    range_type element_range(dof_id dof) const {
        return boost::counting_range(element_ranges[dof].first, element_ranges[dof].second + 1);
    }
};

}

#endif /* ADS_BASIS_DATA_HPP_ */
