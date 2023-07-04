// SPDX-FileCopyrightText: 2015 - 2021 Marcin Łoś <marcin.los.91@gmail.com>
// SPDX-License-Identifier: MIT

#ifndef ADS_EXPERIMENTAL_HORRIBLE_SPARSE_MATRIX_HPP
#define ADS_EXPERIMENTAL_HORRIBLE_SPARSE_MATRIX_HPP

#include <cstddef>
#include <unordered_map>

#include "ads/experimental/all.hpp"
#include "ads/solver/mumps.hpp"

namespace ads {

class horrible_sparse_matrix {
private:
    struct index_pair {
        ads::global_dof i, j;
    };

    struct hasher {
        auto operator()(index_pair const& e) const noexcept -> std::size_t {
            return e.i ^ (e.j + 0x9e3779b9 + (e.i << 6) + (e.i >> 2));
        }
    };

    struct equality {
        auto operator()(index_pair const& a, index_pair const& b) const noexcept -> bool {
            return a.i == b.i && a.j == b.j;
        }
    };

    std::unordered_map<index_pair, double, hasher, equality> storage_;

public:
    auto operator()(ads::global_dof i, ads::global_dof j) -> double& {
        return storage_[index_pair{i, j}];
    }

    auto has_entry(ads::global_dof i, ads::global_dof j) const -> bool {
        return storage_.find(index_pair{i, j}) != storage_.end();
    }

    auto mumpsify(ads::mumps::problem& problem) const -> void {
        for (auto const& [idx, val] : storage_) {
            if (val != 0) {
                problem.add(idx.i + 1, idx.j + 1, val);
            }
        }
    }

    auto begin() const { return std::cbegin(storage_); }

    auto end() const { return std::cend(storage_); }
};

}  // namespace ads

#endif  // ADS_EXPERIMENTAL_HORRIBLE_SPARSE_MATRIX_HPP
