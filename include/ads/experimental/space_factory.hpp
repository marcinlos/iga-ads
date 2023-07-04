// SPDX-FileCopyrightText: 2015 - 2023 Marcin Łoś <marcin.los.91@gmail.com>
// SPDX-License-Identifier: MIT

#ifndef ADS_EXPERIMENTAL_SPACE_FACTORY_HPP
#define ADS_EXPERIMENTAL_SPACE_FACTORY_HPP

#include <utility>

#include "ads/experimental/all.hpp"

namespace ads {

class space_factory {
private:
    ads::global_dof offset_ = 0;

public:
    template <typename Space, typename... Args>
    auto next(Args&&... args) -> Space {
        auto space = Space{std::forward<Args>(args)..., offset_};
        offset_ += space.dof_count();
        return space;
    }

    auto dim() const -> ads::global_dof { return offset_; };
};

}  // namespace ads

#endif  // ADS_EXPERIMENTAL_SPACE_FACTORY_HPP
