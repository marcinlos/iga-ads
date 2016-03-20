#ifndef ADS_PROBLEMS_TUMOR_VASCULATURE_DEFS_HPP_
#define ADS_PROBLEMS_TUMOR_VASCULATURE_DEFS_HPP_

#include <cstddef>
#include "ads/lin/tensor.hpp"
#include "ads/util/math/vec.hpp"


namespace ads {
namespace tumor {
namespace vasc {

constexpr std::size_t Dim = 2;

using vector = math::vec<Dim>;
using val_array = lin::tensor<double, Dim>;


}
}
}


#endif /* ADS_PROBLEMS_TUMOR_VASCULATURE_DEFS_HPP_ */
