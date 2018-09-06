#ifndef PROBLEMS_STOKES_SPACE_SET_HPP_
#define PROBLEMS_STOKES_SPACE_SET_HPP_

#include "ads/simulation.hpp"


namespace ads {

struct space_set {
    dimension U1x, U1y;
    dimension U2x, U2y;
    dimension Px, Py;
};

}

#endif /* PROBLEMS_STOKES_SPACE_SET_HPP_ */
