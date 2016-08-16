#ifndef ADS_SIMULATION_SIMULATION_BASE_HPP_
#define ADS_SIMULATION_SIMULATION_BASE_HPP_

#include "ads/simulation/config.hpp"


namespace ads {


class simulation_base {
protected:
    timesteps_config steps;

private:
    virtual void before() { }

    virtual void after() { }

    virtual void before_step(int /*iter*/, double /*t*/) { }

    virtual void step(int /*iter*/, double /*t*/) { }

    virtual void after_step(int /*iter*/, double /*t*/) { }

public:
    simulation_base(const timesteps_config& steps);

    virtual ~simulation_base() = 0;

    void run();
};

inline simulation_base::~simulation_base() = default;

}


#endif /* ADS_SIMULATION_SIMULATION_BASE_HPP_ */
