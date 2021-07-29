// SPDX-FileCopyrightText: 2015 - 2021 Marcin Łoś <marcin.los.91@gmail.com>
// SPDX-License-Identifier: MIT

#ifndef ADS_SIMULATION_SIMULATION_BASE_HPP
#define ADS_SIMULATION_SIMULATION_BASE_HPP

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
    explicit simulation_base(const timesteps_config& steps);

    simulation_base(const simulation_base&) = delete;
    simulation_base(simulation_base&&) = delete;
    simulation_base& operator=(const simulation_base&) = delete;
    simulation_base& operator=(simulation_base&&) = delete;

    virtual ~simulation_base() = 0;

    void run();
};

inline simulation_base::~simulation_base() = default;

}  // namespace ads

#endif  // ADS_SIMULATION_SIMULATION_BASE_HPP
