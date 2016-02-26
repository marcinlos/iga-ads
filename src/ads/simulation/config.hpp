#ifndef ADS_SIMULATION_CONFIG_HPP_
#define ADS_SIMULATION_CONFIG_HPP_


namespace ads {


struct dim_config {
    int p;
    int elements;
    double a;
    double b;

    dim_config(int p, int elements, double a = 0, double b = 1)
    : p{p}
    , elements{elements}
    , a{a}
    , b{b}
    { }
};


struct timesteps_config {
    int step_count;
    double dt;

    timesteps_config(int step_count, double dt)
    : step_count{step_count}, dt{dt}
    { }
};


struct config_1d {
    dim_config x;
    timesteps_config steps;
    int derivatives;
};

struct config_2d {
    dim_config x, y;
    timesteps_config steps;
    int derivatives;
};


struct config_3d {
    dim_config x, y, z;
    timesteps_config steps;
    int derivatives;
};



}


#endif /* ADS_SIMULATION_CONFIG_HPP_ */
