#ifndef ADS_SIMULATION_CONFIG_HPP_
#define ADS_SIMULATION_CONFIG_HPP_


namespace ads {


struct dim_config {
    int p;
    int elements;
    double a;
    double b;
    int quad_order;

    dim_config(int p, int elements, double a, double b, double quad_order)
    : p{p}
    , elements{elements}
    , a{a}
    , b{b}
    , quad_order{quad_order}
    { }

    dim_config(int p, int elements, double a = 0, double b = 1): dim_config{ p, elements, a, b, p + 1 }
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
