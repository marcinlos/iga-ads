// SPDX-FileCopyrightText: 2015 - 2023 Marcin Łoś <marcin.los.91@gmail.com>
// SPDX-License-Identifier: MIT

#ifndef CO2_SEQUESTRATION_CO2_SEQUESTRATION_2D_TRIM_HPP
#define CO2_SEQUESTRATION_CO2_SEQUESTRATION_2D_TRIM_HPP

#include <galois/Timer.h>
#include <boost/format.hpp>

#include "ads/executor/galois.hpp"
#include "ads/output_manager.hpp"
#include "ads/simulation.hpp"
#include "ads/solver/mumps.hpp"


#include "ads/lin/tensor/tensor.hpp"

namespace ads::problems {

class co2_sequestration_2d_trim : public simulation_2d {
private:
    using Base = simulation_2d;
    vector_type p, s, p_prev, s_prev;

    output_manager<2> output;
    galois_executor executor{4};
    galois::StatTimer integration_timer{"integration"};
    ads::mumps::solver solver;

    // parameter maps
    using map_array = ads::lin::tensor<double, 2>;
    map_array porosity;
    map_array permeability;

    // parameters
    const std::array<int, 2> resolution;  // resolution of the PINN grid / S&P data i/o
    const int iter;        // current iteration
    const double mu_w;     // brine viscosity
    const double mu_g;     // gas viscosity
    const double K;        // permeability tensor
    const double phi;      // porosity
    const double rho_w;    // brine density
    const double rho_g;    // gas density
    const double g;        // gravitational acceleration
    const double qg_x;     // gas injection location - x
    const double qg_y;     // gas injection location - y
    const double qg_rate;  // gas injection rate
    const int mesh_x;
    const int mesh_y;
    std::string porosity_map;
    std::string permeability_map;
    bool verbose;

public:
    explicit co2_sequestration_2d_trim(const config_2d& config, const std::array<int, 2>& resolution, const int iter,
                                  const double mu_w, const double mu_g, const double K, const double phi,
                                  const double rho_w, const double rho_g, const double g,
                                  const double qg_x, const double qg_y, const double qg_rate,
                                  std::string porosity_map, std::string permeability_map, bool verbose)
    : Base{config}
    , p{shape()}
    , s{shape()}
    , p_prev{shape()}
    , s_prev{shape()}
    , resolution{resolution}
    , iter{iter}
    , mu_w{mu_w}
    , mu_g{mu_g}
    , K{K}
    , phi{phi}
    , rho_w{rho_w}
    , rho_g{rho_g}
    , g{g}
    , qg_x{qg_x}
    , qg_y{qg_y}
    , qg_rate{qg_rate}
    , verbose{verbose}
    , porosity_map{porosity_map}
    , permeability_map{permeability_map}
    , mesh_x{config.x.b}
    , mesh_y{config.y.b}
    , output{x.B, y.B, 2 * config.x.elements, 2 * config.y.elements}
    , porosity{{1, 1}}
    , permeability{{1, 1}} { }

    // this function sets the initial state of the gas saturation
    double init_state(double x, double y) {
        /*
        initial state of the gas saturation - zero for now
        */
        return 0;
    };

    double source_g(double x, double y, double t) {
        double dx = x - qg_x;
        double dy = y - qg_y;
        double r2 = std::min(0.5 * (dx * dx + dy * dy), 1.0);
        return qg_rate * ((r2 - 1) * (r2 - 1) * (r2 + 1) * (r2 + 1));
    }

    double source_w(double x, double y, double t) {
        double val = 0;
        return val;
    }

private:
    auto read_map_data(std::string_view filename) -> ads::lin::tensor<double, 2> {
        std::ifstream input{std::string(filename)};

        int nx;
        int ny;
        input >> nx >> ny;
        auto data = ads::lin::tensor<double, 2>{{nx, ny}};

        for (int i = 0; i < nx; i++) {
            for (int j = 0; j < ny; j++) {
                double x;
                double y;
                double val;
                input >> x >> y >> val;
                data(i, j) = val;
            }
        }
        return data;
    }

    static auto point_lookup(double x, int n) -> std::tuple<int, int> {
        int min_1 = std::min(static_cast<int>(x * n), n - 1);
        int min_2 = std::min(static_cast<int>(min_1 + 1), n);

        if (min_1 == min_2) {
            min_2 = min_1 + 1;
            std::cout << "min_1 == min_2" << std::endl;
        }

        return std::make_tuple(min_1, min_2);
    }

    auto approximate_map_data(point_type xy, map_array map) -> double {
        double x = xy[0];
        double y = xy[1];
        double mesh_xf = mesh_x;
        double mesh_yf = mesh_y;
        const auto [ix1, ix2] = point_lookup(x / mesh_xf, map.size(0) - 1);
        const auto [iy1, iy2] = point_lookup(y / mesh_yf, map.size(1) - 1);

        // we will use bilinear interpolation here
        double f_x1y1 = map(ix1, iy1);
        double f_x2y1 = map(ix2, iy1);
        double f_x1y2 = map(ix1, iy2);
        double f_x2y2 = map(ix2, iy2);

        double x1 = ix1 * mesh_xf / (map.size(0) - 1);
        double x2 = ix2 * mesh_xf / (map.size(0) - 1);
        double y1 = iy1 * mesh_yf / (map.size(1) - 1);
        double y2 = iy2 * mesh_yf / (map.size(1) - 1);

        double fx_y1 = f_x1y1 + (f_x2y1 - f_x1y1) * (x - x1) / (x2 - x1);
        double fx_y2 = f_x1y2 + (f_x2y2 - f_x1y2) * (x - x1) / (x2 - x1);
        double fxy = fx_y1 + (fx_y2 - fx_y1) * (y - y1) / (y2 - y1);

        if (std::isnan(fxy)) {
            std::cout << x << " " << y << std::endl;
            std::cout << ix1 << " " << ix2 << " " << iy1 << " " << iy2 << std::endl;
            std::cout << "Interpolation function testing at: x = " << x << ", y = " << y << std::endl;
            std::cout << "x1: " << x1 << ", x2: " << x2 << ", y1: " << y1 << ", y2: " << y2 << std::endl;
            std::cout << "f_x1y1: " << f_x1y1 << ", f_x2y1: " << f_x2y1 << ", f_x1y2: " << f_x1y2
                      << ", f_x2y2: " << f_x2y2 << std::endl;
            std::cout << "fx_y1: " << fx_y1 << ", fx_y2: " << fx_y2 << std::endl;
            std::cout << "fxy: " << fxy << std::endl;
            std::cout << "map size: " << map.size(0) << ", " << map.size(1) << std::endl;
            std::cout << std::endl;
        }

        return fxy;
    }

    void visualize_map_data(map_array map, const std::string& output_file) {

        std::ofstream ofs(output_file);
        if (ofs.is_open()) {
            for (int i = 0; i < mesh_x; i++) {
                for (int j = 0; j < mesh_y * 4; j++) {
                    double x = i * mesh_x / (mesh_x);
                    double y = j * mesh_y / (mesh_y * 4);
                    double val = approximate_map_data({x, y}, map);
                    ofs << x << " " << y << " " << val << std::endl;
                }
            }
        }
    }

    void output_s_to_pinn(const vector_type& v, const std::array<int, 2>& resolution, const int iter) {
        output_manager<2> output_pinn{x.B, y.B, resolution[0], resolution[1]};
        output_pinn.to_file(v, "s_to_pinn_%d.data", iter);
    }

    void read_p_from_pinn(vector_type& v, const int iter) {
        auto name = str(boost::format("p_from_pinn_%d.data") % iter);
        vector_type data = read_map_data(name);


    }

    void compute_rhs_simple(double t, vector_type& v, const vector_type& v_pinn) {

        executor.for_each(elements(), [&](index_type e) {
            auto U = element_rhs();

            double J = jacobian(e);
            for (auto q : quad_points()) {
                double w = weight(q);
                auto x = point(e, q);
                for (auto a : dofs_on_element(e)) {
                    auto aa = dof_global_to_local(e, a);
                    value_type vv = eval_basis(e, q, a);

                    auto fval = approximate_map_data(x, v_pinn);

                    double val = fval * vv.val;
                    U(aa[0], aa[1]) += val * w * J;
                }
            }

            executor.synchronized([&]() { update_global_rhs(v, U, e); });
        });
    }

    void output_s_to_iga(const vector_type& v, const int iter) {
        auto name = str(boost::format("s_to_iga_%d.data") % iter);
        std::ofstream ofs(name, std::ios::out | std::ios::trunc);
        if (ofs.is_open()) {
            for (auto idx : dofs()) {
                double val = v(idx[0], idx[1]);
                ofs << idx[0] << " " << idx[1] << " " << val << "\n";
            }
            ofs.close();
        }
    }

    void read_s_from_iga(vector_type& v, const int iter) {
        auto name = str(boost::format("s_from_iga_%d.data") % iter);
        std::ifstream ifs(name, std::ios::in);
        if (ifs.is_open()) {
            for (auto idx : dofs()) {
                double val;
                ifs >> idx[0] >> idx[1] >> val;
                v(idx[0], idx[1]) = val;
            }
            ifs.close();
        } else {
            std::cerr << "Failed to open s_from_iga file for reading." << std::endl;
        }
    }

    void solve(vector_type& v) {
        // lin::vector buf{{y.dofs()}};
        // compute_projection(buf, y.basis, [](double y) { return std::sin(y * M_PI); });
        // for (int i = 0; i < y.dofs(); ++i) {
        //    v(0, i) = buf(i);
        // }
        Base::solve(v);
    }

    void prepare_matrices() {
        x.fix_left();
        x.fix_right();
        // y.fix_left();
        // y.fix_right();

        // !!all the bcs are Dirichlet for now - this is not correct for the final version!!
        Base::prepare_matrices();
    }

    void before() override {
        prepare_matrices();

        if (porosity_map != "none") {
            porosity = read_map_data(porosity_map);
            if (iter == 0) {
                visualize_map_data(porosity, "porosity_visualization.dat");
                if (verbose) {
                    std::cout << "Porosity map interpolated and visualized." << std::endl;
                }
            }
        }

        // if (permeability_map != "none") {
            // permeability = read_map_data("permeability_k1.data");
        // }

    }

    void before_step(int /*iter*/, double /*t*/) override {
        // load the S and P data from files
        // load the coefficients!!! from the previous step
        if (iter > 0) {
            auto p_from_data =
            auto s_from_data =

            p_prev = p_from_data;
            s_prev = s_from_data;
        }
        else{
            // if this is the first iteration, we use initial values
            auto init = [this](double x, double y) { return init_state(x, y); };
            projection(s_prev, init);
            solve(s_prev);
            p_prev.fill_with_zeros();
        }
    }

    void step(int /*iter*/, double t) override {
        // solve for S after loading data for P
        compute_rhs(t);
        dirichlet_bc(s, boundary::left, x, y, [](double t) { return 0; });
        dirichlet_bc(s, boundary::right, x, y, [](double t) { return 0; });
        solve(s);
    }

    void after_step(int iter, double /*t*/) override {
        // output the S coefficients

        // save the plotting data
        output.to_file(p, "p.out_%d.data", iter);
        output.to_file(s, "s.out_%d.data", iter);
        if (verbose) {
            std::cout << "Iteration " << iter << " passed" << std::endl;
        }
    }

    void compute_rhs(double t) {
        auto& rhs = s;

        zero(rhs);

        executor.for_each(elements(), [&](index_type e) {
            auto U = element_rhs();

            double J = jacobian(e);
            for (auto q : quad_points()) {
                double w = weight(q);
                auto x = point(e, q);
                for (auto a : dofs_on_element(e)) {
                    auto aa = dof_global_to_local(e, a);
                    value_type v = eval_basis(e, q, a);
                    value_type s = eval_fun(s_prev, e, q);
                    value_type p_here = eval_fun(p, e, q);

                    double s_val = std::clamp(s.val, 0.0, 1.0);
                    double term_1 = s_val * grad_dot(p_here, v) * K / mu_g;
                    double term_2 = s_val * v.dy * K * rho_g * g / mu_g;
                    double term_3 = v.val * source_g(x[0], x[1], t);

                    double phi_here;
                    if (porosity_map != "none") {
                        phi_here = approximate_map_data(x, porosity);
                    } else {
                        phi_here = phi;
                    }

                    // temporary porosity increase
                    phi_here += 0.1;
                    double val = (term_2 + term_3 - term_1) * steps.dt / phi_here + s_val * v.val;

                    // NOTE! this term is a temporary enforcement of the upper ceiling on saturation
                    double term_extra = -1 * s.val * v.val * (x[1] >= mesh_y - 1);
                    val += term_extra;

                    U(aa[0], aa[1]) += val * w * J;
                }
            }

            executor.synchronized([&]() { update_global_rhs(rhs, U, e); });
        });
    }

};

}  // namespace ads::problems

#endif
