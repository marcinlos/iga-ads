#ifndef ADS_PROBLEMS_CH_2D_HPP_
#define ADS_PROBLEMS_CH_2D_HPP_

#include "ads/simulation.hpp"
#include "ads/output_manager.hpp"
#include "ads/executor/galois.hpp"
#include "ads/lin/dense_matrix.hpp"
#include "ads/lin/dense_solve.hpp"


namespace ads {

class ch_2d : public simulation_2d {
private:
    using Base = simulation_2d;
    using vector_view = lin::tensor_view<double, 2>;

    std::vector<double> sol, buf;	// Buffers for keeping the combined RHS vectors

    // Two sets of vectors for two coupled equations
    vector_type u, u_prev;
    vector_type u2, u2_prev;

    double theta = 1.5;
    double lambda;	// = 1.0/(config.dim.n*config.dim.n);		// Mesh size parameter

    double eps = 1.0;	//0.033333;
    double P = 0.0;
    double chi = 1.0;	// Gomez, Hughes 2011

    double curt = 0.0;	// Current time

    double GL_free_energy = 0.0;	// current GL free energy

    output_manager<2> output;
    galois_executor executor{4};

    // Large matrices - dense matrices + data for solver
    lin::dense_matrix Ax, Ay;
    lin::solver_ctx Ax_ctx, Ay_ctx;

    // BC-NEW: Mass matrices for the periodic basis
    lin::dense_matrix Mx, My;
    lin::solver_ctx Mx_ctx, My_ctx;

public:
    ch_2d(const config_2d& config)
    : Base{ config }
    , lambda(1.0/(9000.0))		// Gomez et al. 2008
//    , lambda(1.0/(config.x.elements*config.y.elements))	// Some other papers
    , sol(2 * (x.dofs() - 1) * (y.dofs() - 1))
    , buf(2 * (x.dofs() - 1) * (y.dofs() - 1))
    , u{ shape() }
    , u_prev{ shape() }
    , u2{ shape() }
    , u2_prev{ shape() }
    , output{ x.B, y.B, 200 }
    , Ax{ 2 * (x.dofs() - 1), 2 * (x.dofs() - 1) } // 2x2 matrix for x direction
    , Ay{ 2 * (y.dofs() - 1), 2 * (y.dofs() - 1) } // 2x2 matrix for y direction
    , Ax_ctx{ Ax }
    , Ay_ctx{ Ay }
    , Mx{ x.dofs() - 1, x.dofs() - 1 } // BC-NEW: initialization of mass matrices
    , My{ y.dofs() - 1, y.dofs() - 1 }
    , Mx_ctx{ Mx }
    , My_ctx{ My }
    {
        // Fill the large matrices
        matrix(Ax, x.basis, steps.dt);
        matrix(Ay, y.basis, steps.dt);

        // BC-NEW: Fill the mass matrices
        mass_matrix(Mx, x.basis, steps.dt);
        mass_matrix(My, y.basis, steps.dt);

        // Init random
        srand(time(NULL));

/*
	// Check 1D matrices assembly
        std::cout << "Mx " << std::endl;
        std::cout << Mx << std::endl << std::endl;

        std::cout << "My " << std::endl;
        std::cout << My << std::endl << std::endl;

        std::cout << "Ax " << std::endl;
        std::cout << Ax << std::endl << std::endl;

        std::cout << "Ay " << std::endl;
        std::cout << Ay << std::endl << std::endl;

        exit(1);
*/
    }

    double init_c(double x, double y) {
	// Circle in the center
        if (sqrt((x-0.5)*(x-0.5)+(y-0.5)*(y-0.5))<=0.1)
            return 0.88;
        else
            return 0.12;

        // Here we use a random distribution following Gomez et al. (2008)
        double cbar = 0.63; // 0.26;  // 0.7;
        double r = ((rand() % 101) - 50)/1000.;	// random number in [-0.05; 0.05]
        return cbar-r;
    }

/*
    double init_eta(double x, double y) {
        // For the moment we just return F(c)
        return 0; // get_F(init_c(x, y)) + init_c(x, y);
    }
*/

void compute_projection_with_laplacian() {
	// c, eta - components of the solution
        for (auto e : elements()) {
            double J = jacobian(e);
            for (auto q : quad_points()) {
                double w = weigth(q);
                for (auto a : dofs_on_element(e)) {
                    value_type v = eval_basis(e, q, a);
                    value_type cc = eval_fun(u, e, q);
                    double val = get_F(cc.val) * v.val - grad_dot(cc, v); // laplacian integrated by parts
                    u2(a[0], a[1]) += val * w * J;
                }
            }
        }
    }

private:
    // BC-NEW: Mass matrix no longer diagonal
    void mass_matrix(lin::dense_matrix& M, const basis_data& d, double h) {
        auto N = d.basis.dofs() - 1;				// Size of Mx and My
        for (element_id e = 0; e < d.elements; ++ e) {
            for (int q = 0; q < d.quad_order; ++ q) {
                int first = d.first_dof(e);
                int last = d.last_dof(e);
                for (int a = 0; a + first <= last; ++ a) {
                    for (int b = 0; b + first <= last; ++ b) {
                        int ia = a + first;			// ia and ib change from 0 to N (inclusive)
                        int ib = b + first;
                        auto va = d.b[e][q][0][a];
                        auto vb = d.b[e][q][0][b];
//std::cout<<e<<" "<<ia<<" "<<ib<<", "<<va<<" "<<vb<<" ===> "<<ia % N<<" "<<ib % N<<std::endl;
                        M(ia % N, ib % N) += va * vb * d.w[q] * d.J[e]; // upper left
                    }
                }
            }
        }
    }

    void matrix(lin::dense_matrix& K, const basis_data& d, double h) {
        auto N = d.basis.dofs() - 1;
        for (element_id e = 0; e < d.elements; ++ e) {
            for (int q = 0; q < d.quad_order; ++ q) {
                int first = d.first_dof(e);
                int last = d.last_dof(e);
                for (int a = 0; a + first <= last; ++ a) {
                    for (int b = 0; b + first <= last; ++ b) {
                        int ia = a + first;
                        int ib = b + first;
                        auto va = d.b[e][q][0][a];
                        auto vb = d.b[e][q][0][b];
                        auto da = d.b[e][q][1][a];
                        auto db = d.b[e][q][1][b];
                        //K(ia, ib) += (va * vb + 0.5 * h * da * db) * d.w[q] * d.J[e]; // upper left
                        //K(N + ia, N + ib) += (va * vb + 0.5 * h * da * db) * d.w[q] * d.J[e]; // lower right

//			if (ia != 0 && ia != N-1 && ia != N && ia != 2*N-1 &&
//				ib != 0 && ib != N-1 && ib != N && ib != 2*N-1) {
	                        // Cahn-Hilliard Mass-Zero-Stiffness-Mass matrix
	                        K(ia % N, ib % N) += va * vb * d.w[q] * d.J[e];				// upper left
	                        K(N + ia % N, N + ib % N) += va * vb * d.w[q] * d.J[e];			// lower right
	                        K(N + ia % N, ib % N) += - eps * eps * da * db * d.w[q] * d.J[e];	// lower left

	                        K(ia % N, N + ib % N) += 0.5 * steps.dt * get_M(va) * da * db * d.w[q] * d.J[e];	// upper right (new idea)
//			}
                    }
                }
            }
        }

        // Dirichlet BC
//	K(0, 0) = 1;
//	K(N-1, N-1) = 1;
//	K(N, N) = 1;
//	K(2*N-1, 2*N-1) = 1;
//	K(N, 0) = 1;
//	K(2*N-1, N-1) = 1;
    }

    void prepare_matrices() {
        // Dirichlet BC
//        x.fix_left();
//        x.fix_right();
//        y.fix_left();
//        y.fix_right();

        Base::prepare_matrices();
        // Dense factorization (overloaded function)
        lin::factorize(Ax, Ax_ctx);
        lin::factorize(Ay, Ay_ctx);
        // BC-NEW: factorize mass matrices
        lin::factorize(Mx, Mx_ctx);
        lin::factorize(My, My_ctx);

    }

    void before() override {
        prepare_matrices();

        auto init = [this](double x, double y) { return init_c(x, y); };
        projection(u, init);
        solve(u);
        output.to_file(u, "OUT/out1_0.data");


        //auto init2 = [this](double x, double y) { return init_eta(x, y); };
        //projection(u2, init2);
        compute_projection_with_laplacian();
        solve(u2);
        output.to_file(u2, "OUT/out2_0.data");
    }

    void before_step(int /*iter*/, double /*t*/) override {
        using std::swap;
        swap(u, u_prev);
        swap(u2, u2_prev);
    }

/*
    double compute_GL_free_energy(double c, double c_x, double c_y) {
        double alpha = 1.0/(3*lambda); // 3000.0;
        double gradc_norm = c_x*c_x+c_y*c_y;

        return c*log(c) + (1.0-c)*log(1.0-c) + 2.0*theta*c*(1.0-c) + theta/(3.0*alpha)*gradc_norm;
    }
*/

    void step(int /*iter*/, double /*t*/) override {
        compute_rhs_1();

	int boundary = 0;	// 1 for Dirichlet, 0 for Neumann

        // BC-NEW: different size and indices are taken mod Nx, Ny to ensure periodic BC
        int Nx = x.dofs() - 1;
        int Ny = y.dofs() - 1;
        // Copy data from separate RHS vectors to the combined one
        // vector_view view_x{ sol.data(), {2*x.dofs(), y.dofs()}};
        std::fill(begin(sol), end(sol), 0);		// BC-NEW
        vector_view view_x{ sol.data(), {2*Nx, Ny}};

        for (auto i = boundary; i < x.dofs()-boundary; ++ i) {
            for (auto j = boundary; j < y.dofs()-boundary; ++ j) {
                // view_x(i, j) = u(i, j);
                // view_x(i + x.dofs(), j) = u2(i, j);
                view_x(i % Nx, j % Ny) += u(i, j);
                view_x(Nx + i % Nx, j % Ny) += u2(i, j);
            }
        }

        // ads_solve(u, buffer, dim_data{Kx, x.ctx}, y.data());
        // ads_solve(u2, buffer, dim_data{Kx, x.ctx}, y.data());
        // Expanded ADS call:
        lin::solve_with_factorized(Ax, view_x, Ax_ctx);			// 2x2
        auto F = lin::cyclic_transpose(view_x, buf.data());

        // BC-NEW: periodic mass matrix instead of a standard one
        // lin::solve_with_factorized(y.data().M, F, y.data().ctx);
        lin::solve_with_factorized(My, F, My_ctx);
        lin::cyclic_transpose(F, view_x);

        // ... and copy the solutions back to the separate vectors
        for (auto i = boundary; i < x.dofs()-boundary; ++ i) {
            for (auto j = boundary; j < y.dofs()-boundary; ++ j) {
                // u(i, j) = view_x(i, j);
                // u2(i, j) = view_x(i + x.dofs(), j);
                u(i, j) = view_x(i % Nx, j % Ny);
                u2(i, j) = view_x(Nx + i % Nx, j % Ny);
            }
        }

        using std::swap;
        swap(u, u_prev);
        swap(u2, u2_prev);

        compute_rhs_2();

        // Copy data from separate RHS vectors to the combined one
        // vector_view view_y{ sol.data(), {x.dofs(), 2*y.dofs()}};
        std::fill(begin(sol), end(sol), 0);			// BC-NEW
        vector_view view_y{ sol.data(), {Nx, 2*Ny}};
        for (auto i = boundary; i < x.dofs()-boundary; ++ i) {
            for (auto j = boundary; j < y.dofs()-boundary; ++ j) {
                // view_y(i, j) = u(i, j);
                // view_y(i, j + y.dofs()) = u2(i, j);
                view_y(i % Nx, j % Ny) += u(i, j);
                view_y(i % Nx, Ny + j % Ny) += u2(i, j);
            }
        }

        // ads_solve(u, buffer, x.data(), dim_data{Ky, y.ctx});
        // ads_solve(u2, buffer, x.data(), dim_data{Ky, y.ctx});
        // Expanded ADS call:
        // NEW: periodic mass matrix instead of a standard one
        // lin::solve_with_factorized(x.data().M, view_y, x.data().ctx);
        lin::solve_with_factorized(Mx, view_y, Mx_ctx);			// Mx
        auto F2 = lin::cyclic_transpose(view_y, buf.data());
        lin::solve_with_factorized(Ay, F2, Ay_ctx);			// 2x2
        lin::cyclic_transpose(F2, view_y);

        // ... and copy the solutions back to the separate vectors
        for (auto i = boundary; i < x.dofs()-boundary; ++ i) {
            for (auto j = boundary; j < y.dofs()-boundary; ++ j) {
                // u(i, j) = view_y(i, j);
                // u2(i, j) = view_y(i, j + y.dofs());
                u(i, j) = view_y(i % Nx, j % Ny);
                u2(i, j) = view_y(i % Nx, Ny + j % Ny);
            }
        }


    }

    void after_step(int iter, double /*t*/) override {
        curt += steps.dt;

        if ((iter+1) % 1 == 0)
            std::cout << "Iter " << iter+1 << ":\tt = "<< curt << "\tdt = "<< steps.dt << "\t(" << GL_free_energy << ")" << std::endl;
        if ((iter+1) % 10 == 0) {
            output.to_file(u, "OUT/out1_%d.data", iter+1);
            output.to_file(u2, "OUT/out2_%d.data", iter+1);
        }
    }

    // Mobility M(c)
    double get_M(double c) {
        return 1.0; // 200.0;
//        return c*(1-c);
    }

    // Chemical potential F(c)
    double get_F(double c) {
//        double gamma = 0.045;
//        return gamma*(4.0*c*c*c-6.0*c*c+2.0*c);

	double norm = 1.0*1.0/lambda;				// L0^2/lambda

        return norm*(1./(2*theta)*log(c/(1-c)) + 1 - 2*c) /* * (2*theta)*/;	// Cortes1 also multiplies by 2*theta (mistake?)
//        return norm*(c*log(c) + (1-c)*log(1-c) + 2*theta*c*(1-c));		// Gomez, Hughes 2011
////        return 1./(2*theta) + 1 - 2*c;
    }

    void compute_rhs_1() {
        auto& rhs = u;
        auto& rhs2 = u2;

        zero(rhs);
        zero(rhs2);

        GL_free_energy = 0.0;

        executor.for_each(elements(), [&](index_type e) {
            auto U = element_rhs();
            auto U2 = element_rhs();

            double J = jacobian(e);
            for (auto q : quad_points()) {
                double w = weigth(q);
                value_type u;

                for (auto a : dofs_on_element(e)) {
                    auto aa = dof_global_to_local(e, a);
                    value_type v = eval_basis(e, q, a);
                    u = eval_fun(u_prev, e, q);
                    value_type u2 = eval_fun(u2_prev, e, q);

                    // Cahn-Hilliard RHS Step 1
                    double RHS1 = u.val * v.val - 0.5 * steps.dt * get_M(u.val) * (/*u2.dx * v.dx +*/ u2.dy * v.dy);
                    double RHS2 = eps * eps * u.dy * v.dy;
                    double Fval = get_F(u.val) * v.val;
//std::cout<<RHS1<<" "<<RHS2<<" "<<Fval<<" "<<w<<" "<<J<<std::endl;
                    U(aa[0], aa[1]) += RHS1 / (1.0-P) * w * J;
                    U2(aa[0], aa[1]) += (RHS2 + chi * Fval) * w * J;
                }

                double alpha = 1.0/(3*lambda); // 3000.0;
                double gltemp = u.val*log(u.val) + (1.0-u.val)*log(1.0-u.val) + 2.0*theta*u.val*(1.0-u.val) + theta/(3.0*alpha)*grad_dot(u, u); //(u.dx*u.dx+u.dy*u.dy);
                GL_free_energy += gltemp * w * J;

                //double gltemp = compute_GL_free_energy(u.val, u.dx, u.dy) - grad_dot(u, v);
                //GL_free_energy += gltemp * w * J; //* v.val
            }

            executor.synchronized([&]() {
                update_global_rhs(rhs, U, e);
                update_global_rhs(rhs2, U2, e);
            });
        });
    }

    void compute_rhs_2() {
        auto& rhs = u;
        auto& rhs2 = u2;

        zero(rhs);
        zero(rhs2);

        executor.for_each(elements(), [&](index_type e) {
            auto U = element_rhs();
            auto U2 = element_rhs();

            double J = jacobian(e);
            for (auto q : quad_points()) {
                double w = weigth(q);
                for (auto a : dofs_on_element(e)) {
                    auto aa = dof_global_to_local(e, a);
                    value_type v = eval_basis(e, q, a);
                    value_type u = eval_fun(u_prev, e, q);
                    value_type u2 = eval_fun(u2_prev, e, q);

                    // Cahn-Hilliard RHS Step 2
                    double RHS1 = u.val * v.val - 0.5 * steps.dt * get_M(u.val) * (u2.dx * v.dx /*+ u2.dy * v.dy*/);
                    double RHS2 = eps * eps * u.dx * v.dx;
                    double Fval = get_F(u.val) * v.val;
//std::cout<<RHS1<<" "<<RHS2<<" "<<Fval<<" "<<w<<" "<<J<<std::endl;
                    U(aa[0], aa[1]) += RHS1 / (1.0-P) * w * J;
                    U2(aa[0], aa[1]) += (RHS2 + chi * Fval) * w * J;
                }
            }

            executor.synchronized([&]() {
                update_global_rhs(rhs, U, e);
                update_global_rhs(rhs2, U2, e);
            });
        });
    }

};

}

#endif /* ADS_PROBLEMS_CH_2D_HPP_ */
