#ifndef MUMPS_HPP_
#define MUMPS_HPP_

#include <dmumps_c.h>
#include <vector>


namespace mumps {


struct problem {

    problem(double* rhs, int n)
    : rhs_{ rhs }, n{ n }
    { }

    problem(std::vector<double>& rhs)
    : rhs_{ rhs.data() }, n{ rhs.size() }
    { }

    void add(int row, int col, double value) {
        rows_.push_back(row);
        cols_.push_back(col);
        values_.push_back(value);
    }

    int nonzero_entries() {
        return values_.size();
    }

    int dofs() {
        return n;
    }

    int* irn() {
        return rows_.data();
    }

    int* jcn() {
        return cols_.data();
    }

    double* a() {
        return values_.data();
    }

    double* rhs() {
        return rhs_;
    }

private:
    std::vector<int> rows_;
    std::vector<int> cols_;
    std::vector<double> values_;
    double* rhs_;
    int n;
};


class solver {
private:
    DMUMPS_STRUC_C id;

public:

    solver() {
        id.job = -1;
        id.par = 1;
        id.sym = 0;
        id.comm_fortran = 1;
        dmumps_c(&id);

        //streams
        id.icntl[0] = 6;
        id.icntl[1] = 6;
        id.icntl[2] = 6;
        id.icntl[3] = 3;

        id.icntl[0] = 3;
        id.icntl[1] = 3;
        id.icntl[2] = 3;
        id.icntl[3] = 3;


        //ordering metis (5), or pord (4), or AMD (0), AMF (2), QAMD (6)
        id.icntl[6] = 2;
    }

    void solve(problem& problem) {
        id.n = problem.dofs();
        id.nz = problem.nonzero_entries();

        id.irn = problem.irn();
        id.jcn = problem.jcn();
        id.a = problem.a();

        id.rhs = problem.rhs();
        id.nrhs = 1;
        id.lrhs = id.n;

        analyze_();
        factorize_();
        solve_();
    }

    ~solver() {
        id.job = -2;
        dmumps_c(&id);
    }

private:
    void factorize_() {
        id.job = 2;
        dmumps_c(&id);
    }

    void analyze_() {
        id.job = 1;
        dmumps_c(&id);
    }

    void solve_() {
        id.job = 3;
        dmumps_c(&id);
    }
};


}


#endif /* MUMPS_HPP_ */
