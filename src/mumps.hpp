#ifndef MUMPS_HPP_
#define MUMPS_HPP_

#include <dmumps_c.h>
#include <vector>
#include <cstring>
#include <iostream>
#include <mpi.h>


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

    int nonzero_entries() const {
        return values_.size();
    }

    int dofs() const {
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
    DMUMPS_STRUC_C id{};

    int& icntl(int idx) {
        return id.icntl[idx - 1];
    }

    double& cntl(int idx) {
        return id.cntl[idx - 1];
    }

    int info(int idx) const {
        return id.info[idx - 1];
    }

    double rinfo(int idx) const {
        return id.rinfo[idx - 1];
    }

    int infog(int idx) const {
        return id.infog[idx - 1];
    }

    double rinfog(int idx) const {
        return id.rinfog[idx - 1];
    }

public:

    solver() {
        int argc = 0;
        char** argv = nullptr;
        MPI_Init(&argc, &argv);

        id.job = -1;
        id.par = 1;
        id.sym = 0;
        id.comm_fortran = MPI_Comm_c2f(MPI_COMM_SELF);

        dmumps_c(&id);

        // Use ordering scheme: auto (0), sequential (1), parallel (2)
        icntl(28) = 1;

        // Sequential orderings:
        // Metis (5), PORD (4), SCOTCH (3), AMD (0), AMF (2), QAMD (6), automatic choice (7)
        icntl(7) = 7;

        // Parallel orderings:
        // PT-SCOTCH (1), ParMetis (2), automatic choice (0)
        // icntl(29) = 2;

        // extra space factor
        icntl(14) = 20;

        // null pivot detection
        icntl(24) = 1;
        cntl(3) = 1e-5;

        // out-of-core
        // icntl(22) = 1;

        //streams
        icntl(1) = 3;
        icntl(2) = 3;
        icntl(3) = 3;
        icntl(4) = 3;
    }

    void print_state(const char* text) const {
        if (id.info[0] != 0) {
            std::cout << text << ":" << std::endl;
            std::cout << "  INFO(1) = " << id.info[0] << std::endl;
            std::cout << "  INFO(2) = " << id.info[1] << std::endl;
            std::cout << "Error: ";
            report_error(std::cout, id.info[0], id.info[1]);
            std::cout << std::endl;
        }
    }

    void save_to_file(problem& problem, const char* output_path) {
        prepare_(problem);
        strcpy(id.write_problem, output_path);

        analyze_();
    }

    void solve(problem& problem, const char* output_path = nullptr) {
        prepare_(problem);

        if (output_path) {
            strcpy(id.write_problem, output_path);
        }

        analyze_();
        // std::cout << "Analysis type: " << infog(32) << std::endl;
        // std::cout << "Ordering used: " << infog(7) << std::endl;
        // report_after_analysis(std::cout);
        factorize_();
        // std::cout << "Deficiency: " << infog(28) << std::endl;
        solve_();
    }

    double flops_assembly() const {
        return rinfog(2);
    }

    double flops_elimination() const {
        return rinfog(3);
    }

    ~solver() {
        id.job = -2;
        dmumps_c(&id);
        MPI_Finalize();
    }

private:
    void prepare_(problem& problem) {
        id.n = problem.dofs();
        id.nz = problem.nonzero_entries();

        id.irn = problem.irn();
        id.jcn = problem.jcn();
        id.a = problem.a();

        id.rhs = problem.rhs();
        id.nrhs = 1;
        id.lrhs = id.n;
    }

    void factorize_() {
        id.job = 2;
        dmumps_c(&id);
        print_state("After factorize_()");
    }

    void analyze_() {
        id.job = 1;
        dmumps_c(&id);
        print_state("After analyze_()");
    }

    void solve_() {
        id.job = 3;
        dmumps_c(&id);
        print_state("After solve_()");
    }

    void report_error(std::ostream& os, int info1, int info2) const {
        switch (info1) {
        case -6:
            os << "Matrix singular in structure (rank = " << info2 << ")";
            break;
        case -7:
            os << "Problem of integer workspace allocation (size = " << info2 << ")";
            break;
        case -8:
            os << "Main internal integer workarray is too small";
            break;
        case -9:
            os << "Main internal real/complex workarray is too small (missing = " << info2 << ")";
            break;
        case -10:
            os << "Numerically singular matrix";
            break;
        case -13:
            os << "Problem of workspace allocation (size = " << info2 << ")";
            break;
        }
    }

    long handle_neg(int value) const {
        if (value >= 0) {
            return value;
        } else {
            return -value * 1'000'000L;
        }
    }

    double as_MB(long size) const {
        return static_cast<double>(size) / (1024 * 1024) * sizeof(double);
    }

    void report_after_analysis(std::ostream& os) const {
        os << "MUMPS (" << id.version_number << ") after analysis:" << std::endl;
        os << "RINFO:" << std::endl;
        os << "  Estimated FLOPS for the elimination:          " << rinfo(1) << std::endl;
        os << "  Disk space for out-of-core factorization:     " << rinfo(5) << " MB" << std::endl;
        os << "  Size of the file used to save data:           " << rinfo(7) << " MB" << std::endl;
        os << "  Size of the MUMPS structure:                  " << rinfo(8) << " MB" << std::endl;
        os << "INFO:" << std::endl;
        os << "  Success:                                      " << info(1) << std::endl;
        auto real_store = handle_neg(info(3));
        os << "  Size of the real space to store factors:      " << real_store
                                                                 << " (" << as_MB(real_store) << " MB)" << std::endl;
        os << "  Size of the integer space to store factors:   " << info(4) << std::endl;
        os << "  Estimated maximum front size:                 " << info(5) << std::endl;
        os << "  Number of nodes in a tree:                    " << info(6) << std::endl;
        os << "  Size of the integer space to factorize:       " << info(7) << std::endl;
        auto real_factor = handle_neg(info(8));
        os << "  Size of the real space to factorize:          " << real_factor
                                                                 << " (" << as_MB(real_factor) << " MB)" << std::endl;
        os << "  Total memory needed:                          " << info(15) << " MB" << std::endl;
        os << "  Total memory needed (OoC):                    " << info(17) << " MB" << std::endl;
        os << "  Size of the integer space to factorize (OoC): " << info(19) << std::endl;
        auto real_factor_ooc = handle_neg(info(20));
        os << "  Size of the real space to factorize (OoC):    " << real_factor_ooc
                                                                 << " (" << as_MB(real_factor_ooc) << " MB)" << std::endl;
        os << "  Estimated number of entries in factors:       " << handle_neg(info(24)) << std::endl;
        auto low_real_factor = handle_neg(info(29));
        os << "  Size of the real space to factorize (low-r):  " << low_real_factor
                                                                 << " (" << as_MB(low_real_factor) << " MB)" << std::endl;
        os << "  Total memory needed (low-rank):               " << info(30) << " MB" << std::endl;
        os << "  Total memory needed (low-rank, OoC):          " << info(31) << " MB" << std::endl;
    }
};


}


#endif /* MUMPS_HPP_ */
