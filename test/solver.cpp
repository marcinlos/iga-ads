#include <iostream>
#include "ads/lin/band_matrix.hpp"
#include "ads/lin/band_solve.hpp"

using namespace ads::lin;
using namespace std;

int main() {
    int kl = 1;
    int ku = 2;
    int n = 6;
    int d = 4;

    band_matrix m(kl, ku, n);
    matrix b(n, d);

    for (int i = 0; i < n; ++ i) {
        for (int j = max(0, i - kl); j < min(n, i + ku + 1); ++ j) {
            m(i, j) = (i + 1) * 10 + j + 1;
        }
        for (int j = 0; j < d; ++ j) {
            b(i, j) = (j + 1) * (i + 1);
        }
    }

    double solution[] = { 0.230377, -0.126052, -0.001655, -0.001112, 0.203603, -0.109609 };
    matrix x(n, d);
    for (int i = 0; i < n; ++ i) {
        for (int j = 0; j < d; ++ j) {
            x(i, j) = (j + 1) * solution[i];
        }
    }

    solver_ctx ctx(m);

    solve(m, b, ctx);

    cout << "Solution:" << endl << b << endl;
    cout << "Expected:" << endl << x << endl;

    cout << "Equal: " << std::boolalpha << equal(b, x, 1e-5);
}




