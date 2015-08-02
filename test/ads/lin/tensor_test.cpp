#include <libunittest/all.hpp>
#include "ads/lin/tensor.hpp"

using namespace unittest::assertions;

namespace ads {
namespace lin {

struct tensor_test: unittest::testcase<> {

    static void run() {
        UNITTEST_CLASS(tensor_test)
        UNITTEST_RUN(test_tensor_basics)
        UNITTEST_RUN(test_equal_tensors_should_be_equal)
        UNITTEST_RUN(test_cyclic_transpose)
    }

    void test_tensor_basics() {
        tensor<double, 2> t { 5, 3 };
        double x = 2;
        t(2, 1) = x;
        int idx = 5 + 2;
        assert_equal(x, t.data()[idx]);
    }

    void test_equal_tensors_should_be_equal() {
        int p = 5, q = 3;
        tensor<double, 2> a { p, q };
        tensor<double, 2> b { 5, 3 };

        for (int i = 0; i < p; ++ i) {
            for (int j = 0; j < q; ++ j) {
                a(i, j) = b(i, j) = 7;
            }
        }
        assert_true(a == b);
    }

    void test_reshape() {
        int p = 5, q = 3;

        double data[] = {
            1, 2, 3, 4, 5,
            6, 7, 8, 9, 10
        };
        tensor_view<double, 2> tensor2d { data, p, q };
        tensor_view<double, 1> tensor1d { data, p * q };
        tensor_view<double, 1> reshaped = reshape(tensor2d, p * q);
        assert_equal(tensor1d, reshaped);
    }

    void test_cyclic_transpose() {
        int k = 2, n = 3, m = 2;

        tensor<double, 3> a { k, n, m };
        tensor<double, 3> e { n, m, k };

        a(0, 0, 0) = e(0, 0, 0) = 111;
        a(1, 0, 0) = e(0, 0, 1) = 211;
        a(0, 1, 0) = e(1, 0, 0) = 121;
        a(1, 1, 0) = e(1, 0, 1) = 221;
        a(0, 2, 0) = e(2, 0, 0) = 131;
        a(1, 2, 0) = e(2, 0, 1) = 231;

        a(0, 0, 1) = e(0, 1, 0) = 112;
        a(1, 0, 1) = e(0, 1, 1) = 212;
        a(0, 1, 1) = e(1, 1, 0) = 122;
        a(1, 1, 1) = e(1, 1, 1) = 222;
        a(0, 2, 1) = e(2, 1, 0) = 132;
        a(1, 2, 1) = e(2, 1, 1) = 232;

        tensor<double, 3> out { n, m, k };
        cyclic_transpose(a, out);

        assert_true(out == e);

        tensor<double, 3> a2 { m, k, n }, a3 { k, n, m };
        cyclic_transpose(out, a2);
        cyclic_transpose(a2, a3);

        assert_true(a3 == a);
    }

};

REGISTER(tensor_test)

}
}
