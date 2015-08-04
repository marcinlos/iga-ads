#include <libunittest/all.hpp>
#include "ads/util/multi_array.hpp"

using namespace unittest::assertions;
using namespace ads;

struct multi_array_test: unittest::testcase<> {

    static void run() {
        UNITTEST_CLASS(multi_array_test)
        UNITTEST_RUN(test_wrapper_indexing_2d)
        UNITTEST_RUN(test_simple_reshape)
        UNITTEST_RUN(test_constant_access)
        UNITTEST_RUN(test_size_query)
    }

    void test_wrapper_indexing_2d() {
        constexpr int n = 2, m = 5;
        int buffer[n * m] = {0};
        ads::multi_array_wrapper<int, 2, int*> a(buffer, n, m);
        for (int i = 0; i < n; ++ i) {
            for (int j = 0; j < m; ++ j) {
                a(i, j) = 10 * (i + 1) + j + 1;
            }
        }
        int expected[] = {11, 12, 13, 14, 15, 21, 22, 23, 24, 25};
        assert_equal_containers(expected, buffer);
    }

    void test_simple_reshape() {
        constexpr int n = 2, m = 5;
        int buffer[n * m] = {0};
        ads::multi_array_wrapper<int, 1, int*> a(buffer, n * m);
        a(7) = 7;
        ads::multi_array_wrapper<int, 2, int*> a2 = reshape<2>(a, n, m);
        assert_equal(7, a2(1, 2));
    }

    void test_constant_access() {
        constexpr int n = 2, m = 5;
        int buffer[n * m] = {0};
        ads::multi_array_wrapper<int, 2, int*> a(buffer, n, m);
        const auto& ca = a;
        ca(0, 1);
        // ca(0, 1) = 3; // doesn't compile
    }

    void test_size_query() {
        constexpr int p = 2, q = 5, r = 3, s = 2;
        int buffer[p * q * r * s] = {0};
        ads::multi_array_wrapper<int, 4, int*> a(buffer, p, q, r, s);

        assert_equal(p, a.size(0));
        assert_equal(q, a.size(1));
        assert_equal(r, a.size(2));
        assert_equal(s, a.size(3));
    }

};

REGISTER(multi_array_test)
