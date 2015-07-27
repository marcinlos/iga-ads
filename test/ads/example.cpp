#include <libunittest/all.hpp>

TEST(test_value_is_true)
{
    assert_true(true, SPOT);
}

struct fixture
{
    int value;
    fixture(): value(42)
    { }
};

TEST_FIXTURE(fixture, test_with_fixture)
{
    assert_equal(42, value, SPOT);
}
