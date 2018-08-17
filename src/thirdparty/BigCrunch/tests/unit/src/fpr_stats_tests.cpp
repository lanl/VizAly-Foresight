#include <fpr_stats_tests.hpp>

#include <cmath>




void PrintTo(const fpr_stats_test_state &state, std::ostream *out)
{
    *out << "error: " << state.error << std::endl;
    *out << "max_exp: " << state.max_exp << std::endl;
    *out << "offset: " << state.offset << std::endl;
    *out << "expected exp_bits: " << state.exp_bits << std::endl;
    *out << "expected mant_bits: " << state.mant_bits << std::endl;
}


fpr_stats_test::fpr_stats_test()
{
    auto state = GetParam();
    this->stats = new fpr_stats(state.error, state.max_exp, state.offset);
}


fpr_stats_test::~fpr_stats_test()
{
    if(this->stats != nullptr)
    {
        delete this->stats;
        this->stats = nullptr;
    }
}


TEST_P(fpr_stats_test, relative_error_test)
{
    auto state = GetParam();
    ASSERT_TRUE(stats->relative_error < std::pow(10, state.error));
}


TEST_P(fpr_stats_test, exp_bits_test)
{
    auto state = GetParam();
    ASSERT_EQ(state.exp_bits, stats->exp_bits);
}


TEST_P(fpr_stats_test, mant_bits_test)
{
    auto state = GetParam();
    ASSERT_EQ(state.mant_bits, stats->mant_bits);
}


TEST_P(fpr_stats_test, value_bits_test)
{
    auto state = GetParam();
    ASSERT_EQ(1 + state.exp_bits + state.mant_bits, stats->value_bits);
}


INSTANTIATE_TEST_CASE_P(Default, fpr_stats_test, testing::Values(
    fpr_stats_test_state{-3, 2, 0, 2, 9},
    fpr_stats_test_state{-3, 2, -1, 2, 9},
    fpr_stats_test_state{-3, 0, -1, 1, 9},
    fpr_stats_test_state{-4, 2, 0, 2, 13},
    fpr_stats_test_state{-4, 2, -1, 2, 13},
    fpr_stats_test_state{-4, 0, -1, 1, 13},
    fpr_stats_test_state{-5, 2, 0, 2, 16},
    fpr_stats_test_state{-5, 2, -1, 2, 16},
    fpr_stats_test_state{-5, 0, -1, 1, 16},
));