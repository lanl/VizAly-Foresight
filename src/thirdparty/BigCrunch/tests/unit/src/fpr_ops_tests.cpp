#include <fpr_ops_tests.hpp>

#include <omp.h>




void PrintTo(const fpr_ops_test_state &state, std::ostream *out)
{
    *out << "test data: [";
    for(int i = 0; i < state.size; ++i)
    {
        *out << state.data[i] << ", ";
    }
    *out << "]" << std::endl;
    *out << "test data size: " << state.size << std::endl;
    *out << "num_threads: " << state.num_threads << std::endl;
    *out << "expected abs_min: " << state.abs_min << std::endl;
    *out << "expected abs_max: " << state.abs_max << std::endl;
    *out << "tolerance setting: " << state.tolerance << std::endl;
    *out << "expected max_exp: " << state.max_exp << std::endl;
    *out << "expected offset: " << state.offset << std::endl;
}


TEST_P(fpr_ops_test, abs_minmax_test)
{
    auto state = GetParam();
    int num_threads = omp_get_num_threads();
    omp_set_num_threads(state.num_threads);
    auto output = bigcrunch::fpr_ops::abs_minmax(state.data, state.size);
    omp_set_num_threads(num_threads);
    ASSERT_EQ(state.abs_min, output.first);
    ASSERT_EQ(state.abs_max, output.second);
}


TEST_P(fpr_ops_test, exp_range_test)
{
    auto state = GetParam();
    auto output = bigcrunch::fpr_ops::exp_range(state.data, state.size, state.tolerance);
    ASSERT_EQ(state.max_exp, output.max_exp);
    ASSERT_EQ(state.offset, output.offset);
}


INSTANTIATE_TEST_CASE_P(Default, fpr_ops_test, testing::Values(
    abs_minmax_test_state{{-1, 5, -2, 3, -6, 4}, 1, 6, 1, 2, 0},
    abs_minmax_test_state{{-1, 5, -2, 3, -6, 4}, 1, 6, 2, 2, 0},
    abs_minmax_test_state{{-1, 5, -2, 3, -6, 4}, 1, 6, 4, 2, 0},
    abs_minmax_test_state{{-1, 5, -0.25, 3, -6, 4}, 0.25, 6, 1, 2, -1},
    abs_minmax_test_state{{-1, 5, -0.25, 3, -6, 4}, 0.25, 6, 2, 2, -1},
    abs_minmax_test_state{{-1, 5, -0.25, 3, -6, 4}, 0.25, 6, 4, 2, -1},
    abs_minmax_test_state{{-0.1, 0.5, -0.25, 0.3, -0.6, 0.4}, 0.1, 0.6, 1, 0, -1},
    abs_minmax_test_state{{-0.1, 0.5, -0.25, 0.3, -0.6, 0.4}, 0.1, 0.6, 2, 0, -1},
    abs_minmax_test_state{{-0.1, 0.5, -0.25, 0.3, -0.6, 0.4}, 0.1, 0.6, 4, 0, -1}
));