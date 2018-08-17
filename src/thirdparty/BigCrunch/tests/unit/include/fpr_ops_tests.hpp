#pragma once

#include <ostream>

#include <bigcrunch/fpr_ops.hpp>

#include <gtest/gtest.h>




struct fpr_ops_test_state
{
    float data[6], abs_min, abs_max;
    int num_threads, max_exp, offset, size=6, tolerance=0;
};


void PrintTo(const fpr_ops_test_state *state, std::ostream *out);


class fpr_ops_test : public testing::TestWithParam<fpr_ops_test_state>
{
public:
    fpr_ops_test() = default;
    ~fpr_ops_test() = default;
};