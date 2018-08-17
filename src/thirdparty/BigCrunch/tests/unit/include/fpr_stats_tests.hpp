#pragma once

#include <ostream>
#include <cstdint>

#include <bigcrunch/fpr_stats.hpp>

#include <gtest/gtest.h>




struct fpr_stats_test_state
{
    int error, max_exp, offset;
    std::uint64_t exp_bits, mant_bits;
};


void PrintTo(const fpr_stats_test_state &state, std::ostream *out);


class fpr_stats_test : public testing::TestWithParam<fpr_stats_test_state>
{
public:
    fpr_stats *stats;

    fpr_stats_test();
    ~fpr_stats_test();
};