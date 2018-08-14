#include <gtest/gtest.h>

#include <fpr_ops_tests.hpp>
#include <fpr_stats_tests.hpp>
#include <fpr_encoder_tests.hpp>
#include <fpr_decoder_tests.hpp>
#include <bit_ops_tests.hpp>




int main(int argc, char **argv)
{
    testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}