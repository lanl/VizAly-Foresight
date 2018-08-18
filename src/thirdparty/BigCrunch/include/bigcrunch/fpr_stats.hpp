/**
 * @file fpr_stats.hpp
 * @author Max Zeyen
 * @date 2018-05-02
 */

#pragma once

#include <cstdint>




namespace bigcrunch
{


    /**
     * @class fpr_stats
     * 
     * @brief Metadata container for FPR
     * 
     * This class computes properties derived from input/output metadata.
     * It converts the relative error to base 2, computes the new exponent bit size, computes the new 
     * mantissa bit size, and determines the total bit size per value.
     */
    class fpr_stats
    {
    public:
        double rel_error;
        int error, max_exp, offset;
        std::uint64_t exp_bits, mant_bits, value_bits, aligned_chunk_size;

        /**
         * @brief Construct a new fpr stats object
         * 
         * The constructor calls all the compute methods and stores the results in the appropriate 
         * class member variables.
         * 
         * @param error Relative error exponent in base 10
         * @param max_exp Maximum exponent of a data range
         * @param offset Minimum exponent of a data range
         */
        fpr_stats(int error, int max_exp, int offset);

    private:
        /**
         * @brief Compute the relative error in base 2
         * 
         * Converts the relative error to base 2, such that the resulting error is smaller equal to 
         * the corresponding base 10 relative error.
         */
        void relative_error();

        /**
         * @brief Compute the new exponent bit size
         * 
         * Compute the number of bits required for the exponents calculated by FPR.
         */
        void exponent_bits();

        /**
         * @brief Compute the new mantissa bit size
         * 
         * Compute the number of bits required for the mantissas calculated by FPR.
         */
        void mantissa_bits();

        /**
         * @brief Align chunk size
         * 
         * Compute the aligned size of value chunks.
         * This is necessary as most chunk sizes in FPR outputs are misaligned, which causes issues with 
         * Blosc.
         */
        void align_chunk_size();
    };


}