/**
 * @file fpr_ops.hpp
 * @author Max Zeyen
 * @date 2018-05-02
 */

#pragma once

#include <utility>
#include <cstdint>

#include <bigcrunch/types.hpp>




namespace bigcrunch
{


    /**
     * @class fpr_ops
     * 
     * @brief FPR specific operations
     * 
     * This class defines FPR specific operations for computing minmax of the absolute values of a data range and 
     * for computing the maximum and minimum exponents with respect to the tolerance.
     */
    class fpr_ops
    {
    public:
        /**
         * @brief Deleted default constructor
         * 
         * The default constructor is deleted, as the class only contains static methods and therefore does not 
         * require an instance of the class.
         */
        fpr_ops() = delete;

        /**
         * @brief Compute minmax of absolute values
         * 
         * Computes minmax of the absolute values of a data range.
         * The returned result is structured the same way as for the STL minmax method, making the first element 
         * in the pair the minimum and the second the maximum.
         * 
         * @tparam T Type of data range
         * @param data Input data
         * @param size Size of data
         * @return std::pair<T, T> Minmax of absolute values of data range
         */
        template<typename T> static std::pair<T, T> abs_minmax(const T *data, std::uint64_t size);

        /**
         * @brief Compute data exponent range
         * 
         * Computes the minimum and maximum exponents of a data range with respect to the tolerance exponent.
         * The minimum exponent or offset becomes tolerance - 1 if it is smaller than the tolerance.
         * This causes new exponents computed by FPR smaller than the tolerance to be encoded as 0.
         * The maximum exponent becomes 0 if it is smaller than the tolerance.
         * 
         * @tparam T Type of data range
         * @param data Input data
         * @param size Size of data
         * @param tolerance Tolerance exponent
         * @return exp_range_t Aggregate containing offset and maximum exponent
         */
        template<typename T> static exp_range_t exp_range(const T *data, std::uint64_t size, int tolerance);
    };


}


#include <bigcrunch/detail/fpr_ops.inl>