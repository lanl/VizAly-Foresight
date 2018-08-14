#pragma once

#include <cmath>
#include <type_traits>
#include <limits>
#include <algorithm>

#include <omp.h>




namespace bigcrunch
{


    template<typename T> inline
    std::pair<T, T> fpr_ops::abs_minmax(const T *data, std::uint64_t size)
    {
        T min_result = std::numeric_limits<T>::max();
        T max_result = 0;

        #pragma omp parallel
        {
            T min_val = std::numeric_limits<T>::max();
            T max_val = 0;

            #pragma omp for nowait
            for(std::uint64_t i = 0; i < size; ++i)
            {
                min_val = std::min(std::abs(data[i]), min_val);
                max_val = std::max(std::abs(data[i]), max_val);
            }

            #pragma omp critical
            {
                min_result = std::min(min_result, min_val);
                max_result = std::max(max_result, max_val);
            }
        }

        return std::make_pair(min_result, max_result);
    }


    template<typename T> inline
    exp_range_t fpr_ops::exp_range(const T *data, std::uint64_t size, int tolerance)
    {
        // this function can only be used with floating-point data types
        static_assert(std::is_floating_point<T>::value, "FPR does only support floating-point data!");

        // init result object
        exp_range_t result;

        // determine absolute value range
        std::pair<T, T> bounds = fpr_ops::abs_minmax(data, size);
        // compute maximum exponent
        result.max_exp = static_cast<int>(std::floor(std::log2(bounds.second)));
        // compute minimum exponent aka offset
        result.offset = static_cast<int>(std::floor(std::log2(bounds.first)));

        // check for special cases
        if(result.max_exp < tolerance)
        {
            result.max_exp = 0;
        }

        // readjust offset if offset is smaller than tolerance
        // important if minimal value is 0.0 and the offset becomes -inf.
        if(result.offset < tolerance)
        {
            result.offset = tolerance - 1;
        }

        return result;
    }


}