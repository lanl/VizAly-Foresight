#pragma once

#include <typeinfo>
#include <cmath>

#include <bigcrunch/fpr_ops.hpp>
#include <bigcrunch/bit_ops.hpp>

#include <omp.h>




namespace bigcrunch
{


    template<typename T> inline
    fpr_encoder<T>::fpr_encoder(int error, int tolerance)
        : error(error),
        tolerance(tolerance)
    {
    }


    template<typename T> inline
    fpr_data fpr_encoder<T>::encode(const T *data, std::uint64_t size) const
    {
        // init result object
        fpr_data result;
        result.tolerance = this->tolerance;
        result.error = this->error;

        // set type id
        if(typeid(T) == typeid(float))
        {
            result.type = 0;
        }
        else
        {
            result.type = 1;
        }

        // compute maximum exponent and offset
        exp_range_t range = fpr_ops::exp_range(data, size, this->tolerance);
        result.max_exp = range.max_exp;
        result.offset = range.offset;

        // compute fpr stats
        fpr_stats stats(this->error, result.max_exp, result.offset);
        // set type size for blosc compressor (for chunk size 8 this corresponds to value_bits);
        result.type_size = stats.value_bits;

        // allocate output array
        result.num_elems = size;
        std::uint64_t data_size = this->compute_new_size(size, stats);
        result.data = new std::uint8_t[data_size]();

        // recode floating-point values
        #pragma omp parallel for schedule(static, 8)
#ifdef WIN32
		for(long long i=0; i < size; ++i)
#else
        for(std::uint64_t i = 0; i < size; ++i)
#endif
        {
            // set sign bit
            std::uint8_t sign = std::signbit(data[i]);
#ifdef WIN32
            this->write_bits({ sign, (unsigned long long)i, 1, 0, stats.aligned_chunk_size }, result.data);
#else
            this->write_bits({value:sign, value_idx:i, num_bits:1, offset:0, aligned_chunk_size:stats.aligned_chunk_size}, result.data);
#endif

            // set exponent bits
            int exponent = this->compute_exponent(data[i], stats);
#ifdef WIN32
            this->write_bits({ static_cast<std::uint64_t>(exponent), (unsigned long long)i, stats.exp_bits, 1, stats.aligned_chunk_size }, result.data);
#else
            this->write_bits({value:static_cast<std::uint64_t>(exponent), value_idx:i, num_bits:stats.exp_bits, offset:1, aligned_chunk_size:stats.aligned_chunk_size}, result.data);
#endif
            // set mantissa bits
            std::uint64_t mantissa = this->compute_mantissa(data[i], exponent, stats);
#ifdef WIN32
            this->write_bits({mantissa, (unsigned long long)i, stats.mant_bits, 1 + stats.exp_bits, stats.aligned_chunk_size}, result.data);
#else
            this->write_bits({ value:mantissa, value_idx : i, num_bits : stats.mant_bits, offset : 1 + stats.exp_bits, aligned_chunk_size : stats.aligned_chunk_size }, result.data);
#endif
        }

        return result;
    }


    template<typename T> inline
    std::uint64_t fpr_encoder<T>::compute_new_size(std::uint64_t num_elems, const fpr_stats &stats) const
    {
        // return aligned total new array size
        return static_cast<std::uint64_t>(std::ceil(num_elems / 8.0) * stats.aligned_chunk_size);
    }


    template<typename T> inline
    int fpr_encoder<T>::compute_exponent(const T &data_point, const fpr_stats &stats) const
    {
        int result = 0;
        // for large values really close to the interval boundary, the new exponent can jump to the next interval 
        // due to floating-point rounding errors (e.g. 4095.999999f yields the same exponent as 4096.0f).
        // the if-clause below corrects this rounding error by flooring the value before calculating the exponent
        // for values larger 1.
        if(std::abs(data_point) >= 1)
        {
            result = static_cast<int>(std::floor(std::log2(std::floor(std::abs(data_point)))) - stats.offset);
        }
        else
        {
            result = static_cast<int>(std::floor(std::log2(std::abs(data_point))) - stats.offset);
        }

        // catches values smaller than the tolerance and sets the exponent to 0
        if(result < 0)
        {
            result = 0;
        }

        return result;
    }


    template<typename T> inline
    std::uint64_t fpr_encoder<T>::compute_mantissa(const T &data_point, int exponent, const fpr_stats &stats) const
    {
        // compute bin width for new mantissa
        double abs_error = 2.0 * stats.rel_error;
        int no_offset_exp = exponent + stats.offset;

        if(no_offset_exp >= this->tolerance)
        {
            abs_error *= std::pow(2, no_offset_exp);
        }

        // compute bin id of new mantissa
        std::uint64_t result = std::floor(std::abs(data_point) / abs_error);
        if(no_offset_exp >= this->tolerance)
        {
            result -= std::pow(2, no_offset_exp) / abs_error;
        }

        return result;
    }


    template<typename T> inline
    void fpr_encoder<T>::write_bits(const value_t &value, std::uint8_t *output) const
    {        
        // get bit index (same for all bytes within a chunk)
        std::uint8_t bit_idx = bit_ops::bit_idx(value.value_idx);
        // get chunk index (index of first byte in chunk)
        std::uint64_t chunk_idx = bit_ops::chunk_idx(value.value_idx, value.aligned_chunk_size); 

        // iterate over bits in value and write them in output byte array
        for(std::uint64_t j = 0; j < value.num_bits; ++j)
        {
            // extract jth bit
            bool bit = (value.value >> j) & 1;
            // calculate byte index
            std::uint64_t byte_idx = chunk_idx + value.offset + j;
            // set jth bit in output array
            bit_ops::update_bit(bit, bit_idx, output[byte_idx]);
        }
    }



}