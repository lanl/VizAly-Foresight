#pragma once

#include <cmath>

#include <bigcrunch/fpr_stats.hpp>
#include <bigcrunch/bit_ops.hpp>

#include <omp.h>




namespace bigcrunch
{


    template<typename T> inline
    darray fpr_decoder<T>::decode(const fpr_data &data) const
    {
        // compute fpr stats
        fpr_stats stats(data.error, data.max_exp, data.offset);

        // allocate output array
        T *result = new T[data.num_elems];

        // reconstruct values
        #pragma omp parallel for schedule(static, 8)
        for(std::uint64_t i = 0; i < data.num_elems; ++i)
        {
            // get sign bit
            value_t params = {value:0, value_idx:i, num_bits:1, offset:0, aligned_chunk_size:stats.aligned_chunk_size};
            std::uint64_t sign = this->read_bits(params, data.data);

            // get exponent bits
            params.num_bits = stats.exp_bits;
            params.offset += 1;
            int exponent = static_cast<int>(this->read_bits(params, data.data));
            exponent += data.offset;

            // get mantissa bits
            params.num_bits = stats.mant_bits;
            params.offset += stats.exp_bits;
            std::uint64_t mantissa = this->read_bits(params, data.data);

            // reconstruct original value
            result[i] = this->reconstruct_value(data, stats, sign, exponent, mantissa);
        }

        return darray(result, data.num_elems);
    }


    template<typename T> inline
    T fpr_decoder<T>::reconstruct_value(const fpr_data &data, const fpr_stats &stats, std::uint8_t sign, int exponent, std::uint64_t mantissa) const
    {
        // define output
        T value;

        // compute absolute error used to encode this value
        double abs_error = 2.0 * stats.rel_error;
        if(exponent >= data.tolerance)
        {
            abs_error *= std::pow(2, exponent);
        }

        // reconstruct mantissa using the previously determined absolute error
        value = static_cast<T>(static_cast<T>(mantissa) * abs_error + 0.5 * abs_error);
        // reconstruct the value's offset and add it to the value
        if(exponent >= data.tolerance)
        {
            value += static_cast<T>(std::pow(2, exponent));
        }

        // reconstruct the value's sign
        value *= std::pow(-1, sign);

        return value;
    }


    template<typename T> inline
    std::uint64_t fpr_decoder<T>::read_bits(const value_t &value, const std::uint8_t *data) const
    {
        // get bit index (same for all bytes within a chunk)
        std::uint8_t bit_idx = bit_ops::bit_idx(value.value_idx);
        // get chunk index (index of first byte in chunk)
        std::uint64_t chunk_idx = bit_ops::chunk_idx(value.value_idx, value.aligned_chunk_size);
        // init result
        std::uint64_t result = 0;
        // convert result variable into byte array
        std::uint8_t *result_bytes = reinterpret_cast<std::uint8_t *>(&result);

        // iterate over the bits and write them to the output variable
        for(std::uint64_t j = 0; j < value.num_bits; ++j)
        {
            // calculate byte index
            std::uint64_t byte_idx = chunk_idx + value.offset + j;
            // read jth bit
            bool bit = bit_ops::get_bit(bit_idx, data[byte_idx]);
            // set jth bit in output variable
            bit_ops::update_bit(bit, j % 8, result_bytes[j / 8]);
        }

        return result;
    }


}