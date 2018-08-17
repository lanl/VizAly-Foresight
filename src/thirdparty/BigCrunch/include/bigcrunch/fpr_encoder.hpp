/**
 * @file fpr_encoder.hpp
 * @author Max Zeyen
 * @date 2018-05-02
 */

#pragma once

#include <cstdint>

#include <bigcrunch/fpr_data.hpp>
#include <bigcrunch/fpr_stats.hpp>




namespace bigcrunch
{


    /**
     * @class fpr_encoder
     * 
     * @brief Encoder for the FPR format
     * 
     * This class defines an encoder for the FPR format.
     * FPR encoding and decoding is limited to FLOAT32 (i.e. float) and FLOAT64 (i.e. double) data types.
     * 
     * @tparam T Type of the original data
     */
    template<typename T> class fpr_encoder
    {
    public:
        /**
         * @brief Construct a new fpr encoder object
         * 
         * The constructor sets the relative error and tolerance exponents.
         * The error exponent is in base 10 and the tolerance exponent is in base 2.
         * Internally the class methods convert the error exponent into a base 2 representation 
         * which is smaller equal to the base 10 value.
         * The tolerance exponent describes the point where the error modeling switches from relative 
         * error to absolute error (i.e. use the absolute error if value is smaller than the tolerance 
         * else use the relative error).
         * 
         * @param error Relative error exponent in base 10
         * @param tolerance Tolerance exponent in base 2
         */
        fpr_encoder(int error, int tolerance);

        /**
         * @brief Encode data into the FPR format.
         * 
         * Encodes data into the FPR format using the error and tolerance exponents.
         * Based on these criteria, the bits per value are reduced while maintaining the given 
         * target relative error.
         * The function returns an object containing the encoded data and metadata required to 
         * decode the data later on.
         * 
         * @param data Data array to encode
         * @param size Size of the data array
         * @return fpr_data Encoded data and metadata
         */
        fpr_data encode(const T *data, std::uint64_t size) const;

    protected:
        int error, tolerance;

        /**
         * @brief Compute size of FPR encoded data
         * 
         * The output size of FPR is computed with respect to the chunk size of the encoding loop.
         * This means, that if necessary the output size is rounded up to always have full chunks by 
         * appending zero bytes to the last chunk.
         * 
         * @param num_elems Number of values in the original data
         * @param stats Metadata object
         * @return std::uint64_t Output size in bytes
         */
        std::uint64_t compute_new_size(std::uint64_t num_elems, const fpr_stats &stats) const;

        /**
         * @brief Compute new exponent value
         * 
         * Computes the new exponent value under consideration of the offset and tolerance.
         * New exponents smaller than the tolerance are encoded as 0.
         * Values larger 1 are floored prior to the exponent computation to prevent floating-point 
         * rounding errors.
         * 
         * @param data_point Input value
         * @param stats Metadata object
         * @return int New exponent value
         */
        int compute_exponent(const T &data_point, const fpr_stats &stats) const;

        /**
         * @brief Compute new mantissa value
         * 
         * Compute the new mantissa value with respect to the relative error target.
         * This process is very similar to computing a 1-dimensional histogram.
         * The bin width is computed using the new exponent, tolerance, and the relative error target.
         * 
         * @param data_point Input value
         * @param exponent New exponent of input value
         * @param stats Metadata object
         * @return std::uint64_t New mantissa value
         */
        std::uint64_t compute_mantissa(const T &data_point, int exponent, const fpr_stats &stats) const;

        /**
         * @brief Write bits into output buffer
         * 
         * Generic method to write bits of the new exponent and mantissa values into the output buffer.
         * 
         * @param value Structure with value, bit/chunk index, number of bits, and chunk offset
         * @param output Output byte buffer
         */
        void write_bits(const value_t &value, std::uint8_t *output) const;
    };


}


#include <bigcrunch/detail/fpr_encoder.inl>