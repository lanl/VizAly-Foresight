/**
 * @file fpr_decoder.hpp
 * @author Max Zeyen
 * @date 2018-04-19
 */

#pragma once

#include <cstdint>

#include <bigcrunch/fpr_data.hpp>
#include <bigcrunch/darray.hpp>




namespace bigcrunch
{


    /**
     * @class fpr_decoder
     * 
     * @brief Decoder for FPR encoded data
     * 
     * This class defines a decoder for FPR encoded data.
     * FPR encoding and decoding is limited to FLOAT32 (i.e. float) and FLOAT64 (i.e. double) data types.
     * 
     * @tparam T Type of the original data
     */
    template<typename T> class fpr_decoder
    {
    public:
        /**
         * @brief Construct a new fpr_decoder object
         * 
         * Default constructor creating an empty fpr_decoder object.
         */
        fpr_decoder() = default;

        /**
         * @brief Decode FPR encoded data
         * 
         * Decodes FPR encoded data and returns it as a darray.
         * 
         * @param data Data to decode
         * @return darray Decoded data array
         */
        darray decode(const fpr_data &data) const;

    protected:
        /**
         * @brief Reconstruct FPR encoded value
         * 
         * Reconstruct a FPR encoded value from the extracted sign, exponent, and mantissa bits.
         * 
         * @param data Meta data (e.g. tolerance)
         * @param stats Encoding meta data (e.g. relative error)
         * @param sign Extracted sign value
         * @param exponent Extracted exponent value
         * @param mantissa Extracted mantissa value
         * @return T Reconstructed value
         */
        T reconstruct_value(const fpr_data &data, const fpr_stats &stats, std::uint8_t sign, int exponent, std::uint64_t mantissa) const;

        /**
         * @brief Read bits from byte array
         * 
         * Reads a certain number of bits from a byte array from a given offset.
         * The number of bits, the offset and other meta data are defined in \ref value_t.
         * The result is returned as a 64 bit unsigned integer type.
         * 
         * @param value Meta information necessary to extract the correct bits
         * @param data Byte array to read bits from
         * @return std::uint64_t Extracted bits
         */
        std::uint64_t read_bits(const value_t &value, const std::uint8_t *data) const;
    };


}


#include <bigcrunch/detail/fpr_decoder.inl>