/**
 * @file types.hpp
 * @author Max Zeyen
 * @date 2018-05-02
 */

#pragma once

#include <cstdint>
#include <unordered_map>




namespace bigcrunch
{


    /**
     * @brief Supported darray types
     * 
     * This enum contains the supported data types of darray.
     */
    enum dtype_t
    {
        UINT8,
        UINT16,
        UINT32,
        UINT64,
        INT8,
        INT16,
        INT32,
        INT64,
        FLOAT32,
        FLOAT64,
        NONE
    };


    /**
     * @brief Output type of exp_range method
     * 
     * This aggregate structure is used to make the output of the exp_range method more understandable.
     */
    struct exp_range_t
    {
        int max_exp, offset;
    };


    /**
     * @brief Convenience aggregate used by the write_bits method
     * 
     * This aggregate structure is used to group the input of the write_bits method together.
     * It contains the value to be written, the bit/byte index, the number of bits to write, and the chunk offset.
     */
    struct value_t
    {
        std::uint64_t value, value_idx, num_bits, offset, aligned_chunk_size;
    };


    /**
     * @brief Supported blosc filters
     * 
     * This enum contains the supported blosc data filters.
     */
    enum blosc_filter_t
    {
        NOSHUFFLE,
        SHUFFLE,
        BITSHUFFLE
    };


    /**
     * @brief Supported blosc compressors
     * 
     * This enum contains the supported blosc compressors.
     */
    enum blosc_compressor_t
    {
        BLOSCLZ,
        LZ4,
        LZ4HC,
        SNAPPY,
        ZLIB,
        ZSTD
    };


    /**
     * @brief Convenience aggregate used by the blosc wrapper
     * 
     * This aggregate structure is used to group the input of the blosc_wrapper together.
     * It contains the compressor type, compression level, and compression filter.
     */
    struct blosc_config_t
    {
        blosc_compressor_t cname;
        int clevel;
        blosc_filter_t filter;
    };


    /**
     * @brief Endianess type
     * 
     * This enum contains the two possible endianess types of computer systems.
     */
    enum endianess_t
    {
        LITTLE,
        BIG
    };


    /**
     * @brief Supported BigCrunch configuration properties
     * 
     * This enum contains the supported BigCrunch configuration properties used to initialize 
     * the BigCrunch compression pipeline.
     */
    enum config_t
    {
        ERR,
        TOLERANCE,
        BLOSC_NTHREADS,
        BLOSC_COMPRESSOR,
        BLOSC_CLEVEL,
        BLOSC_FILTER
    };


    /**
     * @brief BigCrunch configuration convenience type
     * 
     * This is a convenience type definition for BigCrunch's configuration properties.
     */
    using setting_t = std::unordered_map<int, int>;


}