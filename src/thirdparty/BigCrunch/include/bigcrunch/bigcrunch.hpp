/**
 * @file bigcrunch.hpp
 * @author Max Zeyen
 * @date 2018-04-19
 */

#pragma once

#include <cstdint>

#include <bigcrunch/types.hpp>
#include <bigcrunch/darray.hpp>




namespace bigcrunch
{


    /**
     * @class bigcrunch
     * 
     * @brief Main class for BigCrunch compressor library
     * 
     * This class defines the main access point for using BigCrunch's functionality.
     * It provides a simple intuitive interface to compress and decompress data.
     */
    class bigcrunch
    {
    public:
        /**
         * @brief Construct a new bigcrunch object
         * 
         * Constructs a new bigcrunch object using an unordered_map<int, int> object (= setting_t) for passing configuration settings.
         * The key values for the different settings can be found in the enum \ref config_t.
         * 
         * @param settings An unordered_map<int, int> for settings.
         */
        bigcrunch(const setting_t &settings);

        /**
         * @brief Compress data
         * 
         * This method is the entry point to BigCrunch's compression pipeline.
         * It takes a one-dimensional floating-point data array as input and returns a compressed byte array with its size.
         * 
         * @param data One-dimensional floating-point data array
         * @param output Compressed one-dimensional byte array
         * @return std::uint64_t Size of the byte array
         */
        std::uint64_t compress(const darray &data, std::uint8_t **output) const;

        /**
         * @brief Decompress data
         * 
         * This method is the entry point to BigCrunch's decompression pipeline.
         * It takes a one-dimensional byte array created by the compression method and its size as inputs and returns a darray with the reconstructed data.
         * 
         * @param data Compressed one-dimensional byte array
         * @param size Size of the byte array
         * @return darray Reconstructed floating-point data
         */
        darray decompress(const std::uint8_t *data, std::uint64_t size) const;

    protected:
        int error, tolerance, nthreads;
        blosc_config_t blosc_config;
    };


}