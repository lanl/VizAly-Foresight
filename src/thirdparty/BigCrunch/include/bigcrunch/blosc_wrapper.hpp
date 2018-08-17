/**
 * @file blosc_wrapper.hpp
 * @author Max Zeyen
 * @date 2018-04-19
 */

#pragma once

#include <string>
#include <cstdint>

#include <bigcrunch/types.hpp>
#include <bigcrunch/fpr_data.hpp>




namespace bigcrunch
{


    /**
     * @class blosc_wrapper
     * 
     * @brief C++ wrapper for Blosc
     * 
     * This class defines a C++ wrapper for the Blosc framework written in C.
     * The wrapper pre-configures blosc for use within BigCrunch.
     * It also acts as a convenience layer by abstracting some of Blosc's complexity.
     */
    class blosc_wrapper
    {
    public:
        /**
         * @brief Construct a new blosc_wrapper object
         * 
         * Initializes the blosc compression framework and sets the number of threads to be used.
         * 
         * @param nthreads Number of threads used by blosc
         */
        blosc_wrapper(int nthreads=1);

        /**
         * @brief Destroy the blosc wrapper object
         * 
         * Destroys/deinitializes the blosc framework.
         */
        ~blosc_wrapper();

        /**
         * @brief Compress a byte array
         * 
         * Compresses a given byte array contained in an \ref fpr_data object using the blosc framework.
         * The result is output in the same \ref fpr_data object.
         * The default configuration of blosc can be changed by passing a \ref blosc_config_t object to this method.
         * 
         * @param data \ref fpr_data object with the byte array to be compressed
         * @param config Blosc configuration parameters
         */
        void compress(fpr_data &data, const blosc_config_t &config={blosc_compressor_t::ZSTD, 9, blosc_filter_t::SHUFFLE});

        /**
         * @brief Decompress a byte array
         * 
         * Decompress a given byte array contained in an \ref fpr_data object using the blosc framework.
         * The result is output in the same \ref fpr_data object.
         * 
         * @param data \ref fpr_data object with the byte array to be compressed
         */
        void decompress(fpr_data &data);
    };


}