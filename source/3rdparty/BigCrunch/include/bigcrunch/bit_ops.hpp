/**
 * @file bit_ops.hpp
 * @author Max Zeyen
 * @date 2018-04-19
 */

#pragma once

#include <cstdint>




namespace bigcrunch
{


    /**
     * @class bit_ops
     * 
     * @brief Defines bit operations
     * 
     * This class defines several bit manipulation and indexing operations used by FPR.
     */
    class bit_ops
    {
    public:
        /**
         * @brief Construct a new bit ops object
         * 
         * This classes constructor is deleted as it only contains static methods.
         */
        bit_ops() = delete;

        /**
         * @brief Determine chunk index
         * 
         * This method determines the chunk index based on a value's index and the chunk's aligned byte size 
         * in a memory block (e.g. array).
         * A chunk consists of 8 consequtive values in memory which are processed sequentially by a thread.
         * 
         * @param value_idx A value's index in a memory block
         * @param aligned_chunk_size 
         * @return std::uint64_t Chunk index
         */
        static std::uint64_t chunk_idx(std::uint64_t value_idx, std::uint8_t aligned_chunk_size);

        /**
         * @brief Determine bit index
         * 
         * This method determines the bit index based on a value's index in a memory block (e.g. array).
         * 
         * @param value_idx A value's index in a memory block
         * @return std::uint8_t Bit index
         */
        static std::uint8_t bit_idx(std::uint64_t value_idx);

        /**
         * @brief Update a bit
         * 
         * Updates a bit's value at a given position for a given byte.
         * 
         * @param bit New bit value
         * @param pos Bit position in byte
         * @param byte Byte to be modified
         */
        static void update_bit(bool bit, std::uint8_t pos, std::uint8_t &byte);

        /**
         * @brief Read a bit
         * 
         * Read a bit's value at a given position for a given byte.
         * 
         * @param pos Bit position in byte
         * @param byte Byte to be read from
         * @return true The bit's value is 1
         * @return false Thte bit's value is 0
         */
        static bool get_bit(std::uint8_t pos, const std::uint8_t &byte);
    };


}