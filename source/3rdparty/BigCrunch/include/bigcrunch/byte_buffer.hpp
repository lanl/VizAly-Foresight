/**
 * @file byte_buffer.hpp
 * @author Max Zeyen
 * @date 2018-05-10
 */

#pragma once

#include <cstdint>
#include <vector>

#include <bigcrunch/types.hpp>




namespace bigcrunch
{


    /**
     * @class byte_buffer
     * 
     * @brief Class for byte buffer reading and writing
     * 
     * This class defines methods to read and write data out of and in to a byte buffer with a given target endianess.
     * The byte buffer is used as a helper for data serialization.
     */
    class byte_buffer
    {
    public:
        /**
         * @brief Construct a new byte buffer object
         * 
         * Construct an empty byte buffer with given memory reservations and a given target endianess.
         * 
         * @param size Size of reserved memory for the byte buffer
         * @param target Endianess in which the buffer is written
         */
        byte_buffer(std::uint64_t size, endianess_t target=endianess_t::LITTLE);

        /**
         * @brief Determines the system endianess
         * 
         * This method computes the system endianess.
         * This is required to determine if byte order changes are required or not.
         * 
         * @return endianess_t The system's endianess
         */
        endianess_t endianess() const;

        /**
         * @brief Assign the given data to the byte buffer
         * 
         * Assigns the values of the given data array to the byte buffer.
         * The values are copied to the current position in the buffer.
         * The current position is then updated to the index past the last element in the data array.
         * 
         * @param data 
         * @param size 
         */
        void assign(const std::uint8_t *data, std::uint64_t size);

        /**
         * @brief Write value into byte buffer
         * 
         * Operator to write an arbitrary simple value into the byte buffer.
         * The new value is appended after the current buffer position of the buffer.
         * 
         * @tparam T Type of the value
         * @param value Value to be inserted into the buffer
         * @return byte_buffer& Reference to the modified buffer object
         */
        template<typename T> byte_buffer & operator<<(const T &value);

        /**
         * @brief Read value from byte buffer
         * 
         * Operator to read an arbitrary simple value from the byte buffer.
         * The bytes of the read value are read at the current buffer position of the buffer.
         * 
         * @tparam T Type of the value
         * @param value Value to be extracted from the buffer
         * @return byte_buffer& Reference to the modified buffer object
         */
        template<typename T> byte_buffer & operator>>(T &value);
        
        /**
         * @brief Get current buffer position
         * 
         * Returns the current internal buffer position.
         * The buffer position is an absolute position starting at the first element in the buffer.
         * 
         * @return std::uint64_t Internal buffer position
         */
        std::uint64_t tell() const;

        /**
         * @brief Set current buffer position
         * 
         * Moves the current internal buffer position to the given position.
         * The given buffer position is treated as an absolute buffer position starting at the first element 
         * in the buffer.
         * 
         * @param buffer_pos New buffer position
         */
        void seek(std::uint64_t buffer_pos);

        /**
         * @brief Get byte buffer data
         * 
         * Returns a non-constant version of the byte buffer's internal data representation.
         * 
         * @return std::uint8_t* Byte buffer as non-constant unsigned byte array
         */
        std::uint8_t * buffer();

        /**
         * @brief Get byte buffer data
         * 
         * Returns a constant version of the byte buffer's internal data representation.
         * 
         * @return const std::uint8_t* Byte buffer as constant unsigned byte array
         */
        const std::uint8_t * buffer() const;

        /**
         * @brief Get the byte buffer's size
         * 
         * Returns the size of the byte buffer.
         * 
         * @return std::uint64_t The byte buffer's size
         */
        std::uint64_t size() const;

    protected:
        endianess_t system, target;
        std::uint8_t * buffer_m;
        std::uint64_t size_m;
        std::uint64_t buffer_pos;
    };


}


#include <bigcrunch/detail/byte_buffer.inl>