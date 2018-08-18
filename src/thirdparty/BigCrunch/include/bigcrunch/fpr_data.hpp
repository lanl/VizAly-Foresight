/**
 * @file fpr_data.hpp
 * @author Max Zeyen
 * @date 2018-04-19
 */

#pragma once

#include <cstdint>




namespace bigcrunch
{


    /**
     * @class fpr_data
     * 
     * @brief Internal data exchange object
     * 
     * This class defines the internal memory representation to communicate between pipeline steps.
     * It also provides implementations for de/serialization.
     */
    class fpr_data
    {
    public:
        std::int8_t max_exp;
        std::int8_t offset;
        std::int8_t tolerance;
        std::int8_t error;
        std::uint8_t type;
        std::uint8_t type_size;
        std::uint64_t num_elems;
        std::uint64_t compressed_size;
        std::uint8_t *data;

        /**
         * @brief Construct a new fpr_data object
         * 
         * Default constructor building an empty fpr_data object.
         */
        fpr_data();

        /**
         * @brief Destroy the fpr_data object
         * 
         * Deletes the byte array upon destruction.
         */
        ~fpr_data();

        /**
         * @brief Serialize the contents of this class
         * 
         * Serializes the contents of this class into a byte array.
         * This is used for generating outputs for \ref bigcrunch.
         * 
         * @param output Serialized byte array
         * @param size Size of serialized byte array
         */
        void serialize(std::uint8_t **output, std::uint64_t &size);

        /**
         * @brief Deserialize a byte array to populate this class
         * 
         * Deserializes a byte array to populate the contents of this class.
         * This is used for generating inputs for \ref bigcrunch.
         * 
         * @param input Byte array to deserialize
         * @param size Size of byte array
         */
        void deserialize(const std::uint8_t *input, std::uint64_t size);

    private:
        std::uint64_t header_size() const;
    };


}