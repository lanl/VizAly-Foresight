/**
 * @file darray.hpp
 * @author Max Zeyen
 * @date 2018-04-19
 */

#pragma once

#include <cstdint>
#include <string>

#include <bigcrunch/types.hpp>




namespace bigcrunch
{


    /**
     * @class darray
     * 
     * @brief Data array class for arbitrary simple types
     * 
     * This class defines a data array for arbitrary simple types.
     * The supported types are limited by the \ref dtype_t enum.
     * This class is set up like an any type but for arrays and with more type restrictions.
     */
    class darray
    {
    public:
        /**
         * @brief Construct a new darray object
         * 
         * Construct a darray from a c-style array with a supported type.
         * 
         * @tparam T Type of the data
         * @param data Data to construct the darray from
         * @param size Size of the data
         */
        template<typename T> darray(T *data, std::uint64_t size);

        /**
         * @brief Gets the size
         * 
         * Returns the number of elements in the darray.
         * 
         * @return std::uint64_t Number of elements in the darray.
         */
        std::uint64_t size() const;

        /**
         * @brief Gets the type id
         * 
         * Returns the type id of the elements in the darray.
         * 
         * @return dtype_t Type id of darray data
         */
        dtype_t type() const;

        /**
         * @brief Gets the type name
         * 
         * Returns the type id in a human readible string.
         * 
         * @return std::string Name of the type id
         */
        std::string type_name() const;

        /**
         * @brief Determines the type id of a given type
         * 
         * Returns the type id determined for the template parameter.
         * 
         * @tparam T Type whose id will be determined
         * @return dtype_t Type id of the template parameter
         */
        template<typename T> dtype_t type_id() const;

        /**
         * @brief Verify type correctness
         * 
         * Verifies the equality of the template parameter and the darray's internal type.
         * 
         * @tparam T Type to verify
         * @return true The template parameter is equal to the internal type.
         * @return false The template parameter is not equal to the internal type.
         */
        template<typename T> bool verify_type() const;

        /**
         * @brief Get a constant pointer to the internal data array
         * 
         * Cast the internal data array into the given type if correct and return a constant pointer to the array.
         * 
         * @tparam T Type to cast the data into
         * @return const T* Constant pointer to cast data array
         */
        template<typename T> const T * data() const;

        /**
         * @brief Get a pointer to the internal data array
         * 
         * Cast the internal data array into the given type if correct and return a pointer to the array.
         * 
         * @tparam T Type to cast the data into
         * @return T* Pointer to cast data array
         */
        template<typename T> T * data();

    protected:
        void *data_m;
        std::uint64_t size_m;
        dtype_t type_m;
    };


}


#include <bigcrunch/detail/darray.inl>