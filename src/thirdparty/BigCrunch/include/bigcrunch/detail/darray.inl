#pragma once

#include <typeinfo>
#include <type_traits>
#include <stdexcept>




namespace bigcrunch
{


    template<typename T> inline
    darray::darray(T *data, std::uint64_t size)
        : data_m(data),
        size_m(size)
    {
        this->type_m = this->type_id<T>();
    }


    template<typename T> inline
    dtype_t darray::type_id() const
    {
        if(typeid(T) == typeid(std::uint8_t))
        {
            return dtype_t::UINT8;
        }

        if(typeid(T) == typeid(std::uint16_t))
        {
            return dtype_t::UINT16;
        }

        if(typeid(T) == typeid(std::uint32_t))
        {
            return dtype_t::UINT32;
        }

        if(typeid(T) == typeid(std::uint64_t))
        {
            return dtype_t::UINT64;
        }

        if(typeid(T) == typeid(std::int8_t))
        {
            return dtype_t::INT8;
        }

        if(typeid(T) == typeid(std::int16_t))
        {
            return dtype_t::INT16;
        }

        if(typeid(T) == typeid(std::int32_t))
        {
            return dtype_t::INT32;
        }

        if(typeid(T) == typeid(std::int64_t))
        {
            return dtype_t::INT64;
        }

        if(typeid(T) == typeid(float))
        {
            return dtype_t::FLOAT32;
        }

        if(typeid(T) == typeid(double))
        {
            return dtype_t::FLOAT64;
        }

        throw std::runtime_error("Unsupported data type!");
    }


    template<typename T> inline
    const T * darray::data() const
    {
        if(!this->verify_type<T>())
        {
            throw std::runtime_error("Wrong data type!");
        }

        return static_cast<T *>(this->data_m);
    }


    template<typename T> inline
    T * darray::data()
    {
        if(!this->verify_type<T>())
        {
            throw std::runtime_error("Wrong data type!");
        }

        return static_cast<T *>(this->data_m);
    }


    template<typename T> inline
    bool darray::verify_type() const
    {
        dtype_t cast_type = this->type_id<T>();

        if(cast_type == this->type_m)
        {
            return true;
        }

        return false;
    }


}