#include <bigcrunch/darray.hpp>




namespace bigcrunch
{


    std::uint64_t darray::size() const
    {
        return this->size_m;
    }


    dtype_t darray::type() const
    {
        return this->type_m;
    }
    

    std::string darray::type_name() const
    {
        switch(this->type_m)
        {
        case dtype_t::UINT8:
            return "uint8";
        case dtype_t::UINT16:
            return "uint16";
        case dtype_t::UINT32:
            return "uint32";
        case dtype_t::UINT64:
            return "uint64";
        case dtype_t::INT8:
            return "int8";
        case dtype_t::INT16:
            return "int16";
        case dtype_t::INT32:
            return "int32";
        case dtype_t::INT64:
            return "int64";
        case dtype_t::FLOAT32:
            return "float32";
        case dtype_t::FLOAT64:
            return "float64";
        default:
            return "Unknown";
        }
    }

}