#pragma once

#include <algorithm>
#include <stdexcept>




namespace bigcrunch
{


    template<typename T> inline
    byte_buffer & byte_buffer::operator<<(const T &value)
    {
        // boundary check
        if(this->buffer_pos + sizeof(T) > this->size_m)
        {
            throw std::runtime_error("Buffer out of bounds!");
        }

        // convert value to byte array
        const std::uint8_t *temp = reinterpret_cast<const std::uint8_t *>(&value);

        // check if target and system endianess are matching
        if(this->target != this->system)
        {
            // reverse copy the value bytes into the stream in case of an endianess mismatch
            std::reverse_copy(temp, temp + sizeof(T), this->buffer_m + this->buffer_pos);
        }
        else
        {
            // normal copy the value bytes into the stream in case of an endianess match
            std::copy(temp, temp + sizeof(T), this->buffer_m + this->buffer_pos);
        }

        // adjust current buffer position
        this->buffer_pos += sizeof(T);

        // return this byte_buffer object
        return *this;
    }


    template<typename T> inline
    byte_buffer & byte_buffer::operator>>(T &value)
    {
        // boundary check
        if(this->buffer_pos + sizeof(T) > this->size_m)
        {
            throw std::runtime_error("Buffer out of bounds!");
        }

        // convert value to byte array
        std::uint8_t *temp = reinterpret_cast<std::uint8_t *>(&value);

        // check if target and system endianess are matching
        if(this->target != this->system)
        {
            // reverse copy the buffer bytes into the value bytes in case of an endianess mismatch
            std::reverse_copy(this->buffer_m + this->buffer_pos, this->buffer_m + this->buffer_pos + sizeof(T), temp);
        }
        else
        {
            // normal copy the buffer bytes into the value bytes in case of an endianess match
            std::copy(this->buffer_m + this->buffer_pos, this->buffer_m + this->buffer_pos + sizeof(T), temp);
        }

        // adjust current buffer position
        this->buffer_pos += sizeof(T);

        // return this byte_buffer object
        return *this;
    }


}