#include <bigcrunch/byte_buffer.hpp>




namespace bigcrunch
{


    byte_buffer::byte_buffer(std::uint64_t size, endianess_t target)
        : target(target),
        size_m(size),
        buffer_pos(0)
    {
        // determine system endianess
        this->system = this->endianess();
        // init buffer array
        this->buffer_m = new std::uint8_t[this->size_m];
    }


    endianess_t byte_buffer::endianess() const
    {
        int value = 1;

        if(reinterpret_cast<char *>(&value)[0] == 1)
        {
            return endianess_t::LITTLE;
        }

        return endianess_t::BIG;
    }


    void byte_buffer::assign(const std::uint8_t *data, std::uint64_t size)
    {
        // boundary check
        if(this->buffer_pos + size > this->size_m)
        {
            throw std::runtime_error("Buffer out of bounds!");
        }

        // copy data into buffer array
        std::copy(data, data + size, this->buffer_m + this->buffer_pos);

        // adjust current buffer position
        this->buffer_pos += size;
    }


    std::uint64_t byte_buffer::tell() const
    {
        return this->buffer_pos;
    }


    void byte_buffer::seek(std::uint64_t buffer_pos)
    {
        // boundary check
        if(buffer_pos > this->size_m)
        {
            throw std::runtime_error("Buffer out of bounds!");
        }
        
        this->buffer_pos = buffer_pos;
    }


    std::uint8_t * byte_buffer::buffer()
    {
        return this->buffer_m;
    }


    const std::uint8_t * byte_buffer::buffer() const
    {
        return this->buffer_m;
    }


    std::uint64_t byte_buffer::size() const
    {
        return this->size_m;
    }


}