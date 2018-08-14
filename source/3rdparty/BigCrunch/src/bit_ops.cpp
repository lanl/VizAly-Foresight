#include <bigcrunch/bit_ops.hpp>

#include <cmath>




namespace bigcrunch
{


    std::uint64_t bit_ops::chunk_idx(std::uint64_t value_idx, std::uint8_t aligned_chunk_size)
    {
        // chunk_num_values is fixed to 8
        // (value_idx / chunk_num_values) * aligned_chunk_size
        return (value_idx / 8) * aligned_chunk_size;
    }


    std::uint8_t bit_ops::bit_idx(std::uint64_t value_idx)
    {
        // chunk_num_values is fixed to 8
        // value_idx % chunk_num_values
        return value_idx % 8;
    }
    
    
    void bit_ops::update_bit(bool bit, std::uint8_t pos, std::uint8_t &byte)
    {
        if(bit)
        {
            // set bit to 1
            byte |= 1 << pos;
        }
        else
        {
            // set bit to 0
            byte &= ~(1 << pos);
        }
    }


    bool bit_ops::get_bit(std::uint8_t pos, const std::uint8_t &byte)
    {
        return byte & (1 << pos);
    }


}