#include <bigcrunch/fpr_data.hpp>

#include <bigcrunch/byte_buffer.hpp>




namespace bigcrunch
{


    fpr_data::fpr_data()
        : data(nullptr)
    {
    }


    fpr_data::~fpr_data()
    {
        if(this->data != nullptr)
        {
            delete [] this->data;
            this->data == nullptr;
        }
    }


    void fpr_data::serialize(std::uint8_t **output, std::uint64_t &size)
    {
        // allocate byte buffer
        size = this->header_size() + this->compressed_size;
        byte_buffer serial(size);

        // serialize header
        serial << this->max_exp << this->offset << this->tolerance << this->error << this->type;
        serial << this->type_size << this->num_elems << this->compressed_size;

        // serialize data array
        for(std::uint64_t i = 0; i < this->compressed_size; ++i)
        {
            serial << this->data[i];
        }

        // set output array
        *output = serial.buffer();
    }


    void fpr_data::deserialize(const std::uint8_t *input, std::uint64_t size)
    {
        // init byte buffer
        byte_buffer serial(size);
        serial.assign(input, size);
        // reset buffer position to the first element
        serial.seek(0);

        // deserialize header
        serial >> this->max_exp >> this->offset >> this->tolerance >> this->error >> this->type;
        serial >> this->type_size >> this->num_elems >> this->compressed_size;

        // free previous memory if necessary
        if(this->data != nullptr)
        {
            delete [] this->data;
            this->data = nullptr;
        }
        
        // allocate memory for array data
        this->data = new std::uint8_t[this->compressed_size];

        // deserialize data array
        for(std::uint64_t i = 0; i < this->compressed_size; ++i)
        {
            serial >> this->data[i];
        }
    }


    std::uint64_t fpr_data::header_size() const
    {
        return 4 * sizeof(std::int8_t) + 2 * sizeof(std::uint8_t) + 2 * sizeof(std::uint64_t);
    }


}