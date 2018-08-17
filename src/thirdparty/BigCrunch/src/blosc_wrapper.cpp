#include <bigcrunch/blosc_wrapper.hpp>

#include <cstdint>
#include <string>
#include <stdexcept>
#include <cmath>
#include <algorithm>

extern "C"
{
#include <blosc.h>
#include <shuffle.h>
}




namespace bigcrunch
{


    blosc_wrapper::blosc_wrapper(int nthreads)
    {
        blosc_init();
        blosc_set_nthreads(nthreads);
    }


    blosc_wrapper::~blosc_wrapper()
    {
        blosc_destroy();
    }


    void blosc_wrapper::compress(fpr_data &data, const blosc_config_t &config)
    {
        // set compressor
        std::string compressor = "unknown";
        switch(config.cname)
        {
        case blosc_compressor_t::BLOSCLZ:
            compressor = "blosclz";
            break;
        case blosc_compressor_t::LZ4:
            compressor = "lz4";
            break;
        case blosc_compressor_t::LZ4HC:
            compressor = "lz4hc";
            break;
        case blosc_compressor_t::SNAPPY:
            compressor = "snappy";
            break;
        case blosc_compressor_t::ZLIB:
            compressor = "zlib";
            break;
        case blosc_compressor_t::ZSTD:
            compressor = "zstd";
            break;
        }

        int code = blosc_set_compressor(compressor.c_str());
        if(code < 0)
        {
            throw std::invalid_argument("Error setting " + compressor + " compressor! Does it really exist?");
        }

        // determine aligned size
        std::uint8_t aligned_chunk_size = static_cast<std::uint8_t>(std::ceil(data.type_size / 8.0) * 8.0);

        // determine input and preliminary output sizes
        std::uint64_t isize = static_cast<std::uint64_t>(std::ceil(data.num_elems / 8.0) * aligned_chunk_size);

        // intermediate buffer
        std::uint8_t *temp_data = new std::uint8_t[isize]();

        // shuffle the data without compression first
        // since FPR pre-bitshuffles the chunks, Blosc only needs to perform a normal byte shuffle
        shuffle(aligned_chunk_size, isize, data.data, temp_data);

        // clean temp data
        std::swap(data.data, temp_data);
        if(temp_data != nullptr)
        {
            delete [] temp_data;
            temp_data = nullptr;
        }

        // compute new original size (truncate alignment padding)
        isize = static_cast<std::uint64_t>(std::ceil(data.num_elems / 8.0) * data.type_size);
        std::uint64_t csize = isize + BLOSC_MAX_OVERHEAD;

        // allocate output buffer
        std::uint8_t *cdata = new std::uint8_t[csize]();

        // compress data
        code = blosc_compress(config.clevel, blosc_filter_t::NOSHUFFLE, 1, isize, data.data, cdata, csize);
        if(code < 0)
        {
            throw std::runtime_error("Compression error! Error code: " + std::to_string(code));
        }

        // set output data and resize if necessary
        // a compressed_size of 0 indicates that no compression occured
        if(code > 0)
        {
            std::swap(data.data, cdata);
            if(cdata != nullptr)
            {
                delete [] cdata;
                cdata = nullptr;
            }
            data.compressed_size = code;
        }
        else
        {
            data.compressed_size = 0;
        }
    }


    void blosc_wrapper::decompress(fpr_data &data)
    {
        // check if data was compressed
        if(data.compressed_size == 0)
        {
            return;
        }
        
        // determine aligned size
        std::uint8_t aligned_chunk_size = static_cast<std::uint8_t>(std::ceil(data.type_size / 8.0) * 8.0);

        // determine output size with and without padding
        std::uint64_t temp_osize = static_cast<std::uint64_t>(std::ceil(data.num_elems / 8.0) * data.type_size);
        std::uint64_t osize = static_cast<std::uint64_t>(std::ceil(data.num_elems / 8.0) * aligned_chunk_size);

        // allocate output with padding
        std::uint8_t *temp_data = new std::uint8_t[osize]();

        // first stage decompression: decompress byte stream
        int code = blosc_decompress(data.data, temp_data, temp_osize);
        if(code < 0)
        {
            throw std::runtime_error("Decompression error! Error code: " + std::to_string(code));
        }

        std::swap(data.data, temp_data);
        if(temp_data != nullptr)
        {
            delete [] temp_data;
            temp_data = nullptr;
        }

        // allocate output with padding
        std::uint8_t *udata = new std::uint8_t[osize]();

        // second stage decompression: deshuffle byte stream
        unshuffle(aligned_chunk_size, osize, data.data, udata);

        std::swap(data.data, udata);
        if(udata != nullptr)
        {
            delete [] udata;
            udata = nullptr;
        }
    }


}