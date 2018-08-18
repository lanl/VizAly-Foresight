#include <bigcrunch/bigcrunch.hpp>

#include <iostream>

#include <bigcrunch/fpr.hpp>
#include <bigcrunch/blosc_wrapper.hpp>




namespace bigcrunch
{


    bigcrunch::bigcrunch(const setting_t &settings)
        : error(-4),
        tolerance(0),
        nthreads(1),
        blosc_config({blosc_compressor_t::ZSTD, 9, blosc_filter_t::BITSHUFFLE})
    {
        for(auto &setting : settings)
        {
            switch(setting.first)
            {
            case config_t::ERR:
                this->error = setting.second;
                break;
            case config_t::TOLERANCE:
                this->tolerance = setting.second;
                break;
            case config_t::BLOSC_NTHREADS:
                this->nthreads = setting.second;
                break;
            case config_t::BLOSC_COMPRESSOR:
                this->blosc_config.cname = static_cast<blosc_compressor_t>(setting.second);
                break;
            case config_t::BLOSC_CLEVEL:
                this->blosc_config.clevel = setting.second;
                break;
            case config_t::BLOSC_FILTER:
                this->blosc_config.filter = static_cast<blosc_filter_t>(setting.second);
                break;
            default:
                std::cerr << "Setting " << setting.first << " is unknown and will be ignored!";
            }
        }
    }


    std::uint64_t bigcrunch::compress(const darray &data, std::uint8_t **output) const
    {
        // encode original values with fpr
        fpr_data internal = fpr::encode(data, this->error, this->tolerance);

        // compress encoded values
        blosc_wrapper blosc(this->nthreads);
        blosc.compress(internal, this->blosc_config);

        // serialize result array and return result size
        std::uint64_t result_size;
        internal.serialize(output, result_size);
        return result_size;
    }


    darray bigcrunch::decompress(const std::uint8_t *data, std::uint64_t size) const
    {
        // deserialize binary data
        fpr_data internal;
        internal.deserialize(data, size);

        // uncompress data
        blosc_wrapper blosc(this->nthreads);
        blosc.decompress(internal);

        // decode uncompressed data and return the result
        return fpr::decode(internal);
    }


}