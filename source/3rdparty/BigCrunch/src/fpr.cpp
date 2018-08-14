#include <bigcrunch/fpr.hpp>

#include <stdexcept>

#include <bigcrunch/fpr_encoder.hpp>
#include <bigcrunch/fpr_decoder.hpp>




namespace bigcrunch
{


    fpr_data fpr::encode(const darray &data, int error, int tolerance)
    {
        // distinguish between float and double
        switch(data.type())
        {
        case dtype_t::FLOAT32:
        {
            // call internal encoder for float data
            fpr_encoder<float> encoder(error, tolerance);
            return encoder.encode(data.data<float>(), data.size());
        }
        case dtype_t::FLOAT64:
        {
            // call internal encoder for double data
            fpr_encoder<double> encoder(error, tolerance);
            return encoder.encode(data.data<double>(), data.size());
        }
        default:
            std::runtime_error("FPR does only support floating-point data!");
        }
    }


    darray fpr::decode(const fpr_data &data)
    {
        // distinguish between float and double
        switch(data.type)
        {
        case 0:
        {
            // call internal decoder for float data
            fpr_decoder<float> decoder;
            return decoder.decode(data);
        }
        case 1:
        {
            // call internal decoder for double data
            fpr_decoder<double> decoder;
            return decoder.decode(data);
        }
        default:
            std::runtime_error("FPR does only support floating-point data!");
        }
    }


}