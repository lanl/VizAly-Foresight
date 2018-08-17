#include <bigcrunch/fpr_stats.hpp>

#include <cmath>




namespace bigcrunch
{


    fpr_stats::fpr_stats(int error, int max_exp, int offset)
        : error(error),
        max_exp(max_exp),
        offset(offset)
    {
        this->relative_error();
        this->exponent_bits();
        this->mantissa_bits();
        this->value_bits = 1 + this->exp_bits + this->mant_bits;
        this->align_chunk_size();
    }


    void fpr_stats::relative_error()
    {
        this->rel_error = std::pow(2, std::floor(std::log2(std::pow(10, this->error))));
    }


    void fpr_stats::exponent_bits()
    {
        this->exp_bits = static_cast<std::uint64_t>(std::ceil(std::log2(this->max_exp - this->offset + 1)));
    }


    void fpr_stats::mantissa_bits()
    {
        this->mant_bits = static_cast<std::uint64_t>(std::ceil(std::log2(1.0 / (2.0 * this->rel_error))));
    }


    void fpr_stats::align_chunk_size()
    {
        this->aligned_chunk_size = static_cast<std::uint64_t>(std::ceil(this->value_bits / 8.0) * 8.0);
    }


}