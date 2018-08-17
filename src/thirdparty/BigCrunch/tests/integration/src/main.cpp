#include <iostream>
#include <iomanip>
#include <vector>
#include <random>
#include <algorithm>
#include <functional>
#include <cmath>

#include <bigcrunch/bigcrunch.hpp>




double relative_error(double original_value, double reconstructed_value, double tolerance)
{
    double absolute_error = std::abs(original_value - reconstructed_value);

    if(std::abs(original_value) < tolerance)
    {
        return absolute_error;
    }

    return absolute_error / std::abs(original_value);
}


int main(int argc, char **argv)
{
    // init test data
    std::vector<float> data(1000000, 0);
    std::uint64_t osize = data.size() * sizeof(float);

    // populate with uniform random float numbers
    std::mt19937_64::result_type seed = 123;
    auto uniform = std::bind(std::uniform_real_distribution<float>(-5000.0, 5000.0), std::mt19937_64(seed));
    std::generate(data.begin(), data.end(), uniform);

    // init bigcrunch
    bigcrunch::setting_t settings = {{bigcrunch::config_t::ERROR,                                             -3},
                                     {bigcrunch::config_t::TOLERANCE,                                          0},
                                     {bigcrunch::config_t::BLOSC_NTHREADS,                                     1},
                                     {bigcrunch::config_t::BLOSC_FILTER,      bigcrunch::blosc_filter_t::SHUFFLE},
                                     {bigcrunch::config_t::BLOSC_COMPRESSOR, bigcrunch::blosc_compressor_t::ZSTD}};
    bigcrunch::bigcrunch bc(settings);

    // compress data
    std::uint8_t *cdata = nullptr;
    std::uint64_t csize = bc.compress(bigcrunch::darray(data.data(), data.size()), &cdata);

    // decompress data
    auto rdata_array = bc.decompress(cdata, csize);
    std::vector<float> rdata(rdata_array.data<float>(), rdata_array.data<float>() + rdata_array.size());

    // compute relative error
    std::vector<double> rel_error(data.size(), 0);
    double tolerance = std::pow(2, settings[bigcrunch::config_t::TOLERANCE]);
    for(std::uint64_t i = 0; i < data.size(); ++i)
    {
        rel_error[i] = relative_error(data[i], rdata[i], tolerance);
    }

    double max_rel_error = *std::max_element(rel_error.begin(), rel_error.end());

    // results
    std::cout << "Original size: " << osize << std::endl;
    std::cout << "Compressed size: " << csize << std::endl;
    std::cout << "Compression ratio: " << osize / static_cast<double>(csize) << std::endl;
    std::cout << "Relative error: " << max_rel_error << std::endl;

    return 0;
}