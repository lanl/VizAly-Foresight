/**
 * @file fpr.hpp
 * @author Max Zeyen
 * @date 2018-05-02
 */

#pragma once

#include <bigcrunch/darray.hpp>
#include <bigcrunch/fpr_data.hpp>




namespace bigcrunch
{


    /**
     * @class fpr
     * 
     * @brief Main class for FPR encoder/decoder
     * 
     * This class defines the main access point for using FPR's functionality.
     * It provides a simple intuitive interface to encode and decode data.
     */
    class fpr
    {
    public:
        /**
         * @brief Deleted default constructor
         * 
         * The default constructor is deleted, as the class only contains static methods and therefore does not 
         * require an instance of the class.
         */
        fpr() = delete;

        /**
         * @brief Encode data
         * 
         * This method is the entry point to FPR's encoder.
         * It takes a one-dimensional floating-point data array, target relative error exponent, and tolerance 
         * exponent as input and returns an encoded data object with metadata.
         * 
         * @param data Input data array
         * @param error Relative error exponent
         * @param tolerance Tolerance exponent
         * @return fpr_data Object with encoded data and metadata
         */
        static fpr_data encode(const darray &data, int error, int tolerance=0);

        /**
         * @brief Decode data
         * 
         * This method is the entry point to FPR's decoder.
         * It takes a data object containing encoded data and metadata and returns a one-dimensional floating-point 
         * data array.
         * 
         * @param data Object with encoded data and metadata 
         * @return darray Reconstructed floating-point data
         */
        static darray decode(const fpr_data &data);
    };


}