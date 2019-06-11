/*================================================================================
This software is open source software available under the BSD-3 license.

Copyright (c) 2017, Los Alamos National Security, LLC.
All rights reserved.

Pascal Grosset
================================================================================*/

#ifndef _GIO_PV_GIO_DATA_H_
#define _GIO_PV_GIO_DATA_H_

#include <string>

//
// Stores the genericIO data being read in
class GioData
{
  public:               // coz too lazy to do private!
    int id;
    std::string name;
    int size;           // in bytes
    bool isFloat;
    bool isSigned;
    bool ghost;
    bool xVar, yVar, zVar;
    void* data;

    std::string dataType;
    size_t numElements;

    bool loadData;
    bool doWrite;

  public:
    GioData();
    GioData(int _id, std::string _name, int _size, bool _isFloat, bool _isSigned, bool _xVar, bool _yVar, bool _zVar, void* _data=NULL);
    ~GioData();

    void init(int _id, std::string _name, int _size, bool _isFloat, bool _isSigned, bool _xVar, bool _yVar, bool _zVar);
    void setNumElements(size_t _numElements) { numElements = _numElements; }

    int determineDataType();

  	int allocateMem(int offset = 1);
  	int deAllocateMem();
};


inline GioData::GioData()
{
    dataType = "";
    numElements = 0;
    xVar = yVar = zVar = false;
    loadData = false;
    doWrite = false;
    data = NULL;
}

inline GioData::GioData(int _id, std::string _name, int _size, bool _isFloat, bool _isSigned, bool _xVar, bool _yVar, bool _zVar, void* _data)
{
    init(_id, _name, _size, _isFloat, _isSigned, _xVar, _yVar, _zVar);
    data = _data;
    doWrite = false;
}


inline GioData::~GioData()
{
    dataType = "";
    numElements = 0;
}

inline void GioData::init(int _id, std::string _name, int _size, bool _isFloat, bool _isSigned, bool _xVar, bool _yVar, bool _zVar)
{
    id = _id;
    name = _name;
    size = _size;
    isFloat = _isFloat;
    isSigned = _isSigned;
    xVar = _xVar;
    yVar = _yVar;
    zVar = _zVar;
}

inline int GioData::determineDataType()
{
    if (isFloat)
    {
        if (size == 4)
            dataType = "float";
        else if (size == 8)
            dataType = "double";
        else
            return 0;
    }
    else // not float
    {
        if (isSigned)
        {
            if (size == 1)
                dataType = "int8_t";
            else if (size == 2)
                dataType = "int16_t";
            else if (size == 4)
                dataType = "int32_t";
            else if (size == 8)
                dataType = "int64_t";
            else
                return 0;
        }
        else
        {
            if (size == 1)
                dataType = "uint8_t";
            else if (size == 2)
                dataType = "uint16_t";
            else if (size == 4)
                dataType = "uint32_t";
            else if (size == 8)
                dataType = "uint64_t";
            else
                return 0;
        }
    }

    return 1;
}


inline int GioData::allocateMem(int offset)
{
    determineDataType();

    if (dataType == "float")
        data = new float[numElements + offset];
    else if (dataType == "double")
        data = new double[numElements + offset];
    else if (dataType == "int8_t")
        data = new int8_t[numElements + offset];
    else if (dataType == "int16_t")
        data = new int16_t[numElements + offset];
    else if (dataType == "int32_t")
        data = new int32_t[numElements + offset];
    else if (dataType == "int64_t")
        data = new int64_t[numElements + offset];
    else if (dataType == "uint8_t")
        data = new uint8_t[numElements + offset];
    else if (dataType == "uint16_t")
        data = new uint16_t[numElements + offset];
    else if (dataType == "uint32_t")
        data = new uint32_t[numElements + offset];
    else if (dataType == "uint64_t")
        data = new uint64_t[numElements + offset];
    else
        return 0;

    return 1;
}


inline int GioData::deAllocateMem()
{
    if (data == NULL) // already deallocated!
        return 1;

    if (dataType == "float")
        delete[](float*) data;
    else if (dataType == "double")
        delete[](double*) data;
    else if (dataType == "int8_t")
        delete[](int8_t*) data;
    else if (dataType == "int16_t")
        delete[](int16_t*) data;
    else if (dataType == "int32_t")
        delete[](int32_t*) data;
    else if (dataType == "int64_t")
        delete[](int64_t*) data;
    else if (dataType == "uint8_t")
        delete[](uint8_t*) data;
    else if (dataType == "uint16_t")
        delete[](uint16_t*) data;
    else if (dataType == "uint32_t")
        delete[](uint32_t*) data;
    else if (dataType == "uint64_t")
        delete[](uint64_t*) data;
    else
        return 0;

    data = NULL;

    return 1;
}


#endif