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
  public: // coz too lazy to do private!
	int id;
	std::string name;
	int size; // in bytes
	bool isFloat;
	bool isSigned;
	bool ghost;
	bool xVar, yVar, zVar;

	std::string dataType;
	bool loadData;
	size_t numElements;

  public:
	GioData();
	~GioData();

	void init(int _id, std::string _name, int _size, bool _isFloat, bool _isSigned);
	void setNumElements(size_t _numElements) { numElements = _numElements; }

	int determineDataType();
};

inline GioData::GioData()
{
	dataType = "";
	numElements = 0;
	xVar = yVar = zVar = false;
	loadData = false;
}

inline GioData::~GioData()
{
	dataType = "";
	numElements = 0;
}

inline void GioData::init(int _id, std::string _name, int _size, bool _isFloat, bool _isSigned)
{
	id = _id;
	name = _name;
	size = _size;
	isFloat = _isFloat;
	isSigned = _isSigned;
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

#endif