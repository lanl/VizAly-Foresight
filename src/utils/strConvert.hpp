/*================================================================================
This software is open source software available under the BSD-3 license.

Copyright (c) 2017, Los Alamos National Security, LLC.
All rights reserved.

Authors:
 - Pascal Grosset
================================================================================*/

#ifndef _STR_CONV_H_
#define _STR_CONV_H_

#include <sstream>

namespace strConvert
{
	inline double to_double(std::string value)
	{
		std::stringstream sstr(value);
		double val;
		sstr >> val;
		return val;
	}

	inline float to_float(std::string value)
	{
		std::stringstream sstr(value);
		float val;
		sstr >> val;
		return val;
	}

	inline int to_int(std::string value)
	{
		std::stringstream sstr(value);
		int val;
		sstr >> val;
		return val;
	}


	inline __int64_t to_int64(std::string value)
	{
		std::stringstream sstr(value);
		int64_t val;
		sstr >> val;
		return val;
	}

	inline __int32_t to_int32(std::string value)
	{
		std::stringstream sstr(value);
		int32_t val;
		sstr >> val;
		return val;
	}

	inline __int16_t to_int16(std::string value)
	{
		std::stringstream sstr(value);
		int16_t val;
		sstr >> val;
		return val;
	}

	inline __int8_t to_int8(std::string value)
	{
		std::stringstream sstr(value);
		int16_t val;
		sstr >> val;
		return val;
	}


	inline __uint64_t to_uint64(std::string value)
	{
		std::stringstream sstr(value);
		uint64_t val;
		sstr >> val;
		return val;
	}

	inline __uint32_t to_uint32(std::string value)
	{
		std::stringstream sstr(value);
		uint32_t val;
		sstr >> val;
		return val;
	}

	inline __uint16_t to_uint16(std::string value)
	{
		std::stringstream sstr(value);
		uint16_t val;
		sstr >> val;
		return val;
	}

	inline __uint8_t to_uint8(std::string value)
	{
		std::stringstream sstr(value);
		uint8_t val;
		sstr >> val;
		return val;
	}
}
#endif