#ifndef DATA_HANDLER_HPP
#define DATA_HANDLER_HPP

#include <string>
#include <vector>

#include <data_field.hpp>




class data_handler
{
	public:
		virtual ~data_handler() = default;

		virtual void read(const std::string &in_name, std::vector<data_field> &data) = 0;
		virtual void write(const std::string &in_name, const std::string &out_name, const std::vector<data_field> &data) = 0;
};


#endif // DATA_HANDLER_HPP
