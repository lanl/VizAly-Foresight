#ifndef DATA_FIELD_HPP
#define DATA_FIELD_HPP

#include <cstddef>




class data_field
{
	public:
		enum data_type
		{
			FLOAT,
			DOUBLE,
			BYTE,
			NONE
		};

		data_field();
		~data_field();

		void *data;
		std::size_t size;
		data_field::data_type type;
};


#endif // DATA_FIELD_HPP
