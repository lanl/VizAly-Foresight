#include <data_field.hpp>

#include <cstdlib>




data_field::data_field()
	: data(nullptr),
	size(0),
	type(data_field::data_type::NONE)
{
}


data_field::~data_field()
{
	if(this->data != nullptr)
	{
		std::free(this->data);
		this->data = nullptr;
	}
}
