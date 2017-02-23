#include <fio.hpp>

#include <stdexcept>




fio::fio(std::ios::openmode mode)
	: mode(mode)
{
}


fio::fio(const std::string &filename, std::ios::openmode mode)
	: mode(mode),
	filename_m(filename)
{
	this->open(this->filename_m, this->mode);
}


fio::~fio()
{
	this->close();
}


std::size_t fio::size()
{
	std::size_t current_position = this->file.tellg();
	this->file.seekg(0, std::ios::end);
	std::size_t result = this->file.tellg();
	this->file.seekg(current_position);

	return result;
}


std::string fio::filename() const
{
	return this->filename_m;
}


void fio::open(const std::string &filename, std::ios::openmode mode)
{
	if(this->file.is_open())
	{
		this->file.close();
	}

	this->filename_m = filename;
	this->mode = mode;
	this->file.open(this->filename_m, this->mode);

	if(!this->file.is_open())
	{
		throw std::runtime_error("Opening file failed!\n\tFile: " + this->filename_m);
	}
}


void fio::close()
{
	if(this->file.is_open())
	{
		this->file.close();

		if(this->file.is_open())
		{
			throw std::runtime_error("Closing file failed!\n\tFile: " + this->filename_m);
		}

		this->filename_m = "";
	}
}


bool fio::eof()
{
	return this->file.eof();
}


std::string fio::readline()
{
	std::string line = "";

	if(!this->file.eof())
	{
		std::getline(this->file, line);
	}

	return line;
}


std::vector<std::string> fio::readlines(std::size_t numlines)
{
	std::vector<std::string> result;

	while(!this->file.eof() && result.size() < numlines)
	{
		result.push_back(this->readline());
	}

	if(result.back().empty())
	{
		result.pop_back();
	}

	return result;
}


void fio::writeline(const std::string &line)
{
	this->file << line << "\n";
}


void fio::writelines(const std::vector<std::string> &lines)
{
	for(auto &line : lines)
	{
		this->writeline(line);
	}
}
