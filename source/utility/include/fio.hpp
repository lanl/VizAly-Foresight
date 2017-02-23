#ifndef FIO_HPP
#define FIO_HPP

#include <fstream>
#include <cstddef>
#include <string>
#include <vector>




class fio
{
	public:
		fio(std::ios::openmode mode = std::ios::in | std::ios::out);
		fio(const std::string &filename, std::ios::openmode mode = std::ios::in | std::ios::out);
		~fio();

		std::size_t size();
		std::string filename() const;

		void open(const std::string &filename, std::ios::openmode mode = std::ios::in | std::ios::out);
		void close();
		bool eof();

		std::string readline();
		std::vector<std::string> readlines(std::size_t numlines);
		void writeline(const std::string &line);
		void writelines(const std::vector<std::string> &lines);

	private:
		std::fstream file;
		std::ios::openmode mode;
		std::string filename_m;
};


#endif // FIO_HPP
