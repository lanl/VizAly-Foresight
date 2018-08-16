#ifndef _DATA_LOADER_INTERFACE_H_
#define _DATA_LOADER_INTERFACE_H_

#include <string>
#include <mpi.h>

class DataLoaderInterface
{
  protected:
	std::string filename;
	std::string dataType;
	std::string param;
	//size_t numElements;
	//void *data;

	MPI_Comm comm;

  public:
  	virtual void init(std::string _filename, MPI_Comm _comm) = 0;
  	virtual int loadData(std::string paramName) = 0;
  	virtual std::string getDataInfo() = 0;

    size_t numElements;
    size_t elemSize;
    void *data;
};



#endif
