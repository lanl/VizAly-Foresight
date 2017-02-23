/**
 * @brief A concrete instance of GenericIOReader that implements the GenericIO
 * reader interface using MPI.
 */
#ifndef GENERICIOMPIREADER_H_
#define GENERICIOMPIREADER_H_

// CosmologyTools includes
#include "GenericIOReader.h"

// C/C++ includes
#include <sys/types.h> // for off_t


// MPI includes
#include <mpi.h>

namespace gio
{

class GenericIOMPIReader : public GenericIOReader
{
public:
  GenericIOMPIReader();
  virtual ~GenericIOMPIReader();

  /**
   * @brief Closes the file
   */
  void Close();

protected:
  MPI_File FH;

  /**
   * @brief Open the FileHandle.
   */
  void Open();

  /**
   * @brief Reads data into the user-supplied buffer from the MPI file handle
   * @param buf the buffer where the data will be read into
   * @param count the number of bytes to read
   * @param offset the offset from which to read
   * @param variableName name of the data being read, primarily, used for
   * debugging and error reporting.
   */
  void Read(void *buf, size_t count,
                      off_t offset, const std::string &variableName );

  /**
   * @brief Allocates N the internal readers array.
   * @param N number of readers to allocate.
   * @pre N > 0
   */
  void AllocateInternalReaders(const int N);

  /**
   * @brief Returns a new GenericIOMPIReader instance.
   * @return r pointer to a new GenericIOMPIReader instance.
   * @note caller is responsible for properly de-allocating the return object.
   * @see implements GenericIOReader::GetNewInstance()
   */
  GenericIOReader* GetNewInstance()
    { return( new GenericIOMPIReader() ); };


private:
  GIO_DISABLE_COPY_AND_ASSIGNMENT(GenericIOMPIReader);
};

} /* namespace cosmotk */
#endif /* GENERICIOMPIREADER_H_ */
