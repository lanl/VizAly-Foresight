/**
 * @brief A concrete instance of GenericIOReader that implements the GenericIO
 * interface with POSIX.
 */

#ifndef GENERICIOPOSIXREADER_H_
#define GENERICIOPOSIXREADER_H_

// CosmoTools includes
#include "GenericIOReader.h"

// C/C++ includes
#include <sys/types.h> // for off_t

namespace gio
{

class GenericIOPosixReader : public GenericIOReader
{
public:
  GenericIOPosixReader();
  virtual ~GenericIOPosixReader();

  /**
   * @brief Closes the file
   */
  void Close();

protected:
  int FH;

  /**
   * @brief Opens the file
   */
  void Open();

  /**
   * @brief Reads data into the user-supplied buffer from the MPI file handle
   * @param buf the buffer where the data will be read into
   * @param count the number of bytes to read
   * @param offset the offset from which to read
   * @param name the name of the data being read, primarily, used for
   * debugging and error reporting.
   */
  void Read(void *buf, size_t count, off_t offset, const std::string &name);

  /**
   * @brief Alloacates the internal readers array.
   * @param N number of readersto allocate.
   * @pre N > 0.
   */
  void AllocateInternalReaders(const int N);

  /**
   * @brief Returns a new GenericIOPosixReader instance.
   * @return r pointer to a new GenericIOPosixReader instance.
   * @note caller is responsible for properly de-allocating the return object.
   * @see implements GenericIOReader::GetNewInstance().
   */
  GenericIOReader* GetNewInstance()
    { return(new GenericIOPosixReader()); };

private:
  GIO_DISABLE_COPY_AND_ASSIGNMENT(GenericIOPosixReader);
};

} /* namespace cosmotk */
#endif /* GENERICIOPOSIXREADER_H_ */
