/**
 * @brief DataFilter defines an abstract class for all filters.
 */

#ifndef DATAFILTER_H_
#define DATAFILTER_H_

#include "GenericIODefinitions.hpp"

namespace gio
{

// Forward declarations
struct DataBuffer {
  void* Data;
  size_t TypeSize;
  size_t NumBytes;
};

class DataFilter
{
public:
  DataFilter();
  virtual ~DataFilter();

  /**
   * @brief Returns the name of this DataFilter instance.
   * @return name a string corresponding to the name of this filter instance.
   */
  const char* GetName()
    {return this->Name.c_str();};

  /**
   * @brief Sets the data buffers the filter will operate on.
   * @param src pointer to the source buffer.
   * @param src_typesize the number of bytes of the atomic type in `src`
   * @param src_nbytes the total number of bytes in src.
   */
  void SetInputData(void* src, size_t src_typesize, size_t src_nbytes);


  /**
   * @brief Returns the size of the output after the filter is executed.
   * @return N the number of bytes in the output buffer.
   */
  size_t GetOutputByteSize()
    {return this->Destination.NumBytes;};

  /**
   * @brief Copies the output data into the destination buffer.
   * @param dest pointer to the buffer of the output to the filter (out).
   * @pre dest != NULL
   * @pre dest must be of size this->GetOutputSize()
   */
  void GetOutputData(void* dest);

  /**
   * @brief Instantiates a concrete instance of the data filter given the name.
   * @param name a string representing the name of the filter.
   * @return F a concrete instance of the filter. F == NULL if no filter is
   * found by the given name.
   */
  static DataFilter* InstantiateByName(const char* name);

  /**
   * @brief Encodes the user-supplied data.
   * @note Implemented by concrete implementations of DataFilter.
   * @return rc return code, anything below 0 indicates an error.
   * @throw e exception if an error occurs.
   * @pre this->Destination.Data == NULL.
   */
  virtual int Encode()=0;

  /**
   * @brief Decodes the user-supplied data.
   * @note Implemented by concrete implementations of DataFilter.
   * @return rc return code, anything below 0 indicates an error.
   * @throw e exception if an error occurs.
   * @pre this->Destination.Data == NULL.
   */
  virtual int Decode()=0;

protected:
  std::string Name;       // user-friendly name for this filter instance
  DataBuffer Source;      // the source buffer (supplied by the application)
  DataBuffer Destination; // the destination buffer (computed internally)

  /**
   * @brief Resets the destination buffer iff it is not NULL.
   */
  void Reset();

private:
  GIO_DISABLE_COPY_AND_ASSIGNMENT(DataFilter);
};

} /* namespace gio */

#endif /* DATAFILTER_H_ */
