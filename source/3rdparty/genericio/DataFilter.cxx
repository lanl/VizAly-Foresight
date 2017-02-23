#include "DataFilter.h"

// C/C++ includes
#include <cassert>
#include <cstdlib>
#include <cstring>

// Concrete DataFilter instances
//#include "BloscCompression.h"

namespace gio
{

DataFilter::DataFilter()
{
  this->Source.Data     = NULL;
  this->Source.TypeSize = 0;
  this->Source.NumBytes = 0;

  this->Destination.Data     = NULL;
  this->Destination.TypeSize = 0;
  this->Destination.NumBytes = 0;
}

//-----------------------------------------------------------------------------
DataFilter::~DataFilter()
{
  this->Reset();
}

//-----------------------------------------------------------------------------
void DataFilter::SetInputData(void* src, size_t src_typesize, size_t src_nbytes)
{
  this->Source.Data     = src;
  this->Source.TypeSize = src_typesize;
  this->Source.NumBytes = src_nbytes;
}

//-----------------------------------------------------------------------------
void DataFilter::GetOutputData(void* dest)
{
  memcpy(dest,this->Destination.Data,this->Destination.NumBytes);
}

//-----------------------------------------------------------------------------
void DataFilter::Reset()
{
  if(this->Destination.Data != NULL)
    {
    free(this->Destination.Data);
    this->Destination.Data     = NULL;
    this->Destination.TypeSize = 0;
    this->Destination.NumBytes = 0;
    }
}

//-----------------------------------------------------------------------------
DataFilter* DataFilter::InstantiateByName(const char* name)
{
  DataFilter* f = NULL;

//  if( strcmp(name,"BLOSC")==0)
//    {
//    f = new BloscCompression();
//    }

  return( f );
}

} /* namespace gio */
