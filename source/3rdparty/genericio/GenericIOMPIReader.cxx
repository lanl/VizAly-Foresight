#include "GenericIOMPIReader.h"

#include "CRC64.h"
#include "GenericIOUtilities.h"

// C/C++ includes
#include <cassert>
#include <cstddef>
#include <cstdio>
#include <fstream>
#include <iostream>
#include <set>
#include <sstream>
#include <stdexcept>

namespace gio
{

//------------------------------------------------------------------------------
GenericIOMPIReader::GenericIOMPIReader()
{
  this->FH 			    = MPI_FILE_NULL;
  this->IOStrategy      = FileIOMPI;
}

//------------------------------------------------------------------------------
GenericIOMPIReader::~GenericIOMPIReader()
{
  this->Close();
}

//------------------------------------------------------------------------------
void GenericIOMPIReader::Open()
{
  // sanity checks
  assert("pre: non-null file Handle" && (this->FH == MPI_FILE_NULL) );
  assert("pre: empty filename!" && (!this->FileName.empty()) );

  int rc = MPI_File_open(
                this->Communicator,
                const_cast<char*>(this->FileName.c_str()),
                MPI_MODE_RDONLY,
                MPI_INFO_NULL,
                &this->FH);
  if( rc != MPI_SUCCESS )
    {
    throw std::runtime_error( "Unable to open file: " + this->FileName );
    }

  assert("post: file handle is NULL!" && (this->FH != MPI_FILE_NULL) );
}

//------------------------------------------------------------------------------
void GenericIOMPIReader::AllocateInternalReaders(const int N)
{
  assert("pre: Internal Readers should be NULL!" &&
     (this->InternalReaders==NULL) );
  assert("pre: NumberOfReaders(N) > 0" && (N > 0) );

  this->InternalReaders = new GenericIOReader*[N];
  assert("pre: Could not allocate internal readers array!" &&
      (this->InternalReaders != NULL));

  for( int i=0; i < N; ++i )
    {
  this->InternalReaders[ i ] = new GenericIOMPIReader();
    } // END for all readers
}

//------------------------------------------------------------------------------
void GenericIOMPIReader::Close()
{
  MPI_File_close(&this->FH);
  if( this->SplitMode )
    {
    for(int i=0; i < this->NumberOfFiles; ++i )
      {
      this->InternalReaders[ i ]->Close();
      } // END for all files
    } // END if in split mode
}

//------------------------------------------------------------------------------
void GenericIOMPIReader::Read(
        void *buf, size_t count, off_t offset, const std::string &varName)
{
  assert("pre: file handle is NULL!" && (this->FH != MPI_FILE_NULL) );

  while( count > 0 )
    {
    MPI_Status status;
    int rc = MPI_File_read_at(this->FH,offset,buf,count,MPI_BYTE,&status);
    if( rc != MPI_SUCCESS )
      {
      throw std::runtime_error(
          "Unable to read " + varName + " from file " + this->FileName);
      }

    int scount;
    MPI_Get_count(&status,MPI_BYTE,&scount);

    // Check if we somehow reached EOF
    if( scount == 0 )
      {
      throw std::runtime_error(
          "Unexpected EOF error when reading " + varName +
          " from file " + this->FileName
          );
      } // END if

    count -= scount;
    buf = ((char*) buf) + scount;
    offset += scount;
    } // END while
}

} /* namespace cosmotk */
