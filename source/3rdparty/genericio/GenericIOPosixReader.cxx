#include "GenericIOPosixReader.h"

#include <unistd.h>
#include <cassert>
#include <errno.h>
#include <fcntl.h>
#include <stdexcept>
#include <sys/stat.h>

namespace gio
{

GenericIOPosixReader::GenericIOPosixReader()
{
  this->FH 		   = -1;
  this->IOStrategy = FileIOPOSIX;
}

//------------------------------------------------------------------------------
GenericIOPosixReader::~GenericIOPosixReader()
{
  if(this->FH != -1)
   {
   this->Close();
   }
}

//------------------------------------------------------------------------------
void GenericIOPosixReader::Open()
{
  // sanity checks
  assert("pre: file handle should be invalid!" && (this->FH==-1) );
  assert("pre: empty filename!" && (!this->FileName.empty()) );

  int mode = S_IRUSR | S_IWUSR | S_IRGRP;
  this->FH = open(this->FileName.c_str(),O_RDONLY,mode);
  if(this->FH == -1 )
    {
  throw std::runtime_error("Unable to open file: " + this->FileName);
    }
}

//------------------------------------------------------------------------------
void GenericIOPosixReader::Read(
    void *buf, size_t count, off_t offset, const std::string &name)
{
  while(count > 0)
    {
    ssize_t scount;
    errno = 0;
    scount = pread(this->FH,buf,count,offset);
    if(scount == -1)
      {
      if( errno == EINTR )
        {
        continue;
        } // END if EINTR
      throw std::runtime_error("Unable to read variable " + name);
      } // END if scount == -1

    // Check if we somehow reached EOF
    if( scount==0 )
      {
      throw std::runtime_error(
               "Unexpected EOF error when reading " + name +
               " from file " + this->FileName
               );
      } // END if

    count  -= scount;
    buf     = ((char*) buf)+scount;
    offset += scount;
    } // END while
}

//------------------------------------------------------------------------------
void GenericIOPosixReader::Close()
{
  if( this->FH != -1 )
    {
    close(this->FH);
    }

  if( this->SplitMode )
    {
    for(int i=0; i < this->NumberOfFiles; ++i )
      {
      this->InternalReaders[ i ]->Close();
      } // END for all files
    } // END if split mode
}

//------------------------------------------------------------------------------
void GenericIOPosixReader::AllocateInternalReaders(const int N)
{
  assert("pre: Internal Readers should be NULL!" &&
     (this->InternalReaders==NULL) );
  assert("pre: NumberOfReaders(N) > 0" && (N > 0) );

  this->InternalReaders = new GenericIOReader*[N];
  assert("pre: Could not allocate internal readers array!" &&
        (this->InternalReaders != NULL));

  for(int i=0; i < N; ++i)
    {
  this->InternalReaders[ i ]= new GenericIOPosixReader();
    } // END for all readers
}

} /* namespace cosmotk */
