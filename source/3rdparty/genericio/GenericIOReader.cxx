#include "GenericIOReader.h"

#include "CRC64.h"
#include "DataFilter.h"
#include "GenericIOUtilities.h"

#include <cstddef>
#include <cstdio>
#include <cstring>
#include <iostream>
#include <set>
#include <sstream>
#include <stdexcept>
#include <string>

#ifndef MPI_UINT64_T
#define MPI_UINT64_T (sizeof(long) == 8 ? MPI_LONG : MPI_LONG_LONG)
#endif

namespace gio
{

//------------------------------------------------------------------------------
GenericIOReader::GenericIOReader()
{
  this->Communicator = MPI_COMM_NULL;
  this->NumRanks = 0;
  this->Rank = 0;
  this->SwapEndian = false;
  this->SplitMode = false;
  this->NumberOfFiles = 0;
  this->InternalReaders = NULL;
  this->InternalReader = false;
  this->BlockAssignmentStrategy = RR_BLOCK_ASSIGNMENT;

  this->VH.resize(0);
  this->RH.resize(0);
  this->AssignedBlocks.resize(0);
  this->BlockMetadata.resize(0);
  this->EntireHeader.resize(0);
  this->RankNeighbors.resize(0);

  this->IJKMapToGlobalBlockIdx = NULL;

  this->DoAcquireGlobalDimensions = true;
  this->DoAcquireIJKMapping       = true;
  this->DoBlockAssignment         = true;
  this->DoDetermineFileType       = true;
  this->DoReadBlockHeaders        = true;
  this->DoSetupInternalReaders    = true;

  this->SingleBlockMode = false;

  this->HasBlockMetaData = false;

  this->FunctionDescription = "MASTER";
}

//------------------------------------------------------------------------------
GenericIOReader::~GenericIOReader()
{
   this->ClearInternalReaders();
   this->VH.clear();
   this->RH.clear();
   this->AssignedBlocks.clear();
   this->EntireHeader.clear();
   this->RankNeighbors.clear();

   if( this->IJKMapToGlobalBlockIdx != NULL && !this->InternalReader )
     {
     this->IJKMapToGlobalBlockIdx->clear();
     delete this->IJKMapToGlobalBlockIdx;
     this->IJKMapToGlobalBlockIdx = NULL;
     }
}

//------------------------------------------------------------------------------
void GenericIOReader::GetRankNeighbors(RankNeighbor* neighbors)
{
  if( neighbors == NULL )
    {
    return;
    }

  for(size_t i=0; i < this->RankNeighbors.size(); ++i)
    {
    neighbors[i].RankID = this->RankNeighbors[i].RankID;
    memcpy(neighbors[i].Orient,this->RankNeighbors[i].Orient,3*sizeof(int));
    }
}

//------------------------------------------------------------------------------
void GenericIOReader::DetermineFileType()
{
  if( !this->DoDetermineFileType )
    {
    return;
    }

  if( this->Rank == 0 )
    {
    if(this->HasVariable("$rank") && this->HasVariable("$partition"))
      {
      this->SplitMode = true;
      }

    int split = (this->SplitMode)? 1 : 0;
    MPI_Bcast(&split,1,MPI_INT,0,this->Communicator);
    if( this->SplitMode )
      {
      this->ReadBlockToFileMap();
      } // END if splitmode
    } // END if rank 0
  else
    {
    int split = -1;
    MPI_Bcast(&split,1,MPI_INT,0,this->Communicator);
    this->SplitMode = (split==1)? true : false;
    if( this->SplitMode )
      {
      this->ReadBlockToFileMap();
      } // END if splitmode
    } // END else

}

//------------------------------------------------------------------------------
void GenericIOReader::IndexVariables()
{
  std::string name;
  for(unsigned int i=0; i < this->VH.size(); ++i)
    {
    name = std::string( this->VH[i].Name );
    this->VariableName2IndexMap[ name ]= i;
    }
}

//------------------------------------------------------------------------------
int GenericIOReader::GetVariableIndex(const std::string name)
{
  int idx = -1;
  if( this->VariableName2IndexMap.find(name) !=
       this->VariableName2IndexMap.end())
    {
    idx = this->VariableName2IndexMap[name];
    }
  return( idx );
}

//------------------------------------------------------------------------------
void GenericIOReader::AssignBlock(const int blkIdx)
{

  assert("pre: block index is out of bounds" &&
          (blkIdx >= 0) && (blkIdx < static_cast<int>(this->GH.NRanks)) );
  this->AssignedBlocks.push_back( blkIdx );
  if( this->InternalReader )
    {
    if( this->HasBlock(blkIdx) )
      {
      std::cerr << "WARNING: assigning duplicate block " << blkIdx << std::endl;
      }
    this->GlobalBlockIdx2LocalIdx[blkIdx] = this->AssignedBlocks.size()-1;
    }

}

//------------------------------------------------------------------------------
void GenericIOReader::AssignSingleBlock(const RankHeader blockHeader)
{
  this->RH.clear();
  this->AssignedBlocks.clear();
  this->GlobalBlockIdx2LocalIdx.clear();

  this->RH.push_back( blockHeader );
  this->AssignedBlocks.push_back( blockHeader.GlobalRank );
  this->GlobalBlockIdx2LocalIdx[ blockHeader.GlobalRank ] = 0;
}

//------------------------------------------------------------------------------
VariableInfo GenericIOReader::GetFileVariableInfo(const int i)
{
  assert("pre: variable index is out-of-bounds!" &&
          (i >= 0) && (i < static_cast<int>(this->VH.size())) );

  VariableInfo VI(
      std::string( this->VH[ i ].Name ),
      this->VH[ i ].Size,
      static_cast<bool>(this->VH[i].Flags & FloatValue),
      static_cast<bool>(this->VH[i].Flags & SignedValue),
      static_cast<bool>(this->VH[i].Flags & ValueIsPhysCoordX),
      static_cast<bool>(this->VH[i].Flags & ValueIsPhysCoordY),
      static_cast<bool>(this->VH[i].Flags & ValueIsPhysCoordZ),
      static_cast<bool>(this->VH[i].Flags & ValueMaybePhysGhost)
      );
  return( VI );
}

//------------------------------------------------------------------------------
void GenericIOReader::ClearInternalReaders()
{
 if( !this->SplitMode )
  {
  return;
  }

 for(int i=0; i < this->NumberOfFiles; ++i)
   {
   if( this->InternalReaders[i] != NULL )
     {
     delete this->InternalReaders[i];
     }
   }
 delete [] this->InternalReaders;
}

//------------------------------------------------------------------------------
void GenericIOReader::AcquireGlobalDimensions()
{
  // STEP 0: short-circuit if we are not in split mode or if the
  // caller explicitely turned off this operation.
  if( !this->SplitMode || !this->DoAcquireGlobalDimensions )
    {
    /* the global dimensions are already stored in GH.Dims */
    return;
    }

  // Sanity checks!
  assert("pre: number of files must be at least 1" &&
          (this->NumberOfFiles > 0) );
  assert("pre: FileIDToSymbolicID not properly set!" &&
   (this->NumberOfFiles == static_cast<int>(this->FileIDToSymbolicID.size())));

  // STEP 2:Temporary variables
  std::ostringstream oss;     // used to construct file name
  GenericIOReader *r = NULL;  // temporary internal reader
  uint64_t dims[3];           // buffer to store global dimensions
  double originAndScale[6];   // buffer to store origin and scale

  // STEP 3: Get the global dimensions
  if( this->Rank == 0 )
    {
    /* Rank 0 reads one of the files and broadcasts the dimensions */
    // STEP 3.1: Create the file name for one of the separate files
    oss.clear();
    oss << this->FileName << "#" << this->FileIDToSymbolicID[ 0 ];

    // STEP 3.2: Create the reader
    r = this->GetNewInstance();
    assert("pre: internal reader is NULL!" && (r != NULL) );

    // STEP 3.3: Setup the reader
    r->SetCommunicator(MPI_COMM_SELF);
    r->SetFileName( oss.str() );
    r->SetFunctionDescription("GlobalDimensionsReader");
    r->SetDoReadBlockHeaders(false);

    // STEP 3.4: Read the global header & get the global dimensions
    r->OpenAndReadHeader();
    r->GetGlobalDimensions( dims );
    r->GetPhysOrigin(&originAndScale[0]);
    r->GetPhysScale(&originAndScale[3]);

    // STEP 3.5: Close & cleanup
    r->Close();
    delete r;
    r=NULL;

    // STEP 3.6: Broadcast dimensions & originAndScale info
    MPI_Bcast(dims,3,MPI_UINT64_T,0,this->Communicator);
    MPI_Bcast(originAndScale,6,MPI_DOUBLE,0,this->Communicator);
    }
  else
    {
    MPI_Bcast(dims,3,MPI_UINT64_T,0,this->Communicator);
    MPI_Bcast(originAndScale,6,MPI_DOUBLE,0,this->Communicator);
    }

  assert("pre: temporary reader should be NULL!" && (r==NULL));

  // STEP 4: Copy the global dimensions to the global header
  memcpy(this->GH.Dims,dims,3*sizeof(uint64_t));

  // STEP 5: Copy origin
  memcpy(this->GH.PhysOrigin,originAndScale,3*sizeof(double));

  // STEP 6: Copy scale
  memcpy(this->GH.PhysScale,&originAndScale[3],3*sizeof(double));

  // STEP 5: Barrier synchronization
  this->Barrier();
}

//------------------------------------------------------------------------------
bool GenericIOReader::IsIJKSorted()
{
  bool status = true;
  if( (this->GH.Dims[0]==this->GH.NRanks) &&
      (this->GH.Dims[1]==1) &&
      (this->GH.Dims[2]==1))
    {
    // the data was not written using a cartesian communicator, hence,
    // the blocks are not IJK sorted.
    status = false;
    }
  return( status );
}

//------------------------------------------------------------------------------
void GenericIOReader::AcquireIJKMapping()
{
  // STEP 0: short-circuit here if the blocks are not IJK sorted
  // or if the caller explicitly turned off this operation
  if( !this->IsIJKSorted() || !this->DoAcquireIJKMapping )
    {
    return;
    }

  this->IJKMapToGlobalBlockIdx = new ProcessCoordsToRank;

  // STEP 1: decleare some local variables
  int numBlocks = 0;
  std::vector<uint64_t> flatData;
  GenericIOReader *r = NULL;

  // STEP 2: Read in all the block headers in process 0, extract the ijk &
  // global rank for each block and broadcast to all other processes.
  if(this->Rank == 0)
    {
    // Get a new reader instance to read only the block headers
    r = this->GetNewInstance();
    r->SetFunctionDescription(this->FunctionDescription+".IJKReader");
    r->SetCommunicator(MPI_COMM_SELF);
    r->SetFileName(this->FileName);
    r->SetDoAcquireIJKMapping(false);


    // Read the header (i.e.,global header, variable header and block header)
    r->OpenAndReadHeader();


    // Extract the ijk and global rank of each header
    numBlocks = r->GetNumberOfBlockHeaders();
    flatData.resize(4*numBlocks);

    for(int block=0; block < r->GetNumberOfBlockHeaders(); ++block)
      {
      RankHeader blockHeader = r->GetBlockHeader(block);
      flatData[block*4]      = blockHeader.GlobalRank;
      flatData[block*4+1]    = blockHeader.Coords[0];
      flatData[block*4+2]    = blockHeader.Coords[1];
      flatData[block*4+3]    = blockHeader.Coords[2];
      } // END for all blocks

    // Deallocate temporary reader
    delete r;
    r = NULL;

    // Broadcast number of blocks
    MPI_Bcast(&numBlocks,1,MPI_INT,0,this->Communicator);

    // Broadcast the rank & ijk vectors
    MPI_Bcast(&flatData[0],numBlocks*4,MPI_UINT64_T,0,this->Communicator);
    }
  else
    {
    // Broadcast number of blocks
    MPI_Bcast(&numBlocks,1,MPI_INT,0,this->Communicator);

    // Allocate arrays
    flatData.resize(4*numBlocks);

    // Broadcast the rank & ijk vectors
    MPI_Bcast(&flatData[0],numBlocks*4,MPI_UINT64_T,0,this->Communicator);
    }

  // STEP 3: Construct the IJK mapping
  for(int block=0; block < numBlocks; ++block)
    {
    uint64_t i = flatData[ block*4+1 ];
    uint64_t j = flatData[ block*4+2 ];
    uint64_t k = flatData[ block*4+3 ];

    uint64_t idx = flatData[ block*4 ];

    std::string hashCode = GenericIOUtilities::GetHashCodeFromCoords(i,j,k);
    assert("pre: duplicate element in IJKMap!" &&
        this->IJKMapToGlobalBlockIdx->find(hashCode)==
            this->IJKMapToGlobalBlockIdx->end());

    this->IJKMapToGlobalBlockIdx->insert(
        std::pair<std::string,uint64_t>(hashCode,idx) );
    } // END for all blocks

  // STEP 4: Barrier synchronization
  this->Barrier();
}

//------------------------------------------------------------------------------
void GenericIOReader::SetupInternalReaders()
{
  // STEP 0: Short-circuit if we are reading a single file or if
  // the caller has explicitly turned off this operation.
  if( !this->SplitMode || !this->DoSetupInternalReaders )
    {
    return;
    }

  // STEP 1: Construct all readers, based on number of separate files that
  // we have to read.
  this->AllocateInternalReaders(this->NumberOfFiles);

  std::ostringstream oss;
  for(int i=0; i < this->NumberOfFiles; ++i)
    {
    oss.clear(); oss.str("");
    oss << this->FileName << "#" << this->FileIDToSymbolicID[ i ];
    GenericIOReader *iReader = this->InternalReaders[ i ];
    assert("pre: internal reader should not be NULL!" && (iReader != NULL) );
    iReader->SetCommunicator( this->Communicator );
    iReader->SetFunctionDescription("InternalDataReader");
    iReader->SetFileName(oss.str());
    iReader->SetReaderInternal();
    iReader->SetDoBlockAssignment(false);
    iReader->SetDoReadBlockHeaders(false);
    iReader->SetDoAcquireGlobalDimensions(false);
    iReader->SetDoAcquireIJKMapping(false);
    iReader->OpenAndReadHeader();
    assert("internal reader should not be in SplitMode" &&
           !iReader->IsSplitMode());

    iReader->SetDoReadBlockHeaders(true); /* turn this back on to call later */

    iReader->SetIJKMapToGlobalBlockIdx(this->IJKMapToGlobalBlockIdx);

    /* Manually prescribe block assignement for this reader */
    this->InternalReaders[ i ] = iReader;
    assert("manually assign blocks to internal reader!" &&
            iReader->GetNumberOfAssignedBlocks()==0);
    }

  // STEP 2: Prescribe block assignment. Note, block assignment is done on
  // the total number of blocks, in OpenAndReadHeader. The total number of
  // blocks is implicitly given by the total number of elements read from
  // the block-to-file map, see ReadBlockToFileMap() method.
  for( unsigned int block=0; block < this->AssignedBlocks.size(); ++block)
    {
    int globalBlockIdx = this->AssignedBlocks[ block ];

    // Sanity checks!
    assert("ERROR: Cannot map block to file!" &&
     (this->BlockToFileMap.find(globalBlockIdx)!=this->BlockToFileMap.end()));
    assert("ERROR: Cannot map block to idx within file!" &&
     (this->BlockToIdxWithinFile.find(globalBlockIdx)!=
      this->BlockToIdxWithinFile.end()));

    int fileIdx =
       this->SymbolicIDToFileID[this->BlockToFileMap[ globalBlockIdx ]];

    // Get block Idx within file
    int idxWithinFile = this->BlockToIdxWithinFile[ globalBlockIdx ];

    this->InternalReaders[ fileIdx ]->AssignBlock(idxWithinFile);
    } // END for all assigned blocks

  // STEP 3: Proxy variable headers to proxy reader

  // NOTE: we must first clear VH here since it has been populated with
  // the variable headers of the map file, i.e., the $rank, $partition
  // variables.
  this->VH.clear();
  for(int i=0;i < this->InternalReaders[0]->GetNumberOfVariablesInFile();++i)
    {
    this->VH.push_back(this->InternalReaders[0]->GetVariableHeader(i) );
    } // END for all variables

  // STEP 4: Re-index variables s.t. HasVariable will work properly
  this->IndexVariables();

}

//------------------------------------------------------------------------------
void GenericIOReader::OpenAndReadHeader()
{
  // sanity checks
  assert( "pre: No communicator is supplied" &&
          (this->Communicator != MPI_COMM_NULL) );
  assert("pre: FileName is empty!" && (!this->FileName.empty()) );

  // STEP 0: Get Rank and NumRanks information
  MPI_Comm_rank(this->Communicator, &this->Rank);
  MPI_Comm_size(this->Communicator, &this->NumRanks);

  // STEP 1: Open file, if not successful, throw runtime_error exception
  this->Open();

  // STEP 2: Read Global & Variable header
  this->ReadHeader();

  // STEP 3: Index the variables
  this->IndexVariables();

  // STEP 4: Detect file type, i.e., if the data is written to separate files
  // or if it is in a single monolithic file.
  this->DetermineFileType();

  // STEP 5: Distribute blocks from the file to processes according to
  // the user-supplied assignment policy.
  if( this->SplitMode )
    {
    // STEP 5.1: Acquire the global dimensions, since the map file does not
    // have the global dimensions.
    this->AcquireGlobalDimensions();

    // STEP 5.2: Acquire the global ijk-to-rank mapping
    this->AcquireIJKMapping();

    // STEP 5.3: Distribute blocks
    if( this->DoBlockAssignment )
      {
      GenericIOUtilities::BlockAssignment(
          this->Rank,this->NumRanks,
          this->GH.NRanks,
          this->GH.Dims,
          this->BlockAssignmentStrategy,
          this->IJKMapToGlobalBlockIdx,
          this->AssignedBlocks,
          this->RankNeighbors
          );
      this->ConstructGlobalToLocalBlockMapping();
      } // END if BlockAssignment

    // STEP 5.4: Setup internal readers
    this->SetupInternalReaders();
    } // END if SplitMode
  else
    {

    this->AcquireIJKMapping();

    // STEP 5.1: Distribute blocks
    if( this->DoBlockAssignment )
      {
      GenericIOUtilities::BlockAssignment(
         this->Rank,this->NumRanks,
         this->GH.NRanks,
         this->GH.Dims,
         this->BlockAssignmentStrategy,
         this->IJKMapToGlobalBlockIdx,
         this->AssignedBlocks,
         this->RankNeighbors
         );
      this->ConstructGlobalToLocalBlockMapping();
      } // END if BlockAssignment

    } // END else

  // STEP 6: Read the headers of the assigned blocks
  this->ReadBlockHeaders();

  // STEP 7: Reads the block metadata
  this->ReadBlockMetaData();

  // STEP 7: Barrier synchronization
  this->Barrier();
}

//------------------------------------------------------------------------------
void GenericIOReader::ConstructGlobalToLocalBlockMapping()
{
  for(unsigned int i=0; i < this->AssignedBlocks.size(); ++i)
    {
    int globalIdx = this->AssignedBlocks[ i ];
    this->GlobalBlockIdx2LocalIdx[ globalIdx ] = i;
    } // END for all assigned blocks
}

//------------------------------------------------------------------------------
void GenericIOReader::ReadHeader()
{
 // Internal attributes. Each element in the `attribs` array represents a
 // particular attribute, e.g., whether we should swap endian, etc., as
 // indicated below. The reason for storing these attributes in an array
 // instead of individual ints is so that we can send all these attributes
 // with a single broadcast to all ranks.
 int attribs[3]=
   { 0, // indicates whether to swap or not
     0, // indicates the entire header size,including the CRC checksum
     0  // indicates whether an error occured
   };

 // Integers corresponding to indices in the `attribs` array for
 const int SWAP        = 0;
 const int HEADER_SIZE = 1;
 const int ERROR       = 2;

 // Read in attributes
 if( this->Rank==0 )
   {
   // Read the global header
   this->Read(&this->GH,sizeof(GlobalHeader),0,"GlobalHeader");

   // Byte-swap header if necessary
   if( !GenericIOUtilities::DoesFileEndianMatch(&this->GH) )
    {
    GenericIOUtilities::SwapGlobalHeader(&this->GH);
    this->SwapEndian = true;
    attribs[SWAP]    = 1;
    }
   else
    {
    this->SwapEndian = false;
    attribs[SWAP]    = 0;
    }

   // Read entire header & its checksum
   attribs[HEADER_SIZE] = this->GH.HeaderSize+CRCSize;
   this->EntireHeader.resize(attribs[HEADER_SIZE],0xFE/* poison */);
   this->Read(&this->EntireHeader[0],attribs[HEADER_SIZE],0,"EntireHeader");

   // header checksum -- CRC is endian independent. It must be verified
   // before byte-swapping.
   if(!GenericIOUtilities::VerifyChecksum(
      &this->EntireHeader[0],attribs[HEADER_SIZE]))
     {
     std::cerr << "\nWARNING: header checksum verification failed @"
               << __FILE__ << ":" << __LINE__ << std::endl;
     attribs[ERROR] = 1;
    }

   MPI_Bcast(attribs,3,MPI_INT,0,this->Communicator);
   }
 else
   {
   MPI_Bcast(attribs,3,MPI_INT,0,this->Communicator);
   this->SwapEndian = (attribs[SWAP]==1)? true : false;
   this->EntireHeader.resize( attribs[HEADER_SIZE] );
   }

 // Check for errors
 if(attribs[ERROR]==1)
   {
   std::cerr << "\nWARNING: header checksum verification failed @"
              << __FILE__ << ":" << __LINE__ << std::endl;
   throw std::runtime_error("Header CRC checksum failed!");
   }

 // Broadcast the raw bytes of the entire header
 assert("pre: headers has not been properly allocated!" &&
       (static_cast<int>(this->EntireHeader.size())==attribs[HEADER_SIZE]) );
 MPI_Bcast(
   &this->EntireHeader[0],attribs[HEADER_SIZE],MPI_CHAR,0,this->Communicator);

 // Ensure broadcast of the header was successful
 if(crc64_omp(&this->EntireHeader[0],attribs[HEADER_SIZE])!=(uint64_t)-1)
   {
   std::cerr << "\nWARNING: header checksum verification failed @"
              << __FILE__ << ":" << __LINE__ << std::endl;
   attribs[ERROR] = 1;
   }
 int errors = 0;
 MPI_Allreduce(
     &attribs[ERROR],&errors,1,MPI_INT,MPI_SUM,this->Communicator);
 if(errors > 0 )
   {
   throw std::runtime_error("Error broadcasting header to ranks!");
   }

 // Extract the global header
 this->GH = *(GlobalHeader*)(&this->EntireHeader[0]);
 if( this->SwapEndian )
   {
   GenericIOUtilities::SwapGlobalHeader(&this->GH);
   }

 // Read the variable headers
 this->ReadVariableHeaders();
 this->Barrier();
}

//------------------------------------------------------------------------------
void GenericIOReader::ReadBlockMetaData()
{
  // STEP 0: short-circuit iff indicated to skip reading block headers.
  if( !this->DoReadBlockHeaders )
    {
    return;
    }

  // STEP 1: Check if block metadata is available in the file
  if( (offsetof(GlobalHeader,BlocksStart) < this->GH.GlobalHeaderSize) &&
      (this->GH.BlocksSize > 0) )
    {
    this->HasBlockMetaData = true;
    } // END if there are block metadata
  else
    {
    // There are no block metadata available
    return;
    } // END else no block metadata

  // STEP 2: Allocate storage for block metadata. Each variable within each
  // block has a corresponding block metadata instance.
  size_t nblocks = this->AssignedBlocks.size();
  size_t nvars   = this->GH.NVars;

  // STEP 3: Read the block metadata
  if( this->SplitMode )
    {
    // In splitmode, we proxy the block requests to each corresponding reader.
    // Get the correct value for nvars from an internal reader
    nvars = this->InternalReaders[ 0 ]->GetNumberOfVariablesInFile();
    this->BlockMetadata.resize(nblocks*nvars);

    for(int fileIdx=0; fileIdx < this->NumberOfFiles; ++fileIdx)
      {
      GenericIOReader* r = this->InternalReaders[ fileIdx ];
      r->ReadBlockMetaData();
      int numBlocks = r->GetNumberOfBlockHeaders();

      for(int blockIdx=0; blockIdx < numBlocks; ++blockIdx)
        {
        RankHeader blockInfo = r->GetBlockHeader(blockIdx);
         assert("pre: reading an unassigned block header!" &&
                 this->HasBlock(blockInfo.GlobalRank));
         int localIdx = this->GlobalBlockIdx2LocalIdx[blockInfo.GlobalRank];
         assert("pre: invalid local block index!" &&
            (localIdx >= 0) && (localIdx < static_cast<int>(this->RH.size())));

         for(uint64_t varIdx=0; varIdx < nvars; ++varIdx)
           {
           this->BlockMetadata[localIdx*nvars+varIdx]=
               r->GetBlockMetaData(blockIdx,varIdx);
           } // END for all variables
        } // END for all blocks

      } // END for all files
    } // END if data spans multiple files
  else
    {
    this->BlockMetadata.resize(nblocks*nvars);
    for(unsigned int blk=0; blk < this->AssignedBlocks.size(); ++blk)
      {
      int blkIdx = this->AssignedBlocks[ blk ];
      for(unsigned int varIdx=0; varIdx < nvars; ++varIdx)
        {
        this->ReadBlockVariableMetaData(
            blkIdx,varIdx,this->BlockMetadata[blk*nvars+varIdx]);
        } // END for all variables
      } // END for all blocks

    } // END else single file

}

//------------------------------------------------------------------------------
void GenericIOReader::ReadBlockHeaders()
{
  if( !this->DoReadBlockHeaders )
    {
    return;
    }

  this->RH.resize( this->AssignedBlocks.size() );
  if( this->SplitMode )
    {

    for(int fileIdx=0; fileIdx < this->NumberOfFiles; ++fileIdx)
      {
      this->InternalReaders[ fileIdx ]->ReadBlockHeaders();

      // proxy block headers from internal readers to proxy reader
      int numBlocks =
          this->InternalReaders[ fileIdx ]->GetNumberOfBlockHeaders();

      for(int blockIdx=0; blockIdx < numBlocks; ++blockIdx)
        {
        RankHeader blockInfo =
            this->InternalReaders[fileIdx]->GetBlockHeader(blockIdx);
        assert("pre: reading an unassigned block header!" &&
                this->HasBlock(blockInfo.GlobalRank));
        int localIdx = this->GlobalBlockIdx2LocalIdx[blockInfo.GlobalRank];
        assert("pre: invalid local block index!" &&
           (localIdx >= 0) && (localIdx < static_cast<int>(this->RH.size())) );
        this->RH[ localIdx ] = blockInfo;
        } // END for all block headers of this file
      } // END for all files

    } // END if
  else
    {

    for(unsigned int blk=0; blk < this->AssignedBlocks.size(); ++blk)
      {
      int blkIdx = this->AssignedBlocks[ blk ];
      this->ReadBlockHeader(blkIdx,this->RH[blk]);
      } // END for all assigned blocks

    } // END else

  assert("post: number of block headers must match num assigned blocks!" &&
         ( this->AssignedBlocks.size()==this->RH.size() ) );
}

//------------------------------------------------------------------------------
void GenericIOReader::ReadBlockHeader(
        const int blkIdx, RankHeader& blockHeader)
{
  uint64_t offSet = this->GH.RanksStart + blkIdx*sizeof(RankHeader);
  assert("pre: detected rank header offset out-of-bounds!" &&
             offSet < this->EntireHeader.size()-CRCSize );

   // Copy the bytes of the variable header from the raw header data
   memcpy(&blockHeader,&this->EntireHeader[offSet],sizeof(RankHeader));

  if(this->SwapEndian)
    {
    GenericIOUtilities::SwapRankHeader(&blockHeader);
    }
}

//------------------------------------------------------------------------------
void GenericIOReader::ReadBlockVariableMetaData(
        const int blkIdx, const int varIdx, BlockInfo& blockMetaData)
{
  assert("pre: var index out of bounds!" &&
          (varIdx >= 0) && (varIdx < this->GH.NVars)  );

  uint64_t offSet =
    this->GH.BlocksStart + (blkIdx*this->GH.NVars+varIdx)*sizeof(BlockInfo);

  assert("pre: detected block metadata offset out-of-bounds!" &&
               offSet < this->EntireHeader.size()-CRCSize );

  memcpy(&blockMetaData,&this->EntireHeader[offSet],sizeof(BlockInfo));

  if(this->SwapEndian)
    {
    GenericIOUtilities::SwapBlockInfo(&blockMetaData);
    }
}

//------------------------------------------------------------------------------
bool GenericIOReader::IsSpatiallyDecomposed()
{
  bool status = false;
  uint64_t numRanks = this->GH.Dims[0]*this->GH.Dims[1]*this->GH.Dims[2];
  if( numRanks == this->GH.NRanks )
    {
    status = true;
    }
  return( status );
}

//------------------------------------------------------------------------------
void GenericIOReader::GetBlockBounds(const int i, double min[3], double max[3])
{
  assert("pre: block header index is out of bounds!" &&
           (i >= 0) && (i < static_cast<int>(this->RH.size())) );

  if( !this->IsSpatiallyDecomposed() )
    {
    throw std::runtime_error(
      "Called GetBlockBounds() but dataset is not spatially decomposed!");
    }

  for(int dim=0; dim < 3; ++dim)
    {
    double h = (this->GH.PhysScale[dim]-this->GH.PhysOrigin[dim]) /
                static_cast<double>( this->GH.Dims[dim] );
    min[dim] = this->GH.PhysOrigin[dim] + h*this->RH[i].Coords[dim];
    max[dim] = this->GH.PhysOrigin[dim] + h*(this->RH[i].Coords[dim]+1);
    } // END for all dimensions
}

//------------------------------------------------------------------------------
void GenericIOReader::ReadVariableHeader(
        const int idx, VariableHeader& vh)
{
  uint64_t offSet = this->GH.VarsStart + idx*sizeof(VariableHeader);
  assert("pre: detected variable header offset out-of-bounds!" &&
            offSet < this->EntireHeader.size()-CRCSize );

  // Copy the bytes of the variable header from the raw header data
  memcpy(&vh,&this->EntireHeader[offSet],sizeof(VariableHeader));

  if(this->SwapEndian)
    {
    GenericIOUtilities::SwapVariableHeader(&vh);
    }
}

//------------------------------------------------------------------------------
void GenericIOReader::ReadVariableHeaders()
{
  assert( "pre: file has no variables!" && (this->GH.NVars > 0) );

  this->VH.resize( this->GH.NVars );
  for(unsigned int i=0; i < this->GH.NVars; ++i )
    {
    this->ReadVariableHeader(i, this->VH[i] );
    } // END for all variables
}

//------------------------------------------------------------------------------
void GenericIOReader::ReadBlockToFileMap()
{
  assert("pre: file must be in SplitMode" && this->SplitMode);

  // STEP 0: Rank 0 reads block-to-file mapping and distributes it to all ranks
  int NumElements=0;
  int *rank = NULL;
  int *part = NULL;
  GenericIOReader *r = NULL;

  if(this->Rank==0)
    {
    r = this->GetNewInstance();
    r->SetCommunicator(MPI_COMM_SELF);
    r->SetDoDetermineFileType(false);
    r->SetFileName(this->FileName);
    r->OpenAndReadHeader();

    NumElements = r->GetNumberOfElements();
    MPI_Bcast(&NumElements,1,MPI_INT,0,this->Communicator);
    rank = new int[NumElements];
    part = new int[NumElements];
    r->AddVariable("$rank",rank);
    r->AddVariable("$partition",part);
    r->ReadData();

    MPI_Bcast(rank,NumElements,MPI_INT,0,this->Communicator);
    MPI_Bcast(part,NumElements,MPI_INT,0,this->Communicator);
    delete r;
    r=NULL;
    }
  else
    {
    MPI_Bcast(&NumElements,1,MPI_INT,0,this->Communicator);
    rank = new int[NumElements];
    part = new int[NumElements];
    MPI_Bcast(rank,NumElements,MPI_INT,0,this->Communicator);
    MPI_Bcast(part,NumElements,MPI_INT,0,this->Communicator);
    }

  // Correct global header b/c the map file has only a single block!
  this->GH.NRanks = NumElements;

  // STEP 1: All processes
  std::set<int> parts;
  std::map<int,int> counter;
  int blkIdxInFile = -1;
  for(int i=0; i < NumElements; ++i)
    {
    parts.insert(part[i]);
    this->BlockToFileMap[ rank[i] ] = part[i];

    blkIdxInFile = -1;
    if( counter.find(part[i]) != counter.end() )
      {
      blkIdxInFile    = counter[ part[i] ];
      counter[ part[i] ] += 1;
      } // END if
    else
      {
      blkIdxInFile = 0;
      counter[ part[i] ] = 1;
      } // END else
    this->BlockToIdxWithinFile[ rank[i] ] = blkIdxInFile;
    } // END for all elements
  this->NumberOfFiles = parts.size();

  // STEP 2: Compose FileIDToSymbolicID mapping
  std::set<int>::iterator it = parts.begin();
  for(int fileId=0; it != parts.end(); ++it, ++fileId)
    {
    this->FileIDToSymbolicID[fileId] = *it;
    this->SymbolicIDToFileID[ *it ]  = fileId;
    }

  // STEP 3: Delete dynamically allocated data
  counter.clear();
  parts.clear();
  if( rank != NULL )
    {
    delete [] rank;
    }
  if( part != NULL )
    {
    delete [] part;
    }
}

//------------------------------------------------------------------------------
bool GenericIOReader::ReadBlock(const int globalBlockIdx)
{
  this->SingleBlockMode = true;

  bool status = false;

  this->ClearBlockAssignment();
  assert("pre: assigned blocks must be 0!" &&
          (this->AssignedBlocks.size()==0));

  // Clear the block assignment of all the internal readers
  for(int i=0; i < this->NumberOfFiles; ++i)
    {
    assert("pre: Internal reader should not be NULL!" &&
                 (this->InternalReaders[ i ] != NULL) );
    this->InternalReaders[ i ]->ClearBlockAssignment();
    assert("pre: assigned blocks must be 0!" &&
            (this->InternalReaders[i]->GetNumberOfAssignedBlocks()==0));
    this->InternalReaders[ i ]->SetSingleBlockMode( true );

    } // END for all files

  if( this->HasBlock(globalBlockIdx))
    {
    status = true;
    this->AssignBlock(globalBlockIdx);

    if( this->SplitMode )
      {
      // Assign block to corresponding internal reader
      // Sanity checks!
      assert("ERROR: cannot find global-to-local mapping!" &&
         this->GlobalBlockIdx2LocalIdx.find(globalBlockIdx)!=
             this->GlobalBlockIdx2LocalIdx.end());
      assert("ERROR: Cannot map block to file!" &&
       (this->BlockToFileMap.find(globalBlockIdx)!=this->BlockToFileMap.end()));
      assert("ERROR: Cannot map block to idx within file!" &&
       (this->BlockToIdxWithinFile.find(globalBlockIdx)!=
        this->BlockToIdxWithinFile.end()));

      int localIdx = this->GlobalBlockIdx2LocalIdx[ globalBlockIdx ];
      // Get the file index for this block
      int fileIdx =
         this->SymbolicIDToFileID[this->BlockToFileMap[ globalBlockIdx ]];
//
//      // Get block Idx within file
//      int idxWithinFile = this->BlockToIdxWithinFile[ globalBlockIdx ];

      // Assign block to the internal reader
      this->InternalReaders[fileIdx]->AssignSingleBlock(
          this->GetBlockHeader(localIdx));
      this->InternalReaders[fileIdx]->SetSingleBlockMode(true);
      }
    this->ReadData();
    }
  return( status );
}

//------------------------------------------------------------------------------
void GenericIOReader::ReadData()
{
  if( this->SplitMode )
    {
    this->ReadSplitModeData();
    }
  else
    {
    this->ReadSingleFileData();
    }
}

//------------------------------------------------------------------------------
void GenericIOReader::ReadSplitModeData()
{
  // STEP 0: Propagate variables to internal readers
  int nskip = 0;
  for(int i=0; i < this->NumberOfFiles; ++i)
    {
    assert("pre: Internal reader should not be NULL!" &&
             (this->InternalReaders[ i ] != NULL) );

    this->InternalReaders[ i ]->ClearVariables();

    // Determine offset in to variable array to write
    if( i > 0 )
      {
      nskip += this->InternalReaders[ i-1 ]->GetNumberOfElements();
      }

    for(unsigned int varIdx=0; varIdx < this->Vars.size(); ++varIdx)
      {
      // sanity check!
      assert("pre: Data for variable is NULL!" &&
          (this->Vars[varIdx].Data != NULL) );

      // Get the variable size
      size_t vsize = this->Vars[varIdx].Size;

      // Get pointer where reader "i" will start filling in data for this
      // variable
      void *dataPtr =
          static_cast<char*>(this->Vars[varIdx].Data)+(nskip*vsize);
      assert("pre: dataPtr is NULL!" && (dataPtr != NULL) );

      if(this->Vars[varIdx].HasExtraSpace)
        {
        this->InternalReaders[ i ]->AddVariable(
            this->GetVariableInfo(varIdx),dataPtr,
            GenericIOBase::ValueHasExtraSpace);
        } // END if extra space
      else
        {
        this->InternalReaders[ i ]->AddVariable(
          this->GetVariableInfo(varIdx),dataPtr);
        } // END else extra space

      } // END for all variables

    // Read all the data from this file
    this->InternalReaders[ i ]->ReadData();
    } // END for all files
}

//------------------------------------------------------------------------------
void GenericIOReader::ReadSingleFileData()
{
  // STEP 0: Get the total number of elements in an array for this rank
  // TODO: Need to refactor and clean up the implementation a bit
  int N = 0;
  if( this->SingleBlockMode )
    {
    assert("pre: expecting to have at most a single block assigned!" &&
      ((this->AssignedBlocks.size()==1) || (this->AssignedBlocks.size()==0)));

    if( this->AssignedBlocks.size() == 1 )
      {
      int globalIdx = this->AssignedBlocks[ 0 ];
      if( this->InternalReader )
        {
        N = this->GetNumberOfElementsForBlock( 0 );
        }
      else
        {
        int localIdx = this->GlobalBlockIdx2LocalIdx[ globalIdx ];
        N = this->GetNumberOfElementsForBlock( localIdx );
        }
      }
    }
  else
    {
    N = this->GetNumberOfElements();
    }

  // STEP 1: Loop through all registered variables and read them in
  for(unsigned int varIdx=0; varIdx < this->Vars.size(); ++varIdx)
    {
    // Get the variable size
    size_t vsize = this->Vars[varIdx].Size;

    // pointer to data
    void *dataPtr = this->Vars[varIdx].Data;

    // Get the variable index, used to calculate the offset in the file
    int vidx = this->GetVariableIndex( this->Vars[ varIdx ].Name );
    if(vsize != this->VH[vidx].Size)
      {
      std::cerr << "Variable size mismatch for var: "
                << this->Vars[ varIdx ].Name << std::endl;
      throw std::runtime_error(
          "Variable size mismatch for " + this->Vars[ varIdx ].Name);
      }
    assert("pre: cannot find variable index!" &&
            (vidx >= 0) && (vidx < static_cast<int>(this->VH.size())) );

    // Loop through all blocks and read corresponding variable
    for(unsigned int block=0; block < this->AssignedBlocks.size(); ++block)
      {
      // Get the number of elements in the block
      int globalBlockIdx = this->AssignedBlocks[ block ];
      assert("pre: Reader does not have assigned global block index!" &&
              this->HasBlock(globalBlockIdx));
      int localBlockIdx = this->GlobalBlockIdx2LocalIdx[ globalBlockIdx ];
      int NBlockElements = this->GetNumberOfElementsForBlock(localBlockIdx);

      if(NBlockElements==0)
        {
        // skip empty blocks
        continue;
        }

      // sanity check!
      assert( "pre: Data for variable is NULL!" &&
              (this->Vars[varIdx].Data != NULL) );

      // Compute number of bytes to read (if the data is not compressed)
      size_t bytesize = NBlockElements*vsize;

      // Detect if this variable for the given block is compressed
      DataFilter* f = NULL;
      BlockInfo metadata;
      if( this->HasBlockMetaData )
        {
        metadata = this->GetBlockMetaData(localBlockIdx,vidx);
        // TODO: We assume a single filter was applied to the data(?)
        f = DataFilter::InstantiateByName(metadata.Filters[0]);
        } // END if the block metadata is available

      if( f != NULL )
        {

        // Allocate buffer to store the "filtered", e.g., compressed data
        void* filterBuffer = malloc(metadata.Size+CRCSize);

        // Read encoded data
        std::ostringstream oss;
        oss.str("");
        oss << "Reading encoded data for var=" << this->Vars[varIdx].Name;
        this->Read(filterBuffer,metadata.Size+CRCSize,metadata.Start,oss.str());

        // Verify the CRC of the filtered data
        if(!GenericIOUtilities::VerifyChecksum(
                            filterBuffer,metadata.Size+CRCSize))
          {
          throw std::runtime_error(
              "Variable checksum failed for filtered data variable "+
               this->Vars[varIdx].Name
               );
          } // END if CRC verify

        // Get the compress header, which is appended to the data and stores
        // the CRC of the original data
        CompressHeader CH;
        memcpy(&CH,filterBuffer,sizeof(CompressHeader));
        if(this->SwapEndian)
          {
          GenericIOUtilities::SwapCompressHeader(&CH);
          }

        // Get pointer to the data after the compress header
        void* bufferPtr = filterBuffer;
        bufferPtr       = static_cast<char*>(bufferPtr)+sizeof(CompressHeader);

        // Set the filtered data as an input to the filter
        f->SetInputData(bufferPtr,1,metadata.Size-sizeof(CompressHeader));

        // Decode the filtered data
        int rc = f->Decode();
        if( rc < 0)
          {
          throw std::runtime_error(
               "Failed decoding data for var="+this->Vars[varIdx].Name);
          }

        // Ensure the decoded output bytesize mataches the expected
        // bytesize of the data
        assert("Decoded output size mismatch" &&
                f->GetOutputByteSize()==bytesize);

        // Copy the de-coded data to the data buffer
        f->GetOutputData(dataPtr);

        // Verify Decoded data CRC
        if(!GenericIOUtilities::VerifyChecksum(
              dataPtr,bytesize,CH.OrigCRC,false))
          {
          throw std::runtime_error(
              "Decoded Data CRC failure for var="+this->Vars[varIdx].Name);
          } // END if CRC check

        // free the filtered data buffer
        free(filterBuffer);
        filterBuffer = bufferPtr = NULL;

        // Delete the filter
        delete f;
        f = NULL;

        } // END if a filter was applied to the data
      else
        {

        // Calculate the offset in the file for the given variable
        uint64_t offSet = (this->HasBlockMetaData)?
            metadata.Start : this->GetVariableOffSet(vidx,localBlockIdx);
        if( this->Vars[ varIdx ].HasExtraSpace )
          {
          this->Read(dataPtr,bytesize+CRCSize,offSet,this->Vars[varIdx].Name);
          if(!GenericIOUtilities::VerifyChecksum(dataPtr,bytesize+CRCSize))
            {
            throw std::runtime_error(
             "Variable checksum failed for variable " + this->Vars[varIdx].Name);
            } // END if checksum
          } // END if extra space
        else
          {
          // Read in the block data
          this->Read(dataPtr,bytesize,offSet,this->Vars[varIdx].Name);

          // Read checksum
          uint64_t checksum;
          offSet += bytesize;
          this->Read(&checksum,CRCSize,offSet,"variable checksum");

          // Verify Checksum
          if( !GenericIOUtilities::VerifyChecksum(dataPtr,bytesize,checksum) )
            {
            throw std::runtime_error(
             "Variable checksum failed for variable "+this->Vars[varIdx].Name);
            } // END if checksum
          } // END else no extra space

        } // END else no filter was applied to the data

      // Ensure that we de-allocate the filter instance (if used)
      assert(f==NULL);

      // Move pointer for next block
      dataPtr = static_cast<char*>(dataPtr)+bytesize;
      } // END for all assigned blocks

    // Swap endian if necessary
    if(this->SwapEndian)
      {
      std::vector<char> swapBuffer;
      swapBuffer.resize(vsize);
      void *ptr = this->Vars[varIdx].Data;
      for(int i=0; i < N; ++i, ptr=static_cast<char*>(ptr)+vsize)
        {
        GenericIOUtilities::SwapEndian(ptr,vsize,&swapBuffer[0]);
        }
      } // END if swap endian

    } // END for all variables
}

//------------------------------------------------------------------------------
int GenericIOReader::GetNumberOfElements(int globalBlockIdx)
{
  int NElements = 0;
  if( globalBlockIdx != -1 )
    {
    if( this->HasBlock(globalBlockIdx) )
      {
      int localIdx = this->GlobalBlockIdx2LocalIdx[ globalBlockIdx ];
      NElements = this->GetNumberOfElementsForBlock( localIdx );
      } // END if
    } // END if
  else
    {
    if( this->SplitMode )
      {
      for(int i=0; i < this->NumberOfFiles; ++i)
        {
        assert("pre: internal reader is NULL" &&
            (this->InternalReaders[i] != NULL));

        NElements += this->InternalReaders[ i ]->GetNumberOfElements();
        } // END for all files
      } // END if reading in split mode
      else
      {
      unsigned int blkIdx=0;
      for(; blkIdx < this->AssignedBlocks.size(); ++blkIdx)
        {
        NElements += this->GetNumberOfElementsForBlock(blkIdx);
        } // END for all blocks
      }
    } // END else

  return( NElements );
}

//------------------------------------------------------------------------------
int GenericIOReader::GetNumberOfElementsForBlock(const int blkidx)
{
  assert("pre: blkIdx is out-of-bounds!" &&
          (blkidx >= 0) && (blkidx < static_cast<int>(this->RH.size())) );
  return( this->RH[blkidx].NElems);
}

//------------------------------------------------------------------------------
std::string GenericIOReader::GetBlockHeadersInfo()
{
  std::ostringstream oss;
  oss.clear();

  oss << "NUMBLOCKS: " << this->RH.size() << std::endl;
  for(unsigned int i=0; i < this->RH.size(); ++i)
    {
    oss << "BLOCK[" << i << "]:\n";
    oss << "\t IJK: (" << this->RH[i].Coords[0] << ",";
    oss << this->RH[i].Coords[1] << ", ";
    oss << this->RH[i].Coords[2] << ")\n";
    oss << "\t GlobalRank: "  << this->RH[i].GlobalRank << std::endl;
    oss << "\t NumElements: " << this->RH[i].NElems << std::endl;
    oss << "\t OffSet: " << this->RH[i].Start << std::endl;
    if( this->BlockToFileMap.find(this->RH[i].GlobalRank) !=
        this->BlockToFileMap.end() )
      {
      oss << "\t File: " << this->BlockToFileMap[this->RH[i].GlobalRank];
      oss << std::endl;
      oss << "\t Index within File: ";
      oss << this->BlockToIdxWithinFile[this->RH[i].GlobalRank];
      oss << std::endl;
      }
    oss << "\t ===========\n";
    }
  return( oss.str() );
}

//------------------------------------------------------------------------------
std::string GenericIOReader::GetGlobalHeaderString(int indent)
{
  // STEP 0: setup the prefix to each line
  std::ostringstream prefix;
  prefix.clear();
  prefix.str("");
  for(int i=0; i < indent; ++i, prefix << "\t");

  std::ostringstream oss;
  oss.clear();

  oss << prefix.str() << "===============================" << std::endl;
  oss << prefix.str() << "GENERAL INFO\n";
  oss << prefix.str() << "FILENAME: " << this->FileName << std::endl;
  oss << prefix.str() << "NUMBER OF FILES: "
                      << this->NumberOfFiles << std::endl << std::endl;

  oss << prefix.str() << "GLOBAL HEADER:\n";
  oss << prefix.str() << "HeaderSize: " << this->GH.HeaderSize << std::endl;
  oss << prefix.str() << "NElems: " << this->GH.NElems << std::endl;
  oss << prefix.str() << "Dims: " << this->GH.Dims[0] << " "
                                  << this->GH.Dims[1] << " "
                                  << this->GH.Dims[2] << std::endl;
  oss << prefix.str() << "NVars: "      << this->GH.NVars      << std::endl;
  oss << prefix.str() << "VarsSize: "   << this->GH.VarsSize   << std::endl;
  oss << prefix.str() << "VarsStart: "  << this->GH.VarsStart  << std::endl;
  oss << prefix.str() << "NRanks: "     << this->GH.NRanks     << std::endl;
  oss << prefix.str() << "RanksSize: "  << this->GH.RanksSize  << std::endl;
  oss << prefix.str() << "RanksStart: " << this->GH.RanksStart << std::endl;
  oss << prefix.str() << "GlobalHeaderSize: "
                      << this->GH.GlobalHeaderSize << std::endl;
  oss << prefix.str() << "PhysOrigin: "
                      << this->GH.PhysOrigin[0] << " "
                      << this->GH.PhysOrigin[1] << " "
                      << this->GH.PhysOrigin[2] << std::endl;
  oss << prefix.str() << "PhysScale: "
                      << this->GH.PhysScale[0] << " "
                      << this->GH.PhysScale[1] << " "
                      << this->GH.PhysScale[2] << std::endl;

  for(int file=0; file < this->NumberOfFiles; ++file)
    {
    oss << prefix.str()
        << this->InternalReaders[ file ]->GetGlobalHeaderString(1);
    } // END for all files

  oss << "===============================" << std::endl;
  return( oss.str() );
}

//------------------------------------------------------------------------------
std::string GenericIOReader::GetBlockMapString()
{
  std::ostringstream oss;
  oss.clear();

  if( this->SplitMode )
    {
    std::map<int,int>::iterator iter = this->BlockToFileMap.begin();
    for(;iter!=this->BlockToFileMap.end();++iter)
      {
      oss << iter->first << "\t" << this->FileIDToSymbolicID[iter->second];
      oss << std::endl;
      } // END for
    }
  else
    {
    oss << "File is monolithic, no block mapping information available!";
    }

  return( oss.str() );
}

//------------------------------------------------------------------------------
uint64_t GenericIOReader::GetVariableOffSet(int vidx, int localBlkIdx)
{
  // Sanity checks!
  assert("pre: variable index is out-of-bounds!" &&
       (vidx >= 0) && (vidx < static_cast<int>(this->GH.NVars)) );
  assert("pre: local block index is out-of-bounds!" &&
       (localBlkIdx >= 0) && (localBlkIdx < static_cast<int>(this->RH.size())));


  int N = this->GetNumberOfElementsForBlock(localBlkIdx);
  uint64_t offSet = this->RH[ localBlkIdx ].Start;
  for(int i=0; i < vidx; ++i)
    {
    offSet += N*this->VH[i].Size + CRCSize;
    }
  return( offSet );
}

} /* namespace cosmotk */
