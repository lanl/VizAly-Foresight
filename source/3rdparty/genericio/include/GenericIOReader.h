/**
 * @brief GenericIOReader defines the interface of the GenericIO readers.
 * This is a pure abstract class -- all methods are defined by concrete
 * classes.
 */
#ifndef GENERICIOREADER_H_
#define GENERICIOREADER_H_

#include "GenericIOBase.h" // Base class
#include "GenericIODefinitions.hpp"

#include <cassert>     // for assert()
#include <map>         // for STL map
#include <string>      // for STL string
#include <sys/types.h> // for off_t
#include <vector>      // for STL vector

namespace gio
{

class GenericIOReader : public GenericIOBase
{
public:
  GenericIOReader();
  virtual ~GenericIOReader();

  // In-line Get/Set macros
  GIOGetNSetMacro(BlockAssignmentStrategy,int);
  GIOGetNSetMacro(Communicator,MPI_Comm);
  GIOGetNSetMacro(DoAcquireGlobalDimensions,bool);
  GIOGetNSetMacro(DoAcquireIJKMapping,bool);
  GIOGetNSetMacro(DoBlockAssignment,bool);
  GIOGetNSetMacro(DoDetermineFileType,bool);
  GIOGetNSetMacro(DoReadBlockHeaders,bool);
  GIOGetNSetMacro(DoSetupInternalReaders,bool);
  GIOGetNSetMacro(FunctionDescription,std::string)
  GIOGetNSetMacro(IJKMapToGlobalBlockIdx,ProcessCoordsToRank*);
  GIOGetNSetMacro(NumberOfFiles,int);
  GIOGetNSetMacro(SingleBlockMode,bool);

  /**
   * @brief Barrier synchronization across all MPI ranks.
   */
  void Barrier()
    {MPI_Barrier(this->Communicator);};

  /**
   * @brief Returns the total number of rank neighbors for this process.
   * @return N the number of neighbors.
   * @note This method is applicable only when the data is spatially
   * decomposed.
   */
  int GetNumberOfRankNeighbors()
  {
  return this->RankNeighbors.size();
  }

  /**
   * @brief Returns the rank neighbors of this process
   * @param neighbors user-supplied buffer to fill in the neighbors.
   * @note the user-supplied buffer must be pre-allocated according to
   * the value returned by GetNumberOfRankNeighbors().
   * @see RankNeighbor
   */
  void GetRankNeighbors(RankNeighbor* neighbors);

  /**
   * @brief Returns the number of blocks headers that have been read.
   * @return N the number of block headers.
   */
  int GetNumberOfBlockHeaders()
      { return this->RH.size(); }

  /**
   * @brief Returns the requested block header.
   * @param i the index of the requested block header.
   * @return blockHeader the header of the corresponding block.
   */
  RankHeader GetBlockHeader(const int i)
    {
      assert("pre: block header index is out of bounds!" &&
              (i >= 0) && (i < static_cast<int>(this->RH.size())));
      return(this->RH[i]);
    }

  /**
   * @brief Returns the requested block metadata.
   * @param blockIdx the index of the requested block.
   * @param varIdx the index of the requested variable.
   * @return b the metadata of the corresponding variable in the given block.
   * @pre blockIdx >= 0 && blockIdx < this->GH.NRanks
   * @pre varIdx >= 0 && varIdx < this->GH.NVars
   */
  BlockInfo GetBlockMetaData(const int blockIdx,const int varIdx)
    {
      assert("pre: no block metadata!" && this->HasBlockMetaData );
      assert("pre: blockIdx is out-of-bounds!" &&
              (blockIdx >= 0) && (blockIdx < this->GH.NRanks) );
      assert("pre: varIdx is out-of-bounds!" &&
              (varIdx >= 0) && (varIdx < this->GH.NVars) );

      int idx = blockIdx*this->GH.NVars+varIdx;
      return(this->BlockMetadata[idx]);
    }

  /**
   * @brief Returns the GlobalID of the requested block
   * @param i the index of the requested block.
   * @return idx the global block index.
   */
  uint64_t GetBlockGlobalID(const int i)
    {
    assert("pre: block header index is out of bounds!" &&
             (i >= 0) && (i < static_cast<int>(this->RH.size())) );
    return(this->RH[i].GlobalRank);
    }

  /**
   * @brief Checks if the block with the given global block index is assigned
   * to this reader instance, running at the respective rank.
   * @param globalIdx the global block index.
   * @return status true if the block is assigned to this reader, else false.
   */
  bool HasBlock(const int globalIdx)
    {
      if(this->GlobalBlockIdx2LocalIdx.find(globalIdx)==
          this->GlobalBlockIdx2LocalIdx.end())
        {
        return false;
        }
      return true;
    }

  /**
   * @brief Checks if the underlying data is spatially decomposed.
   * @return status true if the data is spatially decomposed, else, false.
   * @note Data writen in GenericIO are spatially decomposed iff the application
   * that writes the data uses an MPI cartesian communicator.
   */
  bool IsSpatiallyDecomposed();

  /**
   * @brief Returns the cartesian bounds of the ith block.
   * @param i the index of the requested block (in).
   * @param min[3] lower corner cartesian coordinates of the bounding box (out).
   * @param max[3] top corner cartesian coordinates of the bounding box (out).
   * @throw std::runtime_error if this method is called on a non-spatially
   * decomposed dataset.
   * @see GenericIOReader::IsSpatiallyDecomposed()
   */
  void GetBlockBounds(const int i, double min[3], double max[3]);

  /**
   * @brief Returns the variable information object for the ith variable.
   * @param i the index of the variable in query.
   * @return varInfo the variable information object
   * @note In contrast to GetVariableInfo() this method returns the information
   * of the ith variable in the file.
   * @pre (i >= 0) && (i < this->GetNumberOfVariablesInFile())
   * @see VariableInfo in GenericIODefinitions.hpp
   * @see GenericIOBase::GetVariableInfo()
   */
  VariableInfo GetFileVariableInfo(const int i);

  /**
   * @brief Returns the variable header of the ith variable.
   * @param i the index of the variable in query
   * @return varHeader the variable header
   * @see VariableHeader in GenericIODefinitions.hpp
   */
  VariableHeader GetVariableHeader(const int i)
    {
    assert("pre: variable index is out-of-bounds!" &&
            (i >= 0) && (i < static_cast<int>(this->VH.size())) );
    return( this->VH[i] );
    }

  /**
   * @brief Return the variable size of the ith variable
   * @param i the index of the variable in query
   * @return varsize the size of the variable, e.g., 8 if it's a double, etc.
   */
  uint64_t GetVariableSize(const int i)
    {
    assert("pre: variable index is out-of-bounds!" &&
           (i >= 0) && (i < static_cast<int>(this->VH.size())) );
    return( this->VH[ i ].Size );
    };

  /**
   * @brief Return the name of the ith variable
   * @param i the index of the variable in query
   * @return varname the variable name
   */
  std::string GetVariableName(const int i)
    {
    assert("pre: variable index is out-of-bounds!" &&
            (i >= 0) && (i < static_cast<int>(this->VH.size())) );
    return( std::string( this->VH[ i ].Name ) );
    };

  /**
   * @brief Return the number of variables in the file
   * @return nvar the number of variables in the file
   */
  int GetNumberOfVariablesInFile()
    {return this->VH.size();};

  /**
   * @brief Return the number of assigned blocks
   * @return nassigned the number of assigned blocks to this process.
   */
  int GetNumberOfAssignedBlocks()
    {return this->AssignedBlocks.size();};

  /**
   * @brief Return the total number of blocks in the file.
   * @return nblocks the total number of blocks in the file.
   * @note The number of blocks in the file is equivalent to
   * the number of processes that wrote the file.
   */
  int GetTotalNumberOfBlocks()
    {return this->GH.NRanks;};

  /**
   * @brief Return the total number of elements.
   * @return N the total number of elements.
   */
  int GetTotalNumberOfElements()
    {return this->GH.NElems;};

  /**
   * @brief Returns the physical scale read from the Global header.
   * @param physScale user-supplied buffer to store the scale (out).
   */
  void GetPhysScale(double physScale[3])
    {
    physScale[0] = this->GH.PhysScale[0];
    physScale[1] = this->GH.PhysScale[1];
    physScale[2] = this->GH.PhysScale[2];
    }

  /**
   * @brief Returns the physical origin read from the Global header.
   * @param origin user-supplied buffer to store the origin (out).
   */
  void GetPhysOrigin(double origin[3])
    {
    origin[0] = this->GH.PhysOrigin[0];
    origin[1] = this->GH.PhysOrigin[1];
    origin[2] = this->GH.PhysOrigin[2];
    }

  /**
   * @brief Get the global dimensions from the Global header.
   * @param dims user-supplied buffer to store the global dimensions (out).
   */
  void GetGlobalDimensions(uint64_t dims[3])
    {
    dims[0] = this->GH.Dims[0];
    dims[1] = this->GH.Dims[1];
    dims[2] = this->GH.Dims[2];
    };

  /**
   * @brief Checks if the given variable exists.
   * @param varName the variable name.
   * @return status true if the variable exists, else false.
   */
  bool HasVariable(std::string varName)
    { return( (this->GetVariableIndex(varName)!=-1) ); };

  /**
   * @brief Checks if the underlying files has been written in split mode
   * @return status true, if the mode is split, else false.
   */
  bool IsSplitMode() { return( this->SplitMode ); };

  /**
   * @brief Indicates that this reader is an internal reader.
   * @note Internal readers do not compute any block assignment,
   * instead they are assigned blocks by the proxy reader that is
   * used from the calling application.
   */
  void SetReaderInternal() { this->InternalReader=true; };

  /**
   * @brief Assigns the given block to the reader instance running on this rank.
   * @param blkIdx the global index of the block to reader
   */
  void AssignBlock(const int blkIdx);

  /**
   * @brief Assigns a single block to the reader corresponding to the given
   * block header and clears out all other previously read headers.
   * @param blockHeader the block header to read
   * @note Used to read individual blocks in split-mode
   */
  void AssignSingleBlock(const RankHeader blockHeader);

  /**
   * @brief Clears all assigned blocks to this reader.
   * @post this->GetNumberOfAssignedBlocks()==0
   */
  void ClearBlockAssignment()
    {this->AssignedBlocks.clear();};

  /**
   * @brief Returns the number of elements that will be read. If called with
   * an argument for the globalBlockIdx, this method will return the number
   * of elements in the requested block if that block is assigned to this
   * process, otherwise, a zero will be returned.
   * @param globalBlockIdx optional argument,
   * @return N the number of elements to read or the number of elements within
   * a particular block if called with an argument.
   */
  int GetNumberOfElements( int globalBlockIdx=-1 );

  /**
   * @brief Opens and reads the header of the file and initializes internal
   * data-structures. Optionally, the caller can specify to skip reading the
   * block headers of a file. By default, only the blocks that the reader
   * is assigned to will be read.
   * @pre !This->FileName.empty()
   */
  void OpenAndReadHeader();

  /**
   * @brief Read the headers of each assigned block.
   * @see ReadBlockHeader
   */
  void ReadBlockHeaders();

  /**
   * @brief Reads the block metadata if available.
   * @note Block metadata is available, e.g., when the data is compressed
   * and consists of additional information for each variable in each block
   * about what filter (e.g., compression algorithm) was applied etc.
   */
  void ReadBlockMetaData();

  /**
   * @brief Reads the data in the user-supplied registered buffers for the
   * given block.
   * @param globalBlockIdx the global index of the block to read.
   * @return status true if the block was read, else false.
   * @note status==false indicates that this->HasBlock(globalBlockIdx)==false.
   * @note Prior to calling this method the calling code should have:
   * <ol>
   *  <li> Called GetNumberOfElements(globalBlockIdx) </li>
   *  <li> Allocated buffers accordingly </li>
   *  <li> Called AddVariable </li>
   * <ol>
   */
  bool ReadBlock(const int globalBlockIdx);

  /**
   * @brief Reads the data in to the user-supplied registered arrays.
   * @note The user should have registered the arrays to read via calls to
   * the AddVariable method.
   * @see GenericIOBase::AddVariable.
   */
  void ReadData();

  /**
   * @brief Closes the file
   */
  virtual void Close() = 0;

  /**
   * @brief Creates a new instance of GenericIOReader.
   * @return reader pointer to a GenericIOReader instance.
   * @note This method is implemented by concrete classes.
   */
  virtual GenericIOReader* GetNewInstance() = 0;

  /**
   * @brief Returns the block header information as a string.
   * @return str the block header information.
   * @note This method is primarily used for debugging.
   */
  std::string GetBlockHeadersInfo();

  /**
   * @brief Returns the global header associated with this file as a string.
   * @return str the global header information packed in a string.
   * @note This method is used primarily for debugging.
   */
  std::string GetGlobalHeaderString(int indent=0);

  /**
   * @brief Returns the block-to-file mapping information as a string.
   * @return str the block-to-file mapping.
   * @note This method is used primarily for debugging.
   */
  std::string GetBlockMapString();

protected:
  MPI_Comm Communicator;
  int Rank;
  int NumRanks;

  // Property that indicates whether the data must be byte-swapped
  bool SwapEndian;

  // Property that indicates whether the data is written in split-mode
  bool SplitMode;

  // Property that indicates that we are reading a single block
  bool SingleBlockMode;

  // Flag that indicates that in addition to the RankHeader information, each
  // block has some additional metadata.
  bool HasBlockMetaData;

  // Properties for configuring the actions when OpenAndReadHeader is called.
  // By default, all properties are enabled.
  bool DoAcquireGlobalDimensions;
  bool DoAcquireIJKMapping;
  bool DoBlockAssignment;
  bool DoDetermineFileType;
  bool DoReadBlockHeaders;
  bool DoSetupInternalReaders;

  // Adding a function description property to the reader, primarily for
  // debugging purposes.
  std::string FunctionDescription;

  // Flag that indicates that the reader is internal.
  bool InternalReader;

  // List of internal readers used when reading in SplitMode.
  GenericIOReader** InternalReaders;

  // Stores the entire raw bytes of the GenericIO header which consists of the
  // GlobalHeader, the variable headers and the block (rank) headers, including
  // the checksum of the header.
  std::vector< char > EntireHeader;

  // Extracted global header, each process will extract the GlobaHeader and
  // cache it in this ivar.
  GlobalHeader GH;

  // Extracted variable header, each process extracts all of the variable
  // headers in the file and cache them in this vector.
  std::vector< VariableHeader > VH;

  // Extracted rank headers, i.e., blocks headers in the file. Each process
  // will *only* extract and cache the blocks that it is assigned.
  std::vector< RankHeader > RH;

  // Stores metadata for each block if available. BlockMetadata is available
  // if the data is "filtered" e.g., compressed
  std::vector< BlockInfo > BlockMetadata;

  std::map<std::string, int> VariableName2IndexMap;

  // Mapping of a global block index to a local block index.
  std::map<int,int> GlobalBlockIdx2LocalIdx;

  // Stores the mapping of data-blocks to a file ID. Note, this is only
  // applicable iff this->SplitMode==true.
  std::map<int,int> BlockToFileMap;

  // Mapping of the structured coordinates of each block to the global linear
  // block index.
  ProcessCoordsToRank* IJKMapToGlobalBlockIdx;

  // Stores a mapping of the file ID to the actual symbolic ID used to compose
  // the filename. On certain systems, when writing in SplitMode to a number
  // files N, the filenames may not be sequential, e.g., <filename>#0,
  // <filename>#1,...,<filename>#N. The suffix "#x" could be arbitrary and
  // that number is stored in the map file that is read by rank 0. Hence,
  // this mapping stores the logical file ID to the symbolic suffix ID, in
  // order to construct the correct filename for each of the files.
  std::map<int,int> FileIDToSymbolicID;
  std::map<int,int> SymbolicIDToFileID;

  // Maps a global block ID, to the global idx within a given file. Note,
  // this is only applicable iff this->SplitMode==true.
  std::map<int,int> BlockToIdxWithinFile;

  // The number of files the data has been split to. Note, this is only
  // applicable iff this->SplitMode==true.
  int NumberOfFiles;

  // Indicates which strategy is going to be used for distributing blocks
  // to processes. See GenericIOBlockAssignment in GenericIODefinitions for
  // a list of possibilities.
  int BlockAssignmentStrategy;

  // List of blocks that will be read by this process. Recall, the number of
  // blocks in the file may not always match the number of ranks that this
  // reader instance is running. Hence, each rank may be assigned more than
  // one block in the file. The vector of AssignedBlocks holds the list of
  // blocks for this process.
  std::vector< int > AssignedBlocks;

  // List of neighbors for this rank. This list is populated only in the
  // case when the domain is RCB decomposed.
  std::vector< RankNeighbor > RankNeighbors;

  /**
   * @brief Opens the FileHandle from where the data is going to be read.
   * @note This method is implemented by concrete implementations.
   */
  virtual void Open()=0;

  /**
   * @brief Reads data into the user-supplied buffer from the MPI file handle
   * @param buf the buffer where the data will be read into
   * @param count the number of bytes to read
   * @param offset the offset from which to read
   * @param name name of the data being read, primarily, used for
   * debugging and error reporting.
   * @note This method is implemented by concrete implementations.
   */
  virtual void Read(void *buf, size_t count, off_t offset,
                  const std::string &name) = 0;

  /**
   * @brief Allocates internal readers.
   * @param N number of readers to allocate.
   * @note This method is implemented by concrete implementations.
   */
  virtual void AllocateInternalReaders(const int N)=0;

  /**
   * @brief This method constructs a global-to-local block mapping.
   */
  void ConstructGlobalToLocalBlockMapping();


  /**
   * @brief Sets up the internal readers used for SplitMode I/O.
   * @note Applicable iff SplitMode.
   */
  void SetupInternalReaders();

  /**
   * @brief Clears the internal readers.
   * @note Applicable iff SplitMode.
   */
  void ClearInternalReaders();

  /**
   * @brief Based on the global & variable headers, this method determines if
   * the file has been written in split mode or not.
   * @pre The variables must have been indexed prior to calling this method,
   * i.e., this method should be called after IndexVariables() has been called.
   * @note This is a collective operation, all ranks must call this method.
   */
  void DetermineFileType();

  /**
   * @brief Reads the global and variable header of the file and broadcasts
   * them to all ranks.
   * @note This is a collective operation, all ranks must call this method.
   */
  void ReadHeader();

  /**
   * @brief Reads in the BlockToFile mapping.
   * @pre this->SplitMode==true.
   * @note This is a collective operation, all ranks must call this method.
   */
  void ReadBlockToFileMap();

  /**
   * @brief Reads the header of the ith variable
   * @param vh the variable header where the dara will be read in.
   * @pre vh != NULL
   * @see ReadVariableHeaders
   */
  void ReadVariableHeader( const int i, VariableHeader& vh );

  /**
   * @brief Reads the block header corresponding to the given block index.
   * @param blkIdx the index of the block to reader
   * @param blockHeader data-structure where to read in the block header
   * @see ReadBlockHeaders
   */
  void ReadBlockHeader(const int blkIdx, RankHeader& blockHeader);

  /**
   * @brief Reads block metadata corresponding to the given block/variable ID.
   * @param blkIdx the index of the block to read (in).
   * @param varIdx the index of the variable to read within the given block (in).
   * @param blockMetadata data-structure where to read in the block metadata.
   * @note This method is called only if block metadata is available, e.g., if
   * the data is compressed.
   * @see ReadBlockHeaders
   */
  void ReadBlockVariableMetaData(
      const int blkIdx, const int varIdx, BlockInfo& blockMetadata);

  /**
   * @brief Reads in the variable headers for the number of variables.
   * @pre This method assumes that the global header has been read in.
   * @see ReadVariableHeader
   */
  void ReadVariableHeaders();

  /**
   * @brief Builds an index based on variable name.
   */
  void IndexVariables();

  /**
   * @brief Read the data in split mode.
   */
  void ReadSplitModeData();

  /**
   * @brief Reads data from a single file.
   */
  void ReadSingleFileData();

  /**
   * @brief Returns the number of elements at the given block.
   * @param localBlkIdx the local block index of the block in query.
   * @return N the number of elements stored in the requested block.
   * @note Assumes that the rank header information has been read.
   */
  int GetNumberOfElementsForBlock(const int localBlkIdx);

  /**
   * @brief Return the index of the variable with the given name.
   * @param name the name of the variable in query
   * @return idx the corresponding index of the variable
   * @post idx >= 0, idx == -1 iff a variable is not found.
   */
  int GetVariableIndex(const std::string name);

  /**
   * @brief Computes the variable offset within the given global block idx.
   * @param vidx the variable index.
   * @param localBlkIdx the local block index.
   * @return offSet the variable offset.
   * @note This method uses the information in the rank header
   */
  uint64_t GetVariableOffSet(int vidx, int localBlkIdx);

  /**
   * @brief This method is used to get the global dimensions and store them
   * in to the global header when reading in split mode.
   * @note When reading in split-mode, the map file sets the dimensions to
   * [1 1 1] since only rank 0 writes the map file. The global dimensions
   * are in the global header of each of the individual files.
   */
  void AcquireGlobalDimensions();

  /**
   * @brief This method acquires the IJK mapping of each block to its
   * corresponding linear index.
   */
  void AcquireIJKMapping();

  /**
   * @brief Checks if the blocks have been written in IJK-sorted order based
   * on the global dimensions from the header.
   * @return status true if the blocks are IJK-sorted, otherwise, false.
   */
  bool IsIJKSorted();

private:
  GIO_DISABLE_COPY_AND_ASSIGNMENT(GenericIOReader);
};

} /* namespace cosmotk */
#endif /* GENERICIOREADER_H_ */
