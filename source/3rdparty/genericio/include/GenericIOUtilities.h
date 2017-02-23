/**
 * @brief A singleton class that consists of generic functionality used
 * throughout the implementation of the GenericIO readers/writers.
 */
#ifndef GENERICIOUTILITIES_H_
#define GENERICIOUTILITIES_H_

#include "GenericIODefinitions.hpp"

#include <map>
#include <sstream>
#include <string>
#include <vector>

#include <mpi.h>

namespace gio
{

class GenericIOUtilities
{
public:

  /**
   * @brief Checks whether this machine is big endian or not
   * @return status true iff big endian, else false
   */
  static bool IsBigEndian()
    {
    const uint32_t one = 1;
    return !(*((char *)(&one)));
    }

  /**
   * @brief Given structured coordinates, (i,j,k), this method generates
   * a hash-code key that is represented by the string "i.j.k".
   * @note The key is intended to be used with a std::map
   * @param i the ith coordinate (in)
   * @param j the jth coordinate (in)
   * @param k the kth coordinate (in)
   * @return hash the string hashcode corresponding to the given coordinates.
   * @post hash.length() > 0
   * @see GenericIOUtilities::GetCoordsFromHashCode
   */
  static std::string GetHashCodeFromCoords(
      const uint64_t i,
      const uint64_t j,
      const uint64_t k)
    {
    std::ostringstream oss;
    oss << i << '.' << j << '.' << k;
    return( oss.str() );
    }

  /**
   * @brief Extract the structured coordinates from the given hashcode.
   * @param hashCode the string hashcode to be processed.
   * @param i the ith coordinate (out)
   * @param j the jth coordinate (out)
   * @param k the kth coordinate (out)
   * @pre hashCode.length() > 0
   * @see GenericIOUtilities::GetHashCodeFromCoords
   */
  static void GetCoordsFromHashCode(
      std::string &hashCode,
      uint64_t &i,
      uint64_t &j,
      uint64_t &k)
  {
  char dot;
  std::istringstream iss(hashCode);
  iss >> i >> dot
      >> j >> dot
      >> k;
  }

  /**
   * @brief Distributes blocks from the file to processes according to the
   * user-supplied assignment strategy.
   * @param processID the ID of the this process instance (in).
   * @param numProcesses the total number of processes (in).
   * @param totalNumberOfBlocks the total number of blocks in the file (in).
   * @param blockDecomposition the block decomposition dimensions (in).
   * @param assignmentStrategy the strategy to use for distribution (in).
   * @param ijkMap mapping of ijk rank coordinates to the linear index (in).
   * @param assigned list of blocks assinged to this process (out).
   * @param neighbors list of rank neighbors (out).
   * @note the neighbors list is only populated if the data blocks are
   * spatially arranged and RCB is used as the assignment strategy.
   * @see GenericIOBlockAssignement for list of assignment strategies.
   */
  static void BlockAssignment(
      const int processID,
      const int numProcesses,
      const int totalNumberOfBlocks,
      const uint64_t blockDecomposition[3],
      const int assignmentStrategy,
      ProcessCoordsToRank* ijkMap,
      std::vector<int>& assigned,
      std::vector<RankNeighbor>& neighbors
      );

  /**
   * @brief Verifies the checksum supplied with and corresponding to the data.
   * @param data the data buffer which includes the checksum padded at the end.
   * @param nbytes the number of bytes including the checksum
   * @return status true iff the checksum verification passes, else, false.
   * @see GenericIOUtilities::VerifyChecksum
   * @note This version is used primarily when the data+crc are read with a
   * single read, i.e., when the user-supplied buffer has allocated extra
   * space for the crc.
   */
  static bool VerifyChecksum(const void* data, size_t nbytes);

  /**
   * @brief Computes a 64-bit CRC checksum for the supplied data and verifies
   * it against the expected checksum, cs.
   * @param data the data to check.
   * @param nbytes the number of bytes of the data.
   * @param cs the checksum to compare against.
   * @param checkInverse indicates if comparing against the inverse checksum.
   * @return status true iff the checksum verification passes, else, false.
   */
  static bool VerifyChecksum(
      const void* data,size_t nbytes,uint64_t cs,
      bool checkInverse=true);

  /**
   * @brief Checks if the endian of this machine matches the endian of the
   * file encoded in the global header.
   * @param GH pointer to the global header
   * @return status true if the endian matches, else, false.
   */
  static bool DoesFileEndianMatch(GlobalHeader *GH);

  /**
   * @brief Swaps the endian of the data pointed to by Addr.
   * @param Addr user-supplied address of the data to swap (in/out)
   * @param N the number of bytes of the data to swap (in)
   * @param buf buffer space used to swap data (optional)
   * @note When this method is called within a tight loop it is advised that
   * a buffer is provided and allocated externally. Otherwise, the method will
   * do the necessary buffering internally but will have a cost of allocating
   * and de-allocating small memory for potentially a large number of times.
   */
  static void SwapEndian(void *Addr, const int N, char *buf=NULL);

  /**
   * @brief Swaps the endian of the given global header.
   * @param GH pointer to the global header to swap.
   */
  static void SwapGlobalHeader(GlobalHeader *GH);

  /**
   * @brief Swaps the endian of the Variable header.
   * @param VH pointer to the variable header to swap.
   */
  static void SwapVariableHeader(VariableHeader *VH);

  /**
   * @brief Swaps the endian of the rank header.
   * @param RH pointer to the rank header to swap.
   */
  static void SwapRankHeader(RankHeader *RH);

  /**
   * @brief Swaps the endian of the compressed block information.
   * @param blockInfo pointer to the compressed block info object.
   */
  static void SwapBlockInfo(BlockInfo* blockInfo);

  /**
   * @brief Swaps the endian of the compress header.
   * @param CH pointer to the compress header object.
   */
  static void SwapCompressHeader(CompressHeader* CH);

  /**
   * @brief Given the variable information, it detects the corresponding
   * underlying primitive type to represent the data associated with the
   * variable.
   * @param vinfo the variable information of the variable in query.
   * @return t the primitive type.
   * @see GenericIOPrimitiveTypes.
   */
  static int DetectVariablePrimitiveType(const VariableInfo &vinfo);

  /**
   * @brief Allocates and returns a pointer (void*) to a buffer wherein the
   * data of the corresponding variable can be stored.
   * @param vinfo the variable information of the variable in query.
   * @param numElements the number of elements to allocate.
   * @param pad indicated whether the variable is padded with extra space.
   * @return bufPtr pointer to allocated buffer.
   * @post bufPtr != NULL.
   */
  static void* AllocateVariableArray(
            const VariableInfo &vinfo,const int numElements,bool pad=true);

  /**
   * @brief De-Allocates a variable array.
   * @param vinfo the variable information of the variable.
   * @param dataArray pointer to the data array being allocated.
   * @post dataArray == NULL
   */
  static void DeallocateVariableArray(
          const VariableInfo &vinfo,void* dataArray);

  /**
   * @brief Rank 0 prints the formatted message to stdout
   * @param comm the communicator used by the parallel application
   * @param fmt the formatted message to print
   * @note This is a collective call, all processes must call this method.
   * @pre comm != MPI_COMM_NULL
   * @pre fmt != NULL
   */
  static void Printf(MPI_Comm comm, const char *fmt,...);

  /**
   * @brief Each process prints the formatted message in stdout in rank order.
   * That is, process "i" prints its message right after process "i-1".
   * @param comm the communicator used by the parallel application
   * @param fmt the formatted message to print
   * @note This is a collective call, all processes must call this method.
   * @pre comm != MPI_COMM_NULL
   * @pre fmt != NULL
   */
  static void SynchronizedPrintf(MPI_Comm comm, const char *fmt,...);

protected:
  GenericIOUtilities();
  virtual ~GenericIOUtilities();

   /**
    * @brief Round-Robins a user-supplied number of blocks to a bucket with the
    * given index, idx. The resulting blocks that are assigned to bucket in
    * query are inserted in the given STL vector.
    * @param bucketIdx the bucket index (in).
    * @param NumBuckets the total number of buckets (in).
    * @param NumBlocks the total number of blocks (in).
    * @param assigned a list of assigned blocks (in/out).
    * @note This method is used within the BlockAssignment (public) method.
    */
   static void RoundRobin(
       const int bucketIdx, const int NumBuckets,
       const int NumBlocks,
       std::vector<int> &assigned);

   /**
    * @brief Partitions the wholeExtent into the number of processes and
    * returns the subExtent assigned to this process.
    * @param processID the ID of this process instance (in).
    * @param numProcesses the total number of processes (in).
    * @param wholeExtent the wholeExtent to be partitioned (in).
    * @param subExtent the subExtent assigned to this process (out).
    * @note This method is used within the RCBAssignment method.
    */
   static void RCBPartition(
       const int processID,
       const int numProcesses,
       const uint64_t wholeExtent[6],
       uint64_t subExtent[6],
       std::vector<RankNeighbor> &neighbors);

   /**
    * @brief Assigns blocks to processes using Recursive Coordinate Bisection.
    * @param processID the ID of this process instance (in).
    * @param blkdims the dimensions of the blocks in the file (in).
    * @param ijkMap mapping of ijk rank coordinates to the linear index (in).
    * @param assigned list of blocks assigned to this process (out).
    * @note In contrast to round-robin assignment, this will distribute blocks
    * in to spatially contiguous regions.
    */
   static void RCBAssignment(
       const int processID,
       const int numProcesses,
       const uint64_t blkdims[3],
       ProcessCoordsToRank* ijkMap,
       std::vector<int> &assigned,
       std::vector<RankNeighbor> &neighbors);

private:
  GIO_DISABLE_COPY_AND_ASSIGNMENT(GenericIOUtilities);
};


} /* namespace cosmotk */
#endif /* GENERICIOUTILITIES_H_ */
