#include "GenericIOUtilities.h"

#include "CRC64.h"

#include <cassert>
#include <cmath>
#include <cstdarg>
#include <cstdio>
#include <cstring>
#include <fstream>
#include <iostream>
#include <queue>
#include <sstream>
#include <stdexcept>
#include <string>


// Some usefull extent macros
#define IMIN(ext) ext[0]
#define IMAX(ext) ext[1]
#define JMIN(ext) ext[2]
#define JMAX(ext) ext[3]
#define KMIN(ext) ext[4]
#define KMAX(ext) ext[5]

// Some useful IJK macros
#define I(ijk) ijk[0]
#define J(ijk) ijk[1]
#define K(ijk) ijk[2]

//------------------------------------------------------------------------------
//  INTERNAL DATA STRUCTURES DEFINITIONS
//------------------------------------------------------------------------------

/**
 * @brief Elegant hack to do random-access on an STL priority_queue.
 * @note e.g., std::vector<T> &list = Container( q ), where q is priority queue.
 * @see http://stackoverflow.com/questions/4484767/how-to-iterate-over-a-priority-queue
 */
template <class T, class S, class C>
    S& Container(std::priority_queue<T, S, C>& q) {
        struct HackedQueue : private std::priority_queue<T, S, C> {
            static S& Container(std::priority_queue<T, S, C>& q) {
                return q.*&HackedQueue::c;
            }
        };
    return HackedQueue::Container(q);
}

//------------------------------------------------------------------------------

/**
 * @brief A simple extent object struct.
 * @note Used in Recursive Coordinate Bisection.
 */
struct extent_t
{
  uint64_t extent[6];

  bool operator < (const extent_t &t) const
  {
   return( this->size() < t.size() );
  }

  bool operator == (const extent_t &t) const
  {
    for(int i=0; i < 6; ++i)
      {
      if( this->extent[i] != t.extent[i])
        {
        return false;
        }
      }
    return true;
  }

  void GetLength(uint64_t L[3]) const
  {
    L[0] = IMAX(this->extent)-IMIN(this->extent)+1;
    L[1] = JMAX(this->extent)-JMIN(this->extent)+1;
    L[2] = KMAX(this->extent)-KMIN(this->extent)+1;
  }

  int getLargestDimension() const
  {
    uint64_t ilength = IMAX(this->extent)-IMIN(this->extent)+1;
    uint64_t jlength = JMAX(this->extent)-JMIN(this->extent)+1;
    uint64_t klength = KMAX(this->extent)-KMIN(this->extent)+1;

    int longest_dim = -1;
    if( (ilength >= jlength) && (ilength >= klength))
      {
      longest_dim = 0;
      }
    else if( (jlength >= ilength) && (jlength >= klength))
      {
      longest_dim = 1;
      }
    else if( (klength >= ilength) && (klength >= jlength))
      {
      longest_dim = 2;
      }
    return( longest_dim );
  }

  void split(extent_t &s1, extent_t &s2) const
  {
    // STEP 0: temporary variables
    uint64_t numNodes = 0;
    uint64_t midIdx   = -1;
    uint64_t minIdx   = -1;
    uint64_t maxIdx   = -1;

    // STEP 1: Initialize sub-extents to this extent instance.
    for(int i=0; i < 6; ++i)
      {
      s1.extent[ i ] = s2.extent[ i ] = this->extent[ i ];
      }

    // STEP 3: Compute split dimension, always split along the largest dimension
    int dim = this->getLargestDimension();

    // STEP 4: Set minIdx according to dim, w.r.t, the extent array
    switch( dim )
      {
      case 0: // i.e., split along the i-dimension
        minIdx = 0;
        break;
      case 1: // i.e., split along the j-dimension
        minIdx = 2;
        break;
      case 2: // i.e., split along the k-dimension
        minIdx = 4;
        break;
      default:
        assert("pre: undefined dimension! code should not reach here!\n");
        throw std::runtime_error("Code should not reach here!\n");
      } // END switch

    // STEP 5: Since we store the extent [imin,imax,jmin,jmax,kmin,kmax]
    // maxIdx is always minIdx+1, w.r.t. the extent array
    maxIdx = minIdx+1;

    // STEP 6: Split and update sub-extents accordingly
    numNodes = (this->extent[maxIdx]-this->extent[minIdx]+1);
    midIdx   = this->extent[minIdx]+
        static_cast<uint64_t>(std::floor(0.5*numNodes));


//    gio::GenericIOUtilities::SynchronizedPrintf(MPI_COMM_WORLD,
//        "\nmin=%d\nmid=%d\nmax=%d\n",
//        this->extent[minIdx],
//        midIdx,
//        this->extent[maxIdx]
//        );

    if( midIdx == this->extent[maxIdx] )
      {
      midIdx = this->extent[minIdx];
      }

    assert("pre: midIdx out-of-bounds!" &&
            (midIdx >= this->extent[minIdx]) &&
            (midIdx <= this->extent[maxIdx]) );
    s1.extent[maxIdx] = midIdx;
    s2.extent[minIdx] = midIdx+1;
  }

  size_t size() const
  {
    // TODO: deal with different dimensions/orientations
    size_t mySize = 0;
    int dims[3];
    dims[0] = IMAX(this->extent)-IMIN(this->extent)+1;
    dims[1] = JMAX(this->extent)-JMIN(this->extent)+1;
    dims[2] = KMAX(this->extent)-KMIN(this->extent)+1;
    mySize = dims[0]*dims[1]*dims[2];
    return(mySize);
  }

  std::string str() const
  {
    std::ostringstream oss;
    oss << "[ ";
    oss << IMIN(this->extent) << " ";
    oss << IMAX(this->extent) << " ";
    oss << JMIN(this->extent) << " ";
    oss << JMAX(this->extent) << " ";
    oss << KMIN(this->extent) << " ";
    oss << KMAX(this->extent) << " ";
    oss << "] ";
    oss << " size: " << this->size();
    return( oss.str() );
  }

};

//------------------------------------------------------------------------------

namespace ExtentUtilities {

#define MAX(ext,i) ext.extent[i*2+1]
#define MIN(ext,i) ext.extent[i*2]

/**
 * @brief Gets the periodicity information of the given extent along the 6
 * directions imin,imax,jmin,jmax,kmin,kmax.
 * @param ext the extent in query (in).
 * @param whole extent of the entire domain (in).
 * @param periodic the periodic tuple (out).
 * @note For each dimensions `d` \in {1,2,3}, periodic[d*2] and periodic[d*2+1]
 * indicates that the given extent, `ext`, is periodic across the min or max
 * boundary along the `d` dimension respectively.
 */
void GetPeriodicity(
    const extent_t& ext, const extent_t& whole, bool periodic[6])
{
  for(int dim=0; dim < 3; ++dim )
    {
    periodic[dim*2]   = (MIN(ext,dim)==MIN(whole,dim))? true : false;
    periodic[dim*2+1] = (MAX(ext,dim)==MAX(whole,dim))? true : false;
    } // END for all dimensions
}

/**
 * @brief Checks if x \in [a,b]
 * @param x the value in query
 * @param a lower range bound
 * @param b upper range bound
 * @return status true iff a <= x <= b, else, false.
 */
bool WithinRange(const uint64_t x, const uint64_t a, const uint64_t b)
{
  return( ((a <= x)&&(x <= b)) );
}

/**
 * @brief Looks at the orientation tuple of the given neighbor and determines
 * if it is neighboring with respect to a reference extent.
 * @param nei the RankNeighbor object in query.
 * @return
 * @note DetermineOrientation should be called first before calling this method.
 * @note Assumes 3-D topology.
 * @see DetermineOrientation().
 */
bool IsNeighbor(const gio::RankNeighbor& nei)
{
  if( nei.Orient[0]==gio::NOT_CONNECTED ||
      nei.Orient[1]==gio::NOT_CONNECTED ||
      nei.Orient[2]==gio::NOT_CONNECTED )
    {
    return false;
    }
  return true;
}


/**
 * @brief Determines the orientation of neiExtent w.r.t. myExtent.
 * @param wholeExtent extent of the entire block topology.
 * @param myExtent the extent of this process.
 * @param periodic indicates in which direction the given extent is periodic.
 * @param neiExtent the extent of a potentially neighboring process.
 * @param nei the RankNeighbor object.
 * @note Assumes 3-D topology.
 */
void DetermineOrientation(
      const extent_t& wholeExtent,
      const extent_t& myExtent,
      const bool periodic[6],
      const extent_t& neiExtent,
      gio::RankNeighbor& nei)
{
  // STEP 0: Get the length vector of the whole extent to use it when
  // determining periodic neighbors that wrap around the min of a dimension.
  uint64_t L[3];
  wholeExtent.GetLength(L);

  // STEP 1: Get the periodicity of the neighboring extent
  bool neiPeriodic[6];
  GetPeriodicity(neiExtent,wholeExtent,neiPeriodic);

  // STEP 2: Loop through each dimension and determine the orientation of
  // the potentially neighboring extent w.r.t. the given extent, myExtent
  for(int dim=0; dim < 3; ++dim)
    {
    uint64_t myMin  = MIN(myExtent,dim);
    uint64_t myMax  = MAX(myExtent,dim);
    uint64_t neiMin = MIN(neiExtent,dim);
    uint64_t neiMax = MAX(neiExtent,dim);

    // Check if we are neighboring on hi
    if( myMax+1 == neiMin )
      {

      if(neiPeriodic[dim*2+1] && periodic[dim*2] && ((myMin+L[dim])-1==neiMax))
        {
        nei.Orient[dim] = gio::ON_MAX_MIN_PERIODIC;
        } // END check min periodic
      else
        {
        nei.Orient[dim] = gio::ON_MAX;
        } // END else just ON_MAX

      } // END if ON_MAX
    else if( myMin == neiMax+1 )
      {

      if(neiPeriodic[dim*2] && periodic[dim*2+1] && ((myMax+1)%L[dim]==neiMin))
        {
        nei.Orient[dim] = gio::ON_MIN_MAX_PERIODIC;
        } // END check max periodic
      else
        {
        nei.Orient[dim] = gio::ON_MIN;
        }

      } // END if ON_MIN
    else if( WithinRange(myMin,neiMin,neiMax) ||
             WithinRange(myMax,neiMin,neiMax) ||
             WithinRange(neiMin,myMin,myMax)  ||
             WithinRange(neiMax,myMin,myMax))
      {
      nei.Orient[dim] = gio::OVERLAPPING;
      } // END if overlapping
    else if(neiPeriodic[dim*2+1] && periodic[dim*2] && ((myMin+L[dim])-1==neiMax))
      {
      nei.Orient[dim] = gio::MIN_PERIODIC;
      } // END check min periodic
    else if(neiPeriodic[dim*2] && periodic[dim*2+1] && ((myMax+1)%L[dim]==neiMin))
      {
      nei.Orient[dim] = gio::MAX_PERIODIC;
      } // END check max periodic
    else
      {
      nei.Orient[dim] = gio::NOT_CONNECTED;
      } // END else NOT_CONNECTED

    } // END for all dimensions

}

/**
 * @brief Loops through the extents and finds the neighbors for the given rank.
 * @param myRank the rank whose neighbors are in query (in).
 * @param wholeExtent the whole extent of the block topology (in).
 * @param extentList the list of extents, wherein, extentList[i] corresponds
 * to the extent assigned to process, P_i (in).
 * @param neighbors the list of rank neighbors (out).
 * @note Assumes 3-D topology.
 * @see RankNeighbor
 */
void FindRankNeighbors(
      const int myRank,
      const extent_t& wholeExtent,
      std::vector<extent_t>& extentList,
      std::vector<gio::RankNeighbor>& neighbors)
{
  assert("pre: rank is out-of-bounds!" && (myRank >= 0) &&
         (myRank < static_cast<int>(extentList.size())) );

  // STEP 0: Initialize neighbors
  neighbors.resize(0);

  // STEP 1: Get extent for this process/rank
  extent_t myExtent = extentList[myRank];

  // STEP 2: Get the periodicity tuple
  bool periodic[6];
  GetPeriodicity(myExtent, wholeExtent, periodic);

  // STEP 3: Loop through the extents of all other ranks and determine
  // if they are neighbors or not
  gio::RankNeighbor nei;
  for(size_t rank=0; rank < extentList.size(); ++rank)
    {
    if( myRank == static_cast<int>(rank) )
      {
      continue;
      }

    nei.RankID    = rank;
    nei.Orient[0] = nei.Orient[1] = nei.Orient[2] = gio::UNDEFINED;
    DetermineOrientation(wholeExtent,myExtent,periodic,extentList[rank],nei);

    if( IsNeighbor( nei ) )
      {
      neighbors.push_back( nei );
      } // END if IsNeighbor
    } // END for all extent ranks
}


} /* END namespace ExtentUtilities */

//------------------------------------------------------------------------------
//  END INTERNAL DATA STRUCTURES DEFINITIONS
//------------------------------------------------------------------------------

namespace gio
{

GenericIOUtilities::GenericIOUtilities()
{
  // TODO Auto-generated constructor stub
}

//------------------------------------------------------------------------------
GenericIOUtilities::~GenericIOUtilities()
{
  // TODO Auto-generated destructor stub
}

//------------------------------------------------------------------------------
void GenericIOUtilities::Printf(MPI_Comm comm, const char *fmt,...)
{
  // Sanity check
  assert("pre: MPI communicator should not be NULL!" &&
          (comm != MPI_COMM_NULL) );
  assert("pre: NULL formatted message!" && (fmt != NULL) );

  int rank = -1;
  MPI_Comm_rank(comm,&rank);
  if( rank==0 )
    {
    va_list argptr;
    va_start(argptr,fmt);
    vprintf(fmt,argptr);
    fflush(stdout);
    va_end(argptr);
    }
  MPI_Barrier(comm);
}

//------------------------------------------------------------------------------
void GenericIOUtilities::SynchronizedPrintf(MPI_Comm comm, const char *fmt,...)
{
  // Sanity check
  assert("pre: MPI communicator should not be NULL!" &&
          (comm != MPI_COMM_NULL) );
  assert("pre: NULL formatted message!" && (fmt != NULL) );

  int rank     = -1;
  int numRanks = -1;
  MPI_Comm_rank(comm,&rank);
  MPI_Comm_size(comm,&numRanks);
  MPI_Request nullRequest = MPI_REQUEST_NULL;

  if( rank == 0 )
    {
    // This is the first rank

    // STEP 0: Print message
    va_list argptr;
    va_start(argptr,fmt);
    printf("# [%d]: ",rank);
    vprintf(fmt,argptr);
    fflush(stdout);
    va_end(argptr);

    // STEP 1: Tell next process to print, if any
    if( numRanks > 1)
      {
      MPI_Isend(NULL,0,MPI_INT,rank+1,0,comm,&nullRequest);
      }
    }
  else if( rank == numRanks-1 )
    {
    // Last rank

    // STEP 0: Block until previous process completes
    MPI_Recv(NULL,0,MPI_INT,rank-1,MPI_ANY_TAG,comm,MPI_STATUS_IGNORE);

    // STEP 1: Print message
    va_list argptr;
    va_start(argptr,fmt);
    printf("# [%d]: ",rank);
    vprintf(fmt,argptr);
    fflush(stdout);
    va_end(argptr);
    }
  else
    {
    // STEP 0: Block until previous process completes
    MPI_Recv(NULL,0,MPI_INT,rank-1,MPI_ANY_TAG,comm,MPI_STATUS_IGNORE);

    // STEP 1: Print message
    va_list argptr;
    va_start(argptr,fmt);
    printf("# [%d]: ",rank);
    vprintf(fmt,argptr);
    fflush(stdout);
    va_end(argptr);

    // STEP 2: tell next process to print
    MPI_Isend(NULL,0,MPI_INT,rank+1,0,comm,&nullRequest);
    }
}

//------------------------------------------------------------------------------
bool GenericIOUtilities::VerifyChecksum(const void* data, size_t nbytes)
{
  if( data==NULL || nbytes==0)
    {
    return true;
    }

  if(crc64_omp(data,nbytes)!=(uint64_t)(-1))
    {
    return false;
    }
  return true;
}

//------------------------------------------------------------------------------
bool GenericIOUtilities::VerifyChecksum(
      const void* data, size_t nbytes, uint64_t cs, bool checkInverse)
{
  if( data==NULL || nbytes==0)
    {
    return true;
    }

  // Compute the checksum
  uint64_t crc64 = crc64_omp(data,nbytes);

  if( checkInverse )
    {
    // Get the inverse checksum, because, that is what is stored in the file
    // and what we are comparing against
    uint64_t crc64_inv;
    crc64_invert(crc64,&crc64_inv);

    // Compare against the inverse checksum, cs, read from the file.
    if(crc64_inv == cs )
      {
      return true;
      } // END if inverse checksum match

    } // END if check against inverse CRC
  else
    {

    if( crc64 == cs )
      {
      return true;
      } // END If checksum match!

    } // END else check against original CRC

  return false;
}

//------------------------------------------------------------------------------
void GenericIOUtilities::SwapGlobalHeader(GlobalHeader *GH)
{
  assert("pre: GH != NULL" && (GH != NULL) );

  std::vector<char> swapBuffer;
  swapBuffer.resize( sizeof(char) );
  for(unsigned int i=0; i < MagicSize; ++i )
    {
    GenericIOUtilities::SwapEndian(
        &GH->Magic[i],sizeof(char),&swapBuffer[0]);
    }

  swapBuffer.resize( sizeof(uint64_t) );

  GenericIOUtilities::SwapEndian(
      &GH->HeaderSize,sizeof(uint64_t),&swapBuffer[0]);

  GenericIOUtilities::SwapEndian(
      &GH->NElems,sizeof(uint64_t),&swapBuffer[0]);

  GenericIOUtilities::SwapEndian(
      &GH->Dims[0],sizeof(uint64_t),&swapBuffer[0]);
  GenericIOUtilities::SwapEndian(
      &GH->Dims[1],sizeof(uint64_t),&swapBuffer[0]);
  GenericIOUtilities::SwapEndian(
      &GH->Dims[2],sizeof(uint64_t),&swapBuffer[0]);

  GenericIOUtilities::SwapEndian(
      &GH->NVars,sizeof(uint64_t),&swapBuffer[0]);

  GenericIOUtilities::SwapEndian(
      &GH->VarsSize,sizeof(uint64_t),&swapBuffer[0]);

  GenericIOUtilities::SwapEndian(
      &GH->VarsStart,sizeof(uint64_t),&swapBuffer[0]);

  GenericIOUtilities::SwapEndian(
      &GH->NRanks,sizeof(uint64_t),&swapBuffer[0]);

  GenericIOUtilities::SwapEndian(
      &GH->RanksSize,sizeof(uint64_t),&swapBuffer[0]);

  GenericIOUtilities::SwapEndian(
      &GH->RanksStart,sizeof(uint64_t),&swapBuffer[0]);

  GenericIOUtilities::SwapEndian(
      &GH->GlobalHeaderSize,sizeof(uint64_t),&swapBuffer[0]);

  swapBuffer.resize( sizeof(double) );

  GenericIOUtilities::SwapEndian(
      &GH->PhysOrigin[0],sizeof(double),&swapBuffer[0]);
  GenericIOUtilities::SwapEndian(
      &GH->PhysOrigin[1],sizeof(double),&swapBuffer[0]);
  GenericIOUtilities::SwapEndian(
      &GH->PhysOrigin[2],sizeof(double),&swapBuffer[0]);

  GenericIOUtilities::SwapEndian(
      &GH->PhysScale[0],sizeof(double),&swapBuffer[0]);
  GenericIOUtilities::SwapEndian(
      &GH->PhysScale[1],sizeof(double),&swapBuffer[0]);
  GenericIOUtilities::SwapEndian(
      &GH->PhysScale[2],sizeof(double),&swapBuffer[0]);

  swapBuffer.resize( sizeof(uint64_t) );

  GenericIOUtilities::SwapEndian(
      &GH->BlocksSize,sizeof(uint64_t),&swapBuffer[0]);
  GenericIOUtilities::SwapEndian(
      &GH->BlocksStart,sizeof(uint64_t),&swapBuffer[0]);
}

//------------------------------------------------------------------------------
void GenericIOUtilities::SwapVariableHeader(VariableHeader *VH)
{
  assert("pre: VH != NULL" && (VH != NULL) );
  std::vector<char> swapBuffer;
  swapBuffer.resize(sizeof(char));

  for(unsigned int i=0; i < NameSize; ++i)
    {
    GenericIOUtilities::SwapEndian(
        &VH->Name[i],sizeof(char),&swapBuffer[0]);
    }

  swapBuffer.resize(sizeof(uint64_t));
  GenericIOUtilities::SwapEndian(
      &VH->Flags,sizeof(uint64_t),&swapBuffer[0]);

  GenericIOUtilities::SwapEndian(
        &VH->Size,sizeof(uint64_t),&swapBuffer[0]);
}

//------------------------------------------------------------------------------
void GenericIOUtilities::SwapRankHeader(RankHeader *RH)
{
  assert("pre: RH != NULL" && (RH != NULL) );

  std::vector<char> swapBuffer;
  swapBuffer.resize( sizeof(uint64_t) );

  GenericIOUtilities::SwapEndian(
      &RH->Coords[0],sizeof(uint64_t),&swapBuffer[0]);
  GenericIOUtilities::SwapEndian(
      &RH->Coords[1],sizeof(uint64_t),&swapBuffer[0]);
  GenericIOUtilities::SwapEndian(
      &RH->Coords[2],sizeof(uint64_t),&swapBuffer[0]);

  GenericIOUtilities::SwapEndian(
      &RH->NElems,sizeof(uint64_t),&swapBuffer[0]);

  GenericIOUtilities::SwapEndian(
      &RH->Start,sizeof(uint64_t),&swapBuffer[0]);

  GenericIOUtilities::SwapEndian(
      &RH->GlobalRank,sizeof(uint64_t),&swapBuffer[0]);
}

//------------------------------------------------------------------------------
void GenericIOUtilities::SwapBlockInfo(BlockInfo* blockInfo)
{
  assert("pre: blockInfo != NULL" && (blockInfo != NULL) );

  std::vector<char> swapBuffer;

  // Swap the filters name
  for(size_t f=0; f < gio::MaxFilters; ++f)
    {
    for(size_t k=0; k < gio::FilterNameSize; ++k)
      {
      GenericIOUtilities::SwapEndian(
          &blockInfo->Filters[f][k],sizeof(char),&swapBuffer[0]);
      } // END for all characters
    } // END for all filters

  swapBuffer.resize(sizeof(uint64_t));
  GenericIOUtilities::SwapEndian(
      &blockInfo->Start,sizeof(uint64_t),&swapBuffer[0]);
  GenericIOUtilities::SwapEndian(
      &blockInfo->Size,sizeof(uint64_t),&swapBuffer[0]);
}

//------------------------------------------------------------------------------
void GenericIOUtilities::SwapCompressHeader(CompressHeader* CH)
{
  assert("pre: CH != NULL" && (CH != NULL) );

  std::vector<char> swapBuffer;
  swapBuffer.resize(sizeof(uint64_t));
  GenericIOUtilities::SwapEndian(&CH->OrigCRC,sizeof(uint64_t),&swapBuffer[0]);
}

//------------------------------------------------------------------------------
bool GenericIOUtilities::DoesFileEndianMatch(GlobalHeader *GH)
{
  assert("pre: GH != NULL" && (GH != NULL)  );

  const char *Magic =
       GenericIOUtilities::IsBigEndian() ? MagicBE : MagicLE;
  const char *MagicInv =
       GenericIOUtilities::IsBigEndian() ? MagicLE : MagicBE;

  std::string magicString =
     std::string(GH->Magic,GH->Magic+MagicSize-1);

  if( magicString != Magic)
   {
   if( magicString == MagicInv )
     {
     return false;
     } // END if swap
   else
     {
     std::cerr << "ERROR: Could not detect file endian!\n";
     } // END else
   } // END if endian does not match file endian
  return true;
}

//------------------------------------------------------------------------------
void GenericIOUtilities::RCBPartition(
      const int processID,
      const int numProcesses,
      const uint64_t wholeExtent[6],
      uint64_t subExtent[6],
      std::vector<RankNeighbor> &neighbors)
{
  // STEP 0: If there is a single process, all blocks will be assi
  if(numProcesses == 1)
    {
    memcpy(subExtent,wholeExtent,6*sizeof(uint64_t));
    return;
    }

  // STEP 1: Create an extent object for the whole extent
  extent_t globalExtent;
  memcpy(globalExtent.extent,wholeExtent,6*sizeof(uint64_t));

  // STEP 2: Create extent priority queue and vector containers
  std::priority_queue<extent_t> extentQueue;
  std::vector<extent_t> &extentList = Container(extentQueue);

  // STEP 3: Push whole extent into the queue
  extentQueue.push( globalExtent );
  assert("pre: extent queue and list out-of-synch!" &&
         (extentQueue.size()==extentList.size()) );

  // STEP 4: Split extents iteratively to match the number of processes
  extent_t s1,s2;
  while( static_cast<int>(extentQueue.size()) < numProcesses )
    {
    extent_t t = extentQueue.top();
    extentQueue.pop();
    t.split(s1,s2);

    extentQueue.push(s1);
    extentQueue.push(s2);
    } // END while
  assert("pre: extent queue size does not match requested num partitions!" &&
          (static_cast<int>(extentQueue.size())==numProcesses));
  assert("pre: extent queue and list out-of-synch!" &&
           (static_cast<int>(extentQueue.size())==extentList.size()) );

  // STEP 5: get sub-extent for this process
  memcpy(subExtent,extentList[processID].extent,6*sizeof(uint64_t));

  // STEP 6: Find rank neighbors for this process
  ExtentUtilities::FindRankNeighbors(
      processID,globalExtent,extentList,neighbors);
}

//------------------------------------------------------------------------------
void GenericIOUtilities::RCBAssignment(
      const int processID,
      const int numProcesses,
      const uint64_t blkdims[3],
      ProcessCoordsToRank* ijkMap,
      std::vector<int> &assigned,
      std::vector<RankNeighbor> &neighbors)
{
  assert("pre: assigned.size()==0" && (assigned.size()==0) );

  // STEP 0: Form whole-extent
  // TODO: Handle different orientations, IJ,IK, JK, etc.
  uint64_t wholeExtent[6];
  IMIN(wholeExtent) = JMIN(wholeExtent) = KMIN(wholeExtent) = 0;
  IMAX(wholeExtent) = blkdims[0]-1;
  JMAX(wholeExtent) = blkdims[1]-1;
  KMAX(wholeExtent) = blkdims[2]-1;

  // STEP 1: Partition the whole-extent & get the extent assigned to this rank.
  uint64_t subExtent[6];
  GenericIOUtilities::RCBPartition(
      processID, numProcesses, wholeExtent, subExtent, neighbors);

//  GenericIOUtilities::SynchronizedPrintf(MPI_COMM_WORLD,
//      "\nAssigned Extent:[%d %d %d %d %d %d]\n",
//      IMIN(subExtent),IMAX(subExtent),
//      JMIN(subExtent),JMAX(subExtent),
//      KMIN(subExtent),KMAX(subExtent)
//      );

  // STEP 3: Process sub-extent and get corresponding blocks
  uint64_t coords[3];
  for(uint64_t i=IMIN(subExtent); i <= IMAX(subExtent); ++i)
    {
    for(uint64_t j=JMIN(subExtent); j <= JMAX(subExtent); ++j)
      {
      for(uint64_t k=KMIN(subExtent); k <= KMAX(subExtent); ++k)
        {
        I(coords)=i; J(coords)=j; K(coords)=k;
        std::string hashCode =
            GenericIOUtilities::GetHashCodeFromCoords(i,j,k);
        ProcessCoordsToRank::iterator it = ijkMap->find(hashCode);
        assert("pre: cannot find IJK in map!\n" && it != ijkMap->end() );
        uint64_t linearIdx = it->second;
        assigned.push_back(linearIdx);
        } // END for all k
      } // END for all j
    } // END for all i

}

//------------------------------------------------------------------------------
void GenericIOUtilities::RoundRobin(
    const int bucketIdx, const int NumBuckets,
    const int NumBlocks,
    std::vector<int> &assigned)
{
  // Sanity checks!
  assert("pre: NumBuckets < NumBlocks" && (NumBuckets < NumBlocks) );
  assert("pre: assigned.size()==0" && (assigned.size()==0) );

  for(int blkIdx=0; blkIdx < NumBlocks; ++blkIdx)
    {
    if( (blkIdx%NumBuckets) == bucketIdx )
      {
      assigned.push_back(blkIdx);
      } // END if
    } // END for all blocks
}

//------------------------------------------------------------------------------
void GenericIOUtilities::BlockAssignment(
          const int processID,
          const int numProcesses,
          const int totalNumberOfBlocks,
          const uint64_t blockDims[3],
          const int assignmentStrategy,
          ProcessCoordsToRank* ijkMap,
          std::vector<int> &assigned,
          std::vector<RankNeighbor> &neighbors)
{
  // Sanity checks!
  assert("pre: processID is out-of-bounds!" &&
          (processID >= 0) && (processID < numProcesses) );
  assert("pre: totalNumberOfBlocks > 0" && (totalNumberOfBlocks > 0) );

  // TODO: must handle different orientations, i.e., IJ, IK, JK (?)
//  assert("pre: block decomposition does not match number of blocks" &&
//   (totalNumberOfBlocks == (blockDims[0]*blockDims[1]*blockDims[2])));
  assert("pre: invalid assignment strategy" &&
          (assignmentStrategy < NUM_BLOCK_ASSIGNMENTS) );


  assigned.clear();

  if(numProcesses < totalNumberOfBlocks)
    {
    switch(assignmentStrategy)
      {
      case RR_BLOCK_ASSIGNMENT:
        GenericIOUtilities::RoundRobin(
            processID,numProcesses,totalNumberOfBlocks,assigned);
        neighbors.resize(0);
        break;
      case RCB_BLOCK_ASSIGNMENT:
        if((totalNumberOfBlocks==blockDims[0]) &&
            (blockDims[1]==1) && (blockDims[2]==1) )
          {
          // file was not written using a cartesian communicator, hence,
          // we cannot utilize RCB assignment. Resorting to RR_BLOCK_ASSIGNMENT
          GenericIOUtilities::RoundRobin(
             processID,numProcesses,totalNumberOfBlocks,assigned);
          neighbors.resize(0);
          } // END if
        else
          {
          GenericIOUtilities::RCBAssignment(
              processID,numProcesses,blockDims,ijkMap,assigned,neighbors);
          } // END else
        break;
      default:
        throw std::runtime_error("Invalid block assignment strategy!");
      } // END switch
    } // END if more blocks than processes
  else if( numProcesses > totalNumberOfBlocks)
    {
    // TODO: Hmm...should we determine the neighbors in this case?
    //one-to-one, not all processes will have a block
    if( processID < totalNumberOfBlocks )
      {
      assigned.push_back(processID);
      }
    } // END more processes than blocks
  else
    {
    // TODO: Hmm...should we determine the neighbors in this case?
    // one-to-one, numProcess == totalNumberOfBlocks
    assigned.push_back(processID);
    } // END one-to-one
}

//------------------------------------------------------------------------------
void GenericIOUtilities::SwapEndian(
        void *Addr, const int N, char *buf)
{
  assert("pre: Addr is NULL!" && (Addr != NULL) );
  assert("pre: N > 0" && (N > 0) );

  // STEP 0: Setup buffer. Either the user has supplied a swap buffer
  // externally to use, or we will handle the buffering internally.
  bool bufferedExternally = true;
  if( buf == NULL )
    {
    bufferedExternally = false;
    buf = new char[N];
    }

  // STEP 1: swap bytes
  for(int srcOffSet=N-1, idx=0; srcOffSet >= 0; --srcOffSet, ++idx)
    {
    buf[idx] = *( (char*)Addr+srcOffSet);
    } // END for all bytes

  // STEP 2: Copy data to original
  memcpy(Addr, (void*)buf, N);

  // STEP 3: if buffering internally, clean up
  if( !bufferedExternally )
    {
    delete [] buf;
    }
}

//------------------------------------------------------------------------------
int GenericIOUtilities::DetectVariablePrimitiveType(
      const VariableInfo &vinfo)
{
  int type = -1;
  if( vinfo.IsFloat )
    {
    if( vinfo.Size == sizeof(float) )
      {
      type = GENERIC_IO_FLOAT_TYPE;
      }
    else if( vinfo.Size == sizeof(double) )
      {
      type = GENERIC_IO_DOUBLE_TYPE;
      }
    else
      {
      std::cerr << "WARNING: Cannot detect floating point variable type!\n";
      }
    } // END if variable is floating point
  else
    {
    if( vinfo.IsSigned )
      {
//      if( vinfo.Size == sizeof(short) )
//        {
//        type = GENERIC_IO_SHORT_TYPE;
//        }
//      else if( vinfo.Size == sizeof(long) )
//        {
//        type = GENERIC_IO_LONG_TYPE;
//        }
//      else if( vinfo.Size == sizeof(long long) )
//        {
//        type = GENERIC_IO_LONG_LONG_TYPE;
//        }
//      else if( vinfo.Size == sizeof(int32_t) )
      if( vinfo.Size == sizeof(int32_t) )
        {
        type = GENERIC_IO_INT32_TYPE;
        }
      else if( vinfo.Size == sizeof(int64_t) )
        {
        type = GENERIC_IO_INT64_TYPE;
        }
      else if( vinfo.Size == sizeof(float) )
        {
        type = GENERIC_IO_FLOAT_TYPE;
        }
      else if( vinfo.Size == sizeof(double) )
        {
        type = GENERIC_IO_DOUBLE_TYPE;
        }
      else
        {
        std::cerr << "WARNING: Cannot detect signed integer type!";
        }
      } // END if variable is signed integer
    else
      {
      if(vinfo.Size == sizeof(uint16_t))
        {
        type = GENERIC_IO_UINT16_TYPE;
        }
      else if(vinfo.Size == sizeof(uint32_t) )
        {
        type = GENERIC_IO_UINT32_TYPE;
        }
      else if( vinfo.Size == sizeof(uint64_t) )
        {
        type = GENERIC_IO_UINT64_TYPE;
        }
      else if(vinfo.Size == sizeof(unsigned char))
        {
        type = GENERIC_IO_UCHAR_TYPE;
        }
      else
        {
        std::cerr << "WARNING: Cannot detect unsigned integer type!";
        }
      } // END if variable is unsigned inter
    } // END if variable is integer type
  return( type );
}

//------------------------------------------------------------------------------
void GenericIOUtilities::DeallocateVariableArray(
      const VariableInfo &vinfo, void* dataArray)
{
  if( dataArray == NULL )
    {
    return;
    }

  int type = GenericIOUtilities::DetectVariablePrimitiveType( vinfo );
  switch(type)
    {
    case GENERIC_IO_UCHAR_TYPE:
      {
      unsigned char* ptr = static_cast<unsigned char*>(dataArray);
      delete [] ptr;
      ptr = NULL;
      }
      break;
    case GENERIC_IO_INT32_TYPE:
      {
      int32_t* ptr = static_cast<int32_t*>(dataArray);
      delete [] ptr;
      ptr = NULL;
      }
      break;
    case GENERIC_IO_INT64_TYPE:
      {
      int64_t* ptr = static_cast<int64_t*>(dataArray);
      delete [] ptr;
      ptr = NULL;
      }
      break;
    case GENERIC_IO_UINT16_TYPE:
      {
      uint16_t* ptr = static_cast<uint16_t*>(dataArray);
      delete [] ptr;
      ptr = NULL;
      }
      break;
    case GENERIC_IO_UINT32_TYPE:
      {
      uint32_t* ptr = static_cast<uint32_t*>(dataArray);
      delete [] ptr;
      ptr = NULL;
      }
      break;
    case GENERIC_IO_UINT64_TYPE:
      {
      uint64_t* ptr = static_cast<uint64_t*>(dataArray);
      delete [] ptr;
      ptr = NULL;
      }
      break;
    case GENERIC_IO_DOUBLE_TYPE:
      {
      double* ptr = static_cast<double*>(dataArray);
      delete [] ptr;
      ptr = NULL;
      }
      break;
    case GENERIC_IO_FLOAT_TYPE:
      {
      float* ptr = static_cast<float*>(dataArray);
      delete [] ptr;
      ptr = NULL;
      }
      break;
    default:
      {
      throw std::runtime_error("undefined variable datatype!");
      }
    } // END switch

  dataArray = NULL;
}

//------------------------------------------------------------------------------
void* GenericIOUtilities::AllocateVariableArray(
        const VariableInfo &vinfo, const int numElements, bool pad )
{
  int  type = GenericIOUtilities::DetectVariablePrimitiveType( vinfo );
  void *ptr = NULL;
  switch( type )
    {
//    case GENERIC_IO_SHORT_TYPE:
//      {
//      short *data = new short[ numElements ];
//      ptr = static_cast<void*>(data);
//      }
//      break;
//    case GENERIC_IO_LONG_TYPE:
//      {
//      long *data = new long[ numElements ];
//      ptr = static_cast<void*>(data);
//      }
//      break;
//    case GENERIC_IO_LONG_LONG_TYPE:
//      {
//      long long *data = new long long[ numElements];
//      ptr = static_cast<void*>(data);
//      }
//      break;
    case GENERIC_IO_UCHAR_TYPE:
      {
      int padsize = (pad)? gio::CRCSize/sizeof(unsigned char) : 0;
      unsigned char *data = new unsigned char[ numElements+padsize ];
      ptr = static_cast<void*>(data);
      }
      break;
    case GENERIC_IO_INT32_TYPE:
      {
      int padsize = (pad)? gio::CRCSize/sizeof(int32_t) : 0;
      int32_t *data = new int32_t[ numElements+padsize ];
      ptr = static_cast<void*>(data);
      }
      break;
    case GENERIC_IO_INT64_TYPE:
      {
      int padsize = (pad)? gio::CRCSize/sizeof(int64_t) : 0;
      int64_t *data = new int64_t[ numElements+padsize ];
      ptr = static_cast<void*>(data);
      }
      break;
    case GENERIC_IO_UINT16_TYPE:
      {
      int padsize = (pad)? gio::CRCSize/sizeof(uint16_t) : 0;
      uint16_t *data = new uint16_t[ numElements+padsize ];
      ptr = static_cast<void*>(data);
      }
      break;
    case GENERIC_IO_UINT32_TYPE:
      {
      int padsize = (pad)? gio::CRCSize/sizeof(uint32_t) : 0;
      uint32_t *data = new uint32_t[ numElements+padsize ];
      ptr = static_cast<void*>(data);
      }
      break;
    case GENERIC_IO_UINT64_TYPE:
      {
      int padsize = (pad)? gio::CRCSize/sizeof(uint64_t) : 0;
      uint64_t *data = new uint64_t[ numElements+padsize ];
      ptr = static_cast<void*>(data);
      }
      break;
    case GENERIC_IO_DOUBLE_TYPE:
      {
      int padsize = (pad)? gio::CRCSize/sizeof(double) : 0;
      double *data = new double[ numElements+padsize ];
      ptr = static_cast<void*>(data);
      }
      break;
    case GENERIC_IO_FLOAT_TYPE:
      {
      int padsize = (pad)? gio::CRCSize/sizeof(float) : 0;
      float *data = new float[ numElements+padsize ];
      ptr = static_cast<void*>(data);
      }
      break;
    default:
      ptr = NULL;
    }
  return( ptr );
}

} /* namespace cosmotk */
