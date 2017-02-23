/**
 * @brief Definitions of internal datal-structures and globals used by the
 * GenericIO framework.
 */
#ifndef GENERICIODEFINITIONS_HPP_
#define GENERICIODEFINITIONS_HPP_

#include <cstddef>  // For size_t
#include <map>      // For STL map
#include <stdint.h> // For explicit 32-bit/64-bit integer types
#include <string>   // For STL string

#include <mpi.h>

/**
 * @brief A macro that expands to setter & getter routines for the given ivar.
 * @param ivar the name of the ivar
 * @param type the type of the ivar
 */
#define GIOGetNSetMacro(ivar, type)       \
    virtual type Get##ivar() const {      \
      return this->ivar;                  \
    }                                     \
    virtual void Set##ivar( type __arg) { \
        this->ivar = __arg;               \
    }

/**
 * @brief A macro that expands to setter & getter routines for the given ivar.
 * @param ivar the name of the ivar
 * @param type the type of the ivar
 */
#define GIOGetMacro(ivar, type)       \
    virtual type Get##ivar() const { \
      return this->ivar;              \
    }

/**
 * @brief A macro to disable the copy constructor and assignment of a class.
 * @note Must be used within the private block of the class definition.
 */
#define GIO_DISABLE_COPY_AND_ASSIGNMENT(ClassName)    \
  ClassName(const ClassName&);                    \
  void operator=(const ClassName&);

namespace gio {

typedef std::map<std::string,uint64_t> ProcessCoordsToRank;

/**
 * @brief An enum of supported block assignment strategies
 */
enum GenericIOBlockAssignment {
  RR_BLOCK_ASSIGNMENT,  //!< round-robin block assignment
  RCB_BLOCK_ASSIGNMENT, //!< rcb block assignment (spatial aware)
  NUM_BLOCK_ASSIGNMENTS //!< number of block assignments
};

/**
 * @brief An enum defining the different neighbor relationships.
 * @see RankNeighbor
 */
enum NeighborOrientation {

  ON_MIN,              //!< neighbor is abutting at the min end

  ON_MAX,              //!< neighbor is abutting at the max end

  OVERLAPPING,         //!< neighbor extents are overlapping

  ON_MAX_MIN_PERIODIC, //!< neighbor is abutting on the max end and along
                       //!< the min periodic boundary.

  ON_MIN_MAX_PERIODIC, //!< neighbor is abutting on the min end and along
                       //!< the max periodic boundary.

  MIN_PERIODIC,        //!< Abutting along the min periodic boundary

  MAX_PERIODIC,        //!< Abutting along the max periodic boundary

  NOT_CONNECTED,       //!< neighbor is not connected

  UNDEFINED,           //!< neighbor relationship has not been defined
};

static const char* NEIGHBOR_ORIENTATION[] = {
    "ON_MIN",
    "ON_MAX",
    "OVERLAPPING",
    "ON_MAX_MIN_PERIODIC",
    "ON_MIN_MAX_PERIODIC",
    "MIN_PERIODIC",
    "MAX_PERIODIC",
    "NOT_CONNECTED",
    "UNDEFINED"
};

/**
 * @brief A struct to define a RankNeighbor.
 * @note In the case of RCB block distribution, each rank, r_i, is assigned
 * a sub-extent, s_i. Each s_i has a set of neighbors whose extents are
 * abutting and the data is spatially contiguous. A neighbor is identified
 * by the "RankID" of the neighboring process and an orientation tuple, "Orient"
 * that indicates how the neighbor is connected in relationship to given rank
 * along each dimension i,j,k.
 * @see NeighborOrientation
 */
struct RankNeighbor {
  int RankID;
  int Orient[3];
};


/**
 * @brief An enum of supported primitive types of the generic-io infrastructure.
 */
enum GenericIOPrimitiveTypes {
/*  GENERIC_IO_SHORT_TYPE,    //!< GENERIC_IO_SHORT_TYPE
  GENERIC_IO_LONG_TYPE,     //!< GENERIC_IO_LONG_TYPE
  GENERIC_IO_LONG_LONG_TYPE,//!< GENERIC_IO_LONG_LONG_TYPE */
  GENERIC_IO_UCHAR_TYPE,    //!< GENERIC_IO_UCHART_TYPE
  GENERIC_IO_INT32_TYPE,    //!< GENERIC_IO_INT32_TYPE
  GENERIC_IO_INT64_TYPE,    //!< GENERIC_IO_INT64_TYPE
  GENERIC_IO_UINT16_TYPE,   //!< GENERIC_IO_UINT16_TYPE
  GENERIC_IO_UINT32_TYPE,   //!< GENERIC_IO_UINT32_TYPE
  GENERIC_IO_UINT64_TYPE,   //!< GENERIC_IO_UINT64_TYPE
  GENERIC_IO_DOUBLE_TYPE,   //!< GENERIC_IO_DOUBLE_TYPE
  GENERIC_IO_FLOAT_TYPE,    //!< GENERIC_IO_FLOAT_TYPE
  NUM_PRIMITIVE_TYPES       //!< NUM_PRIMITIVE_TYPES
};

/**
 * @brief String representation of generic-io primitive types.
 * @note Used mostly for debugging
 */
static const char* PRIMITIVE_NAME[] = {
/*  "GENERIC_IO_SHORT_TYPE",
  "GENERIC_IO_LONG_TYPE",
  "GENERIC_IO_LONG_LONG_TYPE", */
  "GENERIC_IO_UCHAR_TYPE",
  "GENERIC_IO_INT32_TYPE",
  "GENERIC_IO_INT64_TYPE",
  "GENERIC_IO_UINT16_TYPE",
  "GENERIC_IO_UINT32_TYPE",
  "GENERIC_IO_UINT64_TYPE",
  "GENERIC_IO_DOUBLE_TYPE",
  "GENERIC_IO_FLOAT_TYPE"
};

static const size_t CRCSize   = 8;
static const size_t MagicSize = 8;
static const char *MagicBE    = "HACC01B";
static const char *MagicLE    = "HACC01L";
static const size_t NameSize  = 256;


struct VariableHeader {
  char Name[NameSize];
  uint64_t Flags;
  uint64_t Size;
} __attribute__((packed));

struct RankHeader {
  uint64_t Coords[3];
  uint64_t NElems;
  uint64_t Start;
  uint64_t GlobalRank;
} __attribute__((packed));

struct GlobalHeader {
  char Magic[MagicSize];
  uint64_t HeaderSize;
  uint64_t NElems; // The global total
  uint64_t Dims[3];
  uint64_t NVars;
  uint64_t VarsSize;
  uint64_t VarsStart;
  uint64_t NRanks;
  uint64_t RanksSize;
  uint64_t RanksStart;
  uint64_t GlobalHeaderSize;
  double   PhysOrigin[3];
  double   PhysScale[3];
  uint64_t BlocksSize;
  uint64_t BlocksStart;
} __attribute__((packed));

static const size_t FilterNameSize = 8;
static const size_t MaxFilters     = 4;

struct BlockInfo {
  char Filters[MaxFilters][FilterNameSize];
  uint64_t Start;
  uint64_t Size;
} __attribute__((packed));

struct CompressHeader {
  uint64_t OrigCRC;
} __attribute__((packed));

struct VariableInfo
{
  VariableInfo()
    {
    this->Name = "";
    this->Size = 0;
    this->IsFloat  = true;
    this->IsFloat  = true;
    this->IsSigned = true;
    this->IsPhysCoordX = false;
    this->IsPhysCoordY = false;
    this->IsPhysCoordZ = false;
    this->MaybePhysGhost = false;
    }

  VariableInfo(const std::string &N, std::size_t S, bool IF, bool IS,
               bool PCX, bool PCY, bool PCZ, bool PG)
    : Name(N), Size(S), IsFloat(IF), IsSigned(IS),
      IsPhysCoordX(PCX), IsPhysCoordY(PCY), IsPhysCoordZ(PCZ),
      MaybePhysGhost(PG) {}

  std::string Name;
  std::size_t Size;
  bool IsFloat;
  bool IsSigned;
  bool IsPhysCoordX, IsPhysCoordY, IsPhysCoordZ;
  bool MaybePhysGhost;
};

}


#endif /* GENERICIODEFINITIONS_HPP_ */
