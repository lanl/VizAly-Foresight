/*
 *                    Copyright (C) 2015, UChicago Argonne, LLC
 *                               All Rights Reserved
 * 
 *                               Generic IO (ANL-15-066)
 *                     Hal Finkel, Argonne National Laboratory
 * 
 *                              OPEN SOURCE LICENSE
 * 
 * Under the terms of Contract No. DE-AC02-06CH11357 with UChicago Argonne,
 * LLC, the U.S. Government retains certain rights in this software.
 * 
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions are met:
 * 
 *   1. Redistributions of source code must retain the above copyright notice,
 *      this list of conditions and the following disclaimer.
 * 
 *   2. Redistributions in binary form must reproduce the above copyright
 *      notice, this list of conditions and the following disclaimer in the
 *      documentation and/or other materials provided with the distribution.
 * 
 *   3. Neither the names of UChicago Argonne, LLC or the Department of Energy
 *      nor the names of its contributors may be used to endorse or promote
 *      products derived from this software without specific prior written
 *      permission.
 * 
 * *****************************************************************************
 * 
 *                                  DISCLAIMER
 * THE SOFTWARE IS SUPPLIED “AS IS” WITHOUT WARRANTY OF ANY KIND.  NEITHER THE
 * UNTED STATES GOVERNMENT, NOR THE UNITED STATES DEPARTMENT OF ENERGY, NOR
 * UCHICAGO ARGONNE, LLC, NOR ANY OF THEIR EMPLOYEES, MAKES ANY WARRANTY,
 * EXPRESS OR IMPLIED, OR ASSUMES ANY LEGAL LIABILITY OR RESPONSIBILITY FOR THE
 * ACCURACY, COMPLETENESS, OR USEFULNESS OF ANY INFORMATION, DATA, APPARATUS,
 * PRODUCT, OR PROCESS DISCLOSED, OR REPRESENTS THAT ITS USE WOULD NOT INFRINGE
 * PRIVATELY OWNED RIGHTS.
 * 
 * *****************************************************************************
 */

#ifndef GENERICIO_H
#define GENERICIO_H

#include <cstdlib>
#include <vector>
#include <string>
#include <iostream>
#include <limits>
#include <stdint.h>

#ifndef GENERICIO_NO_MPI
#include <mpi.h>
#else
#include <fstream>
#endif

#include <unistd.h>

namespace gio {

class GenericFileIO {
public:
  virtual ~GenericFileIO() {}

public:
  virtual void open(const std::string &FN, bool ForReading = false, bool MustExist = false) = 0;
  virtual void setSize(size_t sz) = 0;
  virtual void read(void *buf, size_t count, off_t offset,
                    const std::string &D) = 0;
  virtual void write(const void *buf, size_t count, off_t offset,
                     const std::string &D) = 0;

protected:
  std::string FileName;
};

#ifndef GENERICIO_NO_MPI
class GenericFileIO_MPI : public GenericFileIO {
public:
  GenericFileIO_MPI(const MPI_Comm &C) : FH(MPI_FILE_NULL), Comm(C) {}
  virtual ~GenericFileIO_MPI();

public:
  virtual void open(const std::string &FN, bool ForReading = false, bool MustExist = false);
  virtual void setSize(size_t sz);
  virtual void read(void *buf, size_t count, off_t offset, const std::string &D);
  virtual void write(const void *buf, size_t count, off_t offset, const std::string &D);

protected:
  MPI_File FH;
  MPI_Comm Comm;
};

class GenericFileIO_MPICollective : public GenericFileIO_MPI {
public:
  GenericFileIO_MPICollective(const MPI_Comm &C) : GenericFileIO_MPI(C) {}

public:
  void read(void *buf, size_t count, off_t offset, const std::string &D);
  void write(const void *buf, size_t count, off_t offset, const std::string &D);
};
#endif

class GenericFileIO_POSIX : public GenericFileIO {
public:
  GenericFileIO_POSIX() : FH(-1) {}
  ~GenericFileIO_POSIX();

public:
  void open(const std::string &FN, bool ForReading = false, bool MustExist = false);
  void setSize(size_t sz);
  void read(void *buf, size_t count, off_t offset, const std::string &D);
  void write(const void *buf, size_t count, off_t offset, const std::string &D);

protected:
  int FH;
};

namespace detail {
// A standard enable_if idiom (we include our own here for pre-C++11 support).
template <bool B, typename T = void>
struct enable_if {};

template <typename T>
struct enable_if<true, T> { typedef T type; };

// A SFINAE-based trait to detect whether a type has a member named x. This is
// designed to work both with structs/classes and also with OpenCL-style vector
// types.
template <typename T>
class has_x {
  typedef char yes[1];
  typedef char no[2];

  template <typename C>
  static yes &test(char(*)[sizeof((*((C *) 0)).x)]);

  template <typename C>
  static no &test(...);

public:
  enum { value = sizeof(test<T>(0)) == sizeof(yes) };
};

// A SFINAE-based trait to detect whether a type is array-like (i.e. supports
// the [] operator).
template <typename T>
class is_array {
  typedef char yes[1];
  typedef char no[2];

  template <typename C>
  static yes &test(char(*)[sizeof((*((C *) 0))[0])]);

  template <typename C>
  static no &test(...);

public:
  enum { value = sizeof(test<T>(0)) == sizeof(yes) };
};
} // namespace detail

class GenericIO {
public:
  enum VariableFlags {
    VarHasExtraSpace =  (1 << 0), // Note that this flag indicates that the
                                  // extra space is available, but the GenericIO
                                  // implementation is required to
                                  // preserve its contents.
    VarIsPhysCoordX  =  (1 << 1),
    VarIsPhysCoordY  =  (1 << 2),
    VarIsPhysCoordZ  =  (1 << 3),
    VarMaybePhysGhost = (1 << 4)
  };

  struct VariableInfo {
    VariableInfo(const std::string &N, std::size_t S, bool IF, bool IS,
                 bool PCX, bool PCY, bool PCZ, bool PG, std::size_t ES = 0)
      : Name(N), Size(S), IsFloat(IF), IsSigned(IS),
        IsPhysCoordX(PCX), IsPhysCoordY(PCY), IsPhysCoordZ(PCZ),
        MaybePhysGhost(PG), ElementSize(ES ? ES : S) {}

    std::string Name;
    std::size_t Size;
    bool IsFloat;
    bool IsSigned;
    bool IsPhysCoordX, IsPhysCoordY, IsPhysCoordZ;
    bool MaybePhysGhost;
    std::size_t ElementSize;
  };

public:
  struct LossyCompressionInfo {
    enum LCMode {
      LCModeNone,
      LCModeAbs,
      LCModeRel,
      LCModeAbsAndRel,
      LCModeAbsOrRel,
      LCModePSNR
    };

    LCMode Mode;
    double AbsErrThreshold;
    double RelErrThreshold;
    double PSNRThreshold;

    LossyCompressionInfo()
      : Mode(LCModeNone), AbsErrThreshold(0.0),
        RelErrThreshold(0.0), PSNRThreshold(0.0) {}
  };

  class Variable {
  private:
    template <typename ET>
    void deduceTypeInfoFromElement(ET *) {
      ElementSize = sizeof(ET);
      IsFloat = !std::numeric_limits<ET>::is_integer;
      IsSigned = std::numeric_limits<ET>::is_signed;
    }

    // There are specializations here to handle array types
    // (e.g. typedef float float4[4];), struct types
    // (e.g. struct float4 { float x, y, z, w; };), and scalar types.
    // Builtin vector types
    // (e.g. typedef float float4 __attribute__((ext_vector_type(4)));) should
    // also work.
    template <typename T>
    typename detail::enable_if<detail::is_array<T>::value, void>::type
    deduceTypeInfo(T *D) {
      Size = sizeof(T);
      deduceTypeInfoFromElement(&(*D)[0]);
    }

    template <typename T>
    typename detail::enable_if<detail::has_x<T>::value &&
                               !detail::is_array<T>::value, void>::type
    deduceTypeInfo(T *D) {
      Size = sizeof(T);
      deduceTypeInfoFromElement(&(*D).x);
    }

    template <typename T>
    typename detail::enable_if<!detail::has_x<T>::value &&
                               !detail::is_array<T>::value, void>::type
    deduceTypeInfo(T *D) {
      Size = sizeof(T);
      deduceTypeInfoFromElement(D);
    }

  public:
    template <typename T>
    Variable(const std::string &N, T* D, unsigned Flags = 0,
        const LossyCompressionInfo &LCI = LossyCompressionInfo())
      : Name(N), Data((void *) D), HasExtraSpace(Flags & VarHasExtraSpace),
        IsPhysCoordX(Flags & VarIsPhysCoordX),
        IsPhysCoordY(Flags & VarIsPhysCoordY),
        IsPhysCoordZ(Flags & VarIsPhysCoordZ),
        MaybePhysGhost(Flags & VarMaybePhysGhost),
        LCI(LCI) {
      deduceTypeInfo(D);
    }

    template <typename T>
    Variable(const std::string &N, std::size_t NumElements, T* D,
             unsigned Flags = 0,
             const LossyCompressionInfo &LCI = LossyCompressionInfo())
      : Name(N), Data((void *) D), HasExtraSpace(Flags & VarHasExtraSpace),
        IsPhysCoordX(Flags & VarIsPhysCoordX),
        IsPhysCoordY(Flags & VarIsPhysCoordY),
        IsPhysCoordZ(Flags & VarIsPhysCoordZ),
        MaybePhysGhost(Flags & VarMaybePhysGhost),
        LCI(LCI) {
      deduceTypeInfoFromElement(D);
      Size = ElementSize*NumElements;
    }

    Variable(const VariableInfo &VI, void *D, unsigned Flags = 0,
             const LossyCompressionInfo &LCI = LossyCompressionInfo())
      : Name(VI.Name), Size(VI.Size), IsFloat(VI.IsFloat),
        IsSigned(VI.IsSigned), Data(D),
        HasExtraSpace(Flags & VarHasExtraSpace),
        IsPhysCoordX((Flags & VarIsPhysCoordX) || VI.IsPhysCoordX),
        IsPhysCoordY((Flags & VarIsPhysCoordY) || VI.IsPhysCoordY),
        IsPhysCoordZ((Flags & VarIsPhysCoordZ) || VI.IsPhysCoordZ),
        MaybePhysGhost((Flags & VarMaybePhysGhost) || VI.MaybePhysGhost),
        ElementSize(VI.ElementSize), LCI(LCI) {}

    template <typename ET>
    bool hasElementType() {
      if (ElementSize != sizeof(ET))
        return false;
      if (IsFloat != !std::numeric_limits<ET>::is_integer)
        return false;
      if (IsSigned != std::numeric_limits<ET>::is_signed)
        return false;

      return true;
    }

    std::string Name;
    std::size_t Size;
    bool IsFloat;
    bool IsSigned;
    void *Data;
    bool HasExtraSpace;
    bool IsPhysCoordX, IsPhysCoordY, IsPhysCoordZ;
    bool MaybePhysGhost;
    std::size_t ElementSize;

    LossyCompressionInfo LCI;
  };

public:
  enum FileIO {
    FileIOMPI,
    FileIOPOSIX,
    FileIOMPICollective
  };

#ifndef GENERICIO_NO_MPI
  GenericIO(const MPI_Comm &C, const std::string &FN, unsigned FIOT = -1)
    : NElems(0), FileIOType(FIOT == (unsigned) -1 ? DefaultFileIOType : FIOT),
      Partition(DefaultPartition), Comm(C), FileName(FN), Redistributing(false),
      DisableCollErrChecking(false), SplitComm(MPI_COMM_NULL) {
    std::fill(PhysOrigin, PhysOrigin + 3, 0.0);
    std::fill(PhysScale,  PhysScale + 3, 0.0);
  }
#else
  GenericIO(const std::string &FN, unsigned FIOT = -1)
    : NElems(0), FileIOType(FIOT == (unsigned) -1 ? DefaultFileIOType : FIOT),
      Partition(DefaultPartition), FileName(FN), Redistributing(false),
      DisableCollErrChecking(false) {
    std::fill(PhysOrigin, PhysOrigin + 3, 0.0);
    std::fill(PhysScale,  PhysScale + 3, 0.0);
  }
#endif

  ~GenericIO() {
    close();

#ifndef GENERICIO_NO_MPI
    if (SplitComm != MPI_COMM_NULL)
      MPI_Comm_free(&SplitComm);
#endif
  }

public:
  static std::size_t requestedExtraSpace() {
    return 8;
  }

  void setNumElems(std::size_t E) {
    NElems = E;

#ifndef GENERICIO_NO_MPI
    int IsLarge = E >= CollectiveMPIIOThreshold;
    int AllIsLarge;
    MPI_Allreduce(&IsLarge, &AllIsLarge, 1, MPI_INT, MPI_SUM, Comm);
    if (!AllIsLarge)
      FileIOType = FileIOMPICollective;
#endif
  }

  void setPhysOrigin(double O, int Dim = -1) {
    if (Dim >= 0)
      PhysOrigin[Dim] = O;
    else
      std::fill(PhysOrigin, PhysOrigin + 3, O);
  }

  void setPhysScale(double S, int Dim = -1) {
    if (Dim >= 0)
      PhysScale[Dim] = S;
    else
      std::fill(PhysScale,  PhysScale + 3, S);
  }

  template <typename T>
  void addVariable(const std::string &Name, T *Data,
                   unsigned Flags = 0,
                   const LossyCompressionInfo &LCI = LossyCompressionInfo()) {
    Vars.push_back(Variable(Name, Data, Flags, LCI));
  }

  template <typename T, typename A>
  void addVariable(const std::string &Name, std::vector<T, A> &Data,
                   unsigned Flags = 0,
                   const LossyCompressionInfo &LCI = LossyCompressionInfo()) {
    T *D = Data.empty() ? 0 : &Data[0];
    addVariable(Name, D, Flags, LCI);
  }

  void addVariable(const VariableInfo &VI, void *Data,
                   unsigned Flags = 0,
                   const LossyCompressionInfo &LCI = LossyCompressionInfo()) {
    Vars.push_back(Variable(VI, Data, Flags, LCI));
  }

  template <typename T>
  void addScalarizedVariable(const std::string &Name, T *Data,
                             std::size_t NumElements, unsigned Flags = 0,
                             const LossyCompressionInfo &LCI = LossyCompressionInfo()) {
    Vars.push_back(Variable(Name, NumElements, Data, Flags, LCI));
  }

  template <typename T, typename A>
  void addScalarizedVariable(const std::string &Name, std::vector<T, A> &Data,
                             std::size_t NumElements, unsigned Flags = 0,
                             const LossyCompressionInfo &LCI = LossyCompressionInfo()) {
    T *D = Data.empty() ? 0 : &Data[0];
    addScalarizedVariable(Name, D, NumElements, Flags, LCI);
  }

#ifndef GENERICIO_NO_MPI
  // Writing
  void write();
#endif

  enum MismatchBehavior {
    MismatchAllowed,
    MismatchDisallowed,
    MismatchRedistribute
  };

  // Reading
  void openAndReadHeader(MismatchBehavior MB = MismatchDisallowed,
                         int EffRank = -1, bool CheckPartMap = true);

  int readNRanks();
  void readDims(int Dims[3]);

  // Note: For partitioned inputs, this returns -1.
  uint64_t readTotalNumElems();

  void readPhysOrigin(double Origin[3]);
  void readPhysScale(double Scale[3]);

  void clearVariables() { this->Vars.clear(); };

  int getNumberOfVariables() { return this->Vars.size(); };

  void getVariableInfo(std::vector<VariableInfo> &VI);

  std::size_t readNumElems(int EffRank = -1);
  void readCoords(int Coords[3], int EffRank = -1);
  int readGlobalRankNumber(int EffRank = -1);

  void readData(int EffRank = -1, bool PrintStats = true, bool CollStats = true);

  void getSourceRanks(std::vector<int> &SR);

  void close() {
    FH.close();
  }

  void setPartition(int P) {
    Partition = P;
  }

  static void setDefaultFileIOType(unsigned FIOT) {
    DefaultFileIOType = FIOT;
  }

  static void setDefaultPartition(int P) {
    DefaultPartition = P;
  }

  static void setNaturalDefaultPartition();

  static void setDefaultShouldCompress(bool C) {
    DefaultShouldCompress = C;
  }

#ifndef GENERICIO_NO_MPI
  static void setCollectiveMPIIOThreshold(std::size_t T) {
#ifndef GENERICIO_NO_NEVER_USE_COLLECTIVE_IO
    CollectiveMPIIOThreshold = T;
#endif
  }
#endif

private:
  // Implementation functions templated on the Endianness of the underlying
  // data.

#ifndef GENERICIO_NO_MPI
  template <bool IsBigEndian>
  void write();
#endif

  template <bool IsBigEndian>
  void readHeaderLeader(void *GHPtr, MismatchBehavior MB, int Rank, int NRanks,
                        int SplitNRanks, std::string &LocalFileName,
                        uint64_t &HeaderSize, std::vector<char> &Header);

  template <bool IsBigEndian>
  int readNRanks();

  template <bool IsBigEndian>
  void readDims(int Dims[3]);

  template <bool IsBigEndian>
  uint64_t readTotalNumElems();

  template <bool IsBigEndian>
  void readPhysOrigin(double Origin[3]);

  template <bool IsBigEndian>
  void readPhysScale(double Scale[3]);

  template <bool IsBigEndian>
  int readGlobalRankNumber(int EffRank);

  template <bool IsBigEndian>
  size_t readNumElems(int EffRank);

  template <bool IsBigEndian>
  void readCoords(int Coords[3], int EffRank);

  void readData(int EffRank, size_t RowOffset, int Rank,
                uint64_t &TotalReadSize, int NErrs[3]);

  template <bool IsBigEndian>
  void readData(int EffRank, size_t RowOffset,
                int Rank, uint64_t &TotalReadSize, int NErrs[3]);

  template <bool IsBigEndian>
  void getVariableInfo(std::vector<VariableInfo> &VI);

protected:
  std::vector<Variable> Vars;
  std::size_t NElems;

  double PhysOrigin[3], PhysScale[3];

  unsigned FileIOType;
  int Partition;
#ifndef GENERICIO_NO_MPI
  MPI_Comm Comm;
#endif
  std::string FileName;

  static unsigned DefaultFileIOType;
  static int DefaultPartition;
  static bool DefaultShouldCompress;

#ifndef GENERICIO_NO_MPI
  static std::size_t CollectiveMPIIOThreshold;
#endif

  // When redistributing, the rank blocks which this process should read.
  bool Redistributing, DisableCollErrChecking;
  std::vector<int> SourceRanks;

  std::vector<int> RankMap;
#ifndef GENERICIO_NO_MPI
  MPI_Comm SplitComm;
#endif
  std::string OpenFileName;

  // This reference counting mechanism allows the the GenericIO class
  // to be used in a cursor mode. To do this, make a copy of the class
  // after reading the header but prior to adding the variables.
  struct FHManager {
    FHManager() : CountedFH(0) {
      allocate();
    }

    FHManager(const FHManager& F) {
      CountedFH = F.CountedFH;
      CountedFH->Cnt += 1;
    }

    ~FHManager() {
      close();
    }

    GenericFileIO *&get() {
      if (!CountedFH)
        allocate();

      return CountedFH->GFIO;
    }

    std::vector<char> &getHeaderCache() {
      if (!CountedFH)
        allocate();

      return CountedFH->HeaderCache;
    }

    bool isBigEndian() {
      return CountedFH ? CountedFH->IsBigEndian : false;
    }

    void setIsBigEndian(bool isBE) {
      CountedFH->IsBigEndian = isBE;
    }

    void allocate() {
      close();
      CountedFH = new FHWCnt;
    };

    void close() {
      if (CountedFH && CountedFH->Cnt == 1)
        delete CountedFH;
      else if (CountedFH)
        CountedFH->Cnt -= 1;

      CountedFH = 0;
    }

    struct FHWCnt {
      FHWCnt() : GFIO(0), Cnt(1), IsBigEndian(false) {}

      ~FHWCnt() {
        close();
      }

protected:
      void close() {
        delete GFIO;
        GFIO = 0;
      }

public:
      GenericFileIO *GFIO;
      size_t Cnt;

      // Used for reading
      std::vector<char> HeaderCache;
      bool IsBigEndian;
    };

    FHWCnt *CountedFH;
  } FH;

  void setFH(
#ifndef GENERICIO_NO_MPI
    MPI_Comm R
#endif
  );
};

} /* END namespace cosmotk */
#endif // GENERICIO_H

