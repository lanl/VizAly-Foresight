/**
 * @brief Base class for Generic I/O that encapsulates common functionality
 * for readers and writers.
 */
#ifndef GENERICIOBASE_H_
#define GENERICIOBASE_H_

#include "GenericIODefinitions.hpp"

// STL includes
#include <vector>
#include <limits>

namespace gio
{

class GenericIOBase
{
public:

  enum VariableFlags {
    FloatValue          = (1 << 0),
    SignedValue         = (1 << 1),
    ValueIsPhysCoordX   = (1 << 2),
    ValueIsPhysCoordY   = (1 << 3),
    ValueIsPhysCoordZ   = (1 << 4),
    ValueMaybePhysGhost = (1 << 5),
    ValueHasExtraSpace  = (1 << 6)
  };

  enum FileIOStrategy {
    FileIOMPI,
    FileIOPOSIX,
    FileIOMPICollective,
    FileIOUndefined
  };

  GenericIOBase();
  virtual ~GenericIOBase();

  // Get and Set macros
  GIOGetNSetMacro(FileName,std::string);
  GIOGetMacro(IOStrategy,unsigned);

  /**
   * @brief Adds a variable associated with given name and data array.
   * @param Name the name of the variable
   * @param Data the user-supplied array consisting of the data
   * @param Flags bit mask associated with this variable
   */
  template <typename T>
  void AddVariable(const std::string &Name, T *Data, unsigned Flags=0)
    {this->Vars.push_back(Variable(Name,Data,Flags)); }

  /**
   * @brief Adds a variable associated with the given name and STL vector.
   * @param Name the name of the variable
   * @param Data the STL vector with the data associated with the variable.
   * @param Flags bit mask associated with this variable
   */
  template <typename T, typename A>
  void AddVariable(const std::string &Name, std::vector<T,A> &Data,
                      unsigned Flags=0)
    {
    T *D = Data.empty() ? NULL : &Data[0];
    this->Vars.push_back( Variable(Name,D,Flags) );
    }

  /**
   * @brief Constructs and returns a variable info object for variable
   * corresponding to the given, user-supplied variable index.
   * @param variableIdx the index of the variable in query.
   * @return VI the variable info object
   * @see VariableInfo in GenericIODefinitions.hpp
   * @see GenericIOReader::GetFileVariableInfo()
   * @note This returns the information of a variable that the user has
   * registered via calls to AddVariable() -- To get the information of a
   * variable in a given file, use the GetFileVariableInfo() method of
   * a reader.
   * @pre ( variableIdx >= 0 ) && ( variableIdx < this->Vars.size() )
   */
  VariableInfo GetVariableInfo(const int variableIdx);

  /**
   * @brief Adds a variable associated with the given user-supplied
   * VariableInfo object and the data pointed by the void* pointer.
   * @param VI the user-supplied variable information object
   * @param Data pointer to the data corresponding to the data
   * @param Flags bit mask associated with this variable
   * @see VariableInfo
   */
  void AddVariable(const VariableInfo &VI, void *Data, unsigned Flags=0)
    { this->Vars.push_back( Variable(VI,Data,Flags) ); }

  /**
   * @brief Return the number of variables associated with this instance.
   * @return N the number of variable associated with this instance.
   */
  int GetNumberOfVariables()
    { return this->Vars.size(); }

  /**
   * @brief Clears all variables associated with this instance.
   * @post this->GetNumberOfVariables()==0.
   */
  void ClearVariables()
    { this->Vars.clear(); }

protected:

  struct Variable
  {
    template <typename T>
    Variable(const std::string &N, T* D, unsigned Flags = 0)
      : Name(N), Size(sizeof(T)),
        IsFloat(!std::numeric_limits<T>::is_integer),
        IsSigned(std::numeric_limits<T>::is_signed),
        Data((void *) D), HasExtraSpace(Flags & ValueHasExtraSpace),
        IsPhysCoordX(Flags & ValueIsPhysCoordX),
        IsPhysCoordY(Flags & ValueIsPhysCoordY),
        IsPhysCoordZ(Flags & ValueIsPhysCoordZ),
        MaybePhysGhost(Flags & ValueMaybePhysGhost) {}

    Variable(const VariableInfo &VI, void *D, unsigned Flags = 0)
      : Name(VI.Name), Size(VI.Size), IsFloat(VI.IsFloat),
        IsSigned(VI.IsSigned), Data(D),
        HasExtraSpace(Flags & ValueHasExtraSpace),
        IsPhysCoordX((Flags & ValueIsPhysCoordX) || VI.IsPhysCoordX),
        IsPhysCoordY((Flags & ValueIsPhysCoordY) || VI.IsPhysCoordY),
        IsPhysCoordZ((Flags & ValueIsPhysCoordZ) || VI.IsPhysCoordZ),
        MaybePhysGhost((Flags & ValueMaybePhysGhost) || VI.MaybePhysGhost) {}

    std::string Name;
    std::size_t Size;
    bool IsFloat;
    bool IsSigned;
    void *Data;
    bool HasExtraSpace;
    bool IsPhysCoordX, IsPhysCoordY, IsPhysCoordZ;
    bool MaybePhysGhost;
  }; // END definition of Variable

  unsigned IOStrategy;
  std::string FileName;
  std::vector<Variable> Vars;

private:
  GIO_DISABLE_COPY_AND_ASSIGNMENT(GenericIOBase);
}; // END GenericIOBase definition


} /* namespace cosmotk */
#endif /* GENERICIOBASE_H_ */
