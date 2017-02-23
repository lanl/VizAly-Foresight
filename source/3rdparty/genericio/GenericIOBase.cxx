#include "GenericIOBase.h"

#include <cassert>

namespace gio
{

GenericIOBase::GenericIOBase()
{
  this->FileName   = "";
  this->IOStrategy = FileIOUndefined;
}

//------------------------------------------------------------------------------
GenericIOBase::~GenericIOBase()
{
  // TODO Auto-generated destructor stub
}

//------------------------------------------------------------------------------
VariableInfo GenericIOBase::GetVariableInfo(const int variableIdx)
{
  assert("pre: variable index is out-of-bounds!" &&
            (variableIdx >= 0) &&
            (variableIdx < static_cast<int>(this->Vars.size()) ) );

  VariableInfo VI(
      this->Vars[variableIdx].Name,
      this->Vars[variableIdx].Size,
      this->Vars[variableIdx].IsFloat,
      this->Vars[variableIdx].IsSigned,
      this->Vars[variableIdx].IsPhysCoordX,
      this->Vars[variableIdx].IsPhysCoordY,
      this->Vars[variableIdx].IsPhysCoordZ,
      this->Vars[variableIdx].MaybePhysGhost
      );
  return( VI );
}

} /* namespace cosmotk */
