//* This file is part of the MOOSE framework
//* https://mooseframework.inl.gov
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#ifdef MOOSE_MFEM_ENABLED

#include "MFEMScalarVariable.h"
#include "MFEMProblem.h"
#include "MooseVariableBase.h"

registerMooseObject("MooseApp", MFEMScalarVariable);

InputParameters
MFEMScalarVariable::validParams()
{
  InputParameters params = MFEMObject::validParams();
  // Require moose variable parameters (not used!), as they are read by AddVariableAction.
  params += MooseVariableBase::validParams();
  params.addClassDescription(
      "Class for adding a global scalar unknown, holding a single degree of freedom, to an MFEM "
      "problem.");
  params.registerBase("MooseVariableBase");
  params.registerSystemAttributeName("MooseVariableBase");
  return params;
}

MFEMScalarVariable::MFEMScalarVariable(const InputParameters & parameters)
  : MFEMObject(parameters), _owns_dof(getMFEMProblem().getProblemData().myid == 0)
{
}

#endif
