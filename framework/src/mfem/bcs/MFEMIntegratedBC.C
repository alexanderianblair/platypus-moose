//* This file is part of the MOOSE framework
//* https://mooseframework.inl.gov
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#ifdef MOOSE_MFEM_ENABLED

#include "MFEMIntegratedBC.h"

InputParameters
MFEMIntegratedBC::validParams()
{
  InputParameters params = MFEMBoundaryCondition::validParams();
  return params;
}

MFEMIntegratedBC::MFEMIntegratedBC(const InputParameters & parameters)
  : MFEMBoundaryCondition(parameters)
{
}

const std::vector<VariableName> &
MFEMIntegratedBC::getCoupledVariableNames() const
{
  static const std::vector<VariableName> no_coupled_variables;
  return no_coupled_variables;
}

#endif
