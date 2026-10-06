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
  params.addParam<bool>(
      "implicit",
      true,
      "Whether this boundary condition is treated implicitly by implicit-explicit "
      "(IMEX) time integration schemes. Ignored by fully implicit and fully "
      "explicit schemes.");
  params.addParamNamesToGroup("implicit", "Advanced");
  return params;
}

MFEMIntegratedBC::MFEMIntegratedBC(const InputParameters & parameters)
  : MFEMBoundaryCondition(parameters), _is_implicit(getParam<bool>("implicit"))
{
}

#endif
