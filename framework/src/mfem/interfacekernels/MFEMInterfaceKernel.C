//* This file is part of the MOOSE framework
//* https://mooseframework.inl.gov
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#ifdef MOOSE_MFEM_ENABLED

#include "MFEMInterfaceKernel.h"
#include "MFEMProblem.h"

InputParameters
MFEMInterfaceKernel::validParams()
{
  InputParameters params = MFEMObject::validParams();
  params += MFEMBoundaryRestrictable::validParams();
  params.addClassDescription("Base class for MFEM interface kernels, which couple variables on "
                             "either side of an interior boundary.");
  params.registerBase("InterfaceKernel");
  params.registerSystemAttributeName("InterfaceKernel");
  params.addRequiredParam<VariableName>("variable",
                                        "Variable on the primary side of the interface");
  params.addRequiredParam<VariableName>("neighbor_var",
                                        "Variable on the neighbor side of the interface");
  return params;
}

MFEMInterfaceKernel::MFEMInterfaceKernel(const InputParameters & parameters)
  : MFEMObject(parameters),
    MFEMBoundaryRestrictable(
        parameters, getMFEMProblem().getMFEMVariableMesh(getParam<VariableName>("variable"))),
    _var_name(getParam<VariableName>("variable")),
    _neighbor_var_name(getParam<VariableName>("neighbor_var"))
{
  if (!isBoundaryRestricted())
    paramError("boundary", "An interface kernel must be restricted to the interface boundaries.");
  if (_var_name == _neighbor_var_name)
    paramError("neighbor_var",
               "must differ from 'variable': a variable continuous across the interface cannot "
               "jump across it.");
  if (&getMFEMProblem().getMFEMVariableMesh(_neighbor_var_name) != &getMesh())
    paramError("neighbor_var", "must be defined on the same mesh as 'variable'.");
}

#endif
