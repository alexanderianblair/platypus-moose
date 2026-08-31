//* This file is part of the MOOSE framework
//* https://mooseframework.inl.gov
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#ifdef MOOSE_MFEM_ENABLED

#include "MFEMIntegralConstraint.h"
#include "MFEMScalarVariable.h"
#include "MFEMProblem.h"

InputParameters
MFEMIntegralConstraint::validParams()
{
  InputParameters params = MFEMConstraint::validParams();
  params.addClassDescription(
      "Base class for weakly constraining a scalar integral quantity of an MFEM variable.");
  params.addRequiredParam<VariableName>(
      "scalar_variable",
      "Name of the MFEMScalarVariable holding the multiplier of this constraint.");
  return params;
}

MFEMIntegralConstraint::MFEMIntegralConstraint(const InputParameters & parameters)
  : MFEMConstraint(parameters),
    _scalar_var_name(getParam<VariableName>("scalar_variable")),
    _scalar_var(getMFEMProblem().getMFEMScalarVariable(_scalar_var_name))
{
}

#endif
