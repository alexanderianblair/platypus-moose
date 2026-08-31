//* This file is part of the MOOSE framework
//* https://mooseframework.inl.gov
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#ifdef MOOSE_MFEM_ENABLED

#include "MFEMScalarVariableValue.h"
#include "MFEMScalarVariable.h"
#include "MFEMProblem.h"

registerMooseObject("MooseApp", MFEMScalarVariableValue);

InputParameters
MFEMScalarVariableValue::validParams()
{
  InputParameters params = MFEMPostprocessor::validParams();
  params.addClassDescription("Returns the value of an MFEM scalar variable.");
  params.addRequiredParam<VariableName>("variable", "Name of the MFEM scalar variable.");
  return params;
}

MFEMScalarVariableValue::MFEMScalarVariableValue(const InputParameters & parameters)
  : MFEMPostprocessor(parameters),
    _scalar_var(getMFEMProblem().getMFEMScalarVariable(getParam<VariableName>("variable")))
{
}

PostprocessorValue
MFEMScalarVariableValue::getValue() const
{
  return _scalar_var.getValue();
}

#endif
