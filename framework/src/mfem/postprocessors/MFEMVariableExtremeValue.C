//* This file is part of the MOOSE framework
//* https://mooseframework.inl.gov
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#ifdef MOOSE_MFEM_ENABLED

#include "MFEMVariableExtremeValue.h"
#include "MFEMProblem.h"

registerMooseObject("MooseApp", MFEMVariableExtremeValue);

InputParameters
MFEMVariableExtremeValue::validParams()
{
  InputParameters params = MFEMPostprocessor::validParams();
  params += MFEMBlockRestrictable::validParams();
  MFEMExecutedObject::addRequiredDependencyParam<VariableName>(
      params, "variable", "Name of the scalar variable to find the extreme value of.");
  params.addParam<MooseEnum>("value_type",
                             MooseEnum("max min", "max"),
                             "Whether to find the maximum or the minimum value.");
  params.addClassDescription("Finds the extreme value taken by a scalar MFEM variable over its "
                             "degrees of freedom.");
  return params;
}

MFEMVariableExtremeValue::MFEMVariableExtremeValue(const InputParameters & parameters)
  : MFEMPostprocessor(parameters),
    MFEMBlockRestrictable(parameters,
                          getMFEMProblem().getMFEMVariableMesh(getParam<VariableName>("variable"))),
    _var(*getMFEMProblem().getGridFunction(getParam<VariableName>("variable"))),
    _is_max(getParam<MooseEnum>("value_type") == "max")
{
  if (_var.ParFESpace()->GetVDim() != 1)
    paramError("variable", "Only scalar MFEM variables have an extreme value.");
}

void
MFEMVariableExtremeValue::execute()
{
  const auto & fespace = *_var.ParFESpace();

  // Seeded with the neutral element of the reduction so that ranks holding no
  // elements of the restricted subdomains do not contribute.
  mfem::real_t extreme = _is_max ? std::numeric_limits<mfem::real_t>::lowest()
                                 : std::numeric_limits<mfem::real_t>::max();

  mfem::Array<int> dofs;
  mfem::Vector values;
  for (const auto i : make_range(fespace.GetNE()))
  {
    if (isSubdomainRestricted() && !getSubdomainMarkers()[fespace.GetAttribute(i) - 1])
      continue;

    fespace.GetElementDofs(i, dofs);
    _var.GetSubVector(dofs, values);
    for (const auto j : make_range(values.Size()))
      extreme = _is_max ? std::max(extreme, values(j)) : std::min(extreme, values(j));
  }

  if (_is_max)
    _communicator.max(extreme);
  else
    _communicator.min(extreme);

  _value = extreme;
}

PostprocessorValue
MFEMVariableExtremeValue::getValue() const
{
  return _value;
}

#endif
