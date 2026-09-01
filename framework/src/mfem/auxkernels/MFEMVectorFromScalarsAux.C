//* This file is part of the MOOSE framework
//* https://mooseframework.inl.gov
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#ifdef MOOSE_MFEM_ENABLED

#include "MFEMVectorFromScalarsAux.h"
#include "MFEMProblem.h"

registerMooseObject("MooseApp", MFEMVectorFromScalarsAux);

InputParameters
MFEMVectorFromScalarsAux::validParams()
{
  InputParameters params = MFEMAuxKernel::validParams();
  params.addClassDescription(
      "Assembles a vector MFEM auxvariable from one scalar coefficient per component.");
  params.addRequiredParam<std::vector<MFEMScalarCoefficientName>>(
      "component_coefficients",
      "Names of the scalar coefficients giving each component of the vector, in order.");
  return params;
}

MFEMVectorFromScalarsAux::MFEMVectorFromScalarsAux(const InputParameters & parameters)
  : MFEMAuxKernel(parameters), _vec_coef(_result_var.ParFESpace()->GetVDim())
{
  const auto & names = getParam<std::vector<MFEMScalarCoefficientName>>("component_coefficients");
  const auto vdim = _result_var.ParFESpace()->GetVDim();

  if (names.size() != static_cast<std::size_t>(vdim))
    paramError("component_coefficients",
               "Expected ",
               vdim,
               " coefficient names to match the vector dimension of variable '",
               _result_var_name,
               "', but ",
               names.size(),
               " were given.");

  for (const auto i : index_range(names))
    // The coefficients belong to the problem's coefficient manager, so ownership is not taken here
    _vec_coef.Set(i, &getScalarCoefficientByName(names[i]), false);
}

void
MFEMVectorFromScalarsAux::execute()
{
  _result_var.ProjectCoefficient(_vec_coef);
}

#endif // MOOSE_MFEM_ENABLED
