//* This file is part of the MOOSE framework
//* https://mooseframework.inl.gov
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#ifdef MOOSE_MFEM_ENABLED

#include "MFEMScalarVectorProductFunction.h"

registerMooseObject("MooseApp", MFEMScalarVectorProductFunction);

InputParameters
MFEMScalarVectorProductFunction::validParams()
{
  InputParameters params = Function::validParams();
  params.addClassDescription("Parses scalar function of position, time and scalar "
                             "problem coefficients (including scalar variables).");
  params.addRequiredParam<MFEMScalarCoefficientName>("coefficient",
                                                     "Name of the scalar coefficient to scale vector by.");
  params.addRequiredParam<MFEMVectorCoefficientName>("vector_coefficient",
                                                     "Name of the vector coefficient to be scaled.");                                                                         
  return params;
}

MFEMScalarVectorProductFunction::MFEMScalarVectorProductFunction(const InputParameters & parameters)
  : Function(parameters),
    _mfem_problem(static_cast<MFEMProblem &>(*getCheckedPointerParam<FEProblemBase *>("_fe_problem_base"))),
    _scalar_coefficient_name(getParam<MFEMScalarCoefficientName>("coefficient")),
    _vector_coefficient_name(getParam<MFEMVectorCoefficientName>("vector_coefficient")),
    _scaled_vector_coef(_mfem_problem.getCoefficients().declareVector<mfem::ScalarVectorProductCoefficient>(name(), _dummy_scalar_coefficient, _dummy_vector_coefficient))
{
}

void
MFEMScalarVectorProductFunction::initialSetup()
{
  // Set input coefficients. Must be done during initialSetup, after all coefficients are constructed.
  _scaled_vector_coef.SetACoef(_mfem_problem.getCoefficients().getScalarCoefficient(_scalar_coefficient_name));
  _scaled_vector_coef.SetBCoef(_mfem_problem.getCoefficients().getVectorCoefficient(_vector_coefficient_name));
}

#endif
