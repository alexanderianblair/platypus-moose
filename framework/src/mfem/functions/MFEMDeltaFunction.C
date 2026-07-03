//* This file is part of the MOOSE framework
//* https://mooseframework.inl.gov
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#ifdef MOOSE_MFEM_ENABLED

#include "MFEMDeltaFunction.h"

registerMooseObject("MooseApp", MFEMDeltaFunction);

InputParameters
MFEMDeltaFunction::validParams()
{
  InputParameters params = Function::validParams();
  params.addClassDescription("Declare a delta function of position.");
  params.addParam<mfem::real_t>("position_x", 0., "x-coordinate of delta function location.");
  params.addParam<mfem::real_t>("position_y", 0., "y-coordinate of delta function location.");
  params.addParam<mfem::real_t>("position_z", 0., "z-coordinate of delta function location.");
  params.addParam<mfem::real_t>("scale", 1., "Scale factor for delta function.");
  return params;
}

MFEMDeltaFunction::MFEMDeltaFunction(const InputParameters & parameters)
  : Function(parameters),
    _mfem_problem(static_cast<MFEMProblem &>(
        *parameters.getCheckedPointerParam<FEProblemBase *>("_fe_problem_base"))),
    _x_center(getParam<mfem::real_t>("position_x")),
    _y_center(getParam<mfem::real_t>("position_y")),
    _z_center(getParam<mfem::real_t>("position_z")),
    _scale(getParam<mfem::real_t>("scale"))
{
  // create MFEMDeltaCoefficient
  _mfem_problem.getCoefficients().declareScalar<mfem::DeltaCoefficient>(
      name(), _x_center, _y_center, _scale);
}

// Constraint - single real_t (lambda) appended to equation system
// constraint matrix C (one row, [1,...,n]) acts on vector of dofs
// transpose C^T (one column, acts on lambda to control summed contribution of all dofs in
// constraint)

// effectively, adding a single dof to the problem
// operator is a bilinearform but one var is one dof in size
#endif
