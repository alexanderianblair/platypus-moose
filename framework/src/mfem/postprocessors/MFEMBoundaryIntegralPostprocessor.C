//* This file is part of the MOOSE framework
//* https://mooseframework.inl.gov
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#ifdef MOOSE_MFEM_ENABLED

#include "MFEMBoundaryIntegralPostprocessor.h"
#include "MFEMProblem.h"

registerMooseObject("MooseApp", MFEMBoundaryIntegralPostprocessor);

InputParameters
MFEMBoundaryIntegralPostprocessor::validParams()
{
  InputParameters params = MFEMPostprocessor::validParams();
  params += MFEMBoundaryRestrictable::validParams();
  params.addClassDescription("Calculates the integral of a scalar coefficient over a set of "
                             "boundaries, which may lie in the interior of the mesh.");
  params.addRequiredParam<MFEMScalarCoefficientName>("coefficient",
                                                     "Name of the coefficient to integrate.");
  params.addParam<unsigned int>(
      "quadrature_order",
      4,
      "Order of the quadrature rule on the boundary elements, which integrates polynomials of "
      "this degree exactly.");
  return params;
}

MFEMBoundaryIntegralPostprocessor::MFEMBoundaryIntegralPostprocessor(
    const InputParameters & parameters)
  : MFEMPostprocessor(parameters),
    MFEMBoundaryRestrictable(parameters, getMFEMProblem().mesh().getMFEMParMesh()),
    _coef(getScalarCoefficient("coefficient")),
    _h1_fec(1, getMesh().Dimension()),
    _h1_fespace(const_cast<mfem::ParMesh *>(&getMesh()), &_h1_fec),
    _boundary_integrator(&_h1_fespace),
    _unity(&_h1_fespace)
{
  _unity = 1.0;
  // The quadrature order is oa * (order of the H1 space) + ob, so oa = 0 sets it directly.
  auto * integ =
      new mfem::BoundaryLFIntegrator(_coef, 0, int(getParam<unsigned int>("quadrature_order")));
  if (isBoundaryRestricted())
    _boundary_integrator.AddBoundaryIntegrator(integ, getBoundaryMarkers());
  else
    _boundary_integrator.AddBoundaryIntegrator(integ);
}

void
MFEMBoundaryIntegralPostprocessor::execute()
{
  _boundary_integrator.Assemble();
  // The shape functions of the H1 space sum to one, so this is the integral; the parallel
  // inner product sums the contributions of every rank.
  _integral = _boundary_integrator(_unity);
}

PostprocessorValue
MFEMBoundaryIntegralPostprocessor::getValue() const
{
  return _integral;
}

#endif
