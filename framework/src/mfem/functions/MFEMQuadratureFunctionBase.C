//* This file is part of the MOOSE framework
//* https://mooseframework.inl.gov
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#ifdef MOOSE_MFEM_ENABLED

#include "MFEMQuadratureFunctionBase.h"

namespace
{
mfem::QuadratureSpace
makeQuadratureSpace(mfem::ParMesh & mesh, const int order)
{
  // MFEM's device integrators build their quadrature space from an IntegrationRule, recording the
  // rule's exact order (e.g. 3 for the 2-point Gauss rule returned for order 2), and require the
  // order of any quadrature function they read to match. Construct from the same rule so the
  // points and recorded order agree. That constructor rejects mixed-geometry meshes, which device
  // assembly does not support, so those keep the order-based constructor.
  if (mesh.GetNumGeometries(mesh.Dimension()) <= 1)
    return mfem::QuadratureSpace(mesh, mfem::IntRules.Get(mesh.GetTypicalElementGeometry(), order));
  return mfem::QuadratureSpace(&mesh, order);
}
}

InputParameters
MFEMQuadratureFunctionBase::validParams()
{
  InputParameters params = Function::validParams();
  params.addRequiredRangeCheckedParam<int>(
      "order",
      "order>=0",
      "Order of the quadrature rule the projected values are stored on. This must match the "
      "integration rule used by the objects consuming this coefficient.");
  MooseEnum updates(MFEMQuadratureFunctionCoefficientBase::getUpdatePolicyOptions(), "NONLINEAR");
  params.addParam<MooseEnum>(
      "updates",
      updates,
      "When the stored values are re-projected from the source coefficient: 'none' projects "
      "exactly once, 'time' re-projects when the simulation time changes, 'nonlinear' additionally "
      "re-projects whenever solution variables change (i.e. on each nonlinear iteration).");
  return params;
}

MFEMQuadratureFunctionBase::MFEMQuadratureFunctionBase(const InputParameters & parameters)
  : Function(parameters),
    _mfem_problem(
        cast_ref<MFEMProblem &>(*parameters.getCheckedPointerParam<SubProblem *>("_subproblem"))),
    _qspace(makeQuadratureSpace(_mfem_problem.mesh().getMFEMParMesh(), getParam<int>("order")))
{
}

MFEMQuadratureFunctionCoefficientBase::UpdatePolicy
MFEMQuadratureFunctionBase::updatePolicy() const
{
  return getParam<MooseEnum>("updates")
      .getEnum<MFEMQuadratureFunctionCoefficientBase::UpdatePolicy>();
}

#endif
